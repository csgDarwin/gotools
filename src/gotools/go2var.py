#!/usr/bin/env python3
"""
go2var: Optimized integrated Sort + Variant calling pipeline.

-o / --output takes a FULL output path (directory + filename). Parent
directories are created automatically.

Combines custom MAF parsing with NumPy vectorized variant detection
and species filtering/reordering - all in a single multiprocessed pipeline.

"""

import argparse
import gzip
import json
import sys
import time
import threading
import queue
import multiprocessing as mp
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Iterator
from functools import partial
import numpy as np

# Pre-compute uppercase ASCII mapping table
_UPPER_TABLE = np.zeros(256, dtype=np.uint8)
for i in range(256):
    _UPPER_TABLE[i] = ord(chr(i).upper()) if i < 128 else i

_GAP = ord('-')


# =============================================================================
# Configuration
# =============================================================================

def read_json_lists(json_file: str) -> Tuple[List[str], List[str], int]:
    """Read reference and alignment order lists from JSON config."""
    with open(json_file, 'r') as f:
        data = json.load(f)

    reference_list = data.get("ReferenceList", [])
    align_list = data.get("AlignOrderList", [])

    if not reference_list:
        raise ValueError("ReferenceList is empty or missing in JSON file")
    if not align_list:
        raise ValueError("AlignOrderList is empty or missing in JSON file")

    return reference_list, align_list, len(reference_list)


# =============================================================================
# Custom MAF Parser (No BioPython)
# =============================================================================

def parse_maf_block_custom(block_lines: List[str]) -> Optional[Tuple[Dict[str, Tuple[str, int, bytes]], str]]:
    """
    Parse a single MAF block into species -> (id, start, sequence) dictionary.

    Returns:
        Tuple of (species_dict, first_species) or None if invalid
        species_dict maps species_name -> (full_id, start_position, sequence_bytes)
    """
    species_dict: Dict[str, Tuple[str, int, bytes]] = {}
    first_species = None

    for line in block_lines:
        if not line.startswith('s'):
            continue

        parts = line.split()
        if len(parts) < 7:
            continue

        # s src start size strand srcSize text
        src = parts[1]
        start = int(parts[2])
        seq_text = parts[6]

        # Extract species name (before first dot)
        dot_idx = src.find('.')
        if dot_idx > 0:
            species = src[:dot_idx]
        else:
            species = src

        # Check for duplicate species
        if species in species_dict:
            return None

        species_dict[species] = (src, start, seq_text.upper().encode('ascii'))

        if first_species is None:
            first_species = species

    if not species_dict:
        return None

    return species_dict, first_species


def read_maf_blocks_fast(filepath: str, block_size: int = 500) -> Iterator[List[Tuple[int, List[str]]]]:
    """
    Fast MAF block reader that yields batches of (index, block_lines) tuples.
    """
    is_gzipped = filepath.endswith('.gz')
    file_opener = gzip.open if is_gzipped else open
    open_mode = 'rt' if is_gzipped else 'r'

    current_block: List[str] = []
    batch: List[Tuple[int, List[str]]] = []
    global_index = 0
    in_block = False

    with file_opener(filepath, open_mode) as f:
        for line in f:
            if line.startswith('a'):
                # New block - save previous
                if current_block:
                    batch.append((global_index, current_block))
                    global_index += 1
                    if len(batch) >= block_size:
                        yield batch
                        batch = []
                current_block = [line]
                in_block = True

            elif line.startswith('s') and in_block:
                current_block.append(line)

            elif line.strip() == '' and in_block and current_block:
                # End of block
                batch.append((global_index, current_block))
                global_index += 1
                if len(batch) >= block_size:
                    yield batch
                    batch = []
                current_block = []
                in_block = False

        # Handle last block
        if current_block:
            batch.append((global_index, current_block))
        if batch:
            yield batch


# =============================================================================
# NumPy Vectorized Variant Detection
# =============================================================================

def find_variants_numpy(sequences: List[bytes], r_num: int, skip_gaps: bool) -> np.ndarray:
    """
    Find variant positions using vectorized NumPy operations.

    Returns boolean array where True = variant position.
    """
    if len(sequences) < r_num or not sequences:
        return np.array([], dtype=bool)

    seq_len = len(sequences[0])
    n_seqs = len(sequences)

    # Stack sequences into 2D array
    seq_array = np.frombuffer(b''.join(sequences), dtype=np.uint8).reshape(n_seqs, seq_len)

    # Uppercase (should already be, but ensure)
    seq_array = _UPPER_TABLE[seq_array]

    # Reference sequences
    ref_seqs = seq_array[:r_num, :]
    ref_base = ref_seqs[0, :]

    # Check reference consensus (all refs must match first ref)
    ref_consensus = np.all(ref_seqs == ref_base, axis=0)

    # Handle gaps
    first_has_gap = (ref_base == _GAP)
    valid_positions = ~first_has_gap if skip_gaps else np.ones(seq_len, dtype=bool)

    # Check if any non-ref matches reference base
    if n_seqs > r_num:
        non_ref_seqs = seq_array[r_num:, :]
        any_match = np.any(non_ref_seqs == ref_base, axis=0)
        # Variant = refs agree AND no non-ref matches AND valid position
        variants = ref_consensus & ~any_match & valid_positions
    else:
        # No non-ref sequences: "no non-ref matches ref" is vacuously true,
        # so every position where refs agree is a variant (matches v4 behavior)
        variants = ref_consensus & valid_positions

    return variants


# =============================================================================
# Integrated Block Processor
# =============================================================================

def process_block_optimized(
    indexed_block: Tuple[int, List[str]],
    primate_set: frozenset,
    align_order: List[str],
    ref_threshold: int,
    skip_gaps: bool,
    verbose: bool
) -> Tuple[int, Optional[List[str]]]:
    """
    Optimized integrated processing: Parse, Filter, Reorder, Detect variants.

    All done with custom parsing and NumPy - no BioPython.
    """
    block_index, block_lines = indexed_block
    r_num = ref_threshold

    # --- Step 1: Parse block with custom parser ---
    parsed = parse_maf_block_custom(block_lines)
    if parsed is None:
        return (block_index, None)

    species_dict, first_species = parsed

    # --- Step 2: Check minimum reference species ---
    species_in_block = set(species_dict.keys())
    primate_count = len(species_in_block & primate_set)
    if primate_count < ref_threshold:
        return (block_index, None)

    # --- Step 3: Reorder sequences according to align_order ---
    reordered_ids: List[str] = []
    reordered_seqs: List[bytes] = []
    ref_start = 0
    ref_id = ""

    for i, species in enumerate(align_order):
        if species in species_dict:
            full_id, start, seq_bytes = species_dict[species]

            # Use MAF `src` unchanged. It already has the species.chrom form,
            # so we do not prepend species (which would duplicate the prefix
            # for non-hg38 refs, e.g. panTro6.panTro6.chr1).
            new_id = full_id

            # First sequence in reordered list is the reference
            if not reordered_ids:
                ref_id = new_id
                ref_start = start

            reordered_ids.append(new_id)
            reordered_seqs.append(seq_bytes)

    if not reordered_seqs:
        return (block_index, None)

    # --- Step 4: Variant detection using NumPy ---
    variants = find_variants_numpy(reordered_seqs, r_num, skip_gaps)

    if len(variants) == 0:
        return (block_index, [])

    # --- Step 5: Generate output lines ---
    first_seq = reordered_seqs[0]
    is_gap = np.frombuffer(first_seq, dtype=np.uint8) == _GAP

    # Calculate genomic positions
    non_gap_cumsum = np.cumsum(~is_gap)

    variant_lines = []
    variant_indices = np.where(variants)[0]

    num_seqs = len(reordered_seqs)

    for idx in variant_indices:
        # Calculate genomic position
        if is_gap[idx]:
            if idx > 0:
                seq_pos = ref_start + non_gap_cumsum[idx - 1] - 1
            else:
                seq_pos = ref_start - 1
            start_pos = max(0, seq_pos)
            end_pos = seq_pos + 1
        else:
            if idx > 0:
                seq_pos = ref_start + non_gap_cumsum[idx] - 1
            else:
                seq_pos = ref_start
            start_pos = seq_pos
            end_pos = seq_pos + 1

        if verbose:
            column = bytes([s[idx] for s in reordered_seqs]).decode('ascii')
            if num_seqs > r_num:
                subst = f"{chr(reordered_seqs[r_num][idx])}>{chr(reordered_seqs[0][idx])}"
            else:
                subst = f" >{chr(reordered_seqs[0][idx])}"
            variant_lines.append(f"{ref_id}\t{start_pos}\t{end_pos}\t{subst}\t1\t{column}\n")
        else:
            variant_lines.append(f"{ref_id}\t{start_pos}\t{end_pos}\n")

    return (block_index, variant_lines)


# =============================================================================
# Writer Thread
# =============================================================================

def writer_thread_optimized(out_filename: str, write_queue: queue.Queue, merge_output: bool = False) -> None:
    """Optimized writer thread with large buffers and batch writes."""
    result_buffer: Dict[int, List[str]] = {}
    next_expected = 0

    try:
        with open(out_filename, "w", buffering=8*1024*1024) as outfile:  # 8MB buffer
            if merge_output:
                current_chrom = None
                current_start = None
                current_end = None

                while True:
                    result_index, variant_lines = write_queue.get()

                    if result_index is None:
                        write_queue.task_done()
                        break

                    result_buffer[result_index] = variant_lines if variant_lines else []

                    while next_expected in result_buffer:
                        lines = result_buffer.pop(next_expected)
                        for line in lines:
                            if line.startswith("#"):
                                continue
                            fields = line.split("\t", 3)
                            if len(fields) < 3:
                                continue

                            chrom = fields[0]
                            start = int(fields[1])
                            end = int(fields[2])

                            if current_chrom == chrom and current_end == start:
                                current_end = end
                            else:
                                if current_chrom is not None:
                                    outfile.write(f"{current_chrom}\t{current_start}\t{current_end}\n")
                                current_chrom = chrom
                                current_start = start
                                current_end = end

                        next_expected += 1

                    write_queue.task_done()

                # Flush remaining
                while next_expected in result_buffer:
                    lines = result_buffer.pop(next_expected)
                    for line in lines:
                        if line.startswith("#"):
                            continue
                        fields = line.split("\t", 3)
                        if len(fields) < 3:
                            continue
                        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                        if current_chrom == chrom and current_end == start:
                            current_end = end
                        else:
                            if current_chrom is not None:
                                outfile.write(f"{current_chrom}\t{current_start}\t{current_end}\n")
                            current_chrom, current_start, current_end = chrom, start, end
                    next_expected += 1

                if current_chrom is not None:
                    outfile.write(f"{current_chrom}\t{current_start}\t{current_end}\n")

            else:
                # Standard mode with batched writes
                write_buffer: List[str] = []
                buffer_bytes = 0
                max_buffer = 16 * 1024 * 1024  # 16MB

                while True:
                    result_index, variant_lines = write_queue.get()

                    if result_index is None:
                        write_queue.task_done()
                        break

                    result_buffer[result_index] = variant_lines if variant_lines else []

                    while next_expected in result_buffer:
                        lines = result_buffer.pop(next_expected)
                        write_buffer.extend(lines)
                        buffer_bytes += sum(len(l) for l in lines)
                        next_expected += 1

                        if buffer_bytes >= max_buffer:
                            outfile.writelines(write_buffer)
                            write_buffer = []
                            buffer_bytes = 0

                    write_queue.task_done()

                # Flush remaining
                while next_expected in result_buffer:
                    write_buffer.extend(result_buffer.pop(next_expected))
                    next_expected += 1

                if write_buffer:
                    outfile.writelines(write_buffer)

    except Exception as e:
        print(f"Fatal error in writer: {e}", file=sys.stderr)
        raise


# =============================================================================
# Main Pipeline
# =============================================================================

def default_output_filename(in_file: str) -> str:
    """Generate default output filename from input file.

    Removes .maf.gz or .maf extension and adds .bed
    """
    base = in_file
    if base.endswith('.maf.gz'):
        base = base[:-7]
    elif base.endswith('.gz'):
        base = base[:-3]
    if base.endswith('.maf'):
        base = base[:-4]
    return base + ".bed"


def go2var_sorted_optimized(
    in_file: str,
    json_file: str,
    out_file: str = None,
    skip_gaps: bool = True,
    num_workers: int = None,
    block_size: int = 500,
    verbose: bool = False,
    merge_output: bool = False
) -> None:
    """
    Highly optimized integrated Go2Sort + Go2Var pipeline.

    Uses custom MAF parsing and NumPy vectorized operations.
    """
    if num_workers is None:
        num_workers = min(8, max(1, mp.cpu_count() - 1))

    # Default output filename
    if out_file is None:
        out_file = default_output_filename(in_file)

    # Ensure output directory exists (creates parents as needed)
    out_path = Path(out_file)
    if out_path.parent and not out_path.parent.exists():
        print(f"Creating output directory: {out_path.parent}", file=sys.stderr)
        out_path.parent.mkdir(parents=True, exist_ok=True)

    # Read config
    primate_list, align_order, ref_threshold = read_json_lists(json_file)
    primate_set = frozenset(primate_list)

    in_path = Path(in_file)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_file}")

    # Create partial function
    process_func = partial(
        process_block_optimized,
        primate_set=primate_set,
        align_order=align_order,
        ref_threshold=ref_threshold,
        skip_gaps=skip_gaps,
        verbose=verbose
    )

    print(f"=== Go2Var ===", file=sys.stderr)
    print(f"Input: {in_file}", file=sys.stderr)
    print(f"Output: {out_file}", file=sys.stderr)
    print(f"Workers: {num_workers}", file=sys.stderr)
    print(f"Block size: {block_size}", file=sys.stderr)
    print(f"Reference species: {len(primate_list)} (r_num={ref_threshold})", file=sys.stderr)
    print(f"Merge output: {merge_output}", file=sys.stderr)

    start_time = time.time()

    blocks_processed = 0
    blocks_passed = 0

    write_queue: queue.Queue = queue.Queue()
    writer = threading.Thread(
        target=writer_thread_optimized,
        args=(out_file, write_queue, merge_output),
        daemon=False
    )
    writer.start()

    try:
        with mp.Pool(processes=num_workers) as pool:
            for batch in read_maf_blocks_fast(in_file, block_size):
                # Process batch
                results = pool.map(process_func, batch)

                # Send to writer
                for idx, variant_lines in results:
                    if variant_lines is not None:
                        write_queue.put((idx, variant_lines))
                        if variant_lines:
                            blocks_passed += 1
                    else:
                        write_queue.put((idx, []))

                blocks_processed += len(batch)

                if blocks_processed % 10000 == 0:
                    elapsed = time.time() - start_time
                    rate = blocks_processed / elapsed if elapsed > 0 else 0
                    print(f"Processed {blocks_processed} blocks, "
                          f"{blocks_passed} with variants "
                          f"({rate:.0f} blocks/sec)...", file=sys.stderr)

        # Signal completion
        write_queue.put((None, None))
        write_queue.join()
        writer.join()

        elapsed = time.time() - start_time
        print(f"\n=== Complete ===", file=sys.stderr)
        print(f"Total blocks: {blocks_processed}", file=sys.stderr)
        print(f"Blocks with variants: {blocks_passed}", file=sys.stderr)
        print(f"Time: {elapsed:.2f} seconds", file=sys.stderr)
        print(f"Speed: {blocks_processed / elapsed:.0f} blocks/sec", file=sys.stderr)
        print(f"Output: {out_file}", file=sys.stderr)

    except Exception as e:
        try:
            write_queue.put((None, None))
        except:
            pass
        raise IOError(f"Error: {e}") from e


def main():
    """Console script entry point."""
    parser = argparse.ArgumentParser(
        description="Go2Var: High-performance integrated Sort + Variant pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Performance improvements over v4:
  - Custom MAF parser 
  - NumPy vectorized variant detection - processes entire blocks at once
  - Byte-level sequence operations
  - Larger I/O buffers (8-16MB)
  - Maintains all filtering/reordering logic

Examples:
  python go2var.py input.maf config.json
  # Output: input.bed

  python go2var.py input.maf.gz config.json -w 8 --merge
  # Output: input.bed

  python go2var.py input.maf config.json custom_output.bed
        """
    )
    parser.add_argument("input_file", type=str, help="Input MAF file")
    parser.add_argument("json_file", type=str, help="JSON config file")
    parser.add_argument("output_file", type=str, nargs='?', default=None,
                        help="Output BED file (default: input basename + .bed)")
    parser.add_argument("-w", "--workers", type=int, default=None,
                        help="Number of workers (default: 8 or CPU-1)")
    parser.add_argument("-b", "--block-size", type=int, default=500,
                        help="Blocks per batch (default: 500)")
    parser.add_argument("--include-gaps", action="store_true", default=False,
                        help="Include gap positions")
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Include substitution details")
    parser.add_argument("--merge", action="store_true", default=False,
                        help="Merge consecutive positions")

    args = parser.parse_args()

    try:
        go2var_sorted_optimized(
            args.input_file,
            args.json_file,
            args.output_file,
            skip_gaps=not args.include_gaps,
            num_workers=args.workers,
            block_size=args.block_size,
            verbose=args.verbose,
            merge_output=args.merge
        )
    except (FileNotFoundError, ValueError, json.JSONDecodeError, IOError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
