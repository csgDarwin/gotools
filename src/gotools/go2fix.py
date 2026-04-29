#!/usr/bin/env python3
"""go2fix: identify conserved positions in DNA multiple sequence alignments.

The CLI reads a MAF alignment, finds positions where at least ``--max-conserved``
rows agree on the same non-N, non-gap base relative to the first (reference)
sequence, and writes a BED file of those positions (merged into intervals by
default).

The ``-o``/``--output`` option takes a full output path (directory + filename);
parent directories are created automatically.
"""

import argparse
import gzip
import logging
import multiprocessing as mp
import queue
import sys
import threading
import time
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import numpy as np

__all__ = [
    "parse_maf_blocks",
    "find_conserved_vectorized",
    "process_block",
    "writer_thread_fast",
    "default_output_filename",
    "go2fix_optimized",
    "go2fix_single",
    "main",
]

logger = logging.getLogger(__name__)

_UPPER_TABLE = np.zeros(256, dtype=np.uint8)
for _i in range(256):
    _UPPER_TABLE[_i] = ord(chr(_i).upper()) if _i < 128 else _i

_GAP = ord('-')
_N = ord('N')

_ERR = "__ERROR__"


def parse_maf_blocks(file_handle) -> Iterator[Tuple[str, int, List[bytes]]]:
    """Yield ``(ref_id, ref_start, sequences)`` tuples for each MAF block.

    Args:
        file_handle: Iterable of MAF lines (text or bytes).

    Yields:
        Tuples of the reference (first) sequence identifier, its 0-based start
        position, and the list of aligned sequence bytes for the block.
    """
    current_sequences: List[bytes] = []
    ref_id = ""
    ref_start = 0
    in_block = False

    for line in file_handle:
        if isinstance(line, bytes):
            line = line.decode('ascii', errors='ignore')

        if not line or line[0] == '#':
            continue

        if line[0] == 'a':
            if in_block and current_sequences:
                yield (ref_id, ref_start, current_sequences)
            current_sequences = []
            in_block = True

        elif line[0] == 's' and in_block:
            parts = line.split()
            if len(parts) >= 7:
                if not current_sequences:
                    ref_id = parts[1]
                    ref_start = int(parts[2])
                # (#2) no .upper() — _UPPER_TABLE handles case in NumPy.
                current_sequences.append(parts[6].encode('ascii'))

        elif line[0] in ('', '\n', '\r') and in_block and current_sequences:
            pass

    if current_sequences:
        yield (ref_id, ref_start, current_sequences)


def find_conserved_vectorized(sequences: List[bytes], max_conserved: int) -> np.ndarray:
    """Return a boolean mask of conserved alignment columns.

    A column is conserved when every row matches the reference (row 0) base,
    no row contains an ``N``, and the reference base is not a gap. Returns an
    empty array if there are fewer than ``max_conserved`` rows or no input.

    Args:
        sequences: Aligned rows as ASCII byte strings (all the same length).
        max_conserved: Minimum number of rows required for a conserved call.

    Returns:
        A 1-D boolean ``ndarray`` of length ``len(sequences[0])``.
    """
    n_seqs = len(sequences)
    if n_seqs < max_conserved or not sequences:
        return np.array([], dtype=bool)

    seq_len = len(sequences[0])
    seq_array = np.frombuffer(b''.join(sequences), dtype=np.uint8).reshape(n_seqs, seq_len)
    seq_array = _UPPER_TABLE[seq_array]

    ref_base = seq_array[0, :]
    not_gap = ref_base != _GAP
    has_n = np.any(seq_array == _N, axis=0)
    all_identical = np.all(seq_array == ref_base, axis=0)

    return not_gap & ~has_n & all_identical


def process_block(block_data: Tuple) -> Tuple[int, List[Tuple]]:
    """Process a single MAF block into BED tuples.

    Args:
        block_data: ``(block_index, ref_id, ref_start, sequences, max_conserved, verbose)``.

    Returns:
        ``(block_index, [(chrom, start, end, base_or_None), ...])``. On
        per-block failure the second element is a single error tuple keyed by
        ``_ERR`` so the writer can skip it.
    """
    block_index, ref_id, ref_start, sequences, max_conserved, verbose = block_data
    try:
        if not sequences:
            return (block_index, [])

        conserved = find_conserved_vectorized(sequences, max_conserved)
        if len(conserved) == 0:
            return (block_index, [])

        first_seq = np.frombuffer(sequences[0], dtype=np.uint8)
        first_seq = _UPPER_TABLE[first_seq]
        is_gap = (first_seq == _GAP)
        non_gap_cumsum = np.cumsum(~is_gap)
        conserved_indices = np.where(conserved)[0]
        if len(conserved_indices) == 0:
            return (block_index, [])

        bed_tuples: List[Tuple] = []
        for idx in conserved_indices:
            if idx > 0:
                seq_pos = int(ref_start + non_gap_cumsum[idx] - 1)
            else:
                seq_pos = int(ref_start)
            start_pos = seq_pos
            end_pos = seq_pos + 1
            if verbose:
                base = chr(int(first_seq[idx]))
                bed_tuples.append((ref_id, start_pos, end_pos, base))
            else:
                bed_tuples.append((ref_id, start_pos, end_pos, None))
        return (block_index, bed_tuples)

    except Exception as e:
        return (block_index, [(_ERR, 0, 0, f"block {block_index}: {e}")])


def _block_generator(
    parser_iter: Iterable[Tuple[str, int, List[bytes]]],
    max_conserved: int,
    verbose: bool,
) -> Iterator[Tuple]:
    """Yield per-block work units for the worker pool."""
    block_index = 0
    for ref_id, ref_start, sequences in parser_iter:
        if not sequences:
            continue
        yield (block_index, ref_id, ref_start, sequences, max_conserved, verbose)
        block_index += 1


def writer_thread_fast(
    out_filename: str,
    write_queue: "queue.Queue",
    merge_output: bool,
    verbose: bool,
) -> None:
    """Drain ``write_queue`` of per-block result tuples and write a BED file.

    Tuples arrive out of order from ``imap_unordered``; results are buffered
    and emitted in original block order. Merging adjacent positions on the
    same chromosome is enabled by ``merge_output``.
    """
    result_buffer: Dict[int, List[Tuple]] = {}
    next_expected = 0

    def _format_line(t: Tuple) -> str:
        if verbose and t[3] is not None:
            return f"{t[0]}\t{t[1]}\t{t[2]}\t >{t[3]}\n"
        return f"{t[0]}\t{t[1]}\t{t[2]}\n"

    try:
        with open(out_filename, "w", buffering=4 * 1024 * 1024) as outfile:
            if merge_output:
                cur_chrom: Optional[str] = None
                cur_start: Optional[int] = None
                cur_end: Optional[int] = None

                def _consume(tuples: List[Tuple]):
                    nonlocal cur_chrom, cur_start, cur_end
                    for t in tuples:
                        if t[0] == _ERR:
                            continue
                        chrom, start, end = t[0], t[1], t[2]
                        if cur_chrom == chrom and cur_end == start:
                            cur_end = end
                        else:
                            if cur_chrom is not None:
                                outfile.write(
                                    f"{cur_chrom}\t{cur_start}\t{cur_end}\n"
                                )
                            cur_chrom, cur_start, cur_end = chrom, start, end

                while True:
                    item = write_queue.get()
                    if item is None:
                        write_queue.task_done()
                        break
                    block_index, bed_tuples = item
                    result_buffer[block_index] = bed_tuples
                    while next_expected in result_buffer:
                        _consume(result_buffer.pop(next_expected))
                        next_expected += 1
                    write_queue.task_done()

                # Drain out-of-order stragglers (imap_unordered leaves gaps).
                while next_expected in result_buffer:
                    _consume(result_buffer.pop(next_expected))
                    next_expected += 1

                if cur_chrom is not None:
                    outfile.write(f"{cur_chrom}\t{cur_start}\t{cur_end}\n")

            else:
                write_buffer: List[str] = []
                max_lines = 50000

                while True:
                    item = write_queue.get()
                    if item is None:
                        write_queue.task_done()
                        break
                    block_index, bed_tuples = item
                    result_buffer[block_index] = bed_tuples
                    while next_expected in result_buffer:
                        for t in result_buffer.pop(next_expected):
                            if t[0] == _ERR:
                                continue
                            write_buffer.append(_format_line(t))
                        next_expected += 1
                    if len(write_buffer) >= max_lines:
                        outfile.writelines(write_buffer)
                        write_buffer.clear()
                    write_queue.task_done()

                while next_expected in result_buffer:
                    for t in result_buffer.pop(next_expected):
                        if t[0] == _ERR:
                            continue
                        write_buffer.append(_format_line(t))
                    next_expected += 1

                if write_buffer:
                    outfile.writelines(write_buffer)

    except Exception:
        logger.exception("Fatal error in writer thread")
        raise


def default_output_filename(in_file: str, output_dir: str = "fixed_sites/") -> str:
    """Build a default ``<dir>/<basename>_conserved.bed`` path from an input MAF."""
    base = Path(in_file).name
    if base.endswith('.maf.gz'):
        base = base[:-7]
    elif base.endswith('.gz'):
        base = base[:-3]
    if base.endswith('.maf'):
        base = base[:-4]
    return str(Path(output_dir) / f"{base}_conserved.bed")


def _ensure_parent_dir(out_filename: str) -> None:
    parent = Path(out_filename).parent
    if str(parent) not in ("", "."):
        parent.mkdir(parents=True, exist_ok=True)


def go2fix_optimized(
    in_filename: str,
    max_conserved: int = 464,
    num_workers: Optional[int] = None,
    chunksize: int = 8,
    verbose: bool = False,
    merge_output: bool = True,
    out_filename: Optional[str] = None,
    default_output_dir: str = "fixed_sites/",
) -> None:
    """Multiprocess pipeline: parser → workers → writer.

    Args:
        in_filename: Path to MAF or MAF.gz file.
        max_conserved: Minimum number of rows that must agree for a position
            to be called conserved.
        num_workers: Worker process count. Defaults to ``min(8, CPU - 1)``.
        chunksize: ``imap_unordered`` chunk size (blocks per dispatch).
        verbose: If True, the BED output includes the conserved base in a
            fourth column.
        merge_output: Merge adjacent conserved positions on the same chromosome.
        out_filename: Full output path. If None, derived from ``in_filename``
            and ``default_output_dir``.
        default_output_dir: Directory used when ``out_filename`` is None.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If ``max_conserved`` is less than 1.
        IOError: If the pipeline fails during processing.
    """
    max_conserved = int(max_conserved)
    if max_conserved < 1:
        raise ValueError(f"max_conserved must be >= 1, got {max_conserved}")

    in_path = Path(in_filename)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_filename}")

    is_gzipped = in_filename.endswith('.gz')

    if out_filename is None:
        out_filename = default_output_filename(in_filename, default_output_dir)
    _ensure_parent_dir(out_filename)

    if num_workers is None:
        num_workers = min(8, max(1, mp.cpu_count() - 1))

    logger.info("=== Go2Fix Optimized v2 (pipelined, tuple-IO) ===")
    logger.info("Input: %s", in_filename)
    logger.info("Workers: %d", num_workers)
    logger.info("imap_unordered chunksize: %d", chunksize)
    logger.info("Min rows (max_conserved): %d", max_conserved)
    logger.info("Merge output: %s", merge_output)
    logger.info("Output: %s", out_filename)

    start_time = time.time()

    write_queue: queue.Queue = queue.Queue(maxsize=num_workers * 8)

    writer = threading.Thread(
        target=writer_thread_fast,
        args=(out_filename, write_queue, merge_output, verbose),
        daemon=False,
    )
    writer.start()

    file_opener = gzip.open if is_gzipped else open
    open_mode = 'rt' if is_gzipped else 'r'

    blocks_done = 0
    try:
        with mp.Pool(processes=num_workers) as pool, \
             file_opener(in_filename, open_mode) as infile:
            gen = _block_generator(parse_maf_blocks(infile), max_conserved, verbose)
            # (#1) imap_unordered pipelines parser → workers → writer.
            for result in pool.imap_unordered(process_block, gen, chunksize=chunksize):
                write_queue.put(result)
                blocks_done += 1
                if blocks_done % 10000 == 0:
                    logger.info("Processed %d blocks...", blocks_done)

        write_queue.put(None)
        write_queue.join()
        writer.join()

        elapsed = time.time() - start_time
        logger.info("=== Complete ===")
        logger.info("Total blocks: %d", blocks_done)
        logger.info("Time: %.2f seconds", elapsed)
        if elapsed > 0:
            logger.info("Speed: %.1f blocks/sec", blocks_done / elapsed)
        logger.info("Output: %s", out_filename)

    except Exception as e:
        try:
            write_queue.put(None)
        except Exception:
            logger.debug("Failed to enqueue sentinel during error cleanup", exc_info=True)
        raise OSError(f"Error processing file: {e}") from e


def go2fix_single(
    in_filename: str,
    max_conserved: int = 464,
    verbose: bool = False,
    merge_output: bool = True,
    out_filename: Optional[str] = None,
    default_output_dir: str = "fixed_sites/",
) -> None:
    """Single-threaded reference implementation (matches the multiprocess output).

    Useful for debugging, small inputs, and environments where ``multiprocessing``
    is unavailable. Same semantics as :func:`go2fix_optimized` minus the worker
    pool and ``num_workers``/``chunksize`` parameters.
    """
    max_conserved = int(max_conserved)
    if max_conserved < 1:
        raise ValueError(f"max_conserved must be >= 1, got {max_conserved}")

    in_path = Path(in_filename)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_filename}")

    is_gzipped = in_filename.endswith('.gz')

    if out_filename is None:
        out_filename = default_output_filename(in_filename, default_output_dir)
    _ensure_parent_dir(out_filename)

    logger.info("=== Go2Fix (single-threaded) ===")
    start_time = time.time()

    file_opener = gzip.open if is_gzipped else open
    open_mode = 'rt' if is_gzipped else 'r'
    block_count = 0

    with open(out_filename, "w", buffering=4 * 1024 * 1024) as outfile, \
         file_opener(in_filename, open_mode) as infile:

        cur_chrom: Optional[str] = None
        cur_start: Optional[int] = None
        cur_end: Optional[int] = None

        for ref_id, ref_start, sequences in parse_maf_blocks(infile):
            if not sequences:
                continue
            block_count += 1
            if block_count % 5000 == 0:
                logger.info("Processed %d blocks...", block_count)

            _, bed_tuples = process_block(
                (0, ref_id, ref_start, sequences, max_conserved, verbose)
            )

            if merge_output:
                for t in bed_tuples:
                    if t[0] == _ERR:
                        continue
                    chrom, start, end = t[0], t[1], t[2]
                    if cur_chrom == chrom and cur_end == start:
                        cur_end = end
                    else:
                        if cur_chrom is not None:
                            outfile.write(
                                f"{cur_chrom}\t{cur_start}\t{cur_end}\n"
                            )
                        cur_chrom, cur_start, cur_end = chrom, start, end
            else:
                for t in bed_tuples:
                    if t[0] == _ERR:
                        continue
                    if verbose and t[3] is not None:
                        outfile.write(f"{t[0]}\t{t[1]}\t{t[2]}\t >{t[3]}\n")
                    else:
                        outfile.write(f"{t[0]}\t{t[1]}\t{t[2]}\n")

        if merge_output and cur_chrom is not None:
            outfile.write(f"{cur_chrom}\t{cur_start}\t{cur_end}\n")

    elapsed = time.time() - start_time
    logger.info("=== Complete ===")
    logger.info("Total blocks: %d", block_count)
    logger.info("Time: %.2f seconds", elapsed)
    logger.info("Output: %s", out_filename)


def main():
    """Console-script entry point for ``go2fix``."""
    parser = argparse.ArgumentParser(
        description="Go2Fix v2 - pipelined + tuple-IO (explicit output filename)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
-o / --output accepts a FULL output path (directory + filename).
Parent directories are created automatically.

Examples:
  go2fix alignment.maf.gz -o out/mybed.bed
  go2fix alignment.maf -o results/run1/conserved.bed --workers 16 -m 200
  go2fix alignment.maf
  # Output (no -o): fixed_sites/alignment_conserved.bed
        """,
    )
    parser.add_argument("input_file", type=str, help="Input MAF alignment file")
    parser.add_argument(
        "-o", "--output",
        dest="output_file",
        type=str,
        default=None,
        help="Full output file path (directory + filename). "
             "Parent directories are created if needed. "
             "Default: fixed_sites/<basename>_conserved.bed",
    )
    parser.add_argument(
        "-m", "--max-conserved",
        type=int,
        default=464,
        help="Minimum number of rows required for conserved call (default: 464)",
    )
    parser.add_argument(
        "--no-merge",
        action="store_true",
        default=False,
        help="Disable merging of consecutive conserved positions (merging is on by default)",
    )
    parser.add_argument(
        "--single-threaded",
        action="store_true",
        default=False,
        help="Use single-threaded mode",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of worker processes (default: min(8, CPU-1))",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=8,
        help="imap_unordered chunksize — blocks per pool dispatch (default: 8). "
             "Larger values reduce IPC; smaller values spread load more evenly.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Include base identity in output (adds 4th column)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress informational progress messages on stderr.",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(message)s",
        stream=sys.stderr,
        force=True,
    )

    try:
        merge = not args.no_merge

        if args.single_threaded:
            go2fix_single(
                args.input_file,
                max_conserved=args.max_conserved,
                verbose=args.verbose,
                merge_output=merge,
                out_filename=args.output_file,
            )
        else:
            go2fix_optimized(
                args.input_file,
                max_conserved=args.max_conserved,
                num_workers=args.workers,
                chunksize=args.chunksize,
                verbose=args.verbose,
                merge_output=merge,
                out_filename=args.output_file,
            )

    except (OSError, FileNotFoundError, ValueError) as e:
        logger.error("Error: %s", e)
        sys.exit(1)
    except Exception as e:
        logger.error("Unexpected error: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
