#!/usr/bin/env python3
"""addpro: build a gene + 1 kb upstream-promoter BED from a GTF.

The CLI reads a GTF (plain text or gzip) and emits BED rows for each
``gene`` feature. Three output modes are supported:

* ``separate`` (default) — two rows per gene: the gene body, then a
  strand-aware 1 kb upstream promoter. For ``+`` strand genes the
  promoter is ``[max(0, start - 1 - N), start - 1)``; for ``-`` strand
  genes it is ``[end, end + N)``. Default ``N`` is 1000.
* ``merged`` — one row per gene: gene body and the upstream promoter
  joined into a single interval (``[promoter_start, gene_end)`` for
  ``+`` strand, ``[gene_start, promoter_end)`` for ``-`` strand).
* ``both`` — writes the ``separate`` and ``merged`` outputs to two
  files using a shared path prefix.

The BED ``name`` column carries a clean gene symbol resolved from the
GTF attribute string in the order ``gene_name`` -> ``gene`` ->
``gene_id``. When a symbol is shared by multiple ``gene`` records,
every occurrence is suffixed with ``|<gene_id>`` so names are unique
within the file.

RefSeq accessions in the chromosome field (``NC_000001.11``,
``NC_000023.11``, etc.) are renamed to UCSC-style names (``chr1``,
``chrX``, etc.) for GRCh38.p14. Rows whose chromosome name does not
start with ``chr`` after renaming are dropped by default; pass
``--no-filter`` to keep them.

The ``-o``/``--output`` option takes a full output path (directory +
filename); parent directories are created automatically. In ``both``
mode it is treated as a path prefix (any trailing ``.bed`` is stripped)
and two files are written: ``<prefix>.separate.bed`` and
``<prefix>.merged.bed``.
"""

import argparse
import gzip
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Dict, Iterator, List, Optional

__all__ = [
    "REFSEQ_TO_CHR",
    "GeneRecord",
    "open_gtf",
    "parse_gtf_attributes",
    "iter_gene_records",
    "default_output_filename",
    "addpro_gtf_to_bed",
    "main",
]

logger = logging.getLogger(__name__)

REFSEQ_TO_CHR = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9":  "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9":  "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
    "NC_012920.1":  "chrM",
}

DEFAULT_NAME_PRIORITY = ("gene_name", "gene", "gene_id")
ATTR_RE = re.compile(r'\s*([A-Za-z0-9_]+)\s+"([^"]*)"\s*;')
VALID_MODES = ("separate", "merged", "both")


@dataclass
class GeneRecord:
    """Per-gene record carried between parse and write phases."""

    chrom: str
    gene_start0: int  # BED 0-based start (== n_start - 1)
    gene_end: int     # BED half-open end (== n_end)
    strand: str
    gene_id: str
    symbol: str       # resolved name (before any |gene_id disambiguation)


def open_gtf(path: str) -> IO[str]:
    """Open a GTF file as text, transparently handling ``.gz``."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_gtf_attributes(attr_field: str) -> Dict[str, str]:
    """Parse the 9th GTF column into a dict of attribute key -> value."""
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(attr_field)}


def _resolve_symbol(attrs: Dict[str, str], priority) -> str:
    for key in priority:
        val = attrs.get(key, "").strip()
        if val:
            return val
    return ""


def default_output_filename(in_file: str, mode: str = "separate", upstream: int = 1000) -> str:
    """Build a default output path next to the input GTF.

    ``separate`` + 1000 bp -> ``<base>.1kbpromoter.bed`` (legacy).
    ``separate`` + N bp    -> ``<base>.{N}bp_promoter.bed``.
    ``merged`` + 1000 bp   -> ``<base>.gene_1kbpromoter.bed``.
    ``merged`` + N bp      -> ``<base>.gene_{N}bp_promoter.bed``.
    ``both``               -> ``<base>`` (used as a prefix; the writer
    appends ``.separate.bed`` and ``.merged.bed``).
    """
    in_path = Path(in_file).resolve()
    base = in_path.name
    if base.endswith(".gz"):
        base = base[:-3]
    if base.endswith(".gtf"):
        base = base[:-4]
    parent = in_path.parent

    if mode == "both":
        return str(parent / base)

    width_tag = "1kbpromoter" if upstream == 1000 else f"{upstream}bp_promoter"
    if mode == "merged":
        return str(parent / f"{base}.gene_{width_tag}.bed")
    return str(parent / f"{base}.{width_tag}.bed")


def _ensure_parent_dir(out_filename: str) -> None:
    parent = Path(out_filename).parent
    if str(parent) not in ("", "."):
        parent.mkdir(parents=True, exist_ok=True)


def iter_gene_records(path: str) -> Iterator[GeneRecord]:
    """Yield one ``GeneRecord`` per ``gene`` line in the GTF.

    Rows with unrecognised strand or missing ``gene_id`` are skipped;
    counts are logged at INFO level. The BED ``name`` is resolved as
    ``gene_name`` -> ``gene`` -> ``gene_id``.
    """
    n_unknown_strand = 0
    n_missing_id = 0
    with open_gtf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            part = line.rstrip("\n").rstrip("\r").split("\t")
            if len(part) < 9 or part[2] != "gene":
                continue

            strand = part[6]
            if strand not in ("+", "-"):
                n_unknown_strand += 1
                continue

            try:
                n_start = int(part[3])
                n_end = int(part[4])
            except ValueError:
                continue

            attrs = parse_gtf_attributes(part[8])
            gene_id = attrs.get("gene_id", "").strip()
            if not gene_id:
                n_missing_id += 1
                continue

            symbol = _resolve_symbol(attrs, DEFAULT_NAME_PRIORITY) or gene_id

            yield GeneRecord(
                chrom=part[0],
                gene_start0=n_start - 1,
                gene_end=n_end,
                strand=strand,
                gene_id=gene_id,
                symbol=symbol,
            )

    if n_unknown_strand:
        logger.warning(
            "Skipped %d gene features with unrecognised strand", n_unknown_strand
        )
    if n_missing_id:
        logger.warning(
            "Skipped %d gene features with missing gene_id", n_missing_id
        )


def _disambiguate(records: List[GeneRecord]) -> Dict[str, str]:
    """Return ``gene_id -> name`` map; suffix ``|gene_id`` on collisions."""
    counts: Dict[str, int] = {}
    for rec in records:
        counts[rec.symbol] = counts.get(rec.symbol, 0) + 1
    n_dup = 0
    name_by_id: Dict[str, str] = {}
    for rec in records:
        if counts[rec.symbol] > 1:
            name_by_id[rec.gene_id] = f"{rec.symbol}|{rec.gene_id}"
            n_dup += 1
        else:
            name_by_id[rec.gene_id] = rec.symbol
    if n_dup:
        logger.info("Disambiguated %d symbols with |gene_id suffix", n_dup)
    return name_by_id


def _apply_chrom_rename(chrom: str) -> str:
    return REFSEQ_TO_CHR.get(chrom, chrom)


def _row_passes_filter(chrom: str, keep_non_chr: bool) -> bool:
    if keep_non_chr:
        return True
    return chrom.startswith("chr")


def _promoter_interval(rec: GeneRecord, upstream: int):
    """Return ``(prom_start, prom_end)`` in BED 0-based half-open coords."""
    if rec.strand == "+":
        return max(0, rec.gene_start0 - upstream), rec.gene_start0
    # rec.strand == "-"
    return rec.gene_end, rec.gene_end + upstream


def _write_separate(records, name_by_id, upstream, keep_non_chr, fh) -> int:
    n = 0
    for rec in records:
        chrom = _apply_chrom_rename(rec.chrom)
        if not _row_passes_filter(chrom, keep_non_chr):
            continue
        name = name_by_id[rec.gene_id]
        prom_start, prom_end = _promoter_interval(rec, upstream)
        fh.write(f"{chrom}\t{rec.gene_start0}\t{rec.gene_end}\t{name}\t.\t{rec.strand}\n")
        fh.write(f"{chrom}\t{prom_start}\t{prom_end}\t{name}\t.\t{rec.strand}\n")
        n += 2
    return n


def _write_merged(records, name_by_id, upstream, keep_non_chr, fh) -> int:
    n = 0
    for rec in records:
        chrom = _apply_chrom_rename(rec.chrom)
        if not _row_passes_filter(chrom, keep_non_chr):
            continue
        name = name_by_id[rec.gene_id]
        prom_start, prom_end = _promoter_interval(rec, upstream)
        if rec.strand == "+":
            start, end = prom_start, rec.gene_end
        else:
            start, end = rec.gene_start0, prom_end
        fh.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t{rec.strand}\n")
        n += 1
    return n


def _strip_bed_suffix(path: str) -> str:
    return path[:-4] if path.endswith(".bed") else path


def addpro_gtf_to_bed(
    in_filename: str,
    out_filename: Optional[str] = None,
    keep_non_chr: bool = False,
    mode: str = "separate",
    upstream: int = 1000,
) -> None:
    """Write a gene + 1 kb promoter BED file from a GTF.

    Args:
        in_filename: Path to the input GTF (``.gtf`` or ``.gtf.gz``).
        out_filename: For ``mode='separate'`` or ``'merged'``, full output
            BED path. For ``mode='both'``, a path prefix; ``.separate.bed``
            and ``.merged.bed`` are appended. Parent directories are
            created if needed. If ``None``, defaults are derived from the
            input basename.
        keep_non_chr: If True, retain rows whose chromosome name does not
            start with ``chr`` after RefSeq renaming.
        mode: One of ``'separate'`` (default), ``'merged'``, ``'both'``.
        upstream: Promoter window width in bp (default 1000, must be >=0).

    Raises:
        FileNotFoundError: If ``in_filename`` does not exist.
        ValueError: If ``mode`` is invalid, ``upstream`` is negative, or
            the GTF contains no ``gene`` features.
    """
    if mode not in VALID_MODES:
        raise ValueError(f"mode must be one of {VALID_MODES}, got {mode!r}")
    if upstream < 0:
        raise ValueError(f"upstream must be >= 0, got {upstream}")

    in_path = Path(in_filename)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_filename}")

    if out_filename is None:
        out_filename = default_output_filename(in_filename, mode=mode, upstream=upstream)

    logger.info("=== addpro ===")
    logger.info("Input: %s", in_filename)
    logger.info("Mode:  %s", mode)
    logger.info("Upstream: %d bp", upstream)

    records = list(iter_gene_records(in_filename))
    if not records:
        raise ValueError(f"No 'gene' features parsed from {in_filename}")

    name_by_id = _disambiguate(records)

    if mode == "both":
        prefix = _strip_bed_suffix(out_filename)
        sep_path = f"{prefix}.separate.bed"
        mer_path = f"{prefix}.merged.bed"
        _ensure_parent_dir(sep_path)
        _ensure_parent_dir(mer_path)
        with open(sep_path, "w") as fh:
            n_sep = _write_separate(records, name_by_id, upstream, keep_non_chr, fh)
        with open(mer_path, "w") as fh:
            n_mer = _write_merged(records, name_by_id, upstream, keep_non_chr, fh)
        logger.info("Wrote %d lines to %s", n_sep, sep_path)
        logger.info("Wrote %d lines to %s", n_mer, mer_path)
    else:
        _ensure_parent_dir(out_filename)
        with open(out_filename, "w") as fh:
            if mode == "separate":
                n = _write_separate(records, name_by_id, upstream, keep_non_chr, fh)
            else:
                n = _write_merged(records, name_by_id, upstream, keep_non_chr, fh)
        logger.info("Wrote %d lines to %s", n, out_filename)

    if not keep_non_chr:
        logger.info("Filtered rows whose chromosome did not start with 'chr'")


def main():
    """Console-script entry point for ``addpro``."""
    parser = argparse.ArgumentParser(
        description="addpro: GTF -> gene + 1 kb upstream-promoter BED, with "
                    "RefSeq (NC_*.*) -> UCSC chromosome renaming, configurable "
                    "promoter width, and three output modes (separate / merged / both).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
-o / --output accepts a FULL output path (directory + filename) for
'separate' and 'merged' modes. In 'both' mode it is treated as a path
prefix; the writer appends '.separate.bed' and '.merged.bed'. Parent
directories are created automatically.

Examples:
  addpro annotation.gtf
  # default 'separate' mode, output: annotation.1kbpromoter.bed

  addpro annotation.gtf --mode merged -o results/gene_up1kb.bed
  addpro annotation.gtf --mode both -o results/run1
  # writes results/run1.separate.bed and results/run1.merged.bed

  addpro annotation.gtf.gz --upstream 2000 -o promoters_2kb.bed
  addpro annotation.gtf --no-filter            # keep non-chr rows
  addpro annotation.gtf --quiet                # suppress info messages
        """,
    )
    parser.add_argument(
        "input_gtf",
        type=str,
        help="Input GTF annotation (.gtf or .gtf.gz)",
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_file",
        type=str,
        default=None,
        help="Full output BED path (or path prefix in 'both' mode). "
             "Parent directories are created if needed. "
             "Default: <input_basename>.1kbpromoter.bed",
    )
    parser.add_argument(
        "-m", "--mode",
        choices=VALID_MODES,
        default="separate",
        help="Output mode: 'separate' (gene body + promoter as two rows; default), "
             "'merged' (one row per gene = gene body + upstream extension), "
             "or 'both' (write both files using -o as a prefix).",
    )
    parser.add_argument(
        "-u", "--upstream",
        type=int,
        default=1000,
        help="Promoter window width in bp (default: 1000).",
    )
    parser.add_argument(
        "-n", "--no-filter",
        action="store_true",
        default=False,
        help='Keep rows whose chromosome does not start with "chr" '
             "(default: drop them).",
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
        addpro_gtf_to_bed(
            args.input_gtf,
            out_filename=args.output_file,
            keep_non_chr=args.no_filter,
            mode=args.mode,
            upstream=args.upstream,
        )
    except (OSError, FileNotFoundError, ValueError) as e:
        logger.error("Error: %s", e)
        sys.exit(1)
    except Exception as e:
        logger.error("Unexpected error: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
