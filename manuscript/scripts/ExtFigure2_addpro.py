#!/usr/bin/env python3
"""addpro.py — parse a GTF, emit gene-body + 1 kb upstream promoter intervals,
rename RefSeq accessions (NC_*.* -> chrN), and drop alt chromosomes."""

import argparse
import gzip
import os
import sys

import pandas as pd

REFSEQ_TO_CHR = {
    'NC_000001.11': 'chr1',
    'NC_000002.12': 'chr2',
    'NC_000003.12': 'chr3',
    'NC_000004.12': 'chr4',
    'NC_000005.10': 'chr5',
    'NC_000006.12': 'chr6',
    'NC_000007.14': 'chr7',
    'NC_000008.11': 'chr8',
    'NC_000009.12': 'chr9',
    'NC_000010.11': 'chr10',
    'NC_000011.10': 'chr11',
    'NC_000012.12': 'chr12',
    'NC_000013.11': 'chr13',
    'NC_000014.9': 'chr14',
    'NC_000015.10': 'chr15',
    'NC_000016.10': 'chr16',
    'NC_000017.11': 'chr17',
    'NC_000018.10': 'chr18',
    'NC_000019.10': 'chr19',
    'NC_000020.11': 'chr20',
    'NC_000021.9': 'chr21',
    'NC_000022.11': 'chr22',
    'NC_000023.11': 'chrX',
    'NC_000024.10': 'chrY',
    'NC_012920.1': 'chrM',
}


def open_gtf(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def basename_no_ext(path):
    name = os.path.basename(path)
    if name.endswith('.gz'):
        name = name[:-3]
    if '.' in name:
        name = name.rsplit('.', 1)[0]
    return name


def parse_gtf(path):
    """Yield [chrom, start, end, info, '.', strand] for gene + 1 kb promoter."""
    with open_gtf(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            part = line.rstrip('\n').split('\t')
            if len(part) < 7 or part[2] != 'gene':
                continue
            chrom = part[0]
            n_start = int(part[3])
            n_end = int(part[4])
            strand = part[6]
            info = part[8] if len(part) > 8 else '.'

            if strand == '+':
                prom_start = max(0, n_start - 1001)
                prom_end = n_start - 1
            elif strand == '-':
                prom_start = n_end
                prom_end = n_end + 1000
            else:
                print(f"Warning: unknown strand on line: {line.rstrip()}",
                      file=sys.stderr)
                continue

            gene_start = n_start - 1
            gene_end = n_end

            yield [chrom, gene_start, gene_end, info, '.', strand]
            yield [chrom, prom_start, prom_end, info, '.', strand]


def main():
    ap = argparse.ArgumentParser(
        description="GTF -> gene+1kb-promoter BED with RefSeq->chr renaming "
                    "and alt-chromosome filtering.")
    ap.add_argument('input_gtf', help='Input GTF file (.gtf or .gtf.gz)')
    ap.add_argument('-o', '--output',
                    help='Output BED path. '
                         'Default: <input_basename>.1kbpromoter.bed in input dir.')
    ap.add_argument('-n', '--no-filter', action='store_true',
                    help='Keep rows whose chrom does not start with "chr" '
                         '(default: drop them).')
    args = ap.parse_args()

    if args.output:
        out_path = args.output
    else:
        out_path = os.path.join(
            os.path.dirname(os.path.abspath(args.input_gtf)),
            basename_no_ext(args.input_gtf) + '.1kbpromoter.bed',
        )

    rows = list(parse_gtf(args.input_gtf))
    if not rows:
        print("Error: no 'gene' features parsed from input.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(rows)
    df[0] = df[0].map(lambda c: REFSEQ_TO_CHR.get(c, c))

    if not args.no_filter:
        before = len(df)
        df = df[df[0].astype(str).str.startswith('chr')]
        dropped = before - len(df)
        if dropped:
            print(f"Filtered out {dropped} rows with non-standard chromosome names")
    else:
        print("Keeping all rows (no filtering applied)")

    df.to_csv(out_path, sep='\t', header=False, index=False)
    print(f"Wrote {len(df)} lines to {out_path}")


if __name__ == '__main__':
    main()
