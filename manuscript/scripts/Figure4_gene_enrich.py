#!/usr/bin/env python3
"""
Query g:Profiler with gene subsets and download GO enrichment results as CSV.

Usage:
    python gene_enrich.py <gene_list_file> [--output-dir DIR] [--sizes N N ...]
                               [--base-url URL]

Arguments:
    gene_list_file: Path to a tab-separated file with gene names in the first column

Options:
    --output-dir, -o: Directory to save output files (default: same as input file)
    --sizes:          Subset sizes to analyze (default: 250 500 1000 2500)
    --base-url:       Full g:Profiler endpoint URL. Use an archive URL such as
                      https://biit.cs.ut.ee/gprofiler_archive3/e113_eg59_p19/
                      to pin a specific Ensembl release. Default: live/latest.
"""

import argparse
import os
import re
import pandas as pd
import requests
from gprofiler import GProfiler


def normalize_base_url(base_url):
    """
    Normalize a g:Profiler base URL so callers can pass either the root
    (.../e113_eg59_p19/) or a deeper path (.../e113_eg59_p19/gost,
    .../api/gost/profile/) interchangeably. The gprofiler-official package
    and resolve_ambiguous_genes() both append '/api/gost/profile/' to the
    base, so we strip those suffixes if present and trim trailing slashes.
    """
    if base_url is None:
        return None
    url = base_url.strip().rstrip('/')
    # Peel known trailing path segments
    for suffix in ('/api/gost/profile', '/api/gost', '/gost'):
        if url.endswith(suffix):
            url = url[: -len(suffix)]
    return url


def resolve_ambiguous_genes(gene_list, base_url, organism='hsapiens'):
    """
    Query g:Profiler to find ambiguous gene names, then resolve each by
    picking the Ensembl ID with the most GO annotations.
    Returns a new gene list with ambiguous names replaced by Ensembl IDs.
    """
    payload = {
        'organism': organism,
        'query': gene_list,
        'sources': [],
        'user_threshold': 1.0,
        'no_evidences': True,
    }
    r = requests.post(f"{base_url}/api/gost/profile/", json=payload)
    raw = r.json()
    meta = raw.get('meta', {}).get('genes_metadata', {})
    ambiguous = meta.get('ambiguous', {})

    if not ambiguous:
        return gene_list, {}

    # For each ambiguous gene, pick the Ensembl ID with the most GO annotations
    resolutions = {}
    for gene_name, candidates in ambiguous.items():
        best = max(candidates, key=lambda c: c.get('number_of_go_annotations', 0))
        resolutions[gene_name] = best['gene_id']

    # Replace ambiguous gene names with resolved Ensembl IDs
    resolved_list = []
    for g in gene_list:
        if g in resolutions:
            resolved_list.append(resolutions[g])
        else:
            resolved_list.append(g)

    return resolved_list, resolutions


def derive_version_tag(base_url):
    """
    Extract a short version tag from a g:Profiler archive URL, e.g.
    https://biit.cs.ut.ee/gprofiler_archive3/e113_eg59_p19/ -> 'e113_eg59_p19'.
    Returns 'live' for the default endpoint or any URL without an e<n>_eg<n>_p<n> segment.
    """
    if base_url is None:
        return 'live'
    m = re.search(r'(e\d+_eg\d+_p\d+)', base_url)
    if m:
        return m.group(1)
    # Fallback: sanitize any non-alphanumeric characters
    tag = re.sub(r'[^A-Za-z0-9]+', '_', base_url.strip('/').split('/')[-1])
    return tag or 'custom'


def run_gprofiler_analysis(gene_file, subset_sizes=[250, 500, 1000, 2500],
                           output_dir=None, sig_method='g_SCS',
                           resolve_ambiguous=True, base_url=None):
    """
    Run g:Profiler GO enrichment analysis on gene subsets.

    Args:
        gene_file: Path to tab-separated file with gene names in first column
        subset_sizes: List of subset sizes to analyze (default: [250, 500, 1000, 2500])
        output_dir: Directory to save output files (default: same as input file)
        base_url: Optional full g:Profiler endpoint URL (e.g. an archived release).
                  If None, the gprofiler-official package default (live) is used.
    """
    # Get output directory from input file path if not specified
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(gene_file))
    else:
        output_dir = os.path.abspath(output_dir)
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

    # Normalize so users can pass either the root URL or one ending in
    # /gost, /api/gost, /api/gost/profile, or with a trailing slash.
    base_url = normalize_base_url(base_url)

    base_name = os.path.splitext(os.path.basename(gene_file))[0]
    version_tag = derive_version_tag(base_url)

    # Read gene list (genes are in first column, tab-separated)
    genes_df = pd.read_csv(gene_file, sep='\t', header=None, names=['gene', 'score'])
    genes = genes_df['gene'].str.strip().tolist()

    print(f"Input file: {gene_file}")
    print(f"Output directory: {output_dir}")
    print(f"Total genes loaded: {len(genes)}")

    # Initialize g:Profiler (optionally pinned to an archived endpoint)
    gp_kwargs = {'return_dataframe': True}
    if base_url is not None:
        gp_kwargs['base_url'] = base_url
    gp = GProfiler(**gp_kwargs)
    print(f"g:Profiler endpoint: {gp.base_url}")
    print(f"Version tag: {version_tag}")

    # Sources to query (GO and other databases)
    go_sources = ["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP", "TF", "MIRNA", "HPA", "CORUM", "HP"]

    # Process each subset
    for size in subset_sizes:
        if size > len(genes):
            print(f"\nSkipping top {size} genes (only {len(genes)} genes available)")
            continue

        print(f"\nProcessing top {size} genes...")

        gene_subset = genes[:size]
        print(f"  Genes in subset: {len(gene_subset)}")

        # Optionally resolve ambiguous gene names by picking the Ensembl ID
        # with the most GO annotations (matches web interface behavior)
        if resolve_ambiguous:
            resolved_subset, resolutions = resolve_ambiguous_genes(
                gene_subset, gp.base_url)
            if resolutions:
                print(f"  Resolved {len(resolutions)} ambiguous gene(s):")
                for name, ensg in resolutions.items():
                    print(f"    {name} -> {ensg}")
        else:
            resolved_subset = gene_subset
            print(f"  Skipping ambiguous gene resolution (--no-resolve)")

        # Query g:Profiler for GO terms using the resolved gene list
        try:
            results = gp.profile(
                organism='hsapiens',
                query=resolved_subset,
                sources=go_sources,
                user_threshold=0.05,
                significance_threshold_method=sig_method,
                no_evidences=False
            )

            if results is not None and len(results) > 0:
                # Save all results to a single CSV file (versioned filename)
                output_file = os.path.join(
                    output_dir,
                    f"{base_name}_top{size}_gprofiler_{version_tag}.csv"
                )
                results.to_csv(output_file, index=False)
                print(f"  Saved: {output_file} ({len(results)} terms)")

                # Print summary by source
                for source in go_sources:
                    count = len(results[results['source'] == source])
                    if count > 0:
                        print(f"    {source}: {count} terms")
            else:
                print(f"  No significant terms found for top {size} genes")

        except Exception as e:
            print(f"  Error querying g:Profiler for top {size} genes: {e}")

    print("\nDone!")


def main():
    parser = argparse.ArgumentParser(
        description="Query g:Profiler with gene subsets and download GO enrichment "
                    "results as CSV. Supports pinning to an archived Ensembl release "
                    "via --base-url."
    )
    parser.add_argument(
        "gene_file",
        help="Path to tab-separated file with gene names in the first column"
    )
    parser.add_argument(
        "--sizes",
        type=int,
        nargs='+',
        default=[250, 500, 1000, 2500],
        help="Subset sizes to analyze (default: 250 500 1000 2500)"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default=None,
        help="Output directory for CSV files (default: same directory as input file)"
    )
    parser.add_argument(
        "--sig-method", "-s",
        type=str,
        choices=['g_SCS', 'fdr', 'bonferroni'],
        default='g_SCS',
        help="Significance threshold method (default: g_SCS)"
    )
    parser.add_argument(
        "--no-resolve",
        action="store_true",
        default=False,
        help="Skip ambiguous gene resolution (default: resolve ambiguous genes)"
    )
    parser.add_argument(
        "--base-url",
        type=str,
        default=None,
        help="Full g:Profiler endpoint URL. Use an archive URL such as "
             "https://biit.cs.ut.ee/gprofiler_archive3/e112_eg59_p19/ to pin a "
             "specific Ensembl release. Default: package default (live/latest)."
    )

    args = parser.parse_args()

    if not os.path.exists(args.gene_file):
        print(f"Error: File not found: {args.gene_file}")
        return 1

    run_gprofiler_analysis(args.gene_file, args.sizes, args.output_dir,
                           args.sig_method,
                           resolve_ambiguous=not args.no_resolve,
                           base_url=args.base_url)
    return 0


if __name__ == "__main__":
    exit(main())
