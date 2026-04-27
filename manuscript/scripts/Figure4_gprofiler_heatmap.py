#!/usr/bin/env python3
"""
Lineage enrichment heatmap for GO:BP and other sources from g:Profiler CSV results.

Filters:
  1. Source filter (e.g. GO:BP)
  2. Fold enrichment filter (--min-fe): removes terms below cutoff in ALL lineages
  3. Top N most significant terms per lineage (--top)

Generates a combined heatmap — all lineages as columns, union of top terms
as rows, grouped by lineage overlap.

Usage:
    python gprofiler_heatmap.py <csv> --source "GO:BP" --top 20
    python gprofiler_heatmap.py <csv> --source "GO:BP" --top 30 --min-fe 1.0
    python gprofiler_heatmap.py <csv> --source "GO:BP" --top 20 --rank-by fe
"""

import argparse
import os

import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42   # editable text in PDF
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ── helpers ──────────────────────────────────────────────────────────────────

def classify_terms_by_overlap(term_sets, lineage_order, term_pval=None):
    """
    Classify terms by which lineages share them.
    Returns ordered term list and term->lineage_set mapping.
    Groups are sorted: most shared first, then by position in lineage_order.
    Within each group, terms are sorted by min p-value (smallest first)
    if term_pval is provided, otherwise alphabetically.
    """
    lin_rank = {l: i for i, l in enumerate(lineage_order)}

    term_to_lineages = {}
    for l, terms in term_sets.items():
        for term in terms:
            term_to_lineages.setdefault(term, set()).add(l)
    term_to_lineages = {t: frozenset(ls) for t, ls in term_to_lineages.items()}

    groups = {}
    for term, lin_set in term_to_lineages.items():
        groups.setdefault(lin_set, []).append(term)

    sorted_keys = sorted(
        groups.keys(),
        key=lambda s: (-len(s),
                       [lin_rank.get(l, 999)
                        for l in sorted(s, key=lambda x: lin_rank.get(x, 999))]))

    ordered_terms = []
    for gk in sorted_keys:
        if term_pval is not None:
            ordered_terms.extend(
                sorted(groups[gk], key=lambda t: term_pval.get(t, 1.0)))
        else:
            ordered_terms.extend(sorted(groups[gk]))

    return ordered_terms, term_to_lineages


def prepare_data(input_file, source_filter):
    """Read CSV, filter by source, calculate fold enrichment."""
    # encoding='utf-8-sig' strips a BOM if the CSV came from a g:Profiler
    # web-interface export.
    df = pd.read_csv(input_file, encoding='utf-8-sig')

    # Normalize column names so both the Python-API output (p_value, native, name,
    # lineageName) and the web-interface CSV export (adjusted_p_value, term_id,
    # term_name, LineageName) work without further changes downstream.
    df = df.rename(columns={
        'adjusted_p_value': 'p_value',
        'term_id': 'native',
        'term_name': 'name',
        'LineageName': 'lineageName',
    })

    # Build lineage number -> lineageName mapping
    if 'lineageName' in df.columns:
        name_map = df.drop_duplicates('lineage') \
                     .set_index('lineage')['lineageName'].to_dict()
    else:
        name_map = {l: str(l) for l in df['lineage'].unique()}

    lineages = sorted(df['lineage'].unique())
    lineage_names = [name_map[l] for l in lineages]

    # Filter by source
    if source_filter.upper() == 'ALL':
        filtered = df.copy()
        source_label = 'All'
    else:
        filtered = df[df['source'] == source_filter].copy()
        source_label = source_filter

    print(f"{source_label} terms (before filtering): {len(filtered)}")

    # Compute fold enrichment
    filtered['fe'] = (filtered['intersection_size'] / filtered['query_size']) / \
                     (filtered['term_size'] / filtered['effective_domain_size'])
    filtered['log2_fe'] = np.log2(filtered['fe'])
    filtered['neg_log10_pval'] = -np.log10(
        filtered['p_value'].clip(lower=1e-300))
    filtered['lineageLabel'] = filtered['lineage'].map(name_map)

    # Build native ID -> display name mapping
    id_to_name = filtered.drop_duplicates('native').set_index('native')['name'].to_dict()

    for ln in lineage_names:
        count = filtered[filtered['lineageLabel'] == ln]['native'].nunique()
        print(f"  Lineage {ln}: {count} terms")

    return filtered, lineage_names, source_label, id_to_name


# ── combined heatmap ─────────────────────────────────────────────────────────

def plot_heatmap(data, lineages, ordered_terms, term_to_lineages,
                 source_label, top_n, output_file, id_to_name):
    """Combined heatmap: rows = terms (native IDs), columns = lineages."""
    all_terms = [t for t in ordered_terms if t in data['native'].values]

    if len(all_terms) == 0:
        print(f"No terms to plot for '{source_label}'. Skipping heatmap.")
        return

    matrix = pd.DataFrame(index=all_terms, columns=lineages, dtype=float)
    pval_matrix = pd.DataFrame(index=all_terms, columns=lineages, dtype=float)

    for _, row in data.iterrows():
        if row['native'] in matrix.index:
            matrix.loc[row['native'], row['lineageLabel']] = row['log2_fe']
            pval_matrix.loc[row['native'], row['lineageLabel']] = row['p_value']

    # Display labels: use human-readable name for y-axis
    display_labels = [id_to_name.get(t, t) for t in all_terms]

    fig_height = max(6, len(all_terms) * 0.28 + 1.5)
    fig_width = max(5, len(lineages) * 0.9 + 4)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    mat_data = matrix.values.astype(float)
    # Mask both NaN and zero so those cells render as white
    masked = np.ma.masked_where((np.isnan(mat_data)) | (mat_data == 0), mat_data)

    cmap = plt.cm.Reds.copy()
    cmap.set_bad(color='white')

    vmin = 0
    vmax = np.nanmax(mat_data) if not np.all(np.isnan(mat_data)) else 1
    im = ax.imshow(masked, aspect='auto', cmap=cmap, interpolation='nearest',
                   vmin=vmin, vmax=vmax)

    # Significance asterisks
    """  for i, term in enumerate(all_terms):
        for j, lin in enumerate(lineages):
            pv = pval_matrix.loc[term, lin]
            if pd.notna(pv) and pv < 0.05:
                ax.text(j, i, '*', ha='center', va='center',
                        fontsize=7, fontweight='bold', color='black') """

    ax.set_xticks(range(len(lineages)))
    ax.set_xticklabels(lineages, fontsize=10, fontweight='bold')
    ax.set_xlabel('Lineage', fontsize=11)
    ax.set_yticks(range(len(all_terms)))
    ax.set_yticklabels(display_labels, fontsize=7)
    ax.set_title(
        f'{source_label} Enrichment Heatmap ({f"top {top_n}/lineage" if top_n else "all terms"})\n'
        f'color = log$_2$(Fold Enrichment),  * = p < 0.05',
        fontsize=11, fontweight='bold')

    # Group separators
    prev_lins = None
    for i, term in enumerate(all_terms):
        cur_lins = term_to_lineages.get(term)
        if prev_lins is not None and cur_lins != prev_lins:
            ax.axhline(y=i - 0.5, color='black', linewidth=0.6,
                       linestyle='--', alpha=0.4)
        prev_lins = cur_lins

    cbar = plt.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label('log$_2$(Fold Enrichment)', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved heatmap: {output_file}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Publication-quality lineage enrichment heatmap "
                    "from g:Profiler results")
    parser.add_argument("input_file", help="g:Profiler CSV file")
    parser.add_argument("--source", default="GO:BP",
                        help='Source to filter (e.g. GO:BP, GO:CC, GO:MF, '
                             'HP, KEGG, REAC, TF, MIRNA, All). '
                             'Default: GO:BP')
    parser.add_argument("--top", type=int, default=None,
                        help='Top N terms per lineage. '
                             '0 or omit for all terms (default: all)')
    parser.add_argument("--min-fe", type=float, default=None,
                        help='Min log2(fold enrichment) cutoff. '
                             'Removes terms below this in ALL lineages '
                             '(default: no cutoff)')
    parser.add_argument("--rank-by", choices=['pvalue', 'fe'],
                        default='pvalue',
                        help='Rank top N terms by p-value or fold enrichment '
                             '(default: pvalue)')
    parser.add_argument("-o", "--output", default=None,
                        help='Output path for combined heatmap. '
                             'Default: <base>.<source>_heatmap.png')
    parser.add_argument("--all-fe-csv", default=None,
                        help='Output CSV of ALL terms in the input file with '
                             'fe and log2_fe columns. '
                             'Default: <heatmap_base>.allFE.csv')
    args = parser.parse_args()

    input_file = args.input_file
    source_filter = args.source
    top_n = args.top if args.top else None  # treat 0 as None (all terms)
    base = os.path.splitext(os.path.basename(input_file))[0]
    out_dir = os.path.dirname(os.path.abspath(input_file))
    source_tag = source_filter.replace(':', '_')

    # Determine output file path early so allFE CSV can share the same base
    if args.output:
        output_file = args.output
    else:
        top_tag = f"_top{top_n}" if top_n else "_all"
        output_file = os.path.join(
            out_dir, f"{base}.{source_tag}{top_tag}_heatmap.png")

    # ── load & filter ────────────────────────────────────────────────────
    data, lineages, source_label, id_to_name = prepare_data(input_file, source_filter)

    # ── save all-terms FE CSV ────────────────────────────────────────────
    all_fe_file = args.all_fe_csv if args.all_fe_csv else \
        os.path.splitext(output_file)[0] + ".allFE.csv"
    drop_cols = ['neg_log10_pval', 'lineageLabel']
    move_cols = ['fe', 'log2_fe']
    all_base_cols = [c for c in data.columns
                     if c not in drop_cols and c not in move_cols]
    data[all_base_cols + move_cols].to_csv(all_fe_file, index=False)
    print(f"Saved all-terms FE CSV: {all_fe_file} ({len(data)} rows)")

    # Fold enrichment filter (all-lineages: remove term only if below
    # cutoff in EVERY lineage, preserving cross-lineage comparisons)
    if args.min_fe is not None:
        max_fe = data.groupby('native')['log2_fe'].max()
        keep_terms = set(max_fe[max_fe >= args.min_fe].index)
        before = len(data['native'].unique())
        data = data[data['native'].isin(keep_terms)].copy()
        after = len(data['native'].unique())
        print(f"  After log2 FE >= {args.min_fe} filter: {after} terms "
              f"(removed {before - after})")

    # ── select top N per lineage (union) ─────────────────────────────────
    rank_by = args.rank_by
    if top_n is None:
        selected = set(data['native'].unique())
        print(f"  Selected {len(selected)} unique terms (all, no top-N filter)")
    else:
        selected = set()
        for lin in lineages:
            ld = data[data['lineageLabel'] == lin]
            if rank_by == 'fe':
                lin_data = ld.nlargest(top_n, 'log2_fe')
            else:
                lin_data = ld.nsmallest(top_n, 'p_value')
            selected.update(lin_data['native'].tolist())
        print(f"  Selected {len(selected)} unique terms "
              f"(top {top_n} per lineage by {rank_by}, union)")

    # Build term sets for overlap classification
    term_sets = {lin: set(data[data['lineageLabel'] == lin]['native']) & selected
                 for lin in lineages}
    # Min p-value per term (across lineages) for sorting within groups
    term_pval = data[data['native'].isin(selected)].groupby('native')['p_value'].min().to_dict()
    ordered_terms, term_to_lineages = classify_terms_by_overlap(
        term_sets, lineages, term_pval=term_pval)

    # ── plot heatmap ────────────────────────────────────────────────────
    # ── save CSV of selected terms with fe and log2_fe columns ──────────
    csv_file = os.path.splitext(output_file)[0] + ".csv"
    csv_data = data[data['native'].isin(selected)].copy()
    # Move fe and log2_fe to the end, drop internal helper columns
    drop_cols = ['neg_log10_pval', 'lineageLabel']
    move_cols = ['fe', 'log2_fe']
    base_cols = [c for c in csv_data.columns
                 if c not in drop_cols and c not in move_cols]
    csv_data = csv_data[base_cols + move_cols]
    csv_data.to_csv(csv_file, index=False)
    print(f"Saved CSV: {csv_file} ({len(csv_data)} rows)")
    plot_heatmap(data, lineages, ordered_terms, term_to_lineages,
                 source_label, top_n, output_file, id_to_name)

    print("\nDone.")


if __name__ == "__main__":
    main()
