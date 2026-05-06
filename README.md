# GOtools

Fast conserved-site and variant detection from MAF multiple sequence alignments.

`GOtools` provides three command-line tools:

- **`go2fix`** — scans a MAF alignment for conserved positions (identical base across ≥ N rows, no Ns, no gaps) and outputs a BED file of conserved intervals.
- **`go2var`** — scans a MAF alignment for positions where non-reference species diverge from the reference (skipping reference `N` and gap positions), using a JSON config to define species threshold and ordering, and outputs a BED file of variant intervals.
- **`addpro`** — parses a GTF, emits gene-body + 1 kb upstream promoter intervals as a BED file, renames RefSeq accessions (NC_*.* → chrN), and drops non-standard chromosomes. Streams output row-by-row so memory use is O(1) regardless of input size.

---

## System requirements

**Software dependencies:**
- Python ≥ 3.8 (tested on 3.8, 3.10, 3.12)
- `numpy` ≥ 1.20 (tested with numpy 2.4.4)

---

## Install

### Option 1 — install directly from GitHub (recommended)

```bash
pip install git+https://github.com/csgDarwin/gotools.git@main
```

After install, the commands `go2fix`, `go2var`, and `addpro` are on your PATH.

### Option 2 — clone and install editable

```bash
git clone https://github.com/csgDarwin/gotools.git
cd gotools
pip install -e .
```

---

## Usage

### `go2fix` — find conserved positions

```bash
go2fix input.maf.gz -m 464 -o output_dir/ --workers 8
```

- `input.maf` or `input.maf.gz` — MAF-format multiple sequence alignment
- `-m / --max-conserved` — minimum number of rows that must match (default 464)
- `-o / --output` — output path (directory or full filename)
- `--workers` — number of worker processes
- `--no-merge` — writes one record per position instead of merging adjacent conserved bases

See `go2fix --help` for full options.

#### Choosing `-m` (`--max-conserved`)

`-m` is the minimum number of MAF rows that must agree at a position for it to be called conserved. **It must match the number of haplotype rows in your MAF block.** If `-m` is set higher than the row count, the output BED will be **empty**.

| Alignment                         | `-m` |
|-----------------------------------|------|
| Full HPRC pangenome               | 464  |
| HPRCv2 subset                     | 363  |

Verify the row count in your MAF before running:

The default (`-m 464`) is intended for the full HPRC pangenome; tune it for any other alignment.

---

### `go2var` — find variant positions

```bash
go2var input.maf.gz config.json output.bed --workers 8
```

- `input.maf` or `input.maf.gz` — MAF alignment
- `config.json` — JSON with `ReferenceList` (species threshold) and `AlignOrderList` (species ordering)
- `output.bed` — output BED file
- `--workers` — number of worker processes
- `--merge` — merge adjacent variant positions


See `go2var --help` for full options.

**Example `config.json`:**
```json
{
  "ReferenceList": ["hg38", "hs1"],
  "AlignOrderList": ["hg38", "hs1", "panTro6", "gorGor6", "ponAbe3"]
}
```
---

### `addpro` — generate gene + 1 kb promoter BED from a GTF

```bash
addpro annotation.gtf -o annotation.1kbpromoter.bed
addpro annotation.gtf.gz -o annotation.1kbpromoter.bed   # gzip supported
addpro annotation.gtf --mode merged -o annotation.gene1kb.bed
addpro annotation.gtf --mode both   -o annotation          # writes annotation.separate.bed + annotation.merged.bed
addpro annotation.gtf --upstream 500 -o annotation.500bp.bed
addpro annotation.gtf --no-filter                          # keep all chromosomes
```

- `input_gtf` — GTF annotation file (plain text or `.gtf.gz`)
- `-o / --output` — output BED path (default: `<input_basename>.1kbpromoter.bed` alongside the input)
- `--mode {separate,merged,both}` — output layout (default: `separate`)
  - `separate` — two rows per gene: gene body + upstream promoter as distinct intervals
  - `merged` — one row per gene: gene body and upstream promoter as a single merged interval
  - `both` — writes two files using `-o` as a prefix: `<prefix>.separate.bed` and `<prefix>.merged.bed`
- `--upstream N` — promoter width in bp (default: 1000)
- `-n / --no-filter` — keep rows whose chromosome name does not start with `chr` (default: drop them)

The BED `name` column contains the gene symbol (`gene_name` attribute, falling back to `gene` then `gene_id`). RefSeq accessions (`NC_000001.11`, `NC_000023.11`, etc.) are automatically renamed to UCSC-style names (`chr1`, `chrX`, etc.) for GRCh38.p14.


---

## Demo

A tiny MAF fixture is bundled with the repository at `examples/small.maf` (a 2-block, 9-species alignment) along with its config at `examples/small_config.json`. After installing `gotools`, run both tools against this fixture to confirm your install is working.

### `go2fix` demo

```bash
cd examples
go2fix small.maf -m 2 -o demo_conserved.bed --single-threaded
head demo_conserved.bed
wc -l demo_conserved.bed
```

**Expected output** — `demo_conserved.bed` should contain **30 BED intervals** of merged conserved positions. The first lines look like:

```
hg38.chr1	10464	10471
hg38.chr1	10472	10474
hg38.chr1	10475	10478
hg38.chr1	10493	10504
hg38.chr1	10507	10517
```

Columns: `chrom`, `start`, `end` (0-based half-open).

**Expected run time:** < 1 second on a standard laptop.

### `go2var` demo

```bash
cd examples
go2var small.maf small_config.json demo_variants.bed -w 1
head demo_variants.bed
wc -l demo_variants.bed
```

**Expected output** — `demo_variants.bed` should contain **141 variant positions** where the reference species (`hg38`) differs from all aligned non-reference species. The first lines look like:

```
hg38.chr1	10436	10437
hg38.chr1	10437	10438
hg38.chr1	10438	10439
hg38.chr1	10439	10440
hg38.chr1	10440	10441
```

**Expected run time:** < 1 second on a standard laptop.

If both commands complete and produce non-empty BED files matching the formats above, the install is working correctly.

---

## Input / output formats

- **Input (both tools):** MAF (Multiple Alignment Format), plain text or gzip-compressed. The reference (first) sequence in each block defines the coordinate system.
- **Output (both tools):** BED (0-based, half-open intervals). `go2fix` output columns: `chrom`, `start`, `end`, `base`. `go2var` output columns: `chrom`, `start`, `end`, optional substitution details.

---

## Citation

If you use `gotools`, please cite the preprint:

> Goh C-S, Davenport MH, Lee C, Jarvis ED. (2026). Fixed human pangenome sequences reveal origins of common human traits. *bioRxiv*. https://doi.org/10.64898/2026.04.30.722098

Software-level metadata (version, release date) is in [`CITATION.cff`](CITATION.cff).

---

## License

MIT — see [LICENSE](LICENSE).
