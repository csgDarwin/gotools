# gotools

Fast conserved-site and variant detection from MAF multiple sequence alignments.

`gotools` provides two command-line tools:

- **`go2fix`** — scans a MAF alignment for conserved positions (identical base across ≥ N rows, no Ns, no gaps) and emits a BED file of conserved intervals.
- **`go2var`** — scans a MAF alignment for positions where non-reference species diverge from the reference, using a JSON config to define species threshold and ordering, and emits a BED file of variant intervals.

Both tools use a custom NumPy-vectorized MAF parser; the only third-party dependency is `numpy`.

---

## Install

### Option 1 — install directly from GitHub (recommended)

```bash
pip install git+https://github.com/csgDarwin/gotools.git@v0.1.0
```

After install, the commands `go2fix` and `go2var` are on your PATH.

### Option 2 — clone and install editable

```bash
git clone https://github.com/csgDarwin/gotools.git
cd gotools
pip install -e .
```

### Option 3 — standalone single-file download (no install)

Both scripts are fully self-contained — `numpy` is the only runtime dependency.

```bash
# download either script directly
wget https://raw.githubusercontent.com/csgDarwin/gotools/v0.1.0/src/gotools/go2fix.py
pip install numpy
python go2fix.py input.maf -o output.bed
```

PyPI and bioconda install paths are planned for a future release.

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
- `--no-merge` — emit one record per position instead of merging adjacent conserved bases

See `go2fix --help` for full options.

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

## Input / output formats

- **Input (both tools):** MAF (Multiple Alignment Format), plain text or gzip-compressed. The reference (first) sequence in each block defines the coordinate system.
- **Output (both tools):** BED (0-based, half-open intervals). `go2fix` output columns: `chrom`, `start`, `end`, `base`. `go2var` output columns: `chrom`, `start`, `end`, optional substitution details.

---

## Citation

TODO: Add citation once published.

---

## License

MIT — see [LICENSE](LICENSE).
