# Examples

This directory contains a small runnable demo for both `gotools` commands.

## System requirements

**Typical install time:** ~30 seconds on a standard laptop (downloads the `numpy` wheel plus the `gotools` source; faster if `numpy` is already cached in the target environment).

**Operating systems:**
- macOS (tested on Darwin 24.6 / macOS 15)
- Linux (any POSIX distribution with Python ≥ 3.8)
- Windows is not officially supported; may work via WSL.

**Hardware:**
- No non-standard hardware required. A standard desktop or laptop with ≥ 4 GB RAM is sufficient.
- Multi-core CPU recommended — `go2fix` and `go2var` can parallelize across cores.
## Files

| File | Description |
|---|---|
| `small.maf` | 2-block, 9-species MAF fixture (hg38 reference; chimp, bonobo, gorilla, orangutan, gibbon, macaque, marmoset) |
| `small_config.json` | `go2var` config naming `hg38` as the reference species and listing the full species order |

## `go2fix` — find conserved positions

From this `examples/` directory:

```bash
# Require at least 2 matching rows per position; single-threaded for deterministic output.
go2fix small.maf -m 2 -o demo_conserved.bed --single-threaded

head demo_conserved.bed
wc -l demo_conserved.bed
```

**Expected output:** 30 merged BED intervals. First lines:

```
hg38.chr1	10464	10471
hg38.chr1	10472	10474
hg38.chr1	10475	10478
```

**Expected run time:** < 1 second on a standard laptop.

### Other `go2fix` options to try

```bash
# Default parallel run (uses multiple workers).
go2fix small.maf -m 2 -o demo_conserved.bed --workers 4

# Disable adjacent-position merging (emit one record per position).
go2fix small.maf -m 2 -o demo_conserved_unmerged.bed --single-threaded --no-merge
```

## `go2var` — find variant positions

From this `examples/` directory:

```bash
# Single worker for deterministic output.
go2var small.maf small_config.json demo_variants.bed -w 1

head demo_variants.bed
wc -l demo_variants.bed
```

**Expected output:** 141 variant positions. First lines:

```
hg38.chr1	10436	10437
hg38.chr1	10437	10438
hg38.chr1	10438	10439
```

**Expected run time:** < 1 second on a standard laptop.

### Other `go2var` options to try

```bash
# Merge adjacent variants; include gap positions in output.
go2var small.maf small_config.json demo_variants_merged.bed -w 1 --merge --include-gaps

# Verbose mode (adds substitution details and column-wise bases).
go2var small.maf small_config.json demo_variants_verbose.bed -w 1 --verbose
```

## About the `go2var` config

```json
{
  "ReferenceList": ["hg38"],
  "AlignOrderList": ["hg38", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "mSymSyn1", "MFA8", "calJac240"]
}
```

- **`ReferenceList`** — species whose bases must all agree at a position for it to be considered a reference consensus. A block is skipped if fewer than all the species in the reference list are present.
- **`AlignOrderList`** — full species ordering. Only species present in both the MAF block and this list are used for variant calling; they are reordered according to this list before detection.
