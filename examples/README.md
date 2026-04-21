# Examples

## `go2fix` — find conserved positions

```bash
# Require at least 2 matching rows per position; parallel workers.
go2fix input.maf.gz -m 2 -o output.bed --workers 8

# Single-threaded for deterministic ordering (useful for diffs).
go2fix input.maf.gz -m 2 -o output.bed --single-threaded

# Disable adjacent-position merging (emit one record per position).
go2fix input.maf.gz -m 2 -o output.bed --no-merge
```

## `go2var` — find variant positions

```bash
# Default run with JSON config.
go2var input.maf.gz config.json output.bed --workers 8

# Merge adjacent variants; include gap positions in output.
go2var input.maf.gz config.json output.bed --merge --include-gaps
```

Example `config.json`:

```json
{
  "ReferenceList": ["hg38"],
  "AlignOrderList": ["hg38", "panTro6", "gorGor6", "ponAbe3", "macFas5"]
}
```

- `ReferenceList` — species whose bases must all agree at a position for it to be considered a reference consensus.
- `AlignOrderList` — full species ordering; only species present in both the block and this list are used.
