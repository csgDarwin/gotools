# Changelog

All notable changes to `gotools` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-05-05

### Added
- `addpro --mode {separate,merged,both}` (default `separate`).
- `addpro --upstream N` (default 1000) for configurable promoter width.

### Changed
- `addpro` BED `name` column is now a gene symbol (`gene_name` → `gene` →
  `gene_id`); duplicate symbols suffixed with `|<gene_id>`.

## [0.1.1] - 2026-05-03

### Fixed
- `go2var` no longer emits variant calls at reference positions where the
  first reference base is `N`.

## [0.1.0] - 2026-04-28

### Added
- `go2fix` — fast conserved-site detection from MAF multiple sequence alignments.
  Vectorized NumPy core, multiprocess pipeline with pipelined parser → workers → writer.
- `go2var` — variant-position calling from MAF alignments using a JSON config to
  declare reference species and alignment ordering.
- Console entry points `go2fix` and `go2var` installed alongside the package.
- Public Python API: `gotools.go2fix_optimized`, `gotools.go2fix_single`,
  `gotools.go2var_sorted_optimized`.
- Example MAF fixtures (`examples/small.maf`, `examples/hprc.example.maf`) with
  expected BED outputs for end-to-end verification.
- Smoke and correctness tests under `tests/`.
- GitHub Actions CI (pytest + ruff on Python 3.10 and 3.12).
- PyPI release workflow via Trusted Publishing on `v*` tag push.

[Unreleased]: https://github.com/csgDarwin/gotools/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/csgDarwin/gotools/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/csgDarwin/gotools/releases/tag/v0.1.0
