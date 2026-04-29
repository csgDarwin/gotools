# Releasing `gotools`

This document describes how to cut a release. The single source of truth for
the version is `src/gotools/__init__.py` (read by hatch via `dynamic = ["version"]`
in `pyproject.toml`).

## 1. Pre-release checklist

- [ ] CI is green on `main`.
- [ ] `CHANGELOG.md` has an entry for the new version under `## [X.Y.Z]`.
- [ ] `CITATION.cff` `version:` and `date-released:` are updated.
- [ ] README install commands reference the new tag (`@vX.Y.Z`).

## 2. Bump the version

Edit `src/gotools/__init__.py`:

```python
__version__ = "X.Y.Z"
```

Commit:

```bash
git add src/gotools/__init__.py CHANGELOG.md CITATION.cff README.md
git commit -m "Release vX.Y.Z"
```

## 3. Tag and push

```bash
git tag -a vX.Y.Z -m "Release vX.Y.Z"
git push origin main --follow-tags
```

The `release.yml` GitHub workflow runs on every `v*` tag push:

1. Builds an sdist and a wheel with `python -m build`.
2. Uploads them to PyPI via [Trusted Publishing](https://docs.pypi.org/trusted-publishers/)
   (no API token required).
3. Attaches the artifacts to a GitHub Release named `vX.Y.Z`.

## 4. One-time PyPI Trusted Publisher setup

Before the first release, register the project on PyPI and configure a Trusted
Publisher:

1. Reserve the project name on PyPI: <https://pypi.org/manage/projects/>
   (or use TestPyPI first: <https://test.pypi.org/manage/projects/>).
2. Add a Trusted Publisher with these values:
   - **Owner:** `csgDarwin`
   - **Repository:** `gotools`
   - **Workflow:** `release.yml`
   - **Environment:** `pypi`
3. Configure the matching `pypi` environment in the GitHub repo settings:
   `Settings -> Environments -> New environment -> pypi`. No secrets needed —
   OIDC handles authentication.

## 5. conda-forge

`gotools` is pure Python with only NumPy as a runtime dependency, so the
conda-forge path is straightforward and runs after the PyPI release exists:

1. Fork <https://github.com/conda-forge/staged-recipes>.
2. Add `recipes/gotools/meta.yaml` pointing at the PyPI sdist for `vX.Y.Z`.
   Use `grayskull pypi gotools` to scaffold the recipe.
3. Open a PR to `staged-recipes`. Once merged, `conda-forge/gotools-feedstock`
   is created and subsequent PyPI releases are picked up automatically by the
   conda-forge bot.
