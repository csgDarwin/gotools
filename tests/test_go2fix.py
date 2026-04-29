import shutil
import subprocess
import sys
from pathlib import Path

FIXTURE = Path(__file__).parent / "fixtures" / "small.maf"
EXAMPLES = Path(__file__).parent.parent / "examples"
GO2FIX = shutil.which("go2fix") or str(Path(sys.prefix) / "bin" / "go2fix")


def _run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)


def test_go2fix_cli_smoke(tmp_path):
    out = tmp_path / "out.bed"
    _run([GO2FIX, str(FIXTURE), "-m", "2", "-o", str(out), "--single-threaded"], cwd=tmp_path)
    assert out.exists(), "go2fix did not create output file"
    assert out.stat().st_size > 0, "go2fix output is empty"
    assert sum(1 for _ in out.open()) > 0, "go2fix output has zero lines"


def test_go2fix_module_smoke(tmp_path):
    out = tmp_path / "out.bed"
    _run(
        [sys.executable, "-m", "gotools.go2fix", str(FIXTURE), "-m", "2",
         "-o", str(out), "--single-threaded"],
        cwd=tmp_path,
    )
    assert out.exists() and out.stat().st_size > 0


def test_go2fix_hprc_matches_expected(tmp_path):
    """End-to-end: go2fix on the HPRC example must reproduce the committed BED byte-for-byte."""
    maf = EXAMPLES / "hprc.example.maf"
    expected = EXAMPLES / "hprc.example.expected.bed"
    if not maf.exists() or not expected.exists():
        import pytest
        pytest.skip("HPRC example fixture not present")

    out = tmp_path / "hprc.bed"
    _run([GO2FIX, str(maf), "-m", "464", "-o", str(out), "--single-threaded", "--quiet"],
         cwd=tmp_path)

    actual = out.read_text()
    expected_text = expected.read_text()
    assert actual == expected_text, "go2fix HPRC output diverged from examples/hprc.example.expected.bed"
