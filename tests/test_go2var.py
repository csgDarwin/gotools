import shutil
import subprocess
import sys
from pathlib import Path

FIXTURE = Path(__file__).parent / "fixtures" / "small.maf"
CONFIG = Path(__file__).parent / "fixtures" / "small_config.json"
EXAMPLES = Path(__file__).parent.parent / "examples"
GO2VAR = shutil.which("go2var") or str(Path(sys.prefix) / "bin" / "go2var")


def _run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)


def test_go2var_cli_smoke(tmp_path):
    out = tmp_path / "out.bed"
    _run(
        [GO2VAR, str(FIXTURE), str(CONFIG), str(out), "-w", "1"],
        cwd=tmp_path,
    )
    assert out.exists(), "go2var did not create output file"
    assert out.stat().st_size > 0, "go2var output is empty"
    assert sum(1 for _ in out.open()) > 0, "go2var output has zero lines"


def test_go2var_module_smoke(tmp_path):
    out = tmp_path / "out.bed"
    _run(
        [sys.executable, "-m", "gotools.go2var", str(FIXTURE), str(CONFIG), str(out), "-w", "1"],
        cwd=tmp_path,
    )
    assert out.exists() and out.stat().st_size > 0


def test_go2var_small_matches_expected(tmp_path):
    """End-to-end: go2var on small.maf must reproduce the committed BED byte-for-byte."""
    maf = EXAMPLES / "small.maf"
    config = EXAMPLES / "small_config.json"
    expected = EXAMPLES / "small.expected.bed"
    if not (maf.exists() and config.exists() and expected.exists()):
        import pytest
        pytest.skip("small example fixture not present")

    out = tmp_path / "small.bed"
    _run([GO2VAR, str(maf), str(config), str(out), "-w", "1", "--quiet"], cwd=tmp_path)

    actual = out.read_text()
    expected_text = expected.read_text()
    assert actual == expected_text, "go2var small output diverged from examples/small.expected.bed"
