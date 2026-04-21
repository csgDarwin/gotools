import subprocess
import sys
from pathlib import Path

FIXTURE = Path(__file__).parent / "fixtures" / "small.maf"
GO2FIX = Path(sys.prefix) / "bin" / "go2fix"


def _run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)


def test_go2fix_cli_smoke(tmp_path):
    out = tmp_path / "out.bed"
    _run([str(GO2FIX), str(FIXTURE), "-m", "2", "-o", str(out), "--single-threaded"], cwd=tmp_path)
    assert out.exists(), "go2fix did not create output file"
    assert out.stat().st_size > 0, "go2fix output is empty"
    assert sum(1 for _ in out.open()) > 0, "go2fix output has zero lines"


def test_go2fix_module_smoke(tmp_path):
    out = tmp_path / "out.bed"
    _run(
        [sys.executable, "-m", "gotools.go2fix", str(FIXTURE), "-m", "2", "-o", str(out), "--single-threaded"],
        cwd=tmp_path,
    )
    assert out.exists() and out.stat().st_size > 0
