import subprocess
import sys
from pathlib import Path

FIXTURE = Path(__file__).parent / "fixtures" / "small.maf"
CONFIG = Path(__file__).parent / "fixtures" / "small_config.json"
GO2VAR = Path(sys.prefix) / "bin" / "go2var"


def _run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)


def test_go2var_cli_smoke(tmp_path):
    out = tmp_path / "out.bed"
    _run(
        [str(GO2VAR), str(FIXTURE), str(CONFIG), str(out), "-w", "1"],
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
