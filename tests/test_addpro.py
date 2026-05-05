import shutil
import subprocess
import sys
from pathlib import Path

FIXTURE = Path(__file__).parent / "fixtures" / "small.gtf"
ADDPRO = shutil.which("addpro") or str(Path(sys.prefix) / "bin" / "addpro")


def _run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)


def test_addpro_cli_smoke(tmp_path):
    out = tmp_path / "out.bed"
    _run([ADDPRO, str(FIXTURE), "-o", str(out)], cwd=tmp_path)
    assert out.exists(), "addpro did not create output file"
    assert out.stat().st_size > 0, "addpro output is empty"


def test_addpro_module_smoke(tmp_path):
    out = tmp_path / "out.bed"
    _run(
        [sys.executable, "-m", "gotools.addpro", str(FIXTURE), "-o", str(out)],
        cwd=tmp_path,
    )
    assert out.exists() and out.stat().st_size > 0


def test_addpro_refseq_rename(tmp_path):
    """NC_000001.11 accessions must be renamed to chr1 in the output."""
    out = tmp_path / "out.bed"
    _run([ADDPRO, str(FIXTURE), "-o", str(out)], cwd=tmp_path)
    chroms = {line.split("\t")[0] for line in out.read_text().splitlines() if line}
    assert "chr1" in chroms, "NC_000001.11 was not renamed to chr1"
    assert "NC_000001.11" not in chroms, "Raw RefSeq accession still present in output"


def test_addpro_promoter_coords(tmp_path):
    """Verify gene-body and promoter BED intervals for both strands.

    DDX11L1  NC_000001.11 + [11874,14409] -> chr1:11873-14409 (gene), chr1:10873-11873 (promoter)
    WASH7P   NC_000001.11 - [14362,29370] -> chr1:14361-29370 (gene), chr1:29370-30370 (promoter)
    """
    out = tmp_path / "out.bed"
    _run([ADDPRO, str(FIXTURE), "-o", str(out)], cwd=tmp_path)

    rows = [line.split("\t") for line in out.read_text().splitlines() if line]
    intervals = {(r[0], int(r[1]), int(r[2])) for r in rows}

    # DDX11L1 - + strand
    assert ("chr1", 11873, 14409) in intervals, "DDX11L1 gene body wrong"
    assert ("chr1", 10873, 11873) in intervals, "DDX11L1 promoter wrong"

    # WASH7P - - strand
    assert ("chr1", 14361, 29370) in intervals, "WASH7P gene body wrong"
    assert ("chr1", 29370, 30370) in intervals, "WASH7P promoter wrong"


def test_addpro_mode_merged(tmp_path):
    """In 'merged' mode, gene body and 1 kb upstream are emitted as one row per gene.

    DDX11L1  + strand -> chr1:10873-14409  (promoter_start .. gene_end)
    WASH7P   - strand -> chr1:14361-30370  (gene_start .. promoter_end)
    """
    out = tmp_path / "merged.bed"
    _run([ADDPRO, str(FIXTURE), "--mode", "merged", "-o", str(out)], cwd=tmp_path)

    rows = [line.split("\t") for line in out.read_text().splitlines() if line]
    intervals = {(r[0], int(r[1]), int(r[2]), r[5]) for r in rows}

    assert ("chr1", 10873, 14409, "+") in intervals, "DDX11L1 merged interval wrong"
    assert ("chr1", 14361, 30370, "-") in intervals, "WASH7P merged interval wrong"


def test_addpro_mode_both(tmp_path):
    """'both' writes two files at <prefix>.separate.bed and <prefix>.merged.bed.

    Row count: separate == 2 * merged (every gene contributes two rows in
    separate, one in merged).
    """
    prefix = tmp_path / "run1"
    _run([ADDPRO, str(FIXTURE), "--mode", "both", "-o", str(prefix)], cwd=tmp_path)

    sep = prefix.with_name("run1.separate.bed")
    mer = prefix.with_name("run1.merged.bed")
    assert sep.exists() and sep.stat().st_size > 0, "separate.bed missing"
    assert mer.exists() and mer.stat().st_size > 0, "merged.bed missing"

    n_sep = sum(1 for line in sep.read_text().splitlines() if line)
    n_mer = sum(1 for line in mer.read_text().splitlines() if line)
    assert n_sep == 2 * n_mer, f"separate ({n_sep}) should be 2x merged ({n_mer})"


def test_addpro_upstream_param(tmp_path):
    """--upstream N narrows the promoter window to exactly N bp."""
    out = tmp_path / "out_500.bed"
    _run([ADDPRO, str(FIXTURE), "--upstream", "500", "-o", str(out)], cwd=tmp_path)

    rows = [line.split("\t") for line in out.read_text().splitlines() if line]
    intervals = {(r[0], int(r[1]), int(r[2])) for r in rows}

    # Gene bodies unchanged
    assert ("chr1", 11873, 14409) in intervals, "DDX11L1 gene body wrong"
    assert ("chr1", 14361, 29370) in intervals, "WASH7P gene body wrong"

    # 500 bp promoters, strand-aware
    assert ("chr1", 11373, 11873) in intervals, "DDX11L1 500 bp promoter wrong"
    assert ("chr1", 29370, 29870) in intervals, "WASH7P 500 bp promoter wrong"
