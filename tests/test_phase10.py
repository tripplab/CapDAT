import json
import subprocess
import sys
from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")


def _write_dataset(path: Path, capsids=("c1", "c2", "c3", "c4"), include_topology=True, h_values=(1.0,), patch_values=(2.0,)):
    rows = []
    folds = ["2_0", "3_0", "5_0"]
    for h in h_values:
        for patch_R in patch_values:
            for ci, cap in enumerate(capsids):
                for fi, fold in enumerate(folds):
                    row = {
                        "capside": cap,
                        "fold": fold,
                        "h": h,
                        "patch_R": patch_R,
                        "d_mag_max": 1.0 + 0.15 * fi + 0.02 * ci,
                        "t": 10.0 + fi + 0.1 * ci,
                        "Hout": 1.0 + 0.1 * fi + 0.01 * ci,
                        "Hin": 0.5 + 0.05 * fi - 0.01 * ci,
                        "patch_elems": 10,
                        "size_h": 0.5,
                        "D": 10.0,
                        "E": 100.0,
                    }
                    if include_topology:
                        row.update({"H1": 4.0 + fi, "H2": 7.0 + ci, "H0": 2.0})
                    else:
                        row.update({"H1": 4.0 + fi, "H2": 7.0 + ci})
                    rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)


def _run_phase10(tmp_path: Path, csv_path: Path):
    outdir = tmp_path / "out"
    cmd = [
        sys.executable,
        "tools/capa.py",
        "--phase",
        "10",
        "-i",
        str(csv_path),
        "--outdir",
        str(outdir),
        "-o",
        str(tmp_path / "normalized.csv"),
        "-r",
        str(tmp_path / "phase01_report.json"),
        "--perm-n",
        "5",
        "--bootstrap-n",
        "5",
        "--no-plots",
    ]
    return subprocess.run(cmd, cwd=Path(__file__).resolve().parents[1], text=True, capture_output=True, check=False), outdir


def test_phase10_writes_required_outputs_and_excludes_each_capsid(tmp_path):
    csv_path = tmp_path / "stable.csv"
    _write_dataset(csv_path)
    result, outdir = _run_phase10(tmp_path, csv_path)
    assert result.returncode in (0, 1), result.stderr
    required = [
        "phase10_loo_run_index.csv",
        "phase10_centered_association_sensitivity.csv",
        "phase10_ranking_sensitivity.csv",
        "phase10_model_sensitivity.csv",
        "phase10_distance_sensitivity.csv",
        "phase10_overall_stability_summary.csv",
        "phase10_dominated_conclusions.csv",
        "phase10_report.json",
    ]
    for name in required:
        assert (outdir / name).exists(), name
    run_index = pd.read_csv(outdir / "phase10_loo_run_index.csv")
    assert set(run_index["excluded_capside"]) == {"NONE", "c1", "c2", "c3", "c4"}
    assert (run_index["run_label"] == "full").sum() == 1


def test_phase10_keeps_multiple_h_and_patch_r_separate(tmp_path):
    csv_path = tmp_path / "multi_condition.csv"
    _write_dataset(csv_path, h_values=(1.0, 2.0), patch_values=(2.0, 3.0))
    result, outdir = _run_phase10(tmp_path, csv_path)
    assert result.returncode in (0, 1), result.stderr
    report = json.loads((outdir / "phase10_report.json").read_text())
    assert sorted(report["h_values"]) == [1.0, 2.0]
    assert sorted(report["patch_R_values"]) == [2.0, 3.0]
    centered = pd.read_csv(outdir / "phase10_centered_association_sensitivity.csv")
    if not centered.empty:
        assert set(centered["h"]) == {1.0, 2.0}
        assert set(centered["patch_R"]) == {2.0, 3.0}
