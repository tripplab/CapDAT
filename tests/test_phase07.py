import json
import subprocess
from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "capa.py"


def _base_rows(h_values=(0.5,), r_values=(1.0,), capsides=("capA", "capB", "capC", "capD"), folds=("2_0", "3_0", "5_0")):
    rows = []
    for h in h_values:
        for r in r_values:
            for ci, cap in enumerate(capsides):
                for fi, fold in enumerate(folds):
                    rows.append({
                        "capside": cap,
                        "fold": fold,
                        "h": h,
                        "patch_R": r,
                        "d_mag_max": 1.0 + 0.3 * ci + 0.1 * fi + 0.05 * h + 0.02 * r,
                        "D": 10.0,
                        "t": 2.0 + 0.2 * ci + 0.05 * fi,
                        "Hout": 0.10 + 0.03 * ci + 0.01 * fi,
                        "Hin": -0.20 - 0.02 * ci + 0.015 * fi,
                        "H1": 5 + ci + fi,
                        "H2": 2 + ci,
                        "patch_elems": 100,
                        "size_h": 0.1,
                        "E": 1000,
                    })
    return rows


def _run_phase7(tmp_path, rows, extra_args=None):
    input_csv = tmp_path / "input.csv"
    outdir = tmp_path / "results"
    pd.DataFrame(rows).to_csv(input_csv, index=False)
    cmd = ["python3", str(SCRIPT), "--phase", "7", "-i", str(input_csv), "--outdir", str(outdir), "--perm-n", "30", "--seed", "123", "--no-plots"]
    if extra_args:
        cmd.extend(extra_args)
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=120)
    return result, outdir


def test_phase07_valid_dataset_writes_required_outputs(tmp_path):
    result, outdir = _run_phase7(tmp_path, _base_rows())
    assert result.returncode == 0, result.stderr + result.stdout
    summary = pd.read_csv(outdir / "phase07_mantel_summary.csv")
    assert len(summary) == 1
    assert summary.loc[0, "comparison_id"] == "D_mech_vs_D_geom_simple"
    assert pd.notna(summary.loc[0, "p_perm_two_sided"])
    assert (outdir / "phase07_D_mech_h_0p5_R_1p0.csv").exists()
    assert (outdir / "phase07_D_geom_simple_h_0p5_R_1p0.csv").exists()
    assert (outdir / "phase07_distance_pairs_h_0p5_R_1p0.csv").exists()
    report = json.loads((outdir / "phase07_report.json").read_text())
    assert report["condition_cols"] == ["h", "patch_R"]


def test_phase07_separates_multiple_h_and_patch_r_conditions(tmp_path):
    result, outdir = _run_phase7(tmp_path, _base_rows(h_values=(0.5, 1.0), r_values=(1.0, 2.0)))
    assert result.returncode == 0, result.stderr + result.stdout
    summary = pd.read_csv(outdir / "phase07_mantel_summary.csv")
    assert len(summary) == 4
    assert set(summary["h"]) == {0.5, 1.0}
    assert set(summary["patch_R"]) == {1.0, 2.0}
    for htag in ["0p5", "1p0"]:
        for rtag in ["1p0", "2p0"]:
            assert (outdir / f"phase07_D_mech_h_{htag}_R_{rtag}.csv").exists()


def test_phase07_zero_geometry_variance_drops_variable(tmp_path):
    rows = _base_rows()
    for row in rows:
        row["t"] = 3.0
    result, outdir = _run_phase7(tmp_path, rows)
    assert result.returncode == 0, result.stderr + result.stdout
    summary = pd.read_csv(outdir / "phase07_mantel_summary.csv")
    assert summary.loc[0, "variables_used_matrix_b"] == "Hout;Hin"
    assert "geometry_variable_dropped" in str(summary.loc[0, "warning_flags"])


def test_phase07_zero_mechanical_variance_invalid_without_crash(tmp_path):
    rows = _base_rows()
    for row in rows:
        row["d_mag_max"] = 1.0
    result, outdir = _run_phase7(tmp_path, rows)
    assert result.returncode != 0
    summary = pd.read_csv(outdir / "phase07_mantel_summary.csv")
    assert not bool(summary.loc[0, "valid_comparison"])
    assert "zero_mechanical_distance_variance" in str(summary.loc[0, "warning_flags"])


def test_phase07_duplicate_observation_key_fails(tmp_path):
    rows = _base_rows()
    rows.append(dict(rows[0]))
    result, outdir = _run_phase7(tmp_path, rows)
    assert result.returncode != 0
    report = json.loads((outdir / "phase07_report.json").read_text())
    assert report["status"] == "FAIL"
    assert "duplicate_observation_key" in report["severe_flags"]
