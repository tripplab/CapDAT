import json
import subprocess
from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "capa.py"


def _base_rows(h_values=(0.5,), r_values=(1.0,), include_h0=False):
    rows = []
    for h in h_values:
        for r in r_values:
            for ci, cap in enumerate(("capA", "capB", "capC", "capD")):
                for fi, fold in enumerate(("2_0", "3_0", "5_0")):
                    row = {
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
                        "H2": 2 + ci + 2 * fi,
                        "patch_elems": 100,
                        "size_h": 0.1,
                        "E": 1000,
                    }
                    if include_h0:
                        row["H0"] = 20 + ci + fi
                    rows.append(row)
    return rows


def _run_phase9(tmp_path, rows, extra_args=None):
    input_csv = tmp_path / "input.csv"
    outdir = tmp_path / "results"
    pd.DataFrame(rows).to_csv(input_csv, index=False)
    cmd = ["python3", str(SCRIPT), "--phase", "9", "-i", str(input_csv), "--outdir", str(outdir), "--perm-n", "20", "--bootstrap-n", "20", "--seed", "123", "--no-plots"]
    if extra_args:
        cmd.extend(extra_args)
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=120)
    return result, outdir


def test_phase09_valid_dataset_writes_required_outputs(tmp_path):
    result, outdir = _run_phase9(tmp_path, _base_rows())
    assert result.returncode == 0, result.stderr + result.stdout
    topo = pd.read_csv(outdir / "phase09_topology_table.csv")
    assert {"H1_norm_c", "H2_norm_c"}.issubset(topo.columns)
    sums = topo.groupby(["capside", "h", "patch_R"])[["H1_norm_c", "H2_norm_c"]].sum().abs()
    assert (sums < 1e-10).all().all()
    mantel = pd.read_csv(outdir / "phase09_mantel_summary.csv")
    assert set(mantel["comparison_id"]) == {"D_mech_vs_D_geom_simple", "D_mech_vs_D_topo", "D_mech_vs_D_geom_topo", "D_geom_simple_vs_D_topo"}
    assert (outdir / "phase09_D_topo_h_0p5_R_1p0.csv").exists()
    assert (outdir / "phase09_D_geom_topo_h_0p5_R_1p0.csv").exists()
    assert (outdir / "phase09_geom_vs_geom_topo_comparison.csv").exists()
    report = json.loads((outdir / "phase09_report.json").read_text())
    assert report["condition_cols"] == ["h", "patch_R"]
    assert report["topology_variables_primary"] == ["H1_norm", "H2_norm"]


def test_phase09_h0_absent_is_warning_not_failure(tmp_path):
    result, outdir = _run_phase9(tmp_path, _base_rows())
    assert result.returncode == 0, result.stderr + result.stdout
    report = json.loads((outdir / "phase09_report.json").read_text())
    assert report["status"] in {"PASS", "WARN"}
    assert "H0_not_used" in report["warnings"]


def test_phase09_h0_can_be_included_when_requested(tmp_path):
    result, outdir = _run_phase9(tmp_path, _base_rows(include_h0=True), ["--include-H0"])
    assert result.returncode == 0, result.stderr + result.stdout
    report = json.loads((outdir / "phase09_report.json").read_text())
    assert "H0_norm" in report["topology_variables_used"]
    assoc = pd.read_csv(outdir / "phase09_topology_centered_associations.csv")
    assert "H0_norm_c" in set(assoc["predictor"])


def test_phase09_zero_topology_variance_drops_variable(tmp_path):
    rows = _base_rows()
    for row in rows:
        row["H1"] = 5
    result, outdir = _run_phase9(tmp_path, rows)
    assert result.returncode == 0, result.stderr + result.stdout
    mantel = pd.read_csv(outdir / "phase09_mantel_summary.csv")
    topo = mantel[mantel["comparison_id"].eq("D_mech_vs_D_topo")].iloc[0]
    assert topo["variables_used_matrix_b"] == "H2_norm"
    assert "topology_variable_dropped_zero_variance" in str(topo["warning_flags"])
    assert "one_dimensional_topology_distance" in str(topo["warning_flags"])


def test_phase09_separates_multiple_h_and_patch_r_conditions(tmp_path):
    result, outdir = _run_phase9(tmp_path, _base_rows(h_values=(0.5, 1.0), r_values=(1.0, 2.0)))
    assert result.returncode == 0, result.stderr + result.stdout
    comparison = pd.read_csv(outdir / "phase09_geom_vs_geom_topo_comparison.csv")
    assert len(comparison) == 4
    assert set(comparison["h"]) == {0.5, 1.0}
    assert set(comparison["patch_R"]) == {1.0, 2.0}
