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
                    spread = (ci + 1) * fi
                    rows.append({
                        "capside": cap,
                        "fold": fold,
                        "h": h,
                        "patch_R": r,
                        "d_mag_max": 1.0 + 0.2 * ci + 0.12 * spread + 0.03 * h + 0.01 * r,
                        "D": 10.0,
                        "t": 2.0 + 0.1 * ci + 0.04 * spread,
                        "Hout": 0.10 + 0.03 * ci + 0.01 * spread,
                        "Hin": -0.20 - 0.02 * ci + 0.015 * spread,
                        "H1": 5 + ci + fi,
                        "H2": 2 + ci,
                        "patch_elems": 100,
                        "size_h": 0.1,
                        "E": 1000,
                    })
    return rows


def _run_phase8(tmp_path, rows, extra_args=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    input_csv = tmp_path / "input.csv"
    outdir = tmp_path / "results"
    pd.DataFrame(rows).to_csv(input_csv, index=False)
    cmd = ["python3", str(SCRIPT), "--phase", "8", "-i", str(input_csv), "--outdir", str(outdir), "--perm-n", "20", "--bootstrap-n", "20", "--seed", "123", "--no-plots"]
    if extra_args:
        cmd.extend(extra_args)
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=120)
    return result, outdir


def test_phase08_valid_dataset_writes_required_outputs(tmp_path):
    result, outdir = _run_phase8(tmp_path, _base_rows())
    assert result.returncode == 0, result.stderr + result.stdout
    aniso = pd.read_csv(outdir / "phase08_capsid_anisotropy.csv")
    assert len(aniso) == 4
    assert {"A_mech", "A_mech_pct", "A_t", "A_Hout", "A_Hin"}.issubset(aniso.columns)
    assert aniso[["A_mech", "A_t", "A_Hout", "A_Hin"]].notna().all().all()
    comp = pd.read_csv(outdir / "phase08_anisotropy_comparison_summary.csv")
    assert set(comp["comparison_id"]) == {"A_mech_vs_A_t", "A_mech_vs_A_Hout", "A_mech_vs_A_Hin"}
    ranks = pd.read_csv(outdir / "phase08_capsid_anisotropy_rankings.csv")
    assert set(ranks["metric"]) == {"A_mech", "A_t", "A_Hout", "A_Hin"}
    top = pd.read_csv(outdir / "phase08_top_rank_coincidence.csv")
    assert len(top) == 3
    report = json.loads((outdir / "phase08_report.json").read_text())
    assert report["anisotropy_group_cols"] == ["capside", "h", "patch_R"]
    assert report["condition_cols"] == ["h", "patch_R"]


def test_phase08_separates_multiple_h_conditions(tmp_path):
    result, outdir = _run_phase8(tmp_path, _base_rows(h_values=(0.5, 1.0)))
    assert result.returncode == 0, result.stderr + result.stdout
    aniso = pd.read_csv(outdir / "phase08_capsid_anisotropy.csv")
    comp = pd.read_csv(outdir / "phase08_anisotropy_comparison_summary.csv")
    assert set(aniso["h"]) == {0.5, 1.0}
    assert len(aniso) == 8
    assert len(comp) == 6


def test_phase08_separates_multiple_patch_r_conditions(tmp_path):
    result, outdir = _run_phase8(tmp_path, _base_rows(r_values=(1.0, 2.0)))
    assert result.returncode == 0, result.stderr + result.stdout
    aniso = pd.read_csv(outdir / "phase08_capsid_anisotropy.csv")
    comp = pd.read_csv(outdir / "phase08_anisotropy_comparison_summary.csv")
    assert set(aniso["patch_R"]) == {1.0, 2.0}
    assert len(aniso) == 8
    assert len(comp) == 6


def test_phase08_two_folds_warns_but_computes(tmp_path):
    rows = [r for r in _base_rows() if not (r["capside"] == "capD" and r["fold"] == "5_0")]
    result, outdir = _run_phase8(tmp_path, rows)
    assert result.returncode == 0, result.stderr + result.stdout
    aniso = pd.read_csv(outdir / "phase08_capsid_anisotropy.csv")
    capd = aniso[aniso["capside"] == "capD"].iloc[0]
    assert pd.notna(capd["A_mech"])
    assert "low_fold_count_for_anisotropy" in str(capd["phase08_warning_flags"])


def test_phase08_one_valid_hin_fold_only_keeps_other_metrics(tmp_path):
    rows = _base_rows()
    for row in rows:
        if row["capside"] == "capC" and row["fold"] != "2_0":
            row["Hin"] = None
    result, outdir = _run_phase8(tmp_path, rows)
    assert result.returncode == 0, result.stderr + result.stdout
    aniso = pd.read_csv(outdir / "phase08_capsid_anisotropy.csv")
    capc = aniso[aniso["capside"] == "capC"].iloc[0]
    assert pd.isna(capc["A_Hin"])
    assert pd.notna(capc["A_mech"])
    assert "insufficient_folds_for_anisotropy" in str(capc["phase08_warning_flags"])


def test_phase08_zero_curvature_denominator_is_flagged(tmp_path):
    rows = _base_rows()
    for row in rows:
        row["Hout"] = 0.0
    result, outdir = _run_phase8(tmp_path, rows)
    assert result.returncode == 0, result.stderr + result.stdout
    aniso = pd.read_csv(outdir / "phase08_capsid_anisotropy.csv")
    assert aniso["A_Hout"].isna().all()
    assert aniso["phase08_warning_flags"].str.contains("invalid_anisotropy_denominator").any()


def test_phase08_zero_anisotropy_variance_invalidates_comparison(tmp_path):
    rows = _base_rows()
    for row in rows:
        row["t"] = {"2_0": 2.0, "3_0": 3.0, "5_0": 4.0}[row["fold"]]
    result, outdir = _run_phase8(tmp_path, rows)
    assert result.returncode == 0, result.stderr + result.stdout
    comp = pd.read_csv(outdir / "phase08_anisotropy_comparison_summary.csv")
    row = comp[comp["comparison_id"] == "A_mech_vs_A_t"].iloc[0]
    assert not bool(row["valid_comparison"])
    assert "zero_anisotropy_variance" in str(row["warning_flags"])


def test_phase08_duplicate_observation_key_fails(tmp_path):
    rows = _base_rows()
    rows.append(dict(rows[0]))
    result, outdir = _run_phase8(tmp_path, rows)
    assert result.returncode != 0
    report = json.loads((outdir / "phase08_report.json").read_text())
    assert report["status"] == "FAIL"
    assert "duplicate_observation_key" in report["severe_flags"]


def test_phase08_missing_geometry_values_are_variable_wise(tmp_path):
    rows = _base_rows()
    for row in rows:
        if row["capside"] == "capB" and row["fold"] == "3_0":
            row["Hout"] = None
        if row["capside"] == "capD" and row["fold"] == "5_0":
            row["Hin"] = None
    result, outdir = _run_phase8(tmp_path, rows)
    assert result.returncode == 0, result.stderr + result.stdout
    aniso = pd.read_csv(outdir / "phase08_capsid_anisotropy.csv")
    assert aniso["A_mech"].notna().all()
    assert aniso["A_t"].notna().all()
    assert "nonfinite_Hout" in ";".join(aniso["phase08_warning_flags"].fillna("").astype(str))
    assert "nonfinite_Hin" in ";".join(aniso["phase08_warning_flags"].fillna("").astype(str))


def test_phase08_reproducible_tables_across_repeated_runs(tmp_path):
    rows = _base_rows(h_values=(0.5, 1.0), r_values=(1.0,))
    result1, outdir1 = _run_phase8(tmp_path / "run1", rows)
    result2, outdir2 = _run_phase8(tmp_path / "run2", rows)
    assert result1.returncode == 0, result1.stderr + result1.stdout
    assert result2.returncode == 0, result2.stderr + result2.stdout
    for filename in [
        "phase08_capsid_anisotropy.csv",
        "phase08_anisotropy_comparison_summary.csv",
        "phase08_capsid_anisotropy_rankings.csv",
        "phase08_top_rank_coincidence.csv",
    ]:
        assert (outdir1 / filename).read_text() == (outdir2 / filename).read_text()
