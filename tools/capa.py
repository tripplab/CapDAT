#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
capa.py

Phase 1:
Read, validate, normalize, and audit the master CSV for the
mechanical-geometrical-topological capsid analysis pipeline.

Expected observational unit:
    (capside, h, patch_R, fold)

Required raw columns, conceptually:
    capside, fold, d_mag_max, <t>, Hout, Hin,
    H1, H2, patch_elems, size_h, h, patch_R, D, E

Main derived variables:
    d_max_pct     = 100 * d_mag_max / D
    patch_volume  = patch_elems * size_h^3
    H1_norm       = H1 / patch_volume
    H2_norm       = H2 / patch_volume

Optional:
    H0_norm       = H0 / patch_volume, if H0 exists
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime, timezone
from typing import Any, Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# User-visible configuration
# ---------------------------------------------------------------------

REQUIRED_CANONICAL_COLUMNS = [
    "capside",
    "fold",
    "d_mag_max",
    "t",
    "Hout",
    "Hin",
    "H1",
    "H2",
    "patch_elems",
    "size_h",
    "h",
    "patch_R",
    "D",
    "E",
]

OBSERVATION_KEY = ["capside", "h", "patch_R", "fold"]
PHASE02_GROUP_COLS = ["capside", "h", "patch_R"]
PHASE02_FOLD_SUMMARY_COLS = ["h", "patch_R", "fold"]
PHASE02_REQUIRED_COLUMNS = [
    "capside",
    "fold",
    "h",
    "patch_R",
    "d_mag_max",
    "D",
    "d_max_pct",
]
PHASE02_PRESERVED_COLUMNS = [
    "t",
    "Hout",
    "Hin",
    "H1",
    "H2",
    "H1_norm",
    "H2_norm",
    "H0",
    "H0_norm",
    "patch_elems",
    "size_h",
    "E",
]
DEFAULT_EXPECTED_FOLDS = ["2_0", "2_1", "3_0", "3_1", "5_0"]
PHASE02_RANKING_METHOD = {"method": "dense", "descending": True}

NUMERIC_COLUMNS = [
    "d_mag_max",
    "t",
    "Hout",
    "Hin",
    "H1",
    "H2",
    "patch_elems",
    "size_h",
    "h",
    "patch_R",
    "D",
    "E",
]

OPTIONAL_NUMERIC_COLUMNS = [
    "H0",
]

# Accepted aliases.
# This allows the raw CSV to contain "<t>" but the internal pipeline
# will use the cleaner canonical name "t".
COLUMN_ALIASES = {
    "capside": ["capside", "capsid", "pdbid", "PDBID", "structure", "id"],
    "fold": ["fold", "fold_id", "axis", "symmetry_axis"],
    "d_mag_max": ["d_mag_max", "dmax", "d_max", "max_displacement", "disp_max"],
    "t": ["<t>", "t", "thickness", "mean_t", "t_mean"],
    "Hout": ["Hout", "H_out", "Houter", "H_outer", "mean_Hout"],
    "Hin": ["Hin", "H_in", "Hinner", "H_inner", "mean_Hin"],
    "H0": ["H0", "h0"],
    "H1": ["H1", "h1"],
    "H2": ["H2", "h2"],
    "patch_elems": ["patch_elems", "voxels", "n_voxels", "elements", "n_elements"],
    "size_h": ["size_h", "voxel_size", "voxel_h", "element_size"],
    "h": ["h", "mesh_h", "mesh_resolution", "resolution"],
    "patch_R": ["patch_R", "R", "radius", "patch_radius"],
    "D": ["D", "diameter", "capsid_D", "capside_D"],
    "E": ["E", "young", "young_modulus", "youngs_modulus"],
}


# Existing derived columns that may be present in the CSV and should be
# compared against recomputed values for verification.
DERIVED_COLUMN_ALIASES = {
    "d_max_pct": ["d_max_pct", "dmax_pct", "d_max_percent", "dmax%", "d_max%"],
    "patch_volume": ["patch_volume", "volume_patch", "patch_vol"],
    "H1_norm": ["H1_norm", "H1norm", "H1_norm_volume", "H1_per_volume"],
    "H2_norm": ["H2_norm", "H2norm", "H2_norm_volume", "H2_per_volume"],
    "H0_norm": ["H0_norm", "H0norm", "H0_norm_volume", "H0_per_volume"],
}


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def _clean_name(name: str) -> str:
    """Normalize a column name for matching, without changing meaning."""
    return str(name).strip()


def _match_columns(
    df: pd.DataFrame,
    aliases: Dict[str, List[str]]
) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """
    Match raw CSV columns to canonical column names.

    Returns
    -------
    mapping:
        Dictionary raw_name -> canonical_name.
    collisions:
        Dictionary canonical_name -> list of raw columns if more than one
        raw column maps to the same canonical name.
    """
    raw_cols = list(df.columns)
    raw_lookup = {_clean_name(c).lower(): c for c in raw_cols}

    canonical_to_raw: Dict[str, List[str]] = {}

    for canonical, alias_list in aliases.items():
        hits = []
        for alias in alias_list:
            key = _clean_name(alias).lower()
            if key in raw_lookup:
                hits.append(raw_lookup[key])
        if hits:
            canonical_to_raw[canonical] = hits

    mapping: Dict[str, str] = {}
    collisions: Dict[str, List[str]] = {}

    for canonical, hits in canonical_to_raw.items():
        unique_hits = list(dict.fromkeys(hits))
        if len(unique_hits) > 1:
            collisions[canonical] = unique_hits
        else:
            mapping[unique_hits[0]] = canonical

    return mapping, collisions


def rename_to_canonical(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict]:
    """
    Rename raw columns to canonical names using COLUMN_ALIASES.

    If ambiguous aliases are found, the function reports them but does
    not guess.
    """
    mapping, collisions = _match_columns(df, COLUMN_ALIASES)

    report = {
        "column_mapping": mapping,
        "alias_collisions": collisions,
    }

    if collisions:
        collision_text = json.dumps(collisions, indent=2, ensure_ascii=False)
        raise ValueError(
            "Ambiguous column aliases detected. Fix the CSV headers or alias table.\n"
            f"{collision_text}"
        )

    renamed = df.rename(columns=mapping).copy()

    return renamed, report


def find_existing_derived_columns(df: pd.DataFrame) -> Dict[str, str]:
    """
    Find existing derived columns in the input CSV, if any, so they can
    be compared against recomputed values.
    """
    mapping, collisions = _match_columns(df, DERIVED_COLUMN_ALIASES)

    if collisions:
        raise ValueError(
            "Ambiguous derived-column aliases detected:\n"
            + json.dumps(collisions, indent=2, ensure_ascii=False)
        )

    # mapping is raw -> canonical; invert to canonical -> raw
    return {canonical: raw for raw, canonical in mapping.items()}


def coerce_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Convert expected numeric columns to numeric dtype."""
    out = df.copy()

    for col in NUMERIC_COLUMNS + OPTIONAL_NUMERIC_COLUMNS:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    return out


def validate_required_columns(df: pd.DataFrame) -> List[str]:
    """Return missing required canonical columns."""
    return [col for col in REQUIRED_CANONICAL_COLUMNS if col not in df.columns]


def validate_observation_key(df: pd.DataFrame) -> Dict:
    """
    Check whether the expected observational key uniquely identifies rows:
        (capside, h, patch_R, fold)
    """
    key_report = {
        "observation_key": OBSERVATION_KEY,
        "n_rows": int(len(df)),
        "n_unique_observations": None,
        "n_duplicate_rows": None,
        "duplicate_examples": [],
    }

    missing_key_cols = [c for c in OBSERVATION_KEY if c not in df.columns]
    if missing_key_cols:
        key_report["missing_key_columns"] = missing_key_cols
        return key_report

    duplicated_mask = df.duplicated(subset=OBSERVATION_KEY, keep=False)
    duplicated = df.loc[duplicated_mask, OBSERVATION_KEY]

    key_report["n_unique_observations"] = int(
        df[OBSERVATION_KEY].drop_duplicates().shape[0]
    )
    key_report["n_duplicate_rows"] = int(duplicated_mask.sum())

    if not duplicated.empty:
        key_report["duplicate_examples"] = (
            duplicated.drop_duplicates()
            .head(20)
            .to_dict(orient="records")
        )

    return key_report


def validate_values(df: pd.DataFrame) -> Dict:
    """
    Basic sanity checks:
    - missing values
    - non-finite values
    - non-positive denominators
    - negative counts
    """
    report = {}

    cols_to_check = [c for c in REQUIRED_CANONICAL_COLUMNS if c in df.columns]

    missing_counts = df[cols_to_check].isna().sum().to_dict()
    report["missing_values"] = {
        col: int(n) for col, n in missing_counts.items() if n > 0
    }

    numeric_present = [c for c in NUMERIC_COLUMNS if c in df.columns]
    non_finite = {}
    for col in numeric_present:
        values = df[col].to_numpy(dtype=float)
        bad = np.sum(~np.isfinite(values))
        if bad > 0:
            non_finite[col] = int(bad)
    report["non_finite_values"] = non_finite

    non_positive_checks = {}
    for col in ["D", "patch_elems", "size_h", "t"]:
        if col in df.columns:
            n_bad = int((df[col] <= 0).sum())
            if n_bad > 0:
                non_positive_checks[col] = n_bad
    report["non_positive_values"] = non_positive_checks

    negative_count_checks = {}
    for col in ["H0", "H1", "H2"]:
        if col in df.columns:
            n_bad = int((df[col] < 0).sum())
            if n_bad > 0:
                negative_count_checks[col] = n_bad
    report["negative_topology_counts"] = negative_count_checks

    return report


def compute_normalized_values(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute the phase-1 derived variables.

    Notes
    -----
    d_max_pct uses D as the characteristic length:
        d_max_pct = 100 * d_mag_max / D

    patch_volume assumes size_h is a length scale:
        patch_volume = patch_elems * size_h^3

    H1_norm and H2_norm are density-like topological descriptors:
        Hn_norm = Hn / patch_volume
    """
    out = df.copy()

    out["d_max_pct"] = 100.0 * out["d_mag_max"] / out["D"]

    out["patch_volume"] = out["patch_elems"] * (out["size_h"] ** 3)

    out["H1_norm"] = out["H1"] / out["patch_volume"]
    out["H2_norm"] = out["H2"] / out["patch_volume"]

    if "H0" in out.columns:
        out["H0_norm"] = out["H0"] / out["patch_volume"]

    return out


def compare_existing_derived_values(
    original_df: pd.DataFrame,
    normalized_df: pd.DataFrame,
    tolerance: float = 1e-8
) -> Dict:
    """
    Compare existing derived columns, if present in the original CSV,
    against recomputed values.

    This is only a verification step. The recomputed values are treated
    as authoritative for the downstream pipeline.
    """
    existing = find_existing_derived_columns(original_df)
    comparison = {}

    for canonical, raw_col in existing.items():
        if canonical not in normalized_df.columns:
            continue

        old = pd.to_numeric(original_df[raw_col], errors="coerce")
        new = normalized_df[canonical]

        diff = old - new
        abs_diff = diff.abs()

        comparison[canonical] = {
            "input_column": raw_col,
            "n_compared": int(diff.notna().sum()),
            "max_abs_diff": float(abs_diff.max(skipna=True)),
            "mean_abs_diff": float(abs_diff.mean(skipna=True)),
            "passes_tolerance": bool(abs_diff.max(skipna=True) <= tolerance),
            "tolerance": tolerance,
        }

    return comparison


def summarize_dataset(df: pd.DataFrame) -> Dict:
    """Create a compact dataset summary for the validation report."""
    summary = {
        "n_rows": int(len(df)),
        "n_columns": int(df.shape[1]),
        "capsides": sorted(df["capside"].dropna().astype(str).unique().tolist())
        if "capside" in df.columns else [],
        "folds": sorted(df["fold"].dropna().astype(str).unique().tolist())
        if "fold" in df.columns else [],
        "h_values": sorted(df["h"].dropna().unique().tolist())
        if "h" in df.columns else [],
        "patch_R_values": sorted(df["patch_R"].dropna().unique().tolist())
        if "patch_R" in df.columns else [],
    }

    if all(c in df.columns for c in ["capside", "fold"]):
        summary["n_capsides"] = int(df["capside"].nunique(dropna=True))
        summary["n_folds"] = int(df["fold"].nunique(dropna=True))

    return summary


def build_validation_report(
    original_df: pd.DataFrame,
    canonical_df: pd.DataFrame,
    normalized_df: pd.DataFrame,
    rename_report: Dict,
    tolerance: float
) -> Dict:
    """Assemble full phase-1 validation report."""
    missing_required = validate_required_columns(canonical_df)

    report = {
        "phase": "01_read_validate_normalize",
        "status": "PASS" if not missing_required else "FAIL",
        "missing_required_columns": missing_required,
        "rename_report": rename_report,
        "dataset_summary": summarize_dataset(normalized_df),
        "observation_key_report": validate_observation_key(normalized_df),
        "value_validation": validate_values(canonical_df),
        "derived_value_verification": compare_existing_derived_values(
            original_df=original_df,
            normalized_df=normalized_df,
            tolerance=tolerance,
        ),
        "derived_columns_created": [
            c for c in ["d_max_pct", "patch_volume", "H1_norm", "H2_norm", "H0_norm"]
            if c in normalized_df.columns
        ],
    }

    severe_flags = []

    if missing_required:
        severe_flags.append("Missing required columns.")

    value_validation = report["value_validation"]

    if value_validation["non_positive_values"].get("D", 0) > 0:
        severe_flags.append("D contains non-positive values; d_max_pct is invalid.")

    if value_validation["non_positive_values"].get("patch_elems", 0) > 0:
        severe_flags.append("patch_elems contains non-positive values; topology normalization is invalid.")

    if value_validation["non_positive_values"].get("size_h", 0) > 0:
        severe_flags.append("size_h contains non-positive values; topology normalization is invalid.")

    if report["observation_key_report"].get("n_duplicate_rows", 0) > 0:
        severe_flags.append("Duplicate observational keys detected.")

    if severe_flags:
        report["status"] = "FAIL"

    report["severe_flags"] = severe_flags

    return report


def write_outputs(
    normalized_df: pd.DataFrame,
    report: Dict,
    output_csv: Path,
    report_json: Path
) -> None:
    """Write normalized CSV and validation report."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    report_json.parent.mkdir(parents=True, exist_ok=True)

    normalized_df.to_csv(output_csv, index=False)

    with report_json.open("w", encoding="utf-8") as f:
        json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)




def _json_scalar(value: Any) -> Any:
    """Convert numpy/pandas scalar values to JSON-compatible Python values."""
    if pd.isna(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def _sorted_json_values(series: pd.Series) -> List[Any]:
    """Return deterministic, JSON-compatible unique values for report fields."""
    values = [_json_scalar(v) for v in series.dropna().unique().tolist()]
    return sorted(values, key=lambda x: str(x))


def _nan_to_none(obj: Any) -> Any:
    """Recursively replace NaN/Inf values with None for strict JSON output."""
    if isinstance(obj, dict):
        return {k: _nan_to_none(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_nan_to_none(v) for v in obj]
    if isinstance(obj, tuple):
        return [_nan_to_none(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def _semicolon_join(values: pd.Series | List[Any]) -> str:
    """Join categorical values deterministically without changing fold labels."""
    if isinstance(values, pd.Series):
        items = values.dropna().astype(str).unique().tolist()
    else:
        items = [str(v) for v in values if not pd.isna(v)]
    return ";".join(sorted(dict.fromkeys(items)))


def _append_flag(existing: Any, flag: str) -> str:
    """Append a semicolon-separated warning flag exactly once."""
    flags = [] if pd.isna(existing) or existing == "" else str(existing).split(";")
    if flag not in flags:
        flags.append(flag)
    return ";".join(flags)


def run_phase01(args: argparse.Namespace) -> Tuple[pd.DataFrame, Dict]:
    """Run Phase 1 in memory and write its normalized CSV/report outputs."""
    original_df = pd.read_csv(args.input)
    canonical_df, rename_report = rename_to_canonical(original_df)
    canonical_df = coerce_numeric_columns(canonical_df)

    missing_required = validate_required_columns(canonical_df)
    if missing_required:
        report = {
            "phase": "01_read_validate_normalize",
            "status": "FAIL",
            "missing_required_columns": missing_required,
            "rename_report": rename_report,
            "dataset_summary": summarize_dataset(canonical_df),
            "severe_flags": ["Missing required columns."],
        }
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open("w", encoding="utf-8") as f:
            json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)
        return canonical_df, report

    normalized_df = compute_normalized_values(canonical_df)
    report = build_validation_report(
        original_df=original_df,
        canonical_df=canonical_df,
        normalized_df=normalized_df,
        rename_report=rename_report,
        tolerance=args.derived_tolerance,
    )
    write_outputs(
        normalized_df=normalized_df,
        report=report,
        output_csv=args.output,
        report_json=args.report,
    )
    return normalized_df, report


def validate_phase02_input(df: pd.DataFrame, expected_folds: List[str]) -> Dict:
    """Validate Phase 2 inputs without silently aggregating problematic rows."""
    missing_columns = [col for col in PHASE02_REQUIRED_COLUMNS if col not in df.columns]
    duplicate_examples: List[Dict[str, Any]] = []
    n_duplicate_rows = 0
    if not missing_columns:
        duplicated_mask = df.duplicated(subset=OBSERVATION_KEY, keep=False)
        n_duplicate_rows = int(duplicated_mask.sum())
        if n_duplicate_rows:
            duplicate_examples = (
                df.loc[duplicated_mask, OBSERVATION_KEY]
                .drop_duplicates()
                .head(20)
                .to_dict(orient="records")
            )
    finite_mask = pd.Series(False, index=df.index)
    nonfinite_d_max_pct = 0
    if "d_max_pct" in df.columns:
        d_values = pd.to_numeric(df["d_max_pct"], errors="coerce")
        finite_mask = np.isfinite(d_values)
        nonfinite_d_max_pct = int((~finite_mask).sum())

    valid_rows_for_counts = df.loc[finite_mask].copy() if not missing_columns else pd.DataFrame()
    stratum_counts: Dict[Tuple[Any, Any, Any], int] = {}
    single_fold_strata: List[Dict[str, Any]] = []
    missing_expected_fold_strata: List[Dict[str, Any]] = []
    if not valid_rows_for_counts.empty:
        for key, group in valid_rows_for_counts.groupby(PHASE02_GROUP_COLS, sort=True, dropna=False):
            n_folds = int(group["fold"].astype(str).nunique(dropna=True))
            stratum_counts[key] = n_folds
            if n_folds < 2:
                single_fold_strata.append({
                    "capside": _json_scalar(key[0]),
                    "h": _json_scalar(key[1]),
                    "patch_R": _json_scalar(key[2]),
                    "n_valid_folds": n_folds,
                })
            observed = set(group["fold"].dropna().astype(str).unique().tolist())
            missing_expected = [f for f in expected_folds if f not in observed]
            if missing_expected:
                missing_expected_fold_strata.append({
                    "capside": _json_scalar(key[0]),
                    "h": _json_scalar(key[1]),
                    "patch_R": _json_scalar(key[2]),
                    "missing_expected_folds": missing_expected,
                })

    return {
        "missing_columns": missing_columns,
        "duplicate_key_status": {
            "observation_key": OBSERVATION_KEY,
            "is_unique": n_duplicate_rows == 0,
            "n_duplicate_rows": n_duplicate_rows,
            "duplicate_examples": duplicate_examples,
        },
        "nonfinite_d_max_pct": nonfinite_d_max_pct,
        "n_valid_d_max_pct_rows": int(finite_mask.sum()) if "d_max_pct" in df.columns else 0,
        "single_fold_strata": single_fold_strata,
        "missing_expected_fold_strata": missing_expected_fold_strata,
        "n_strata_total": int(len(stratum_counts)),
        "n_valid_strata": int(sum(1 for n in stratum_counts.values() if n >= 2)),
        "n_invalid_single_fold_strata": int(sum(1 for n in stratum_counts.values() if n < 2)),
    }


def compute_fold_centered_table(df: pd.DataFrame, validation: Dict) -> pd.DataFrame:
    """Add within-stratum means, centered deviations, ranks, labels, and flags."""
    out = df.copy()
    out["fold"] = out["fold"].astype(str)
    out["phase02_warning_flags"] = ""
    out["d_max_pct"] = pd.to_numeric(out["d_max_pct"], errors="coerce")
    out["d_mag_max"] = pd.to_numeric(out["d_mag_max"], errors="coerce")
    finite_mask = np.isfinite(out["d_max_pct"])
    out.loc[~finite_mask, "phase02_warning_flags"] = out.loc[~finite_mask, "phase02_warning_flags"].apply(
        lambda x: _append_flag(x, "nonfinite_d_max_pct")
    )

    valid = out.loc[finite_mask].copy()
    fold_counts = valid.groupby(PHASE02_GROUP_COLS, sort=False, dropna=False)["fold"].transform("nunique")
    valid_calc = valid.loc[fold_counts >= 2].copy()
    invalid_single_index = valid.index[fold_counts < 2]
    out.loc[invalid_single_index, "phase02_warning_flags"] = out.loc[invalid_single_index, "phase02_warning_flags"].apply(
        lambda x: _append_flag(x, "invalid_stratum_single_fold")
    )

    for col in [
        "mean_d_max_pct_s", "delta_d_max_pct", "rel_delta_d_max_pct",
        "mean_d_mag_max_s", "delta_d_mag_max", "rank_d_max_pct_within_capsid",
    ]:
        out[col] = np.nan
    out["is_most_deformable_fold"] = False
    out["is_least_deformable_fold"] = False

    if not valid_calc.empty:
        means_pct = valid_calc.groupby(PHASE02_GROUP_COLS, sort=False, dropna=False)["d_max_pct"].transform("mean")
        means_abs = valid_calc.groupby(PHASE02_GROUP_COLS, sort=False, dropna=False)["d_mag_max"].transform("mean")
        ranks = valid_calc.groupby(PHASE02_GROUP_COLS, sort=False, dropna=False)["d_max_pct"].rank(method="dense", ascending=False)
        max_values = valid_calc.groupby(PHASE02_GROUP_COLS, sort=False, dropna=False)["d_max_pct"].transform("max")
        min_values = valid_calc.groupby(PHASE02_GROUP_COLS, sort=False, dropna=False)["d_max_pct"].transform("min")

        out.loc[valid_calc.index, "mean_d_max_pct_s"] = means_pct
        out.loc[valid_calc.index, "delta_d_max_pct"] = valid_calc["d_max_pct"] - means_pct
        rel = (valid_calc["d_max_pct"] - means_pct) / means_pct.where(means_pct > 0)
        out.loc[valid_calc.index, "rel_delta_d_max_pct"] = rel
        out.loc[valid_calc.index, "mean_d_mag_max_s"] = means_abs
        out.loc[valid_calc.index, "delta_d_mag_max"] = valid_calc["d_mag_max"] - means_abs
        out.loc[valid_calc.index, "rank_d_max_pct_within_capsid"] = ranks.astype("Int64")
        out.loc[valid_calc.index, "is_most_deformable_fold"] = valid_calc["d_max_pct"].eq(max_values)
        out.loc[valid_calc.index, "is_least_deformable_fold"] = valid_calc["d_max_pct"].eq(min_values)

        invalid_mean_index = valid_calc.index[~(means_pct > 0)]
        out.loc[invalid_mean_index, "phase02_warning_flags"] = out.loc[invalid_mean_index, "phase02_warning_flags"].apply(
            lambda x: _append_flag(x, "invalid_mean_d_max_pct")
        )

    required = [
        "capside", "h", "patch_R", "fold", "d_mag_max", "D", "d_max_pct",
        "mean_d_max_pct_s", "delta_d_max_pct", "rel_delta_d_max_pct",
        "mean_d_mag_max_s", "delta_d_mag_max", "rank_d_max_pct_within_capsid",
        "is_most_deformable_fold", "is_least_deformable_fold",
    ]
    preserved = [c for c in PHASE02_PRESERVED_COLUMNS if c in out.columns]
    extras = ["phase02_warning_flags"]
    ordered = required + preserved + extras + [c for c in out.columns if c not in required + preserved + extras]
    return out[ordered].sort_values(["capside", "h", "patch_R", "fold"], kind="mergesort").reset_index(drop=True)


def compute_capsid_anisotropy(centered_df: pd.DataFrame, expected_folds: List[str]) -> pd.DataFrame:
    """Produce one mechanical anisotropy row per capsid/h/patch_R stratum."""
    rows = []
    calc = centered_df[np.isfinite(centered_df["d_max_pct"])].copy()
    for key, group in calc.groupby(PHASE02_GROUP_COLS, sort=True, dropna=False):
        valid_group = group[group["rank_d_max_pct_within_capsid"].notna()].copy()
        flags: List[str] = []
        observed = group["fold"].dropna().astype(str).unique().tolist()
        missing_expected = [f for f in expected_folds if f not in set(observed)]
        if missing_expected:
            flags.append("missing_expected_folds")
        if valid_group.empty or valid_group["fold"].nunique(dropna=True) < 2:
            flags.append("invalid_stratum_single_fold")
            rows.append({
                "capside": key[0], "h": key[1], "patch_R": key[2],
                "n_folds": int(group["fold"].nunique(dropna=True)),
                "folds_present": _semicolon_join(group["fold"]),
                "mean_d_max_pct": np.nan, "max_d_max_pct": np.nan, "min_d_max_pct": np.nan,
                "range_d_max_pct": np.nan, "A_mech_pct": np.nan,
                "mean_d_mag_max": np.nan, "max_d_mag_max": np.nan, "min_d_mag_max": np.nan,
                "range_d_mag_max": np.nan, "A_mech_abs": np.nan,
                "most_deformable_fold": "", "least_deformable_fold": "",
                "anisotropy_valid": False, "warning_flags": ";".join(flags),
            })
            continue
        mean_pct = float(valid_group["d_max_pct"].mean())
        max_pct = float(valid_group["d_max_pct"].max())
        min_pct = float(valid_group["d_max_pct"].min())
        mean_abs = float(valid_group["d_mag_max"].mean())
        max_abs = float(valid_group["d_mag_max"].max())
        min_abs = float(valid_group["d_mag_max"].min())
        most = valid_group.loc[valid_group["d_max_pct"].eq(max_pct), "fold"].astype(str).tolist()
        least = valid_group.loc[valid_group["d_max_pct"].eq(min_pct), "fold"].astype(str).tolist()
        if len(most) > 1:
            flags.append("tied_most_deformable_fold")
        if len(least) > 1:
            flags.append("tied_least_deformable_fold")
        anisotropy_valid = bool(mean_pct > 0)
        if not anisotropy_valid:
            flags.append("invalid_mean_d_max_pct")
        rows.append({
            "capside": key[0], "h": key[1], "patch_R": key[2],
            "n_folds": int(valid_group["fold"].nunique(dropna=True)),
            "folds_present": _semicolon_join(valid_group["fold"]),
            "mean_d_max_pct": mean_pct,
            "max_d_max_pct": max_pct,
            "min_d_max_pct": min_pct,
            "range_d_max_pct": max_pct - min_pct,
            "A_mech_pct": (max_pct - min_pct) / mean_pct if mean_pct > 0 else np.nan,
            "mean_d_mag_max": mean_abs,
            "max_d_mag_max": max_abs,
            "min_d_mag_max": min_abs,
            "range_d_mag_max": max_abs - min_abs,
            "A_mech_abs": (max_abs - min_abs) / mean_abs if mean_abs > 0 else np.nan,
            "most_deformable_fold": ";".join(sorted(most)),
            "least_deformable_fold": ";".join(sorted(least)),
            "anisotropy_valid": anisotropy_valid,
            "warning_flags": ";".join(dict.fromkeys(flags)),
        })
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=[
            "capside", "h", "patch_R", "n_folds", "folds_present", "mean_d_max_pct",
            "max_d_max_pct", "min_d_max_pct", "range_d_max_pct", "A_mech_pct",
            "mean_d_mag_max", "max_d_mag_max", "min_d_mag_max", "range_d_mag_max",
            "A_mech_abs", "most_deformable_fold", "least_deformable_fold",
            "anisotropy_valid", "warning_flags",
        ])
    return result.sort_values(["h", "patch_R", "A_mech_pct", "capside"], ascending=[True, True, False, True], na_position="last", kind="mergesort").reset_index(drop=True)


def compute_fold_tendency_summary(centered_df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Summarize recurrent centered deformation tendency by h, patch_R, and fold."""
    calc = centered_df[centered_df["delta_d_max_pct"].notna()].copy()
    rows = []
    for key, group in calc.groupby(PHASE02_FOLD_SUMMARY_COLS, sort=True, dropna=False):
        deltas = group["delta_d_max_pct"]
        rel = group["rel_delta_d_max_pct"]
        n = int(len(group))
        n_capsides = int(group["capside"].nunique(dropna=True))
        frac_pos = float((deltas > 0).sum() / n) if n else np.nan
        frac_neg = float((deltas < 0).sum() / n) if n else np.nan
        frac_zero = float((deltas == 0).sum() / n) if n else np.nan
        if frac_pos >= threshold:
            label = "recurrently_above_capsid_mean"
        elif frac_neg >= threshold:
            label = "recurrently_below_capsid_mean"
        else:
            label = "inconsistent"
        sd = float(deltas.std(ddof=1)) if n > 1 else np.nan
        rows.append({
            "h": key[0], "patch_R": key[1], "fold": key[2],
            "n_observations": n,
            "n_capsides": n_capsides,
            "capsides_present": _semicolon_join(group["capside"]),
            "mean_delta_d_max_pct": float(deltas.mean()),
            "median_delta_d_max_pct": float(deltas.median()),
            "sd_delta_d_max_pct": sd,
            "sem_delta_d_max_pct": sd / np.sqrt(n) if n > 1 else np.nan,
            "min_delta_d_max_pct": float(deltas.min()),
            "max_delta_d_max_pct": float(deltas.max()),
            "fraction_positive": frac_pos,
            "fraction_negative": frac_neg,
            "fraction_zero": frac_zero,
            "mean_rel_delta_d_max_pct": float(rel.mean(skipna=True)),
            "median_rel_delta_d_max_pct": float(rel.median(skipna=True)),
            "tendency_label": label,
            "low_capsid_count": bool(n_capsides < 3),
        })
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=[
            "h", "patch_R", "fold", "n_observations", "n_capsides", "capsides_present",
            "mean_delta_d_max_pct", "median_delta_d_max_pct", "sd_delta_d_max_pct",
            "sem_delta_d_max_pct", "min_delta_d_max_pct", "max_delta_d_max_pct",
            "fraction_positive", "fraction_negative", "fraction_zero",
            "mean_rel_delta_d_max_pct", "median_rel_delta_d_max_pct",
            "tendency_label", "low_capsid_count",
        ])
    return result.sort_values(["h", "patch_R", "mean_delta_d_max_pct", "fold"], ascending=[True, True, False, True], na_position="last", kind="mergesort").reset_index(drop=True)


def build_phase02_report(
    df: pd.DataFrame,
    centered_df: pd.DataFrame,
    anisotropy_df: pd.DataFrame,
    tendency_df: pd.DataFrame,
    validation: Dict,
    args: argparse.Namespace,
    output_files: Dict[str, str],
    warnings: List[str],
    severe_flags: List[str],
) -> Dict:
    """Create the Phase 2 JSON audit report."""
    status = "PASS"
    if severe_flags:
        status = "FAIL"
    elif warnings:
        status = "WARN"
    valid_strata = int(anisotropy_df["anisotropy_valid"].sum()) if "anisotropy_valid" in anisotropy_df.columns else 0
    if valid_strata == 0 and len(df) > 0:
        status = "FAIL"
        if "all_strata_invalid" not in severe_flags:
            severe_flags.append("all_strata_invalid")

    return {
        "phase": "02_fold_level_mechanical_comparison",
        "status": status,
        "input_file": str(args.input),
        "output_directory": str(args.outdir),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "required_columns_checked": PHASE02_REQUIRED_COLUMNS,
        "missing_columns": validation["missing_columns"],
        "n_rows_input": int(len(df)),
        "n_rows_output": int(len(centered_df)),
        "n_capsides": int(df["capside"].nunique(dropna=True)) if "capside" in df.columns else 0,
        "n_folds": int(df["fold"].astype(str).nunique(dropna=True)) if "fold" in df.columns else 0,
        "h_values": _sorted_json_values(df["h"]) if "h" in df.columns else [],
        "patch_R_values": _sorted_json_values(df["patch_R"]) if "patch_R" in df.columns else [],
        "n_strata": int(validation["n_strata_total"]),
        "n_valid_strata": valid_strata,
        "n_invalid_strata": int(len(anisotropy_df) - valid_strata),
        "observation_key": OBSERVATION_KEY,
        "duplicate_key_status": validation["duplicate_key_status"],
        "ranking_method": PHASE02_RANKING_METHOD,
        "tendency_threshold": float(args.tendency_threshold),
        "expected_folds": args.expected_folds,
        "nonfinite_d_max_pct_rows": validation["nonfinite_d_max_pct"],
        "low_capsid_count_summaries": int(tendency_df["low_capsid_count"].sum()) if "low_capsid_count" in tendency_df.columns else 0,
        "output_files": output_files,
        "warnings": warnings,
        "severe_flags": severe_flags,
        "validation_details": validation,
    }


def write_phase02_summary(
    path: Path,
    report: Dict,
    anisotropy_df: pd.DataFrame,
    tendency_df: pd.DataFrame,
) -> None:
    """Write a compact human-readable Phase 2 summary."""
    lines = [
        f"Phase 2 fold-level mechanical comparison: {report['status']}",
        f"Number of capsids: {report['n_capsides']}",
        f"Number of folds: {report['n_folds']}",
        f"h values: {', '.join(map(str, report['h_values']))}",
        f"patch_R values: {', '.join(map(str, report['patch_R_values']))}",
        "",
        "Top fold tendencies by mean_delta_d_max_pct:",
    ]
    if tendency_df.empty:
        lines.append("  none")
    else:
        for _, row in tendency_df.groupby(["h", "patch_R"], sort=True).head(5).iterrows():
            lines.append(
                "  h={h}, patch_R={r}, fold={fold}: mean_delta={delta:.6g}, "
                "fraction_positive={fp:.3f}, label={label}".format(
                    h=row["h"], r=row["patch_R"], fold=row["fold"],
                    delta=row["mean_delta_d_max_pct"], fp=row["fraction_positive"],
                    label=row["tendency_label"],
                )
            )
    lines.extend(["", "Capsids with highest A_mech_pct:"])
    valid_aniso = anisotropy_df[anisotropy_df["anisotropy_valid"]] if not anisotropy_df.empty else anisotropy_df
    if valid_aniso.empty:
        lines.append("  none")
    else:
        for _, row in valid_aniso.groupby(["h", "patch_R"], sort=True).head(5).iterrows():
            lines.append(
                "  h={h}, patch_R={r}, capside={c}: A_mech_pct={a:.6g}, most={most}, least={least}".format(
                    h=row["h"], r=row["patch_R"], c=row["capside"], a=row["A_mech_pct"],
                    most=row["most_deformable_fold"], least=row["least_deformable_fold"],
                )
            )
    lines.extend(["", "Warnings:"])
    if report["warnings"]:
        lines.extend(f"  - {w}" for w in report["warnings"])
    else:
        lines.append("  none")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_phase02_outputs(
    centered_df: pd.DataFrame,
    anisotropy_df: pd.DataFrame,
    tendency_df: pd.DataFrame,
    report: Dict,
    outdir: Path,
) -> None:
    """Write all Phase 2 CSV, JSON, and text summary artifacts."""
    outdir.mkdir(parents=True, exist_ok=True)
    centered_df.to_csv(outdir / "phase02_fold_level_centered.csv", index=False)
    anisotropy_df.to_csv(outdir / "phase02_capsid_anisotropy.csv", index=False)
    tendency_df.to_csv(outdir / "phase02_fold_tendency_summary.csv", index=False)
    with (outdir / "phase02_report.json").open("w", encoding="utf-8") as f:
        json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)
    write_phase02_summary(outdir / "phase02_summary.txt", report, anisotropy_df, tendency_df)


def run_phase02(normalized_df: pd.DataFrame, phase01_report: Dict, args: argparse.Namespace) -> Dict:
    """Run Phase 2 fold-level mechanical comparison from the Phase 1 DataFrame."""
    args.outdir.mkdir(parents=True, exist_ok=True)
    output_files = {
        "fold_level_centered": str(args.outdir / "phase02_fold_level_centered.csv"),
        "capsid_anisotropy": str(args.outdir / "phase02_capsid_anisotropy.csv"),
        "fold_tendency_summary": str(args.outdir / "phase02_fold_tendency_summary.csv"),
        "report": str(args.outdir / "phase02_report.json"),
        "summary": str(args.outdir / "phase02_summary.txt"),
    }
    if phase01_report.get("status") != "PASS":
        report = {
            "phase": "02_fold_level_mechanical_comparison",
            "status": "FAIL",
            "input_file": str(args.input),
            "output_directory": str(args.outdir),
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "required_columns_checked": PHASE02_REQUIRED_COLUMNS,
            "missing_columns": [],
            "n_rows_input": int(len(normalized_df)),
            "n_rows_output": 0,
            "n_capsides": int(normalized_df["capside"].nunique(dropna=True)) if "capside" in normalized_df.columns else 0,
            "n_folds": int(normalized_df["fold"].astype(str).nunique(dropna=True)) if "fold" in normalized_df.columns else 0,
            "h_values": _sorted_json_values(normalized_df["h"]) if "h" in normalized_df.columns else [],
            "patch_R_values": _sorted_json_values(normalized_df["patch_R"]) if "patch_R" in normalized_df.columns else [],
            "n_strata": 0,
            "observation_key": OBSERVATION_KEY,
            "duplicate_key_status": phase01_report.get("observation_key_report", {}),
            "ranking_method": PHASE02_RANKING_METHOD,
            "tendency_threshold": float(args.tendency_threshold),
            "output_files": output_files,
            "warnings": [],
            "severe_flags": [
                "phase01_validation_failed",
                "fold-level analysis is invalid until the normalized input table passes validation",
            ],
        }
        with (args.outdir / "phase02_report.json").open("w", encoding="utf-8") as f:
            json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)
        (args.outdir / "phase02_summary.txt").write_text(
            "Phase 2 fold-level mechanical comparison: FAIL\n"
            "Phase 1 validation failed; Phase 2 was not run.\n",
            encoding="utf-8",
        )
        return report

    validation = validate_phase02_input(normalized_df, args.expected_folds)
    severe_flags: List[str] = []
    warnings: List[str] = []
    if validation["missing_columns"]:
        severe_flags.append("missing_required_column")
    if not validation["duplicate_key_status"]["is_unique"]:
        severe_flags.append("duplicate_observation_key")
    if validation["nonfinite_d_max_pct"]:
        warnings.append(f"nonfinite_d_max_pct_rows={validation['nonfinite_d_max_pct']}")
    if validation["n_invalid_single_fold_strata"]:
        warnings.append(f"invalid_stratum_single_fold_count={validation['n_invalid_single_fold_strata']}")
    if validation["missing_expected_fold_strata"]:
        warnings.append(f"missing_expected_folds_strata={len(validation['missing_expected_fold_strata'])}")

    if severe_flags:
        empty = pd.DataFrame()
        report = build_phase02_report(
            df=normalized_df,
            centered_df=empty,
            anisotropy_df=pd.DataFrame({"anisotropy_valid": []}),
            tendency_df=empty,
            validation=validation,
            args=args,
            output_files=output_files,
            warnings=warnings,
            severe_flags=severe_flags,
        )
        with (args.outdir / "phase02_report.json").open("w", encoding="utf-8") as f:
            json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)
        (args.outdir / "phase02_summary.txt").write_text(
            f"Phase 2 fold-level mechanical comparison: {report['status']}\n"
            f"Severe flags: {', '.join(severe_flags)}\n",
            encoding="utf-8",
        )
        return report

    centered_df = compute_fold_centered_table(normalized_df, validation)
    anisotropy_df = compute_capsid_anisotropy(centered_df, args.expected_folds)
    tendency_df = compute_fold_tendency_summary(centered_df, args.tendency_threshold)

    if not anisotropy_df.empty and anisotropy_df["warning_flags"].astype(bool).any():
        flagged = sorted(set(";".join(anisotropy_df["warning_flags"].dropna().astype(str)).split(";")) - {""})
        warnings.extend(flagged)
    if not tendency_df.empty and tendency_df["low_capsid_count"].any():
        warnings.append("low_capsid_count")
    if not centered_df.empty:
        sums = centered_df.dropna(subset=["delta_d_max_pct"]).groupby(PHASE02_GROUP_COLS, dropna=False)["delta_d_max_pct"].sum().abs()
        if not sums.empty and bool((sums > 1e-10).any()):
            warnings.append("centered_delta_sum_tolerance_exceeded")

    report = build_phase02_report(
        df=normalized_df,
        centered_df=centered_df,
        anisotropy_df=anisotropy_df,
        tendency_df=tendency_df,
        validation=validation,
        args=args,
        output_files=output_files,
        warnings=sorted(dict.fromkeys(warnings)),
        severe_flags=severe_flags,
    )
    write_phase02_outputs(centered_df, anisotropy_df, tendency_df, report, args.outdir)
    return report



# ---------------------------------------------------------------------
# Phase 3: local centered intra-capsid geometry
# ---------------------------------------------------------------------

PHASE03_CENTER_GROUP_COLS = ["capside", "h", "patch_R"]
PHASE03_CONDITION_COLS = ["h", "patch_R"]
PHASE03_PREDICTORS = ["t_c", "Hout_c", "Hin_c"]
PHASE03_RAW_BY_CENTERED = {"t_c": "t", "Hout_c": "Hout", "Hin_c": "Hin"}
PHASE03_RESPONSE = "delta_d_max_pct"
PHASE03_REQUIRED_COLUMNS = [
    "capside", "fold", "h", "patch_R", "d_max_pct", "delta_d_max_pct", "t", "Hout", "Hin",
]
PHASE03_PRESERVE_COLUMNS = [
    "d_mag_max", "D", "mean_d_max_pct_s", "rel_delta_d_max_pct", "rank_d_max_pct_within_capsid",
    "is_most_deformable_fold", "is_least_deformable_fold", "H1", "H2", "H1_norm",
    "H2_norm", "H0", "H0_norm", "patch_elems", "size_h", "E",
]
PHASE03_SLOPE_UNITS = {
    "t_c": "percent deformation per unit thickness",
    "Hout_c": "percent deformation per unit Hout",
    "Hin_c": "percent deformation per unit Hin",
}


def _flag_string(flags: List[str]) -> str:
    return ";".join(dict.fromkeys([f for f in flags if f]))


def _format_condition_value(value: Any) -> str:
    text = str(value).replace(".", "p").replace("-", "m")
    return "".join(ch if ch.isalnum() or ch in ["_", "p", "m"] else "_" for ch in text)


def validate_phase03_input(df: pd.DataFrame) -> Dict:
    missing = [col for col in PHASE03_REQUIRED_COLUMNS if col not in df.columns]
    duplicate_examples: List[Dict[str, Any]] = []
    n_duplicate_rows = 0
    if all(c in df.columns for c in OBSERVATION_KEY):
        dup_mask = df.duplicated(subset=OBSERVATION_KEY, keep=False)
        n_duplicate_rows = int(dup_mask.sum())
        if n_duplicate_rows:
            duplicate_examples = (
                df.loc[dup_mask, OBSERVATION_KEY].drop_duplicates().head(20).to_dict(orient="records")
            )
    nonfinite = {}
    for col in ["delta_d_max_pct", "t", "Hout", "Hin", "h", "patch_R"]:
        if col in df.columns:
            vals = pd.to_numeric(df[col], errors="coerce") if col not in ["capside", "fold"] else df[col]
            if col not in ["capside", "fold"]:
                n_bad = int((~np.isfinite(vals)).sum())
                if n_bad:
                    nonfinite[col] = n_bad
    return {
        "missing_columns": missing,
        "duplicate_key_status": {
            "observation_key": OBSERVATION_KEY,
            "is_unique": n_duplicate_rows == 0,
            "n_duplicate_rows": n_duplicate_rows,
            "duplicate_examples": duplicate_examples,
        },
        "nonfinite_counts": nonfinite,
    }


def compute_centered_geometry(df: pd.DataFrame, upstream_phase02_warn: bool = False) -> pd.DataFrame:
    out = df.copy()
    out["fold"] = out["fold"].astype(str)
    out["phase03_warning_flags"] = ""
    if upstream_phase02_warn:
        out["phase03_warning_flags"] = out["phase03_warning_flags"].apply(lambda x: _append_flag(x, "upstream_phase02_warn"))
    for col in ["t", "Hout", "Hin", "delta_d_max_pct", "d_max_pct", "h", "patch_R"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    for raw in ["t", "Hout", "Hin"]:
        bad = ~np.isfinite(out[raw])
        if bad.any():
            out.loc[bad, "phase03_warning_flags"] = out.loc[bad, "phase03_warning_flags"].apply(
                lambda x, r=raw: _append_flag(x, f"nonfinite_{r}")
            )
    bad_delta = ~np.isfinite(out["delta_d_max_pct"])
    if bad_delta.any():
        out.loc[bad_delta, "phase03_warning_flags"] = out.loc[bad_delta, "phase03_warning_flags"].apply(
            lambda x: _append_flag(x, "nonfinite_delta_d_max_pct")
        )
    for raw in ["t", "Hout", "Hin"]:
        mean_col = f"mean_{raw}_s"
        centered_col = f"{raw}_c"
        out[mean_col] = out.groupby(PHASE03_CENTER_GROUP_COLS, sort=False, dropna=False)[raw].transform("mean")
        out[centered_col] = out[raw] - out[mean_col]
    ordered = [
        "capside", "h", "patch_R", "fold", "d_max_pct", "delta_d_max_pct",
        "rel_delta_d_max_pct", "t", "mean_t_s", "t_c", "Hout", "mean_Hout_s",
        "Hout_c", "Hin", "mean_Hin_s", "Hin_c", "rank_d_max_pct_within_capsid",
        "is_most_deformable_fold", "is_least_deformable_fold", "phase03_warning_flags",
    ]
    preserve = [c for c in PHASE03_PRESERVE_COLUMNS if c in out.columns and c not in ordered]
    ordered = [c for c in ordered if c in out.columns] + preserve + [c for c in out.columns if c not in ordered + preserve]
    return out[ordered].sort_values(OBSERVATION_KEY, kind="mergesort").reset_index(drop=True)


def check_phase03_centering_integrity(centered_df: pd.DataFrame, tol: float = 1e-10) -> Dict:
    checks: Dict[str, Any] = {"tolerance_abs": tol, "by_variable": {}, "failed_delta_strata": [], "geometry_warning_strata": []}
    variables = ["delta_d_max_pct", "t_c", "Hout_c", "Hin_c"]
    for var in variables:
        var_checks = []
        for key, group in centered_df.groupby(PHASE03_CENTER_GROUP_COLS, sort=True, dropna=False):
            finite_vals = pd.to_numeric(group[var], errors="coerce")
            finite_vals = finite_vals[np.isfinite(finite_vals)]
            centered_sum = float(finite_vals.sum()) if len(finite_vals) else np.nan
            passed = bool(np.isfinite(centered_sum) and np.isclose(centered_sum, 0.0, atol=tol, rtol=tol))
            entry = {
                "capside": _json_scalar(key[0]), "h": _json_scalar(key[1]), "patch_R": _json_scalar(key[2]),
                "n_finite": int(len(finite_vals)), "sum": centered_sum, "passes": passed,
            }
            var_checks.append(entry)
            if not passed:
                if var == "delta_d_max_pct":
                    checks["failed_delta_strata"].append(entry)
                else:
                    raw = PHASE03_RAW_BY_CENTERED[var]
                    missing_raw = int((~np.isfinite(pd.to_numeric(group[raw], errors="coerce"))).sum())
                    entry["raw_missing_count"] = missing_raw
                    checks["geometry_warning_strata"].append({"variable": var, **entry})
        checks["by_variable"][var] = var_checks
    checks["delta_passes"] = len(checks["failed_delta_strata"]) == 0
    checks["geometry_passes_or_missing_only"] = True
    return checks


def _ols_fit(df: pd.DataFrame, predictor: str, response: str = PHASE03_RESPONSE) -> Optional[Dict[str, float]]:
    x = pd.to_numeric(df[predictor], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[response], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = len(x)
    if n < 2:
        return None
    x_var = float(np.var(x, ddof=0))
    y_var = float(np.var(y, ddof=0))
    if x_var <= 0 or y_var <= 0:
        return None
    slope, intercept = np.polyfit(x, y, 1)
    y_hat = slope * x + intercept
    residuals = y - y_hat
    r_value = float(np.corrcoef(x, y)[0, 1]) if n > 1 else np.nan
    r_squared = r_value ** 2 if np.isfinite(r_value) else np.nan
    residual_std_error = float(np.sqrt(np.sum(residuals ** 2) / (n - 2))) if n > 2 else np.nan
    return {
        "slope": float(slope), "intercept": float(intercept), "r_value": r_value,
        "r_squared": float(r_squared), "residual_std_error": residual_std_error,
    }


def _valid_model_df(condition_df: pd.DataFrame, predictor: str) -> pd.DataFrame:
    x = pd.to_numeric(condition_df[predictor], errors="coerce")
    y = pd.to_numeric(condition_df[PHASE03_RESPONSE], errors="coerce")
    h = pd.to_numeric(condition_df["h"], errors="coerce")
    r = pd.to_numeric(condition_df["patch_R"], errors="coerce")
    return condition_df[np.isfinite(x) & np.isfinite(y) & np.isfinite(h) & np.isfinite(r)].copy().reset_index(drop=True)


def bootstrap_centered_slope(model_df: pd.DataFrame, predictor: str, args: argparse.Namespace, rng: np.random.Generator) -> Tuple[Dict, List[Dict]]:
    capsids = sorted(model_df["capside"].dropna().astype(str).unique().tolist())
    slopes: List[float] = []
    distribution: List[Dict] = []
    if not capsids or args.bootstrap_n <= 0:
        return {"slopes": slopes, "successful": 0, "failed": int(args.bootstrap_n)}, distribution
    groups = {str(k): g.copy() for k, g in model_df.groupby(model_df["capside"].astype(str), sort=False)}
    for i in range(int(args.bootstrap_n)):
        sampled = rng.choice(capsids, size=len(capsids), replace=True)
        sample_df = pd.concat([groups[c] for c in sampled], ignore_index=True)
        fit = _ols_fit(sample_df, predictor)
        if fit is None:
            continue
        slopes.append(fit["slope"])
        distribution.append({"bootstrap_index": i, "slope_boot": fit["slope"]})
    return {"slopes": slopes, "successful": len(slopes), "failed": int(args.bootstrap_n) - len(slopes)}, distribution


def permute_centered_slope(model_df: pd.DataFrame, predictor: str, args: argparse.Namespace, rng: np.random.Generator) -> Tuple[Dict, List[Dict]]:
    slopes: List[float] = []
    distribution: List[Dict] = []
    if args.perm_n <= 0:
        return {"slopes": slopes, "successful": 0, "failed": 0}, distribution
    base = model_df.copy()
    group_indices = [idx.to_numpy() for _, idx in base.groupby(PHASE03_CENTER_GROUP_COLS, sort=False, dropna=False).groups.items()]
    y = pd.to_numeric(base[PHASE03_RESPONSE], errors="coerce").to_numpy(dtype=float)
    for i in range(int(args.perm_n)):
        permuted_y = y.copy()
        for idx in group_indices:
            if len(idx) > 1:
                permuted_y[idx] = rng.permutation(permuted_y[idx])
        permuted = base.copy()
        permuted[PHASE03_RESPONSE] = permuted_y
        fit = _ols_fit(permuted, predictor)
        if fit is None:
            continue
        slopes.append(fit["slope"])
        distribution.append({"permutation_index": i, "slope_perm": fit["slope"]})
    return {"slopes": slopes, "successful": len(slopes), "failed": int(args.perm_n) - len(slopes)}, distribution


def fit_centered_geometry_models(centered_df: pd.DataFrame, args: argparse.Namespace) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(int(args.seed))
    rows: List[Dict[str, Any]] = []
    boot_rows: List[Dict[str, Any]] = []
    perm_rows: List[Dict[str, Any]] = []
    ci = float(args.ci)
    alpha = (1.0 - ci) / 2.0
    for (h_value, r_value), condition_df in centered_df.groupby(PHASE03_CONDITION_COLS, sort=True, dropna=False):
        for predictor in PHASE03_PREDICTORS:
            flags: List[str] = []
            response_values = pd.to_numeric(condition_df[PHASE03_RESPONSE], errors="coerce")
            pred_values = pd.to_numeric(condition_df[predictor], errors="coerce")
            n_excl_resp = int((~np.isfinite(response_values)).sum())
            n_excl_pred = int((~np.isfinite(pred_values)).sum())
            if n_excl_resp:
                flags.append("nonfinite_delta_d_max_pct")
            raw = PHASE03_RAW_BY_CENTERED[predictor]
            if n_excl_pred:
                flags.append(f"nonfinite_{raw}")
            model_df = _valid_model_df(condition_df, predictor)
            n_obs = int(len(model_df))
            n_caps = int(model_df["capside"].nunique(dropna=True)) if n_obs else 0
            n_strata = int(model_df[PHASE03_CENTER_GROUP_COLS].drop_duplicates().shape[0]) if n_obs else 0
            if n_obs < 4:
                flags.append("insufficient_observations")
            if n_caps < 2:
                flags.append("low_capsid_count")
            x_var = float(np.var(pd.to_numeric(model_df[predictor], errors="coerce"), ddof=0)) if n_obs else 0.0
            y_var = float(np.var(pd.to_numeric(model_df[PHASE03_RESPONSE], errors="coerce"), ddof=0)) if n_obs else 0.0
            if x_var <= 0:
                flags.append("zero_predictor_variance")
            if y_var <= 0:
                flags.append("zero_response_variance")
            valid = n_obs >= 4 and n_caps >= 2 and x_var > 0 and y_var > 0
            fit = _ols_fit(model_df, predictor) if valid else None
            base = {
                "h": h_value, "patch_R": r_value, "predictor": predictor, "response": PHASE03_RESPONSE,
                "n_observations": n_obs, "n_capsides": n_caps, "n_strata": n_strata,
                "slope": np.nan, "intercept": np.nan, "r_value": np.nan, "r_squared": np.nan,
                "residual_std_error": np.nan, "slope_units": PHASE03_SLOPE_UNITS[predictor],
                "slope_bootstrap_mean": np.nan, "slope_bootstrap_sd": np.nan, "slope_ci_low": np.nan,
                "slope_ci_high": np.nan, "ci_level": ci, "bootstrap_n_requested": int(args.bootstrap_n),
                "bootstrap_n_successful": 0, "bootstrap_method": "cluster_by_capside_within_h_patch_R",
                "p_perm_two_sided": np.nan, "perm_n_requested": int(args.perm_n), "perm_n_successful": 0,
                "permutation_method": "shuffle_delta_d_max_pct_within_capside_h_patch_R",
                "perm_slope_mean": np.nan, "perm_slope_sd": np.nan, "perm_slope_min": np.nan,
                "perm_slope_max": np.nan, "n_rows_excluded_nonfinite_response": n_excl_resp,
                "n_rows_excluded_nonfinite_predictor": n_excl_pred, "n_rows_used": n_obs,
                "valid_model": False, "warning_flags": "",
            }
            if fit is not None:
                base.update(fit)
                if n_caps < 3:
                    flags.append("low_capsid_count_for_bootstrap")
                boot, bdist = bootstrap_centered_slope(model_df, predictor, args, rng)
                bslopes = np.asarray(boot["slopes"], dtype=float)
                base["bootstrap_n_successful"] = int(boot["successful"])
                if len(bslopes):
                    base["slope_bootstrap_mean"] = float(np.mean(bslopes))
                    base["slope_bootstrap_sd"] = float(np.std(bslopes, ddof=1)) if len(bslopes) > 1 else 0.0
                    base["slope_ci_low"] = float(np.quantile(bslopes, alpha))
                    base["slope_ci_high"] = float(np.quantile(bslopes, 1.0 - alpha))
                    for row in bdist:
                        boot_rows.append({"h": h_value, "patch_R": r_value, "predictor": predictor, **row})
                else:
                    flags.append("bootstrap_failed")
                if int(args.bootstrap_n) > 0 and boot["successful"] < 0.7 * int(args.bootstrap_n):
                    flags.append("unstable_bootstrap")
                perm, pdist = permute_centered_slope(model_df, predictor, args, rng)
                pslopes = np.asarray(perm["slopes"], dtype=float)
                base["perm_n_successful"] = int(perm["successful"])
                if len(pslopes):
                    base["p_perm_two_sided"] = float((1 + np.sum(np.abs(pslopes) >= abs(base["slope"]))) / (1 + len(pslopes)))
                    base["perm_slope_mean"] = float(np.mean(pslopes))
                    base["perm_slope_sd"] = float(np.std(pslopes, ddof=1)) if len(pslopes) > 1 else 0.0
                    base["perm_slope_min"] = float(np.min(pslopes))
                    base["perm_slope_max"] = float(np.max(pslopes))
                    for row in pdist:
                        perm_rows.append({"h": h_value, "patch_R": r_value, "predictor": predictor, **row})
                else:
                    flags.append("permutation_failed")
                if int(args.perm_n) > 0 and perm["successful"] < 0.7 * int(args.perm_n):
                    flags.append("unstable_permutation")
                base["valid_model"] = not any(f in flags for f in [
                    "zero_predictor_variance", "zero_response_variance", "insufficient_observations", "low_capsid_count",
                    "bootstrap_failed", "permutation_failed", "unstable_bootstrap", "unstable_permutation",
                ])
            base["warning_flags"] = _flag_string(flags)
            rows.append(base)
    return pd.DataFrame(rows), pd.DataFrame(boot_rows), pd.DataFrame(perm_rows)


def make_phase03_scatterplots(centered_df: pd.DataFrame, slopes_df: pd.DataFrame, args: argparse.Namespace) -> Tuple[List[str], List[str]]:
    output_files: List[str] = []
    warnings: List[str] = []
    if args.no_plots:
        return output_files, warnings
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        return output_files, [f"plot_failed: matplotlib unavailable ({exc})"]
    one_condition = centered_df[PHASE03_CONDITION_COLS].drop_duplicates().shape[0] == 1
    for _, slope_row in slopes_df.iterrows():
        predictor = slope_row["predictor"]
        h_value = slope_row["h"]
        r_value = slope_row["patch_R"]
        subset = _valid_model_df(centered_df[(centered_df["h"] == h_value) & (centered_df["patch_R"] == r_value)], predictor)
        if subset.empty:
            continue
        fig, ax = plt.subplots(figsize=(7, 5))
        capsides = sorted(subset["capside"].astype(str).unique().tolist())
        cmap = plt.get_cmap("tab10")
        for idx, capside in enumerate(capsides):
            g = subset[subset["capside"].astype(str) == capside]
            ax.scatter(g[predictor], g[PHASE03_RESPONSE], label=capside, color=cmap(idx % 10))
            for _, point in g.iterrows():
                ax.annotate(str(point["fold"]), (point[predictor], point[PHASE03_RESPONSE]), fontsize=7, alpha=0.75)
        if bool(slope_row.get("valid_model", False)) and np.isfinite(slope_row.get("slope", np.nan)):
            x_vals = pd.to_numeric(subset[predictor], errors="coerce")
            x_line = np.linspace(float(x_vals.min()), float(x_vals.max()), 100)
            y_line = float(slope_row["slope"]) * x_line + float(slope_row["intercept"])
            ax.plot(x_line, y_line, color="black", linewidth=1.5)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax.axvline(0, color="gray", linewidth=0.8, linestyle="--")
        ax.set_xlabel(predictor)
        ax.set_ylabel(PHASE03_RESPONSE)
        ci_text = f"CI=[{slope_row['slope_ci_low']:.3g}, {slope_row['slope_ci_high']:.3g}]" if np.isfinite(slope_row.get("slope_ci_low", np.nan)) else "CI=NA"
        p_text = f"p_perm={slope_row['p_perm_two_sided']:.3g}" if np.isfinite(slope_row.get("p_perm_two_sided", np.nan)) else "p_perm=NA"
        ax.set_title(f"{PHASE03_RESPONSE} ~ {predictor}; h={h_value}, R={r_value}; slope={slope_row['slope']:.3g}; {ci_text}; {p_text}")
        ax.legend(fontsize=7, title="capside")
        fig.tight_layout()
        if one_condition:
            filename = f"phase03_scatter_centered_{predictor}.png"
        else:
            filename = f"phase03_scatter_h_{_format_condition_value(h_value)}_R_{_format_condition_value(r_value)}_{predictor}.png"
        path = args.outdir / filename
        fig.savefig(path, dpi=150)
        plt.close(fig)
        output_files.append(str(path))
    return output_files, warnings


def write_phase03_summary(path: Path, report: Dict, slopes_df: pd.DataFrame) -> None:
    lines = [
        f"Phase 3 centered geometry-mechanics summary: {report['status']}",
        "This summary describes centered associations only; it does not make causal claims.",
        f"Capsids: {report['n_capsides']}",
        f"Folds: {report['n_folds']}",
        f"h values: {report['h_values']}",
        f"patch_R values: {report['patch_R_values']}",
        f"Predictors analyzed: {', '.join(report['predictors'])}",
        "",
    ]
    for (h_value, r_value), group in slopes_df.groupby(["h", "patch_R"], sort=True, dropna=False):
        lines.append(f"Condition h={h_value}, patch_R={r_value}:")
        for _, row in group.iterrows():
            lines.append(
                "  - {pred}: slope={slope:.6g}, CI=[{lo:.6g}, {hi:.6g}], p_perm={p}, valid={valid}, warnings={warn}".format(
                    pred=row["predictor"], slope=row["slope"] if np.isfinite(row["slope"]) else np.nan,
                    lo=row["slope_ci_low"] if np.isfinite(row["slope_ci_low"]) else np.nan,
                    hi=row["slope_ci_high"] if np.isfinite(row["slope_ci_high"]) else np.nan,
                    p=(f"{row['p_perm_two_sided']:.6g}" if np.isfinite(row["p_perm_two_sided"]) else "NA"),
                    valid=bool(row["valid_model"]), warn=row["warning_flags"] or "none",
                )
            )
        lines.append("")
    lines.append("Warnings:")
    if report["warnings"]:
        lines.extend(f"  - {w}" for w in report["warnings"])
    else:
        lines.append("  none")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase03(centered_phase02_df: pd.DataFrame, phase01_report: Dict, phase02_report: Dict, args: argparse.Namespace) -> Dict:
    args.outdir.mkdir(parents=True, exist_ok=True)
    output_files = {
        "centered_geometry": str(args.outdir / "phase03_centered_geometry.csv"),
        "slopes": str(args.outdir / "phase03_centered_geometry_slopes.csv"),
        "report": str(args.outdir / "phase03_report.json"),
        "summary": str(args.outdir / "phase03_summary.txt"),
    }
    warnings: List[str] = []
    severe_flags: List[str] = []
    if phase01_report.get("status") != "PASS":
        severe_flags.append("phase01_validation_failed")
    if phase02_report.get("status") == "FAIL":
        severe_flags.append("phase02_failed")
    upstream_warn = phase02_report.get("status") == "WARN"
    if upstream_warn:
        warnings.append("upstream_phase02_warn")
    validation = validate_phase03_input(centered_phase02_df)
    if validation["missing_columns"]:
        severe_flags.append("missing_required_column")
    if not validation["duplicate_key_status"]["is_unique"]:
        severe_flags.append("duplicate_observation_key")
    if severe_flags:
        status = "FAIL"
        report = {
            "phase": "03_centered_local_geometry", "status": status,
            "timestamp": datetime.now(timezone.utc).isoformat(), "input_file": str(args.input),
            "output_directory": str(args.outdir), "upstream_phase01_status": phase01_report.get("status"),
            "upstream_phase02_status": phase02_report.get("status"),
            "required_columns_checked": PHASE03_REQUIRED_COLUMNS, "missing_columns": validation["missing_columns"],
            "n_rows_input": int(len(centered_phase02_df)), "n_rows_output": 0,
            "n_capsides": 0, "n_folds": 0, "h_values": [], "patch_R_values": [], "n_centering_strata": 0,
            "centering_group_cols": PHASE03_CENTER_GROUP_COLS, "model_conditions": [],
            "predictors": PHASE03_PREDICTORS, "response": PHASE03_RESPONSE,
            "bootstrap_n": int(args.bootstrap_n), "bootstrap_method": "cluster_by_capside_within_h_patch_R",
            "ci_level": float(args.ci), "perm_n": int(args.perm_n),
            "permutation_method": "shuffle_delta_d_max_pct_within_capside_h_patch_R", "seed": int(args.seed),
            "centered_integrity_checks": {}, "output_files": output_files, "warnings": warnings,
            "severe_flags": severe_flags,
        }
        with (args.outdir / "phase03_report.json").open("w", encoding="utf-8") as f:
            json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)
        return report
    centered_df = compute_centered_geometry(centered_phase02_df, upstream_phase02_warn=upstream_warn)
    integrity = check_phase03_centering_integrity(centered_df)
    if not integrity["delta_passes"]:
        severe_flags.append("centered_delta_d_max_pct_integrity_failed")
    if integrity["geometry_warning_strata"]:
        warnings.append("geometry_centering_warnings_due_to_missing_values")
    slopes_df, boot_df, perm_df = fit_centered_geometry_models(centered_df, args)
    if slopes_df.empty or not slopes_df["valid_model"].any():
        severe_flags.append("no_valid_models")
    if not slopes_df.empty:
        for flag in sorted(set(";".join(slopes_df["warning_flags"].fillna("").astype(str)).split(";")) - {""}):
            warnings.append(flag)
    plot_files, plot_warnings = make_phase03_scatterplots(centered_df, slopes_df, args)
    warnings.extend(plot_warnings)
    for w in plot_warnings:
        if "plot_failed" in w:
            if "plot_failed" not in warnings:
                warnings.append("plot_failed")
    if severe_flags:
        status = "FAIL"
    elif warnings or (not slopes_df.empty and not bool(slopes_df["valid_model"].all())):
        status = "WARN"
    else:
        status = "PASS"
    output_files["scatterplots"] = plot_files
    if args.save_bootstrap_distributions:
        output_files["bootstrap_slopes"] = str(args.outdir / "phase03_bootstrap_slopes.csv")
    if args.save_permutation_distributions:
        output_files["permutation_slopes"] = str(args.outdir / "phase03_permutation_slopes.csv")
    model_conditions = centered_df[PHASE03_CONDITION_COLS].drop_duplicates().sort_values(PHASE03_CONDITION_COLS, kind="mergesort").to_dict(orient="records")
    report = {
        "phase": "03_centered_local_geometry", "status": status,
        "timestamp": datetime.now(timezone.utc).isoformat(), "input_file": str(args.input),
        "output_directory": str(args.outdir), "upstream_phase01_status": phase01_report.get("status"),
        "upstream_phase02_status": phase02_report.get("status"),
        "required_columns_checked": PHASE03_REQUIRED_COLUMNS, "missing_columns": validation["missing_columns"],
        "n_rows_input": int(len(centered_phase02_df)), "n_rows_output": int(len(centered_df)),
        "n_capsides": int(centered_df["capside"].nunique(dropna=True)),
        "n_folds": int(centered_df["fold"].astype(str).nunique(dropna=True)),
        "h_values": _sorted_json_values(centered_df["h"]), "patch_R_values": _sorted_json_values(centered_df["patch_R"]),
        "n_centering_strata": int(centered_df[PHASE03_CENTER_GROUP_COLS].drop_duplicates().shape[0]),
        "centering_group_cols": PHASE03_CENTER_GROUP_COLS, "model_conditions": model_conditions,
        "predictors": PHASE03_PREDICTORS, "response": PHASE03_RESPONSE,
        "bootstrap_n": int(args.bootstrap_n), "bootstrap_method": "cluster_by_capside_within_h_patch_R",
        "ci_level": float(args.ci), "perm_n": int(args.perm_n),
        "permutation_method": "shuffle_delta_d_max_pct_within_capside_h_patch_R", "seed": int(args.seed),
        "centered_integrity_checks": integrity, "output_files": output_files,
        "warnings": sorted(dict.fromkeys(warnings)), "severe_flags": severe_flags,
    }
    centered_df.to_csv(args.outdir / "phase03_centered_geometry.csv", index=False)
    slopes_df.to_csv(args.outdir / "phase03_centered_geometry_slopes.csv", index=False)
    if args.save_bootstrap_distributions:
        boot_df.to_csv(args.outdir / "phase03_bootstrap_slopes.csv", index=False)
    if args.save_permutation_distributions:
        perm_df.to_csv(args.outdir / "phase03_permutation_slopes.csv", index=False)
    with (args.outdir / "phase03_report.json").open("w", encoding="utf-8") as f:
        json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)
    write_phase03_summary(args.outdir / "phase03_summary.txt", report, slopes_df)
    return report


def _print_phase03_summary(report: Dict) -> None:
    print(f"[{report['status']}] Phase 3 completed.")
    files = report.get("output_files", {})
    print(f"[INFO] Centered geometry table: {files.get('centered_geometry', 'not written')}")
    print(f"[INFO] Slope summary table: {files.get('slopes', 'not written')}")
    print(f"[INFO] Phase 3 report: {files.get('report', 'not written')}")
    slopes_path = Path(files.get("slopes", ""))
    if slopes_path.exists():
        slopes = pd.read_csv(slopes_path)
        valid = slopes[slopes["valid_model"] == True] if "valid_model" in slopes.columns else pd.DataFrame()
        for _, row in valid.iterrows():
            ci_low = row["slope_ci_low"] if np.isfinite(row["slope_ci_low"]) else np.nan
            ci_high = row["slope_ci_high"] if np.isfinite(row["slope_ci_high"]) else np.nan
            p_value = row["p_perm_two_sided"] if np.isfinite(row["p_perm_two_sided"]) else np.nan
            print(
                f"[MODEL] h={row['h']} patch_R={row['patch_R']} predictor={row['predictor']} "
                f"slope={row['slope']:.6g} CI=[{ci_low:.6g}, {ci_high:.6g}] p_perm={p_value:.6g}"
            )
    if report.get("warnings"):
        print("[WARN] Phase 3 completed with warnings. See phase03_report.json.")


# ---------------------------------------------------------------------
# Main CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Capsid mechanical-geometrical-topological analysis. Phase 1 reads, "
            "validates, and normalizes the master CSV; Phase 2 performs fold-level "
            "mechanical comparisons within each (capside, h, patch_R) stratum; "
            "Phase 3 evaluates centered local geometry-mechanics associations."
        )
    )

    parser.add_argument(
        "--phase",
        choices=["1", "2", "3"],
        default="1",
        help="Analysis phase to run. Phase 1 is the default. Phase 2 runs Phase 1 first; Phase 3 runs Phases 1 and 2 first."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        type=Path,
        help="Input master CSV file."
    )

    parser.add_argument(
        "-o", "--output",
        default=Path("phase01_normalized.csv"),
        type=Path,
        help="Output normalized CSV file."
    )

    parser.add_argument(
        "-r", "--report",
        default=Path("phase01_validation_report.json"),
        type=Path,
        help="Output JSON validation report."
    )

    parser.add_argument(
        "--outdir",
        default=Path("results"),
        type=Path,
        help="Output directory for Phase 2 artifacts."
    )

    parser.add_argument(
        "--tendency-threshold",
        default=0.75,
        type=float,
        help="Fraction threshold for recurrent above/below capsid-mean tendency labels."
    )

    parser.add_argument(
        "--expected-folds",
        default=DEFAULT_EXPECTED_FOLDS,
        nargs="*",
        help="Expected fold labels for warning-only reporting of missing folds."
    )

    parser.add_argument(
        "--derived-tolerance",
        default=1e-8,
        type=float,
        help="Tolerance for comparing existing derived columns against recomputed values."
    )

    parser.add_argument(
        "--bootstrap-n",
        default=10000,
        type=int,
        help="Number of Phase 3 cluster-bootstrap replicates."
    )

    parser.add_argument(
        "--perm-n",
        default=10000,
        type=int,
        help="Number of Phase 3 within-stratum permutation replicates."
    )

    parser.add_argument(
        "--ci",
        default=0.95,
        type=float,
        help="Phase 3 bootstrap percentile confidence interval level."
    )

    parser.add_argument(
        "--seed",
        default=12345,
        type=int,
        help="Random seed for Phase 3 bootstrap and permutation reproducibility."
    )

    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip Phase 3 centered scatterplot generation."
    )

    parser.add_argument(
        "--save-bootstrap-distributions",
        action="store_true",
        help="Save full Phase 3 bootstrap slope distributions."
    )

    parser.add_argument(
        "--save-permutation-distributions",
        action="store_true",
        help="Save full Phase 3 permutation slope distributions."
    )

    parser.add_argument(
        "--fail-on-warning",
        action="store_true",
        help="Exit with non-zero status if the requested phase status is not PASS."
    )

    return parser.parse_args()


def _print_phase01_summary(report: Dict, args: argparse.Namespace) -> None:
    print(f"[{report['status']}] Phase 1 completed.")
    print(f"[INFO] Normalized CSV written to: {args.output}")
    print(f"[INFO] Validation report written to: {args.report}")
    if report.get("severe_flags"):
        print("[WARN] Severe validation flags:")
        for flag in report["severe_flags"]:
            print(f"  - {flag}")


def _print_phase02_summary(report: Dict) -> None:
    print(f"[{report['status']}] Phase 2 completed.")
    files = report.get("output_files", {})
    print(f"[INFO] Row-level centered table: {files.get('fold_level_centered', 'not written')}")
    print(f"[INFO] Capsid anisotropy table: {files.get('capsid_anisotropy', 'not written')}")
    print(f"[INFO] Fold tendency summary: {files.get('fold_tendency_summary', 'not written')}")
    print(f"[INFO] Phase 2 report: {files.get('report', 'not written')}")
    print(f"[INFO] Valid strata: {report.get('n_valid_strata', 0)}")
    print(f"[INFO] Invalid strata: {report.get('n_invalid_strata', 0)}")

    outdir = Path(report.get("output_directory", "."))
    tendency_path = outdir / "phase02_fold_tendency_summary.csv"
    anisotropy_path = outdir / "phase02_capsid_anisotropy.csv"
    if tendency_path.exists():
        tendency = pd.read_csv(tendency_path)
        if not tendency.empty:
            print("[INFO] Top recurrently above-mean fold per h and patch_R:")
            above = tendency[tendency["mean_delta_d_max_pct"] > 0]
            for _, row in above.groupby(["h", "patch_R"], sort=True).head(1).iterrows():
                print(
                    f"  - h={row['h']}, patch_R={row['patch_R']}: "
                    f"fold {row['fold']} (mean_delta={row['mean_delta_d_max_pct']:.6g}, "
                    f"fraction_positive={row['fraction_positive']:.3f})"
                )
    if anisotropy_path.exists():
        aniso = pd.read_csv(anisotropy_path)
        if not aniso.empty and "anisotropy_valid" in aniso.columns:
            valid = aniso[aniso["anisotropy_valid"] == True]
            if not valid.empty:
                print("[INFO] Capsid with highest A_mech_pct per h and patch_R:")
                for _, row in valid.groupby(["h", "patch_R"], sort=True).head(1).iterrows():
                    print(
                        f"  - h={row['h']}, patch_R={row['patch_R']}: "
                        f"{row['capside']} (A_mech_pct={row['A_mech_pct']:.6g})"
                    )


def main() -> int:
    args = parse_args()

    if not 0 <= args.tendency_threshold <= 1:
        print("[ERROR] --tendency-threshold must be between 0 and 1.", file=sys.stderr)
        return 2

    if not args.input.exists():
        print(f"[ERROR] Input file does not exist: {args.input}", file=sys.stderr)
        return 2

    if args.bootstrap_n < 0 or args.perm_n < 0:
        print("[ERROR] --bootstrap-n and --perm-n must be non-negative.", file=sys.stderr)
        return 2

    if not 0 < args.ci < 1:
        print("[ERROR] --ci must be between 0 and 1.", file=sys.stderr)
        return 2

    try:
        normalized_df, phase01_report = run_phase01(args)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2

    if args.phase == "1":
        _print_phase01_summary(phase01_report, args)
        if args.fail_on_warning and phase01_report["status"] != "PASS":
            return 1
        return 0

    phase02_report = run_phase02(normalized_df, phase01_report, args)
    _print_phase02_summary(phase02_report)

    if args.phase == "2":
        if phase02_report["status"] == "FAIL":
            return 1
        if args.fail_on_warning and phase02_report["status"] != "PASS":
            return 1
        return 0

    centered_path = args.outdir / "phase02_fold_level_centered.csv"
    if phase02_report["status"] == "FAIL" or not centered_path.exists():
        phase03_report = run_phase03(pd.DataFrame(), phase01_report, phase02_report, args)
        _print_phase03_summary(phase03_report)
        return 1

    centered_phase02_df = pd.read_csv(centered_path)
    phase03_report = run_phase03(centered_phase02_df, phase01_report, phase02_report, args)
    _print_phase03_summary(phase03_report)

    if phase03_report["status"] == "FAIL":
        return 1
    if args.fail_on_warning and phase03_report["status"] != "PASS":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
