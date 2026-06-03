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
# Phase 4: ordinal fold rankings and rank-agreement summaries
# ---------------------------------------------------------------------

PHASE04_RANK_GROUP_COLS = ["capside", "h", "patch_R"]
PHASE04_KEY_COLUMNS = ["capside", "fold", "h", "patch_R"]
PHASE04_REQUIRED_COLUMNS = PHASE04_KEY_COLUMNS + ["d_max_pct"]
PHASE04_GEOMETRY_COLUMNS = ["t", "Hout", "Hin"]
PHASE04_RANKED_VARIABLES = ["d_max_pct", "t", "t_thinnest", "Hout", "Hin"]
PHASE04_CONSISTENCY_THRESHOLD = 0.75
PHASE04_PRESERVE_COLUMNS = [
    "d_mag_max", "D", "delta_d_max_pct", "rel_delta_d_max_pct", "t_c", "Hout_c", "Hin_c",
    "H1", "H2", "H1_norm", "H2_norm", "H0", "H0_norm", "patch_elems", "size_h", "E",
    "rank_d_max_pct_within_capsid", "is_most_deformable_fold", "is_least_deformable_fold",
]
PHASE04_RANK_COLUMNS = {
    "d_max_pct": "rank_d_max_pct",
    "t": "rank_t",
    "t_thinnest": "rank_t_thinnest",
    "Hout": "rank_Hout",
    "Hin": "rank_Hin",
}
PHASE04_TOP_FLAG_COLUMNS = {
    "d_max_pct": "is_top_d_max_pct",
    "t": "is_top_t",
    "t_thinnest": "is_thinnest_t",
    "Hout": "is_top_Hout",
    "Hin": "is_top_Hin",
}
PHASE04_TIE_FLAG_COLUMNS = {
    "d_max_pct": "tied_d_max_pct_rank",
    "t": "tied_t_rank",
    "t_thinnest": "tied_t_thinnest_rank",
    "Hout": "tied_Hout_rank",
    "Hin": "tied_Hin_rank",
}
PHASE04_COMPARISONS = [
    ("d_max_pct_vs_t", "rank_d_max_pct", "rank_t"),
    ("d_max_pct_vs_t_thinnest", "rank_d_max_pct", "rank_t_thinnest"),
    ("d_max_pct_vs_Hout", "rank_d_max_pct", "rank_Hout"),
    ("d_max_pct_vs_Hin", "rank_d_max_pct", "rank_Hin"),
]
PHASE04_COINCIDENCE_COMPARISONS = [
    ("top_d_max_pct_vs_top_t", "is_top_d_max_pct", "is_top_t"),
    ("top_d_max_pct_vs_thinnest_t", "is_top_d_max_pct", "is_thinnest_t"),
    ("top_d_max_pct_vs_top_Hout", "is_top_d_max_pct", "is_top_Hout"),
    ("top_d_max_pct_vs_top_Hin", "is_top_d_max_pct", "is_top_Hin"),
]


def validate_phase04_input(df: pd.DataFrame) -> Dict:
    missing_keys = [col for col in PHASE04_KEY_COLUMNS if col not in df.columns]
    missing_required = [col for col in PHASE04_REQUIRED_COLUMNS if col not in df.columns]
    missing_geometry = [col for col in PHASE04_GEOMETRY_COLUMNS if col not in df.columns]
    duplicate_examples: List[Dict[str, Any]] = []
    n_duplicate_rows = 0
    if all(c in df.columns for c in PHASE04_KEY_COLUMNS):
        dup_mask = df.duplicated(subset=PHASE04_KEY_COLUMNS, keep=False)
        n_duplicate_rows = int(dup_mask.sum())
        if n_duplicate_rows:
            duplicate_examples = df.loc[dup_mask, PHASE04_KEY_COLUMNS].drop_duplicates().head(20).to_dict(orient="records")
    nonfinite_counts: Dict[str, int] = {}
    for col in ["d_max_pct", "t", "Hout", "Hin"]:
        if col in df.columns:
            vals = pd.to_numeric(df[col], errors="coerce")
            bad = int((~np.isfinite(vals)).sum())
            if bad:
                nonfinite_counts[col] = bad
    n_groups = 0
    insufficient_groups: List[Dict[str, Any]] = []
    if all(c in df.columns for c in PHASE04_RANK_GROUP_COLS + ["fold", "d_max_pct"]):
        finite_d = df[np.isfinite(pd.to_numeric(df["d_max_pct"], errors="coerce"))]
        n_groups = int(finite_d[PHASE04_RANK_GROUP_COLS].drop_duplicates().shape[0])
        for key, group in finite_d.groupby(PHASE04_RANK_GROUP_COLS, sort=True, dropna=False):
            n_folds = int(group["fold"].astype(str).nunique(dropna=True))
            if n_folds < 2:
                insufficient_groups.append({"capside": _json_scalar(key[0]), "h": _json_scalar(key[1]), "patch_R": _json_scalar(key[2]), "n_folds": n_folds})
    return {
        "missing_key_columns": missing_keys,
        "missing_required_columns": missing_required,
        "missing_geometry_columns": missing_geometry,
        "duplicate_key_status": {"observation_key": PHASE04_KEY_COLUMNS, "is_unique": n_duplicate_rows == 0, "n_duplicate_rows": n_duplicate_rows, "duplicate_examples": duplicate_examples},
        "nonfinite_counts": nonfinite_counts,
        "n_ranking_groups": n_groups,
        "insufficient_ranking_groups": insufficient_groups,
    }


def _phase04_rank_series(series: pd.Series, method: str, ascending: bool) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce")
    return vals.rank(method=method, ascending=ascending, na_option="keep")


def compute_fold_rankings(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    out = df.copy()
    out["fold"] = out["fold"].astype(str)
    out["phase04_warning_flags"] = ""
    for col in ["d_max_pct", "t", "Hout", "Hin", "h", "patch_R"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    rank_method = getattr(args, "rank_method", "average")
    if "d_max_pct" in out.columns:
        out["rank_d_max_pct"] = out.groupby(PHASE04_RANK_GROUP_COLS, sort=False, dropna=False)["d_max_pct"].transform(lambda s: _phase04_rank_series(s, rank_method, ascending=False))
    for raw, rank_col in [("t", "rank_t"), ("Hout", "rank_Hout"), ("Hin", "rank_Hin")]:
        if raw in out.columns:
            out[rank_col] = out.groupby(PHASE04_RANK_GROUP_COLS, sort=False, dropna=False)[raw].transform(lambda s: _phase04_rank_series(s, rank_method, ascending=False))
        else:
            out[rank_col] = np.nan
            out["phase04_warning_flags"] = out["phase04_warning_flags"].apply(lambda x, r=raw: _append_flag(x, f"missing_{r}"))
    if "t" in out.columns:
        out["rank_t_thinnest"] = out.groupby(PHASE04_RANK_GROUP_COLS, sort=False, dropna=False)["t"].transform(lambda s: _phase04_rank_series(s, rank_method, ascending=True))
    else:
        out["rank_t_thinnest"] = np.nan
    # Tie flags are based on duplicate finite raw values inside each ranking group.
    raw_for_var = {"d_max_pct": "d_max_pct", "t": "t", "t_thinnest": "t", "Hout": "Hout", "Hin": "Hin"}
    for var, tie_col in PHASE04_TIE_FLAG_COLUMNS.items():
        raw = raw_for_var[var]
        if raw in out.columns:
            out[tie_col] = out.groupby(PHASE04_RANK_GROUP_COLS, sort=False, dropna=False)[raw].transform(lambda s: pd.to_numeric(s, errors="coerce").duplicated(keep=False) & pd.to_numeric(s, errors="coerce").notna()).astype(bool)
        else:
            out[tie_col] = False
    for var, rank_col in PHASE04_RANK_COLUMNS.items():
        flag_col = PHASE04_TOP_FLAG_COLUMNS[var]
        out[flag_col] = False
        if rank_col in out.columns:
            mins = out.groupby(PHASE04_RANK_GROUP_COLS, sort=False, dropna=False)[rank_col].transform("min")
            out[flag_col] = (out[rank_col].notna() & (out[rank_col] == mins))
    for rank_col in PHASE04_RANK_COLUMNS.values():
        if rank_col in out.columns:
            missing = out[rank_col].isna()
            if missing.any():
                out.loc[missing, "phase04_warning_flags"] = out.loc[missing, "phase04_warning_flags"].apply(lambda x: _append_flag(x, "missing_rank_values"))
    for tie_col in PHASE04_TIE_FLAG_COLUMNS.values():
        if tie_col in out.columns and out[tie_col].any():
            out.loc[out[tie_col], "phase04_warning_flags"] = out.loc[out[tie_col], "phase04_warning_flags"].apply(lambda x: _append_flag(x, "ties_present"))
    ordered = PHASE04_KEY_COLUMNS + ["d_max_pct", "t", "Hout", "Hin", "rank_d_max_pct", "rank_t", "rank_t_thinnest", "rank_Hout", "rank_Hin", "tied_d_max_pct_rank", "tied_t_rank", "tied_t_thinnest_rank", "tied_Hout_rank", "tied_Hin_rank", "is_top_d_max_pct", "is_top_t", "is_thinnest_t", "is_top_Hout", "is_top_Hin", "phase04_warning_flags"]
    preserve = [c for c in PHASE04_PRESERVE_COLUMNS if c in out.columns and c not in ordered]
    ordered = [c for c in ordered if c in out.columns] + preserve + [c for c in out.columns if c not in ordered + preserve]
    return out[ordered].sort_values(["h", "patch_R", "capside", "fold"], kind="mergesort").reset_index(drop=True)


def _rank_correlations(x: pd.Series, y: pd.Series) -> Tuple[float, float, float, float]:
    try:
        from scipy import stats  # type: ignore
        s = stats.spearmanr(x, y, nan_policy="omit")
        k = stats.kendalltau(x, y, nan_policy="omit")
        return float(s.statistic), float(s.pvalue), float(k.statistic), float(k.pvalue)
    except Exception:
        # Pandas supplies tie-aware Spearman/Kendall coefficients; p-values are unavailable without scipy.
        return float(x.corr(y, method="spearman")), np.nan, float(x.corr(y, method="kendall")), np.nan


def compute_rank_correlations_by_capsid(rankings_df: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for key, group in rankings_df.groupby(PHASE04_RANK_GROUP_COLS, sort=True, dropna=False):
        capside, h_value, r_value = key
        for comparison, mech_rank, geom_rank in PHASE04_COMPARISONS:
            flags: List[str] = []
            required = [mech_rank, geom_rank]
            valid_rows = group.dropna(subset=required).copy() if all(c in group.columns for c in required) else pd.DataFrame()
            n = int(len(valid_rows))
            spearman_rho = spearman_p = kendall_tau = kendall_p = np.nan
            valid = True
            if n < 3:
                valid = False
                flags.append("insufficient_folds_for_correlation")
            elif valid_rows[mech_rank].nunique(dropna=True) <= 1 or valid_rows[geom_rank].nunique(dropna=True) <= 1:
                valid = False
                flags.append("zero_rank_variance")
            else:
                if n in [3, 4]:
                    flags.append("low_fold_count_rank_correlation")
                tie_cols = []
                if mech_rank == "rank_d_max_pct":
                    tie_cols.append("tied_d_max_pct_rank")
                tie_cols.append({"rank_t": "tied_t_rank", "rank_t_thinnest": "tied_t_thinnest_rank", "rank_Hout": "tied_Hout_rank", "rank_Hin": "tied_Hin_rank"}.get(geom_rank, ""))
                tie_cols = [c for c in tie_cols if c and c in valid_rows.columns]
                if tie_cols and valid_rows[tie_cols].any(axis=None):
                    flags.append("ties_present")
                spearman_rho, spearman_p, kendall_tau, kendall_p = _rank_correlations(valid_rows[mech_rank], valid_rows[geom_rank])
                if not np.isfinite(spearman_rho) or not np.isfinite(kendall_tau):
                    valid = False
                    flags.append("zero_rank_variance")
            rows.append({
                "capside": _json_scalar(capside), "h": _json_scalar(h_value), "patch_R": _json_scalar(r_value),
                "comparison": comparison, "mechanical_rank": mech_rank, "geometry_rank": geom_rank,
                "n_folds_used": n, "folds_used": _semicolon_join(valid_rows["fold"]) if not valid_rows.empty else "",
                "spearman_rho": spearman_rho, "spearman_p": spearman_p, "kendall_tau": kendall_tau, "kendall_p": kendall_p,
                "valid_comparison": bool(valid), "warning_flags": _flag_string(flags),
            })
    return pd.DataFrame(rows).sort_values(["h", "patch_R", "capside", "comparison"], kind="mergesort").reset_index(drop=True)


def _fraction_positive(series: pd.Series) -> float:
    vals = pd.to_numeric(series, errors="coerce").dropna()
    return float((vals > 0).mean()) if len(vals) else np.nan


def _fraction_negative(series: pd.Series) -> float:
    vals = pd.to_numeric(series, errors="coerce").dropna()
    return float((vals < 0).mean()) if len(vals) else np.nan


def _consistency_label(pos: float, neg: float, threshold: float) -> str:
    if np.isfinite(pos) and pos >= threshold:
        return "mostly_positive"
    if np.isfinite(neg) and neg >= threshold:
        return "mostly_negative"
    return "inconsistent"


def compute_rank_correlation_summary(correlations_df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for key, group in correlations_df.groupby(["h", "patch_R", "comparison"], sort=True, dropna=False):
        h_value, r_value, comparison = key
        valid = group[group["valid_comparison"] == True].copy()
        s = pd.to_numeric(valid["spearman_rho"], errors="coerce") if not valid.empty else pd.Series(dtype=float)
        k = pd.to_numeric(valid["kendall_tau"], errors="coerce") if not valid.empty else pd.Series(dtype=float)
        pos_s = _fraction_positive(s)
        neg_s = _fraction_negative(s)
        flags = sorted(set(";".join(group["warning_flags"].dropna().astype(str)).split(";")) - {""})
        if valid.empty:
            flags.append("no_valid_rank_comparisons")
        rows.append({
            "h": _json_scalar(h_value), "patch_R": _json_scalar(r_value), "comparison": comparison,
            "n_capsids_total": int(group["capside"].nunique(dropna=True)), "n_valid_capsids": int(valid["capside"].nunique(dropna=True)),
            "capsids_used": _semicolon_join(valid["capside"]) if not valid.empty else "",
            "mean_spearman_rho": float(s.mean()) if len(s) else np.nan, "median_spearman_rho": float(s.median()) if len(s) else np.nan,
            "sd_spearman_rho": float(s.std(ddof=1)) if len(s) > 1 else np.nan, "min_spearman_rho": float(s.min()) if len(s) else np.nan, "max_spearman_rho": float(s.max()) if len(s) else np.nan,
            "fraction_positive_spearman": pos_s, "fraction_negative_spearman": neg_s,
            "mean_kendall_tau": float(k.mean()) if len(k) else np.nan, "median_kendall_tau": float(k.median()) if len(k) else np.nan,
            "sd_kendall_tau": float(k.std(ddof=1)) if len(k) > 1 else np.nan, "min_kendall_tau": float(k.min()) if len(k) else np.nan, "max_kendall_tau": float(k.max()) if len(k) else np.nan,
            "fraction_positive_kendall": _fraction_positive(k), "fraction_negative_kendall": _fraction_negative(k),
            "consistency_label": _consistency_label(pos_s, neg_s, threshold), "warning_flags": _flag_string(flags),
        })
    return pd.DataFrame(rows).sort_values(["h", "patch_R", "comparison"], kind="mergesort").reset_index(drop=True)


def _group_identifier(capside: Any, h_value: Any, r_value: Any) -> str:
    return f"{capside}|h={h_value}|R={r_value}"


def compute_top_rank_coincidence(rankings_df: pd.DataFrame) -> pd.DataFrame:
    group_rows: List[Dict[str, Any]] = []
    for key, group in rankings_df.groupby(PHASE04_RANK_GROUP_COLS, sort=True, dropna=False):
        capside, h_value, r_value = key
        gid = _group_identifier(capside, h_value, r_value)
        for comparison, a_flag, b_flag in PHASE04_COINCIDENCE_COMPARISONS:
            a = set(group.loc[group[a_flag] == True, "fold"].astype(str).tolist()) if a_flag in group.columns else set()
            b = set(group.loc[group[b_flag] == True, "fold"].astype(str).tolist()) if b_flag in group.columns else set()
            union = a | b
            inter = a & b
            valid = bool(a and b)
            group_rows.append({"h": h_value, "patch_R": r_value, "comparison": comparison, "group_id": gid, "valid": valid, "match": bool(inter) if valid else False, "jaccard": float(len(inter) / len(union)) if valid and union else np.nan})
    detail = pd.DataFrame(group_rows)
    rows: List[Dict[str, Any]] = []
    if detail.empty:
        return pd.DataFrame(columns=["h", "patch_R", "comparison", "n_groups_total", "n_groups_valid", "n_top1_matches", "fraction_top1_match", "mean_top1_jaccard", "median_top1_jaccard", "groups_with_match", "groups_without_match", "warning_flags"])
    for key, group in detail.groupby(["h", "patch_R", "comparison"], sort=True, dropna=False):
        h_value, r_value, comparison = key
        valid = group[group["valid"] == True]
        matched = valid[valid["match"] == True]
        unmatched = valid[valid["match"] == False]
        flags = [] if len(valid) else ["no_valid_rank_comparisons"]
        rows.append({
            "h": _json_scalar(h_value), "patch_R": _json_scalar(r_value), "comparison": comparison,
            "n_groups_total": int(len(group)), "n_groups_valid": int(len(valid)), "n_top1_matches": int(len(matched)),
            "fraction_top1_match": float(len(matched) / len(valid)) if len(valid) else np.nan,
            "mean_top1_jaccard": float(valid["jaccard"].mean()) if len(valid) else np.nan,
            "median_top1_jaccard": float(valid["jaccard"].median()) if len(valid) else np.nan,
            "groups_with_match": ";".join(sorted(matched["group_id"].astype(str).tolist())),
            "groups_without_match": ";".join(sorted(unmatched["group_id"].astype(str).tolist())),
            "warning_flags": _flag_string(flags),
        })
    return pd.DataFrame(rows).sort_values(["h", "patch_R", "comparison"], kind="mergesort").reset_index(drop=True)


def compute_fold_rank_consistency(rankings_df: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for variable, rank_col in PHASE04_RANK_COLUMNS.items():
        if rank_col not in rankings_df.columns:
            continue
        last_rank = rankings_df.groupby(PHASE04_RANK_GROUP_COLS, sort=False, dropna=False)[rank_col].transform("max")
        tmp = rankings_df[["h", "patch_R", "capside", "fold", rank_col]].copy()
        tmp["is_rank1"] = tmp[rank_col] == 1
        tmp["is_rank_last"] = tmp[rank_col] == last_rank
        tmp = tmp.dropna(subset=[rank_col])
        for key, group in tmp.groupby(["h", "patch_R", "fold"], sort=True, dropna=False):
            h_value, r_value, fold = key
            ranks = pd.to_numeric(group[rank_col], errors="coerce").dropna()
            rows.append({
                "h": _json_scalar(h_value), "patch_R": _json_scalar(r_value), "variable": variable, "fold": str(fold),
                "n_capsides": int(group["capside"].nunique(dropna=True)), "capsides_present": _semicolon_join(group["capside"]),
                "mean_rank": float(ranks.mean()) if len(ranks) else np.nan, "median_rank": float(ranks.median()) if len(ranks) else np.nan,
                "sd_rank": float(ranks.std(ddof=1)) if len(ranks) > 1 else np.nan, "min_rank": float(ranks.min()) if len(ranks) else np.nan, "max_rank": float(ranks.max()) if len(ranks) else np.nan,
                "fraction_rank1": float(group["is_rank1"].mean()) if len(group) else np.nan, "fraction_rank_last": float(group["is_rank_last"].mean()) if len(group) else np.nan,
                "warning_flags": "",
            })
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["h", "patch_R", "variable", "mean_rank", "fold"], kind="mergesort").reset_index(drop=True)


def make_phase04_ordinal_figures(rankings_df: pd.DataFrame, args: argparse.Namespace) -> Tuple[List[str], List[str]]:
    output_files: List[str] = []
    warnings: List[str] = []
    if args.no_plots:
        return output_files, warnings
    try:
        import matplotlib.pyplot as plt
        import numpy.ma as ma
    except Exception as exc:
        return output_files, [f"plot_failed: matplotlib unavailable ({exc})"]
    variables = [("d_max_pct", "rank_d_max_pct"), ("t_thinnest", "rank_t_thinnest"), ("Hout", "rank_Hout"), ("Hin", "rank_Hin")]
    one_condition = rankings_df[["h", "patch_R"]].drop_duplicates().shape[0] == 1
    try:
        for (h_value, r_value), subset in rankings_df.groupby(["h", "patch_R"], sort=True, dropna=False):
            folds = sorted(subset["fold"].dropna().astype(str).unique().tolist())
            capsides = sorted(subset["capside"].dropna().astype(str).unique().tolist())
            fig, axes = plt.subplots(2, 4, figsize=(18, 8), constrained_layout=True)
            cmap = plt.get_cmap("tab10")
            for col_idx, (label, rank_col) in enumerate(variables):
                ax = axes[0, col_idx]
                for idx, capside in enumerate(capsides):
                    g = subset[subset["capside"].astype(str) == capside].set_index("fold").reindex(folds)
                    ax.plot(folds, g[rank_col], marker="o", linewidth=1.2, color=cmap(idx % 10), label=capside)
                ax.invert_yaxis()
                ax.set_title(f"{label} ranks; rank 1 at top")
                ax.set_xlabel("fold")
                ax.set_ylabel("rank")
                ax.tick_params(axis="x", rotation=45)
                if col_idx == 0:
                    ax.legend(fontsize=7, title="capside")
                heat = np.full((len(folds), len(capsides)), np.nan)
                for j, capside in enumerate(capsides):
                    g = subset[subset["capside"].astype(str) == capside].set_index("fold")
                    for i, fold in enumerate(folds):
                        if fold in g.index:
                            heat[i, j] = pd.to_numeric(g.loc[fold, rank_col], errors="coerce")
                hax = axes[1, col_idx]
                masked = ma.masked_invalid(heat)
                im = hax.imshow(masked, aspect="auto", cmap="viridis_r")
                hax.set_title(f"{label} rank heatmap")
                hax.set_xticks(range(len(capsides)))
                hax.set_xticklabels(capsides, rotation=45, ha="right")
                hax.set_yticks(range(len(folds)))
                hax.set_yticklabels(folds)
                fig.colorbar(im, ax=hax, fraction=0.046, pad=0.04, label="rank (1 = top)")
            fig.suptitle(f"Phase 4 ordinal rankings (rank 1 = highest value; t_thinnest rank 1 = thinnest); h={h_value}, patch_R={r_value}")
            filename = "phase04_ordinal_rankings.png" if one_condition else f"phase04_ordinal_rankings_h_{_format_condition_value(h_value)}_R_{_format_condition_value(r_value)}.png"
            path = args.outdir / filename
            fig.savefig(path, dpi=150)
            plt.close(fig)
            output_files.append(str(path))
    except Exception as exc:
        warnings.append(f"plot_failed: {exc}")
    return output_files, warnings


def write_phase04_summary(path: Path, report: Dict, fold_consistency_df: pd.DataFrame, coincidence_df: pd.DataFrame, corr_summary_df: pd.DataFrame) -> None:
    lines = [
        f"Phase 4 ordinal fold-rank summary: {report['status']}",
        "This summary describes ordinal rank agreement only; it does not make causal or significance claims.",
        f"Capsids: {report['n_capsides']}",
        f"Folds: {report['n_folds']}",
        f"h values: {report['h_values']}",
        f"patch_R values: {report['patch_R_values']}",
        f"Ranking convention: {report['ranking_direction']}; rank_t_thinnest = 1 means thinnest fold.",
        "",
    ]
    for h_value in report["h_values"]:
        for r_value in report["patch_R_values"]:
            sub_fold = fold_consistency_df[(fold_consistency_df["h"].astype(str) == str(h_value)) & (fold_consistency_df["patch_R"].astype(str) == str(r_value)) & (fold_consistency_df["variable"] == "d_max_pct")]
            if sub_fold.empty:
                continue
            top = sub_fold.sort_values(["fraction_rank1", "mean_rank", "fold"], ascending=[False, True, True], kind="mergesort").iloc[0]
            lines.append(f"Condition h={h_value}, patch_R={r_value}:")
            lines.append(f"  - Fold most often rank 1 mechanically: {top['fold']} (fraction_rank1={top['fraction_rank1']:.3f})")
            match = coincidence_df[(coincidence_df["h"].astype(str) == str(h_value)) & (coincidence_df["patch_R"].astype(str) == str(r_value)) & (coincidence_df["comparison"] == "top_d_max_pct_vs_thinnest_t")]
            if not match.empty and np.isfinite(match.iloc[0]["fraction_top1_match"]):
                lines.append(f"  - top_d_max_pct_vs_thinnest_t fraction: {match.iloc[0]['fraction_top1_match']:.3f}")
            sub_corr = corr_summary_df[(corr_summary_df["h"].astype(str) == str(h_value)) & (corr_summary_df["patch_R"].astype(str) == str(r_value))]
            for _, row in sub_corr.iterrows():
                rho = row["median_spearman_rho"] if np.isfinite(row["median_spearman_rho"]) else np.nan
                lines.append(f"  - {row['comparison']}: median Spearman rho={rho:.6g}, consistency={row['consistency_label']}")
            lines.append("")
    lines.append("Warnings:")
    if report["warnings"]:
        lines.extend(f"  - {w}" for w in report["warnings"])
    else:
        lines.append("  none")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase04(input_df: pd.DataFrame, phase01_report: Dict, phase02_report: Dict, phase03_report: Dict, args: argparse.Namespace) -> Dict:
    args.outdir.mkdir(parents=True, exist_ok=True)
    output_files = {
        "fold_rankings": str(args.outdir / "phase04_fold_rankings.csv"),
        "rank_correlations_by_capsid": str(args.outdir / "phase04_rank_correlations_by_capsid.csv"),
        "rank_correlation_summary": str(args.outdir / "phase04_rank_correlation_summary.csv"),
        "top_rank_coincidence": str(args.outdir / "phase04_top_rank_coincidence.csv"),
        "fold_rank_consistency": str(args.outdir / "phase04_fold_rank_consistency.csv"),
        "report": str(args.outdir / "phase04_report.json"),
        "summary": str(args.outdir / "phase04_summary.txt"),
        "ordinal_figures": [],
    }
    warnings: List[str] = []
    severe_flags: List[str] = []
    if phase02_report.get("status") == "WARN":
        warnings.append("upstream_phase02_warn")
    if phase03_report.get("status") == "WARN":
        warnings.append("upstream_phase03_warn")
    if phase01_report.get("status") == "FAIL":
        severe_flags.append("phase01_validation_failed")
    if phase02_report.get("status") == "FAIL":
        severe_flags.append("phase02_failed")
    validation = validate_phase04_input(input_df) if not input_df.empty else {
        "missing_key_columns": PHASE04_KEY_COLUMNS, "missing_required_columns": PHASE04_REQUIRED_COLUMNS,
        "missing_geometry_columns": PHASE04_GEOMETRY_COLUMNS, "duplicate_key_status": {"is_unique": True, "n_duplicate_rows": 0, "duplicate_examples": []},
        "nonfinite_counts": {}, "n_ranking_groups": 0, "insufficient_ranking_groups": [],
    }
    for col in validation["missing_key_columns"]:
        severe_flags.append("missing_required_column")
    if "d_max_pct" in validation["missing_required_columns"]:
        severe_flags.append("missing_required_column")
    if not validation["duplicate_key_status"].get("is_unique", True):
        severe_flags.append("duplicate_observation_key")
    if validation["n_ranking_groups"] == 0:
        severe_flags.append("no_valid_ranking_group")
    for col in validation["missing_geometry_columns"]:
        warnings.append(f"missing_{col}")
    if validation["insufficient_ranking_groups"]:
        warnings.append("insufficient_folds_for_ranking")
    if validation["nonfinite_counts"].get("d_max_pct", 0):
        severe_flags.append("missing_d_max_pct")
    for col in ["t", "Hout", "Hin"]:
        if validation["nonfinite_counts"].get(col, 0):
            warnings.append(f"missing_{col}")
    if severe_flags:
        rankings_df = pd.DataFrame()
        corr_df = pd.DataFrame()
        corr_summary_df = pd.DataFrame()
        coincidence_df = pd.DataFrame()
        fold_consistency_df = pd.DataFrame()
        status = "FAIL"
    else:
        rankings_df = compute_fold_rankings(input_df, args)
        corr_df = compute_rank_correlations_by_capsid(rankings_df)
        corr_summary_df = compute_rank_correlation_summary(corr_df, float(args.rank_consistency_threshold))
        coincidence_df = compute_top_rank_coincidence(rankings_df)
        fold_consistency_df = compute_fold_rank_consistency(rankings_df)
        all_conditions_have_valid = True
        for _, group in corr_summary_df.groupby(["h", "patch_R"], sort=True, dropna=False):
            if int(group["n_valid_capsids"].sum()) == 0:
                all_conditions_have_valid = False
        if not all_conditions_have_valid:
            warnings.append("no_valid_rank_comparisons")
        if corr_df["warning_flags"].astype(str).str.contains("low_fold_count_rank_correlation", regex=False).any():
            warnings.append("low_fold_count_rank_correlation")
        if corr_df["warning_flags"].astype(str).str.contains("ties_present", regex=False).any() or rankings_df[[c for c in PHASE04_TIE_FLAG_COLUMNS.values() if c in rankings_df.columns]].any(axis=None):
            warnings.append("ties_present")
        if corr_df["warning_flags"].astype(str).str.contains("zero_rank_variance", regex=False).any():
            warnings.append("zero_rank_variance")
        plot_files, plot_warnings = make_phase04_ordinal_figures(rankings_df, args)
        output_files["ordinal_figures"] = plot_files
        if plot_warnings:
            warnings.extend(["plot_failed" if "plot_failed" in w else w for w in plot_warnings])
        status = "WARN" if warnings else "PASS"
        if not all_conditions_have_valid:
            status = "WARN"
        rankings_df.to_csv(output_files["fold_rankings"], index=False)
        corr_df.to_csv(output_files["rank_correlations_by_capsid"], index=False)
        corr_summary_df.to_csv(output_files["rank_correlation_summary"], index=False)
        coincidence_df.to_csv(output_files["top_rank_coincidence"], index=False)
        fold_consistency_df.to_csv(output_files["fold_rank_consistency"], index=False)
    report = {
        "phase": "04_ordinal_fold_rankings", "status": status, "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_file": str(args.input), "output_directory": str(args.outdir),
        "upstream_phase01_status": phase01_report.get("status"), "upstream_phase02_status": phase02_report.get("status"), "upstream_phase03_status": phase03_report.get("status"),
        "required_columns_checked": PHASE04_REQUIRED_COLUMNS + PHASE04_GEOMETRY_COLUMNS, "missing_columns": validation["missing_key_columns"] + [c for c in validation["missing_required_columns"] if c not in validation["missing_key_columns"]] + validation["missing_geometry_columns"],
        "n_rows_input": int(len(input_df)), "n_rows_output": int(len(rankings_df)),
        "n_capsides": int(input_df["capside"].nunique(dropna=True)) if "capside" in input_df.columns else 0,
        "n_folds": int(input_df["fold"].astype(str).nunique(dropna=True)) if "fold" in input_df.columns else 0,
        "h_values": _sorted_json_values(input_df["h"]) if "h" in input_df.columns else [], "patch_R_values": _sorted_json_values(input_df["patch_R"]) if "patch_R" in input_df.columns else [],
        "n_ranking_groups": validation["n_ranking_groups"], "rank_group_cols": PHASE04_RANK_GROUP_COLS,
        "ranked_variables": PHASE04_RANKED_VARIABLES, "ranking_method": getattr(args, "rank_method", "average"),
        "ranking_direction": "descending; rank 1 = highest value (rank_t_thinnest is ascending; rank 1 = thinnest fold)",
        "tie_handling": "average ranks; tied top folds are retained as sets",
        "rank_correlation_comparisons": [c[0] for c in PHASE04_COMPARISONS], "consistency_threshold": float(args.rank_consistency_threshold),
        "duplicate_key_status": validation["duplicate_key_status"], "output_files": output_files,
        "warnings": sorted(dict.fromkeys(warnings)), "severe_flags": sorted(dict.fromkeys(severe_flags)),
    }
    if status == "FAIL":
        # Write empty placeholders so failure diagnostics are deterministic, but only after no upstream commit decision depends on them.
        for key in ["fold_rankings", "rank_correlations_by_capsid", "rank_correlation_summary", "top_rank_coincidence", "fold_rank_consistency"]:
            pd.DataFrame().to_csv(output_files[key], index=False)
        (args.outdir / "phase04_summary.txt").write_text(f"Phase 4 ordinal fold-rank summary: FAIL\nSevere flags: {', '.join(report['severe_flags'])}\n", encoding="utf-8")
    else:
        write_phase04_summary(args.outdir / "phase04_summary.txt", report, fold_consistency_df, coincidence_df, corr_summary_df)
    with (args.outdir / "phase04_report.json").open("w", encoding="utf-8") as f:
        json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)
    return report


def _print_phase04_summary(report: Dict) -> None:
    print(f"[{report['status']}] Phase 4 completed.")
    files = report.get("output_files", {})
    print(f"[INFO] Fold rankings table: {files.get('fold_rankings', 'not written')}")
    print(f"[INFO] Rank correlations by capsid: {files.get('rank_correlations_by_capsid', 'not written')}")
    print(f"[INFO] Rank correlation summary: {files.get('rank_correlation_summary', 'not written')}")
    print(f"[INFO] Top-rank coincidence summary: {files.get('top_rank_coincidence', 'not written')}")
    print(f"[INFO] Fold-rank consistency summary: {files.get('fold_rank_consistency', 'not written')}")
    print(f"[INFO] Phase 4 report: {files.get('report', 'not written')}")
    fold_path = Path(files.get("fold_rank_consistency", ""))
    if fold_path.exists() and fold_path.stat().st_size > 0:
        try:
            fold_df = pd.read_csv(fold_path)
            mech = fold_df[fold_df["variable"] == "d_max_pct"] if "variable" in fold_df.columns else pd.DataFrame()
            for _, row in mech.sort_values(["h", "patch_R", "fraction_rank1", "mean_rank", "fold"], ascending=[True, True, False, True, True], kind="mergesort").groupby(["h", "patch_R"], sort=True).head(1).iterrows():
                print(f"[RANK] h={row['h']} patch_R={row['patch_R']} top_mechanical_fold={row['fold']} fraction_rank1={row['fraction_rank1']:.3f}")
        except Exception:
            pass
    match_path = Path(files.get("top_rank_coincidence", ""))
    if match_path.exists() and match_path.stat().st_size > 0:
        try:
            match_df = pd.read_csv(match_path)
            sub = match_df[match_df["comparison"] == "top_d_max_pct_vs_thinnest_t"] if "comparison" in match_df.columns else pd.DataFrame()
            for _, row in sub.iterrows():
                frac = row["fraction_top1_match"] if np.isfinite(row["fraction_top1_match"]) else np.nan
                print(f"[MATCH] h={row['h']} patch_R={row['patch_R']} comparison=top_d_vs_thinnest_t fraction={frac:.3f}")
        except Exception:
            pass
    if report.get("warnings"):
        print("[WARN] Phase 4 completed with warnings. See phase04_report.json.")


# ---------------------------------------------------------------------
# Main CLI
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------
# Phase 5: simple additive variance partitioning models
# ---------------------------------------------------------------------

PHASE05_CONDITION_COLS = ["h", "patch_R"]
PHASE05_REQUIRED_COLUMNS = ["capside", "fold", "h", "patch_R", "d_max_pct", "t", "Hout", "Hin"]
PHASE05_CORE_COLUMNS = ["d_max_pct", "capside", "fold", "h", "patch_R"]
PHASE05_GEOMETRY_COLUMNS = ["t", "Hout", "Hin"]
PHASE05_MODEL_ORDER = [
    "M0",
    "M1_capside",
    "M2_fold",
    "M3_capside_fold",
    "M3_geometry_base",
    "M4_capside_fold_geometry",
]
PHASE05_COMPARISON_ORDER = [
    "capside_vs_null",
    "fold_vs_null",
    "fold_added_after_capside",
    "capside_added_after_fold",
    "capside_plus_fold_vs_null",
    "geometry_added_after_capside_fold",
]
PHASE05_MODEL_FORMULAS = {
    "M0": {"label": "d_max_pct ~ 1", "categorical": [], "numeric": [], "dataset_type": "core_complete_case"},
    "M1_capside": {"label": "d_max_pct ~ capside", "categorical": ["capside"], "numeric": [], "dataset_type": "core_complete_case"},
    "M2_fold": {"label": "d_max_pct ~ fold", "categorical": ["fold"], "numeric": [], "dataset_type": "core_complete_case"},
    "M3_capside_fold": {"label": "d_max_pct ~ capside + fold", "categorical": ["capside", "fold"], "numeric": [], "dataset_type": "core_complete_case"},
    "M3_geometry_base": {"label": "d_max_pct ~ capside + fold", "categorical": ["capside", "fold"], "numeric": [], "dataset_type": "geometry_complete_case"},
    "M4_capside_fold_geometry": {"label": "d_max_pct ~ capside + fold + t + Hout + Hin", "categorical": ["capside", "fold"], "numeric": ["t", "Hout", "Hin"], "dataset_type": "geometry_complete_case"},
}
PHASE05_PRESERVE_COLUMNS = [
    "d_mag_max", "D", "delta_d_max_pct", "rel_delta_d_max_pct", "t_c", "Hout_c", "Hin_c",
    "H1", "H2", "H1_norm", "H2_norm", "H0", "H0_norm", "patch_elems", "size_h", "E",
    "rank_d_max_pct", "rank_t", "rank_t_thinnest", "rank_Hout", "rank_Hin",
]


def _phase05_float(value: Any) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out if np.isfinite(out) else float("nan")


def _phase05_condition_label(value: Any) -> str:
    text = str(value).strip().replace("-", "m").replace(".", "p")
    return "NA" if text == "" or text.lower() == "nan" else "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in text)


def _phase05_join_flags(flags: List[str]) -> str:
    return ";".join(dict.fromkeys([f for f in flags if f]))


def validate_phase05_input(df: pd.DataFrame) -> Dict[str, Any]:
    key_required = ["capside", "fold", "h", "patch_R", "d_max_pct"]
    missing_key = [c for c in key_required if c not in df.columns]
    missing_geometry = [c for c in PHASE05_GEOMETRY_COLUMNS if c not in df.columns]
    dup_rows = 0
    dup_examples: List[Dict[str, Any]] = []
    if all(c in df.columns for c in OBSERVATION_KEY):
        dup_mask = df.duplicated(subset=OBSERVATION_KEY, keep=False)
        dup_rows = int(dup_mask.sum())
        if dup_rows:
            dup_examples = df.loc[dup_mask, OBSERVATION_KEY].drop_duplicates().head(20).to_dict(orient="records")
    nonfinite: Dict[str, int] = {}
    for col in ["d_max_pct"] + [c for c in PHASE05_GEOMETRY_COLUMNS if c in df.columns]:
        vals = pd.to_numeric(df[col], errors="coerce")
        bad = int((~np.isfinite(vals)).sum())
        if bad:
            nonfinite[col] = bad
    return {
        "required_columns_checked": PHASE05_REQUIRED_COLUMNS,
        "missing_key_columns": missing_key,
        "missing_geometry_columns": missing_geometry,
        "missing_columns": missing_key + missing_geometry,
        "duplicate_key_status": {"observation_key": OBSERVATION_KEY, "is_unique": dup_rows == 0, "n_duplicate_rows": dup_rows, "duplicate_examples": dup_examples},
        "nonfinite_counts": nonfinite,
    }


def _phase05_clean_model_df(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    out = df.copy()
    for col in ["capside", "fold"]:
        if col in out.columns:
            out[col] = out[col].astype(str)
    for col in ["d_max_pct", "h", "patch_R", "t", "Hout", "Hin"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    mask = pd.Series(True, index=out.index)
    for col in columns:
        if col in ["capside", "fold"]:
            mask &= out[col].notna()
        else:
            mask &= np.isfinite(pd.to_numeric(out[col], errors="coerce"))
    return out.loc[mask].sort_values(["h", "patch_R", "capside", "fold"], kind="mergesort").reset_index(drop=True)


def prepare_phase05_model_datasets(df: pd.DataFrame) -> Dict[Tuple[Any, Any], Dict[str, pd.DataFrame]]:
    core = _phase05_clean_model_df(df, PHASE05_CORE_COLUMNS)
    geom_cols = PHASE05_CORE_COLUMNS + PHASE05_GEOMETRY_COLUMNS
    geometry = _phase05_clean_model_df(df, geom_cols) if all(c in df.columns for c in geom_cols) else pd.DataFrame(columns=df.columns)
    conditions = sorted(
        set(map(tuple, core[PHASE05_CONDITION_COLS].drop_duplicates().to_numpy().tolist())) |
        (set(map(tuple, geometry[PHASE05_CONDITION_COLS].drop_duplicates().to_numpy().tolist())) if not geometry.empty else set()),
        key=lambda x: (str(x[0]), str(x[1])),
    )
    out: Dict[Tuple[Any, Any], Dict[str, pd.DataFrame]] = {}
    for h, r in conditions:
        out[(h, r)] = {
            "core_complete_case": core[(core["h"] == h) & (core["patch_R"] == r)].reset_index(drop=True),
            "geometry_complete_case": geometry[(geometry["h"] == h) & (geometry["patch_R"] == r)].reset_index(drop=True) if not geometry.empty else pd.DataFrame(columns=df.columns),
        }
    return out


def _phase05_reference_levels(df: pd.DataFrame, cat_cols: List[str]) -> Dict[str, str]:
    refs = {}
    for col in cat_cols:
        levels = sorted(df[col].dropna().astype(str).unique().tolist())
        refs[col] = levels[0] if levels else ""
    return refs


def _phase05_build_design_matrix(df: pd.DataFrame, spec: Dict[str, Any], levels: Optional[Dict[str, List[str]]] = None) -> Tuple[np.ndarray, List[str], Dict[str, List[str]]]:
    n = len(df)
    cols = [np.ones(n, dtype=float)]
    names = ["Intercept"]
    used_levels: Dict[str, List[str]] = {}
    for cat in spec["categorical"]:
        if levels and cat in levels:
            cat_levels = levels[cat]
        else:
            cat_levels = sorted(df[cat].dropna().astype(str).unique().tolist())
        used_levels[cat] = cat_levels
        for level in cat_levels[1:]:
            cols.append((df[cat].astype(str).to_numpy() == level).astype(float))
            names.append(f"{cat}[T.{level}]")
    for num in spec["numeric"]:
        vals = pd.to_numeric(df[num], errors="coerce").to_numpy(dtype=float)
        cols.append(vals)
        names.append(num)
    X = np.column_stack(cols) if cols else np.empty((n, 0))
    return X, names, used_levels


def _phase05_invalid_fit(h: Any, patch_R: Any, model_id: str, dataset_type: str, df: pd.DataFrame, flags: List[str]) -> Dict[str, Any]:
    y = pd.to_numeric(df["d_max_pct"], errors="coerce") if "d_max_pct" in df.columns and len(df) else pd.Series(dtype=float)
    return {
        "h": h, "patch_R": patch_R, "model_id": model_id, "model_formula_label": PHASE05_MODEL_FORMULAS[model_id]["label"], "dataset_type": dataset_type,
        "n_observations": int(len(df)), "n_capsides": int(df["capside"].nunique()) if "capside" in df.columns and len(df) else 0,
        "n_folds": int(df["fold"].nunique()) if "fold" in df.columns and len(df) else 0, "n_predictor_columns": np.nan,
        "n_effective_parameters": np.nan, "df_model": np.nan, "df_residual": np.nan,
        "response_mean": float(y.mean()) if len(y) else np.nan, "response_sd": float(y.std(ddof=1)) if len(y) > 1 else np.nan,
        "rss": np.nan, "tss": np.nan, "mse": np.nan, "rmse": np.nan, "mae": np.nan, "r_squared": np.nan, "r_squared_adjusted": np.nan,
        "loocv_mse": np.nan, "loocv_rmse": np.nan, "loocv_mae": np.nan, "loocv_n_predictions": 0, "loocv_n_failed_predictions": int(len(df)),
        "aic": np.nan, "bic": np.nan, "design_matrix_rank": np.nan, "design_matrix_condition_number": np.nan,
        "rank_deficient": False, "overfit_warning": _phase05_join_flags([f for f in flags if "overfit" in f or "saturated" in f]),
        "valid_model": False, "warning_flags": _phase05_join_flags(flags),
        "fitted_values": np.full(len(df), np.nan), "residuals": np.full(len(df), np.nan), "coefficients": [], "loocv_predictions": np.full(len(df), np.nan),
        "loocv_valid": np.zeros(len(df), dtype=bool), "loocv_flags": ["model_invalid"] * len(df), "design_levels": {}, "coef_names": [], "params": np.array([]),
    }


def fit_phase05_model(df: pd.DataFrame, model_id: str, h: Any, patch_R: Any, args: argparse.Namespace) -> Dict[str, Any]:
    spec = PHASE05_MODEL_FORMULAS[model_id]
    dataset_type = spec["dataset_type"]
    flags: List[str] = []
    n = int(len(df))
    n_caps = int(df["capside"].nunique()) if "capside" in df.columns and n else 0
    n_folds = int(df["fold"].nunique()) if "fold" in df.columns and n else 0
    if n < args.model_min_n:
        flags.append("insufficient_observations")
    if model_id in ["M1_capside", "M3_capside_fold", "M3_geometry_base", "M4_capside_fold_geometry"] and n_caps < 2:
        flags.append("insufficient_capsides")
    if model_id in ["M2_fold", "M3_capside_fold", "M3_geometry_base", "M4_capside_fold_geometry"] and n_folds < 2:
        flags.append("insufficient_folds")
    if model_id == "M4_capside_fold_geometry":
        if any(c not in df.columns for c in PHASE05_GEOMETRY_COLUMNS):
            flags.append("m4_not_available")
        elif not any(float(pd.to_numeric(df[c], errors="coerce").var(ddof=0)) > 0 for c in PHASE05_GEOMETRY_COLUMNS if c in df.columns):
            flags.append("zero_geometry_variance")
    y = pd.to_numeric(df["d_max_pct"], errors="coerce").to_numpy(dtype=float) if "d_max_pct" in df.columns else np.array([])
    if n and (not np.isfinite(y).all()):
        flags.append("nonfinite_d_max_pct")
    tss = float(np.sum((y - np.mean(y)) ** 2)) if n else np.nan
    if n and (not np.isfinite(tss) or tss <= 0):
        flags.append("zero_response_variance")
    fatal = ["insufficient_observations", "insufficient_capsides", "insufficient_folds", "zero_response_variance", "nonfinite_d_max_pct", "zero_geometry_variance", "m4_not_available"]
    if any(f in flags for f in fatal):
        return _phase05_invalid_fit(h, patch_R, model_id, dataset_type, df, flags)
    try:
        X, coef_names, design_levels = _phase05_build_design_matrix(df, spec)
        if not np.isfinite(X).all():
            return _phase05_invalid_fit(h, patch_R, model_id, dataset_type, df, flags + ["design_matrix_nonfinite"])
        params, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        fitted = X @ params
        resid = y - fitted
        rss = float(np.sum(resid ** 2))
        rank = int(np.linalg.matrix_rank(X))
        n_cols = int(X.shape[1])
        rank_def = rank < n_cols
        if rank_def:
            flags.append("rank_deficient_design")
        cond = float(np.linalg.cond(X)) if n_cols and rank else np.nan
        if np.isfinite(cond) and cond > args.condition_number_threshold:
            flags.append("high_condition_number")
        p_eff = max(rank - 1, 0)
        df_model = p_eff
        df_resid = n - rank
        if n <= rank or n <= p_eff:
            flags.append("saturated_or_overfit_model")
        if df_resid <= 0:
            flags.append("no_residual_degrees_of_freedom")
        r2 = np.nan if tss <= 0 else 1.0 - rss / tss
        if np.isfinite(r2) and np.isclose(r2, 1.0) and ("saturated_or_overfit_model" in flags or df_resid <= 1):
            flags.append("saturated_r2_not_interpretable")
        if n - p_eff - 1 <= 0 or not np.isfinite(r2):
            adj_r2 = np.nan
            flags.append("adjusted_r2_undefined")
        else:
            adj_r2 = 1.0 - (1.0 - r2) * (n - 1) / (n - p_eff - 1)
        mse = rss / n
        rmse = float(np.sqrt(mse))
        mae = float(np.mean(np.abs(resid)))
        k = rank
        aic = float(n * np.log(rss / n) + 2 * k) if rss > 0 else np.nan
        bic = float(n * np.log(rss / n) + np.log(n) * k) if rss > 0 else np.nan
        loocv = compute_loocv_metrics(df, model_id, args)
        flags.extend(loocv["model_flags"])
        coefs = []
        se = np.full(len(params), np.nan)
        tvals = np.full(len(params), np.nan)
        pvals = np.full(len(params), np.nan)
        if df_resid > 0 and not rank_def and rss >= 0:
            try:
                sigma2 = rss / df_resid
                cov = sigma2 * np.linalg.pinv(X.T @ X)
                se = np.sqrt(np.diag(cov))
                tvals = params / se
            except Exception:
                flags.append("coefficient_standard_errors_unavailable")
        for name, val, sev, tv, pv in zip(coef_names, params, se, tvals, pvals):
            coefs.append({
                "h": h, "patch_R": patch_R, "model_id": model_id, "coefficient_name": name,
                "coefficient_value": float(val), "coefficient_standard_error": float(sev) if np.isfinite(sev) else np.nan,
                "coefficient_t_value": float(tv) if np.isfinite(tv) else np.nan,
                "coefficient_p_value": float(pv) if np.isfinite(pv) else np.nan,
                "coefficient_warning_flags": "" if np.isfinite(sev) else "standard_error_unavailable",
            })
        overfit_flags = [f for f in flags if f in {"saturated_or_overfit_model", "no_residual_degrees_of_freedom", "saturated_r2_not_interpretable", "high_condition_number", "unstable_loocv", "loocv_failed"}]
        return {
            "h": h, "patch_R": patch_R, "model_id": model_id, "model_formula_label": spec["label"], "dataset_type": dataset_type,
            "n_observations": n, "n_capsides": n_caps, "n_folds": n_folds, "n_predictor_columns": n_cols - 1,
            "n_effective_parameters": p_eff, "df_model": df_model, "df_residual": df_resid,
            "response_mean": float(np.mean(y)), "response_sd": float(np.std(y, ddof=1)) if n > 1 else np.nan,
            "rss": rss, "tss": tss, "mse": mse, "rmse": rmse, "mae": mae, "r_squared": float(r2) if np.isfinite(r2) else np.nan,
            "r_squared_adjusted": float(adj_r2) if np.isfinite(adj_r2) else np.nan, "loocv_mse": loocv["mse"], "loocv_rmse": loocv["rmse"], "loocv_mae": loocv["mae"],
            "loocv_n_predictions": loocv["n_predictions"], "loocv_n_failed_predictions": loocv["n_failed"], "aic": aic, "bic": bic,
            "design_matrix_rank": rank, "design_matrix_condition_number": cond, "rank_deficient": bool(rank_def),
            "overfit_warning": _phase05_join_flags(overfit_flags), "valid_model": bool(np.isfinite(rss) and np.isfinite(fitted).all()), "warning_flags": _phase05_join_flags(flags),
            "fitted_values": fitted, "residuals": resid, "coefficients": coefs, "loocv_predictions": loocv["predictions"], "loocv_valid": loocv["valid"],
            "loocv_flags": loocv["row_flags"], "design_levels": design_levels, "coef_names": coef_names, "params": params,
        }
    except Exception as exc:
        return _phase05_invalid_fit(h, patch_R, model_id, dataset_type, df, flags + [f"fit_failed:{exc}"])


def compute_loocv_metrics(df: pd.DataFrame, model_id: str, args: argparse.Namespace) -> Dict[str, Any]:
    spec = PHASE05_MODEL_FORMULAS[model_id]
    n = int(len(df))
    y = pd.to_numeric(df["d_max_pct"], errors="coerce").to_numpy(dtype=float)
    preds = np.full(n, np.nan)
    valid = np.zeros(n, dtype=bool)
    row_flags: List[str] = [""] * n
    failed = 0
    model_flags: List[str] = []
    for i in range(n):
        train = df.drop(index=df.index[i]).reset_index(drop=True)
        test = df.iloc[[i]].reset_index(drop=True)
        split_flags: List[str] = []
        if len(train) < args.model_min_n:
            split_flags.append("insufficient_observations")
        for cat in spec["categorical"]:
            if str(test.iloc[0][cat]) not in set(train[cat].astype(str).unique().tolist()):
                split_flags.append("loocv_unseen_category")
        if split_flags:
            failed += 1
            row_flags[i] = _phase05_join_flags(split_flags)
            continue
        try:
            levels = {cat: sorted(train[cat].dropna().astype(str).unique().tolist()) for cat in spec["categorical"]}
            Xtr, _, _ = _phase05_build_design_matrix(train, spec, levels=levels)
            ytr = pd.to_numeric(train["d_max_pct"], errors="coerce").to_numpy(dtype=float)
            if not np.isfinite(Xtr).all() or not np.isfinite(ytr).all() or float(np.sum((ytr - np.mean(ytr)) ** 2)) <= 0:
                raise ValueError("invalid_training_split")
            beta, _, _, _ = np.linalg.lstsq(Xtr, ytr, rcond=None)
            Xte, _, _ = _phase05_build_design_matrix(test, spec, levels=levels)
            preds[i] = float(Xte @ beta)
            valid[i] = np.isfinite(preds[i])
            if not valid[i]:
                row_flags[i] = "loocv_failed"
                failed += 1
        except Exception:
            failed += 1
            row_flags[i] = "loocv_failed"
    errors = y[valid] - preds[valid]
    n_pred = int(valid.sum())
    if n_pred == 0:
        model_flags.append("loocv_failed")
        mse = rmse = mae = np.nan
    else:
        mse = float(np.mean(errors ** 2))
        rmse = float(np.sqrt(mse))
        mae = float(np.mean(np.abs(errors)))
    min_success = max(3, float(args.loocv_min_success_fraction) * n)
    if n_pred < min_success:
        model_flags.append("unstable_loocv")
    if any("loocv_unseen_category" in f for f in row_flags):
        model_flags.append("loocv_unseen_category")
    return {"mse": mse, "rmse": rmse, "mae": mae, "n_predictions": n_pred, "n_failed": int(n - n_pred), "predictions": preds, "valid": valid, "row_flags": row_flags, "model_flags": model_flags}


def _phase05_metric_rows(fits: List[Dict[str, Any]]) -> pd.DataFrame:
    cols = [
        "h", "patch_R", "model_id", "model_formula_label", "dataset_type", "n_observations", "n_capsides", "n_folds", "n_predictor_columns", "n_effective_parameters", "df_model", "df_residual", "response_mean", "response_sd",
        "rss", "tss", "mse", "rmse", "mae", "r_squared", "r_squared_adjusted", "loocv_mse", "loocv_rmse", "loocv_mae", "loocv_n_predictions", "loocv_n_failed_predictions", "aic", "bic",
        "design_matrix_rank", "design_matrix_condition_number", "rank_deficient", "overfit_warning", "valid_model", "warning_flags",
    ]
    df = pd.DataFrame([{c: f.get(c, np.nan) for c in cols} for f in fits])
    df["_order"] = df["model_id"].map({m: i for i, m in enumerate(PHASE05_MODEL_ORDER)})
    return df.sort_values(["h", "patch_R", "_order"], kind="mergesort").drop(columns="_order")


def build_phase05_predictions_table(condition_datasets: Dict[Tuple[Any, Any], Dict[str, pd.DataFrame]], fits: List[Dict[str, Any]]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for fit in fits:
        df = condition_datasets[(fit["h"], fit["patch_R"])][fit["dataset_type"]]
        y = pd.to_numeric(df["d_max_pct"], errors="coerce").to_numpy(dtype=float) if len(df) else np.array([])
        fitted = fit.get("fitted_values", np.full(len(df), np.nan))
        loocv = fit.get("loocv_predictions", np.full(len(df), np.nan))
        valid = fit.get("loocv_valid", np.zeros(len(df), dtype=bool))
        row_flags = fit.get("loocv_flags", [""] * len(df))
        for i, row in df.reset_index(drop=True).iterrows():
            resid = y[i] - fitted[i] if i < len(fitted) and np.isfinite(fitted[i]) else np.nan
            cv_resid = y[i] - loocv[i] if i < len(loocv) and np.isfinite(loocv[i]) else np.nan
            rows.append({
                "h": fit["h"], "patch_R": fit["patch_R"], "model_id": fit["model_id"], "capside": row["capside"], "fold": row["fold"],
                "observed_d_max_pct": y[i], "fitted_d_max_pct": fitted[i] if i < len(fitted) else np.nan, "residual": resid,
                "abs_residual": abs(resid) if np.isfinite(resid) else np.nan, "squared_residual": resid ** 2 if np.isfinite(resid) else np.nan,
                "loocv_predicted_d_max_pct": loocv[i] if i < len(loocv) else np.nan, "loocv_residual": cv_resid,
                "loocv_abs_residual": abs(cv_resid) if np.isfinite(cv_resid) else np.nan, "loocv_squared_residual": cv_resid ** 2 if np.isfinite(cv_resid) else np.nan,
                "loocv_prediction_valid": bool(valid[i]) if i < len(valid) else False, "loocv_warning_flags": row_flags[i] if i < len(row_flags) else "",
            })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["_order"] = out["model_id"].map({m: i for i, m in enumerate(PHASE05_MODEL_ORDER)})
        out = out.sort_values(["h", "patch_R", "_order", "capside", "fold"], kind="mergesort").drop(columns="_order")
    return out


def build_phase05_coefficients_table(fits: List[Dict[str, Any]]) -> pd.DataFrame:
    rows = [coef for fit in fits for coef in fit.get("coefficients", [])]
    cols = ["h", "patch_R", "model_id", "coefficient_name", "coefficient_value", "coefficient_standard_error", "coefficient_t_value", "coefficient_p_value", "coefficient_warning_flags"]
    out = pd.DataFrame(rows, columns=cols)
    if not out.empty:
        out["_order"] = out["model_id"].map({m: i for i, m in enumerate(PHASE05_MODEL_ORDER)})
        out = out.sort_values(["h", "patch_R", "_order", "coefficient_name"], kind="mergesort").drop(columns="_order")
    return out


def compute_phase05_contributions(metrics_df: pd.DataFrame) -> pd.DataFrame:
    comp_specs = [
        ("capside_vs_null", "capside_vs_null", "M0", "M1_capside"),
        ("fold_vs_null", "fold_vs_null", "M0", "M2_fold"),
        ("fold_added_after_capside", "fold_added_after_capside", "M1_capside", "M3_capside_fold"),
        ("capside_added_after_fold", "capside_added_after_fold", "M2_fold", "M3_capside_fold"),
        ("capside_plus_fold_vs_null", "capside_plus_fold_vs_null", "M0", "M3_capside_fold"),
        ("geometry_added_after_capside_fold", "geometry_added_after_capside_fold", "M3_geometry_base", "M4_capside_fold_geometry"),
    ]
    rows: List[Dict[str, Any]] = []
    for (h, r), group in metrics_df.groupby(["h", "patch_R"], sort=True, dropna=False):
        by_model = {row["model_id"]: row for _, row in group.iterrows()}
        for comp_id, label, base_id, ext_id in comp_specs:
            flags: List[str] = []
            base = by_model.get(base_id)
            ext = by_model.get(ext_id)
            if base is None or ext is None or not bool(base.get("valid_model", False)) or not bool(ext.get("valid_model", False)):
                flags.append("comparison_invalid")
                interp = "comparison_invalid"
            else:
                interp = "comparison_invalid"
            same_rows = bool(base is not None and ext is not None and base["dataset_type"] == ext["dataset_type"] and int(base["n_observations"]) == int(ext["n_observations"]))
            if not same_rows:
                flags.append("different_rows_for_comparison")
            def val(row, col):
                return row[col] if row is not None else np.nan
            d_r2 = val(ext, "r_squared") - val(base, "r_squared") if base is not None and ext is not None else np.nan
            d_adj = val(ext, "r_squared_adjusted") - val(base, "r_squared_adjusted") if base is not None and ext is not None else np.nan
            d_rmse = val(ext, "loocv_rmse") - val(base, "loocv_rmse") if base is not None and ext is not None else np.nan
            d_mae = val(ext, "loocv_mae") - val(base, "loocv_mae") if base is not None and ext is not None else np.nan
            if "comparison_invalid" not in flags:
                if np.isfinite(d_r2) and d_r2 > 0 and np.isfinite(d_rmse) and d_rmse < 0:
                    interp = "added_variance_and_improved_loocv"
                elif np.isfinite(d_r2) and d_r2 > 0:
                    interp = "added_in_sample_variance_only"
                elif np.isfinite(d_r2):
                    interp = "no_in_sample_gain"
                else:
                    interp = "comparison_invalid"; flags.append("comparison_invalid")
            rows.append({
                "h": h, "patch_R": r, "comparison_id": comp_id, "contribution_label": label, "base_model": base_id, "extended_model": ext_id,
                "base_dataset_type": val(base, "dataset_type"), "extended_dataset_type": val(ext, "dataset_type"), "same_rows_for_comparison": same_rows,
                "base_n_observations": val(base, "n_observations"), "extended_n_observations": val(ext, "n_observations"),
                "base_r_squared": val(base, "r_squared"), "extended_r_squared": val(ext, "r_squared"), "delta_r_squared": d_r2,
                "base_r_squared_adjusted": val(base, "r_squared_adjusted"), "extended_r_squared_adjusted": val(ext, "r_squared_adjusted"), "delta_r_squared_adjusted": d_adj,
                "base_loocv_rmse": val(base, "loocv_rmse"), "extended_loocv_rmse": val(ext, "loocv_rmse"), "delta_loocv_rmse": d_rmse,
                "base_loocv_mae": val(base, "loocv_mae"), "extended_loocv_mae": val(ext, "loocv_mae"), "delta_loocv_mae": d_mae,
                "interpretation_label": interp, "warning_flags": _phase05_join_flags(flags),
            })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["_order"] = out["comparison_id"].map({m: i for i, m in enumerate(PHASE05_COMPARISON_ORDER)})
        out = out.sort_values(["h", "patch_R", "_order"], kind="mergesort").drop(columns="_order")
    return out


def build_phase05_overfit_diagnostics(metrics_df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    for _, row in metrics_df.iterrows():
        n_params = _phase05_float(row["n_effective_parameters"])
        n_obs = _phase05_float(row["n_observations"])
        denom = n_params + 1 if np.isfinite(n_params) else np.nan
        opp = n_obs / denom if np.isfinite(denom) and denom > 0 else np.nan
        cond = _phase05_float(row["design_matrix_condition_number"])
        flags = [] if pd.isna(row["warning_flags"]) or row["warning_flags"] == "" else str(row["warning_flags"]).split(";")
        if np.isfinite(opp) and opp < args.min_observations_per_parameter:
            flags.append("low_observations_per_parameter")
        if np.isfinite(cond) and cond > args.condition_number_threshold:
            flags.append("high_condition_number")
        if _phase05_float(row["df_residual"]) <= 1:
            flags.append("low_residual_degrees_of_freedom")
        rows.append({
            "h": row["h"], "patch_R": row["patch_R"], "model_id": row["model_id"], "n_observations": row["n_observations"],
            "n_effective_parameters": row["n_effective_parameters"], "df_residual": row["df_residual"], "observations_per_parameter": opp,
            "design_matrix_rank": row["design_matrix_rank"], "rank_deficient": row["rank_deficient"], "condition_number": row["design_matrix_condition_number"],
            "r_squared": row["r_squared"], "r_squared_adjusted": row["r_squared_adjusted"], "loocv_rmse": row["loocv_rmse"],
            "saturated_or_near_saturated": bool("saturated_or_overfit_model" in flags or "low_residual_degrees_of_freedom" in flags),
            "high_condition_number": bool("high_condition_number" in flags), "loocv_unstable": bool("unstable_loocv" in flags or "loocv_failed" in flags),
            "overfit_warning": row["overfit_warning"], "warning_flags": _phase05_join_flags(flags),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["_order"] = out["model_id"].map({m: i for i, m in enumerate(PHASE05_MODEL_ORDER)})
        out = out.sort_values(["h", "patch_R", "_order"], kind="mergesort").drop(columns="_order")
    return out


def make_phase05_model_comparison_figures(metrics_df: pd.DataFrame, contrib_df: pd.DataFrame, outdir: Path) -> Tuple[List[str], List[str]]:
    files: List[str] = []
    warnings: List[str] = []
    try:
        import matplotlib.pyplot as plt
        groups = list(metrics_df.groupby(["h", "patch_R"], sort=True, dropna=False))
        for (h, r), group in groups:
            fig, axes = plt.subplots(2, 2, figsize=(12, 8))
            x = np.arange(len(group))
            labels = group["model_id"].tolist()
            for ax, col, title in [(axes[0,0], "r_squared", "R2"), (axes[0,1], "r_squared_adjusted", "Adjusted R2"), (axes[1,0], "loocv_rmse", "LOOCV RMSE")]:
                vals = pd.to_numeric(group[col], errors="coerce").to_numpy(dtype=float)
                ax.bar(x, vals)
                ax.set_xticks(x); ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
                ax.set_title(title)
                for j, (_, rr) in enumerate(group.iterrows()):
                    if not bool(rr["valid_model"]):
                        ax.text(j, 0, "invalid", rotation=90, ha="center", va="bottom", fontsize=7, color="red")
                    elif "saturated" in str(rr["warning_flags"]):
                        ax.text(j, vals[j] if np.isfinite(vals[j]) else 0, "sat", ha="center", va="bottom", fontsize=7, color="orange")
            cg = contrib_df[(contrib_df["h"] == h) & (contrib_df["patch_R"] == r)]
            axes[1,1].bar(np.arange(len(cg)), pd.to_numeric(cg["delta_r_squared"], errors="coerce"))
            axes[1,1].set_xticks(np.arange(len(cg))); axes[1,1].set_xticklabels(cg["contribution_label"].tolist(), rotation=35, ha="right", fontsize=8)
            axes[1,1].set_title("Delta R2")
            fig.suptitle(f"Phase 5 additive model comparison: h={h}, patch_R={r}")
            fig.tight_layout()
            name = "phase05_model_comparison.png" if len(groups) == 1 else f"phase05_model_comparison_h_{_phase05_condition_label(h)}_R_{_phase05_condition_label(r)}.png"
            path = outdir / name
            fig.savefig(path, dpi=150)
            plt.close(fig)
            files.append(str(path))
    except Exception as exc:
        warnings.append(f"plot_failed:{exc}")
    return files, warnings


def write_phase05_summary(path: Path, report: Dict[str, Any], metrics_df: pd.DataFrame, contrib_df: pd.DataFrame, overfit_df: pd.DataFrame) -> None:
    lines = [
        f"Phase 5 additive variance partitioning: {report['status']}",
        f"Capsids: {report['n_capsides']}",
        f"Folds: {report['n_folds']}",
        f"h values: {report['h_values']}",
        f"patch_R values: {report['patch_R_values']}",
        "",
    ]
    for (h, r), group in metrics_df.groupby(["h", "patch_R"], sort=True, dropna=False):
        lines.append(f"Condition h={h}, patch_R={r}")
        for _, row in group.iterrows():
            lines.append(f"  {row['model_id']}: R2={row['r_squared']:.6g} adjusted_R2={row['r_squared_adjusted']:.6g} LOOCV_RMSE={row['loocv_rmse']:.6g} warnings={row['warning_flags'] or 'none'}")
        cg = contrib_df[(contrib_df["h"] == h) & (contrib_df["patch_R"] == r)]
        for _, row in cg.iterrows():
            lines.append(f"  Delta {row['contribution_label']}: delta_R2={row['delta_r_squared']:.6g} delta_LOOCV_RMSE={row['delta_loocv_rmse']:.6g} label={row['interpretation_label']}")
        geom = cg[cg["comparison_id"] == "geometry_added_after_capside_fold"]
        if not geom.empty:
            grow = geom.iloc[0]
            if np.isfinite(_phase05_float(grow["delta_loocv_rmse"])) and grow["delta_loocv_rmse"] > 0:
                lines.append("  Conservative interpretation: M4 increases in-sample R2 relative to M3, but LOOCV RMSE does not improve, suggesting that the added geometry terms may not provide stable predictive information under this dataset.")
        od = overfit_df[(overfit_df["h"] == h) & (overfit_df["patch_R"] == r)]
        warn = sorted(set(";".join(od["warning_flags"].dropna().astype(str)).split(";")) - {""})
        lines.append(f"  strongest overfit warning: {warn[0] if warn else 'none'}")
        lines.append("")
    lines.append("Warnings:")
    lines.extend([f"  - {w}" for w in report.get("warnings", [])] or ["  none"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase05(input_df: pd.DataFrame, phase01_report: Dict, phase02_report: Dict, phase03_report: Dict, phase04_report: Dict, args: argparse.Namespace) -> Dict:
    args.outdir.mkdir(parents=True, exist_ok=True)
    output_files = {
        "model_metrics": str(args.outdir / "phase05_model_metrics.csv"),
        "model_contributions": str(args.outdir / "phase05_model_contributions.csv"),
        "model_predictions": str(args.outdir / "phase05_model_predictions.csv"),
        "model_coefficients": str(args.outdir / "phase05_model_coefficients.csv"),
        "overfit_diagnostics": str(args.outdir / "phase05_overfit_diagnostics.csv"),
        "report": str(args.outdir / "phase05_report.json"),
        "summary": str(args.outdir / "phase05_summary.txt"),
        "figures": [],
    }
    validation = validate_phase05_input(input_df)
    warnings: List[str] = []
    severe_flags: List[str] = []
    for phase_name, rep in [("phase02", phase02_report), ("phase03", phase03_report), ("phase04", phase04_report)]:
        if rep.get("status") == "WARN":
            warnings.append(f"upstream_{phase_name}_warn")
    if phase01_report.get("status") == "FAIL":
        severe_flags.append("phase01_validation_failed")
    if phase02_report.get("status") == "FAIL":
        severe_flags.append("phase02_validation_failed")
    if validation["missing_key_columns"]:
        severe_flags.append("missing_required_column")
    if not validation["duplicate_key_status"]["is_unique"]:
        severe_flags.append("duplicate_observation_key")
    if validation["nonfinite_counts"].get("d_max_pct", 0):
        warnings.append("nonfinite_d_max_pct")
    if validation["missing_geometry_columns"]:
        warnings.append("missing_geometry_column")
    for geom in PHASE05_GEOMETRY_COLUMNS:
        if validation["nonfinite_counts"].get(geom, 0):
            warnings.append(f"nonfinite_{geom}")
    empty_cols = {
        "metrics": ["h", "patch_R", "model_id"],
        "contrib": ["h", "patch_R", "comparison_id"],
    }
    if severe_flags:
        pd.DataFrame(columns=empty_cols["metrics"]).to_csv(output_files["model_metrics"], index=False)
        pd.DataFrame(columns=empty_cols["contrib"]).to_csv(output_files["model_contributions"], index=False)
        pd.DataFrame().to_csv(output_files["model_predictions"], index=False)
        pd.DataFrame().to_csv(output_files["model_coefficients"], index=False)
        pd.DataFrame().to_csv(output_files["overfit_diagnostics"], index=False)
        report = _phase05_report(input_df, validation, output_files, warnings, severe_flags, "FAIL", phase01_report, phase02_report, phase03_report, phase04_report, args)
        with (args.outdir / "phase05_report.json").open("w", encoding="utf-8") as f:
            json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)
        write_phase05_summary(args.outdir / "phase05_summary.txt", report, pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
        return report
    condition_datasets = prepare_phase05_model_datasets(input_df)
    fits: List[Dict[str, Any]] = []
    reference_levels: Dict[str, Dict[str, str]] = {}
    for (h, r), datasets in condition_datasets.items():
        for model_id in PHASE05_MODEL_ORDER:
            ds_type = PHASE05_MODEL_FORMULAS[model_id]["dataset_type"]
            fit = fit_phase05_model(datasets[ds_type], model_id, h, r, args)
            fits.append(fit)
            refs = _phase05_reference_levels(datasets[ds_type], PHASE05_MODEL_FORMULAS[model_id]["categorical"])
            reference_levels[f"h={h};patch_R={r};{model_id}"] = refs
    metrics_df = _phase05_metric_rows(fits)
    contrib_df = compute_phase05_contributions(metrics_df)
    predictions_df = build_phase05_predictions_table(condition_datasets, fits)
    coef_df = build_phase05_coefficients_table(fits)
    overfit_df = build_phase05_overfit_diagnostics(metrics_df, args)
    warnings.extend(sorted(set(";".join(metrics_df["warning_flags"].dropna().astype(str)).split(";")) - {""}))
    warnings.extend(sorted(set(";".join(contrib_df["warning_flags"].dropna().astype(str)).split(";")) - {""}))
    if not bool(metrics_df["valid_model"].any()):
        severe_flags.append("no_valid_model")
    if not condition_datasets:
        severe_flags.append("no_valid_model_condition")
    if not args.no_plots:
        figs, plot_warnings = make_phase05_model_comparison_figures(metrics_df, contrib_df, args.outdir)
        output_files["figures"] = figs
        warnings.extend(plot_warnings)
    status = "FAIL" if severe_flags else ("WARN" if warnings or not metrics_df["valid_model"].all() else "PASS")
    report = _phase05_report(input_df, validation, output_files, sorted(dict.fromkeys(warnings)), severe_flags, status, phase01_report, phase02_report, phase03_report, phase04_report, args)
    report["categorical_encoding"]["reference_levels_by_condition_model"] = reference_levels
    metrics_df.to_csv(output_files["model_metrics"], index=False)
    contrib_df.to_csv(output_files["model_contributions"], index=False)
    predictions_df.to_csv(output_files["model_predictions"], index=False)
    coef_df.to_csv(output_files["model_coefficients"], index=False)
    overfit_df.to_csv(output_files["overfit_diagnostics"], index=False)
    with (args.outdir / "phase05_report.json").open("w", encoding="utf-8") as f:
        json.dump(_nan_to_none(report), f, indent=2, ensure_ascii=False, allow_nan=False)
    write_phase05_summary(args.outdir / "phase05_summary.txt", report, metrics_df, contrib_df, overfit_df)
    return report


def _phase05_report(input_df: pd.DataFrame, validation: Dict[str, Any], output_files: Dict[str, Any], warnings: List[str], severe_flags: List[str], status: str, phase01_report: Dict, phase02_report: Dict, phase03_report: Dict, phase04_report: Dict, args: argparse.Namespace) -> Dict[str, Any]:
    return {
        "phase": "05_simple_additive_variance_partitioning", "status": status, "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_file": str(args.input), "output_directory": str(args.outdir),
        "upstream_phase01_status": phase01_report.get("status"), "upstream_phase02_status": phase02_report.get("status"),
        "upstream_phase03_status": phase03_report.get("status"), "upstream_phase04_status": phase04_report.get("status"),
        "required_columns_checked": validation["required_columns_checked"], "missing_columns": validation["missing_columns"],
        "n_rows_input": int(len(input_df)), "n_capsides": int(input_df["capside"].nunique(dropna=True)) if "capside" in input_df.columns else 0,
        "n_folds": int(input_df["fold"].astype(str).nunique(dropna=True)) if "fold" in input_df.columns else 0,
        "h_values": _sorted_json_values(input_df["h"]) if "h" in input_df.columns else [],
        "patch_R_values": _sorted_json_values(input_df["patch_R"]) if "patch_R" in input_df.columns else [],
        "n_model_conditions": int(input_df[PHASE05_CONDITION_COLS].drop_duplicates().shape[0]) if all(c in input_df.columns for c in PHASE05_CONDITION_COLS) else 0,
        "response_variable": "d_max_pct", "model_condition_cols": PHASE05_CONDITION_COLS,
        "model_definitions": {k: v["label"] for k, v in PHASE05_MODEL_FORMULAS.items()},
        "categorical_encoding": {"coding": "treatment/reference", "reference_level_rule": "lexicographically first category per condition/model"},
        "complete_case_policy": {"core_complete_case": PHASE05_CORE_COLUMNS, "geometry_complete_case": PHASE05_CORE_COLUMNS + PHASE05_GEOMETRY_COLUMNS, "m3_vs_m4": "M3_geometry_base is refit on the same geometry complete-case rows as M4"},
        "loocv_method": "leave-one-out cross-validation", "loocv_unseen_category_policy": "held-out rows with categorical levels absent from training are marked failed and excluded from LOOCV error metrics",
        "overfit_thresholds": {"model_min_n": int(args.model_min_n), "loocv_min_success_fraction": float(args.loocv_min_success_fraction), "condition_number_threshold": float(args.condition_number_threshold), "min_observations_per_parameter": float(args.min_observations_per_parameter)},
        "duplicate_key_status": validation["duplicate_key_status"], "nonfinite_counts": validation["nonfinite_counts"],
        "output_files": output_files, "warnings": warnings, "severe_flags": severe_flags,
    }


def _print_phase05_summary(report: Dict) -> None:
    print(f"[{report['status']}] Phase 5 completed.")
    files = report.get("output_files", {})
    print(f"[INFO] Model metrics table: {files.get('model_metrics', 'not written')}")
    print(f"[INFO] Model contributions table: {files.get('model_contributions', 'not written')}")
    print(f"[INFO] Model predictions table: {files.get('model_predictions', 'not written')}")
    print(f"[INFO] Overfit diagnostics table: {files.get('overfit_diagnostics', 'not written')}")
    print(f"[INFO] Phase 5 report: {files.get('report', 'not written')}")
    metrics_path = Path(files.get("model_metrics", ""))
    contrib_path = Path(files.get("model_contributions", ""))
    if metrics_path.exists() and metrics_path.stat().st_size > 0:
        try:
            metrics = pd.read_csv(metrics_path)
            for (h, r), group in metrics.groupby(["h", "patch_R"], sort=True, dropna=False):
                ok = group[np.isfinite(pd.to_numeric(group["loocv_rmse"], errors="coerce"))]
                if not ok.empty:
                    best = ok.sort_values("loocv_rmse", kind="mergesort").iloc[0]
                    print(f"[MODEL] h={h} patch_R={r} best_LOOCV_model={best['model_id']} loocv_rmse={best['loocv_rmse']:.6g}")
        except Exception:
            pass
    if contrib_path.exists() and contrib_path.stat().st_size > 0:
        try:
            contrib = pd.read_csv(contrib_path)
            geom = contrib[contrib["comparison_id"] == "geometry_added_after_capside_fold"]
            for _, row in geom.iterrows():
                print(f"[GEOM] h={row['h']} patch_R={row['patch_R']} delta_R2_M4_vs_M3={row['delta_r_squared']:.6g} delta_LOOCV_RMSE={row['delta_loocv_rmse']:.6g} label={row['interpretation_label']}")
        except Exception:
            pass
    if report.get("warnings"):
        print("[WARN] Phase 5 completed with warnings. See phase05_report.json.")

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Capsid mechanical-geometrical-topological analysis. Phase 1 reads, "
            "validates, and normalizes the master CSV; Phase 2 performs fold-level "
            "mechanical comparisons within each (capside, h, patch_R) stratum; "
            "Phase 3 evaluates centered local geometry-mechanics associations; Phase 4 compares ordinal fold rankings."
        )
    )

    parser.add_argument(
        "--phase",
        choices=["1", "2", "3", "4", "5"],
        default="1",
        help="Analysis phase to run. Phase 1 is the default. Phase 5 runs Phases 1 through 5, then fits simple additive variance-partitioning models."
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
        help="Output directory for Phase 2, Phase 3, Phase 4, and Phase 5 artifacts."
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
        help="Skip Phase 3 centered scatterplots, Phase 4 ordinal ranking plots, and Phase 5 model comparison figures."
    )

    parser.add_argument(
        "--save-bootstrap-distributions",
        action="store_true",
        help="Save full Phase 3 bootstrap slope distributions."
    )

    parser.add_argument(
        "--rank-method",
        default="average",
        choices=["average", "min", "max", "dense", "first"],
        help="Phase 4 ranking method for ties. Default: average."
    )

    parser.add_argument(
        "--rank-consistency-threshold",
        default=PHASE04_CONSISTENCY_THRESHOLD,
        type=float,
        help="Phase 4 threshold for mostly_positive / mostly_negative consistency labels. Default: 0.75."
    )

    parser.add_argument(
        "--save-permutation-distributions",
        action="store_true",
        help="Save full Phase 3 permutation slope distributions."
    )

    parser.add_argument(
        "--model-min-n",
        default=3,
        type=int,
        help="Phase 5 minimum number of observations required to attempt a model. Default: 3."
    )

    parser.add_argument(
        "--loocv-min-success-fraction",
        default=0.70,
        type=float,
        help="Phase 5 minimum fraction of successful LOOCV predictions before LOOCV is considered stable. Default: 0.70."
    )

    parser.add_argument(
        "--condition-number-threshold",
        default=1000.0,
        type=float,
        help="Phase 5 high condition number warning threshold. Default: 1000."
    )

    parser.add_argument(
        "--min-observations-per-parameter",
        default=3.0,
        type=float,
        help="Phase 5 observations-per-parameter overfit warning threshold. Default: 3."
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

    if not 0 <= args.rank_consistency_threshold <= 1:
        print("[ERROR] --rank-consistency-threshold must be between 0 and 1.", file=sys.stderr)
        return 2

    if args.model_min_n < 2:
        print("[ERROR] --model-min-n must be at least 2.", file=sys.stderr)
        return 2

    if not 0 < args.loocv_min_success_fraction <= 1:
        print("[ERROR] --loocv-min-success-fraction must be between 0 and 1.", file=sys.stderr)
        return 2

    if args.condition_number_threshold <= 0 or args.min_observations_per_parameter <= 0:
        print("[ERROR] Phase 5 threshold options must be positive.", file=sys.stderr)
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
    centered_phase02_df = pd.DataFrame()
    if phase02_report["status"] == "FAIL" or not centered_path.exists():
        phase03_report = run_phase03(pd.DataFrame(), phase01_report, phase02_report, args)
        _print_phase03_summary(phase03_report)
        if args.phase == "3":
            return 1
    else:
        centered_phase02_df = pd.read_csv(centered_path)
        phase03_report = run_phase03(centered_phase02_df, phase01_report, phase02_report, args)
        _print_phase03_summary(phase03_report)

    if args.phase == "3":
        if phase03_report["status"] == "FAIL":
            return 1
        if args.fail_on_warning and phase03_report["status"] != "PASS":
            return 1
        return 0

    phase04_input_df = centered_phase02_df
    phase03_centered_path = args.outdir / "phase03_centered_geometry.csv"
    if phase03_centered_path.exists() and phase03_centered_path.stat().st_size > 0:
        try:
            phase04_input_df = pd.read_csv(phase03_centered_path)
        except Exception:
            phase04_input_df = centered_phase02_df
    phase04_report = run_phase04(phase04_input_df, phase01_report, phase02_report, phase03_report, args)
    _print_phase04_summary(phase04_report)

    if args.phase == "4":
        if phase04_report["status"] == "FAIL":
            return 1
        if args.fail_on_warning and phase04_report["status"] != "PASS":
            return 1
        return 0

    phase05_input_df = phase04_input_df
    if phase05_input_df.empty or not all(c in phase05_input_df.columns for c in ["capside", "fold", "h", "patch_R", "d_max_pct"]):
        phase05_input_df = centered_phase02_df if not centered_phase02_df.empty else normalized_df
    phase05_report = run_phase05(phase05_input_df, phase01_report, phase02_report, phase03_report, phase04_report, args)
    _print_phase05_summary(phase05_report)

    if phase05_report["status"] == "FAIL":
        return 1
    if args.fail_on_warning and phase05_report["status"] != "PASS":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
