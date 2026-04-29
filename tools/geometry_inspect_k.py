#!/usr/bin/env python3
"""
CapDAT geometry K diagnostics — Step 1
File discovery and CSV loading foundation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd


@dataclass(frozen=True)
class GeometryCsvFiles:
    result_dir: Path
    prefix: str

    outer_curvature: Optional[Path] = None
    inner_curvature: Optional[Path] = None

    outer_derivatives: Optional[Path] = None
    inner_derivatives: Optional[Path] = None

    curvature_valid_mask: Optional[Path] = None
    derivative_valid_mask: Optional[Path] = None
    derivative_failure_reason: Optional[Path] = None

    metric_domain_mask: Optional[Path] = None
    smooth_valid_mask: Optional[Path] = None


@dataclass
class LoadedGeometryCsvs:
    files: GeometryCsvFiles
    outer_curvature: Optional[pd.DataFrame] = None
    inner_curvature: Optional[pd.DataFrame] = None
    outer_derivatives: Optional[pd.DataFrame] = None
    inner_derivatives: Optional[pd.DataFrame] = None
    curvature_valid_mask: Optional[pd.DataFrame] = None
    derivative_valid_mask: Optional[pd.DataFrame] = None
    derivative_failure_reason: Optional[pd.DataFrame] = None
    metric_domain_mask: Optional[pd.DataFrame] = None
    smooth_valid_mask: Optional[pd.DataFrame] = None


COLUMN_ALIAS_MAP: dict[str, tuple[str, ...]] = {
    "point_index": ("point_index", "point_idx", "idx", "index", "vertex_index", "vertex_idx"),
    "x": ("x", "coord_x", "px", "pos_x"),
    "y": ("y", "coord_y", "py", "pos_y"),
    "z": ("z", "coord_z", "pz", "pos_z"),
    "k": ("k", "curvature_k", "principal_k", "k_value"),
    "valid": ("valid", "is_valid", "mask", "valid_mask"),
    "reason": ("reason", "failure_reason", "error_reason", "status_reason"),
}

REQUIRED_COLUMNS_BY_LABEL: dict[str, tuple[str, ...]] = {
    "outer curvature": ("point_index", "k"),
    "inner curvature": ("point_index", "k"),
    "outer derivatives": ("point_index",),
    "inner derivatives": ("point_index",),
    "curvature valid mask": ("point_index", "valid"),
    "derivative valid mask": ("point_index", "valid"),
    "derivative failure reason": ("point_index", "reason"),
    "metric domain mask": ("point_index", "valid"),
    "smooth valid mask": ("point_index", "valid"),
}


def _normalize_column_name(value: str) -> str:
    return "".join(ch for ch in value.lower() if ch.isalnum())


def resolve_column_aliases(df: pd.DataFrame) -> dict[str, str]:
    normalized_to_actual = {_normalize_column_name(col): col for col in df.columns}
    resolved: dict[str, str] = {}
    for canonical, aliases in COLUMN_ALIAS_MAP.items():
        for alias in aliases:
            actual = normalized_to_actual.get(_normalize_column_name(alias))
            if actual is not None:
                resolved[canonical] = actual
                break
    return resolved


def validate_required_columns(
    df: pd.DataFrame,
    *,
    label: str,
) -> dict[str, str]:
    required_columns = REQUIRED_COLUMNS_BY_LABEL.get(label, ())
    resolved_aliases = resolve_column_aliases(df)
    missing = [name for name in required_columns if name not in resolved_aliases]

    print(f"[INFO] {label} discovered columns: {list(df.columns)}")

    if resolved_aliases:
        mapping_str = ", ".join(
            f"{canonical}->{actual}" for canonical, actual in sorted(resolved_aliases.items())
        )
        print(f"[INFO] {label} alias resolution: {mapping_str}")
    else:
        print(f"[WARN] {label} alias resolution: no known aliases found")

    if missing:
        raise ValueError(
            f"Missing required columns for {label}: {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    return resolved_aliases


def _glob_one(
    result_dir: Path,
    patterns: list[str],
    *,
    label: str,
    required: bool = False,
) -> Optional[Path]:
    matches: list[Path] = []

    for pattern in patterns:
        matches.extend(sorted(result_dir.glob(pattern)))

    unique_matches = sorted(set(matches))

    if not unique_matches:
        if required:
            raise FileNotFoundError(
                f"Required CSV not found for '{label}'. "
                f"Searched patterns: {patterns}"
            )
        print(f"[WARN] Missing optional CSV: {label}")
        return None

    if len(unique_matches) > 1:
        print(f"[WARN] Multiple matches for '{label}'. Using first:")
        for path in unique_matches:
            print(f"       - {path}")

    selected = unique_matches[0]
    print(f"[INFO] {label}: {selected}")
    return selected


def discover_geometry_csv_files(
    result_dir: Path,
    prefix: str,
    *,
    surface: str = "both",
    strict: bool = False,
) -> GeometryCsvFiles:
    if not result_dir.exists():
        raise FileNotFoundError(f"Result directory does not exist: {result_dir}")

    if not result_dir.is_dir():
        raise NotADirectoryError(f"Not a directory: {result_dir}")

    p = f"{prefix}*" if prefix else "*"

    need_outer = surface in {"outer", "both"}
    need_inner = surface in {"inner", "both"}

    outer_curvature = _glob_one(
        result_dir,
        [
            f"{p}outer*curv*.csv",
            f"{p}outer*curvature*.csv",
        ],
        label="outer curvature CSV",
        required=strict and need_outer,
    )

    inner_curvature = _glob_one(
        result_dir,
        [
            f"{p}inner*curv*.csv",
            f"{p}inner*curvature*.csv",
        ],
        label="inner curvature CSV",
        required=strict and need_inner,
    )

    outer_derivatives = _glob_one(
        result_dir,
        [
            f"{p}outer*deriv*.csv",
            f"{p}outer*derivative*.csv",
        ],
        label="outer derivatives CSV",
        required=strict and need_outer,
    )

    inner_derivatives = _glob_one(
        result_dir,
        [
            f"{p}inner*deriv*.csv",
            f"{p}inner*derivative*.csv",
        ],
        label="inner derivatives CSV",
        required=strict and need_inner,
    )

    curvature_valid_mask = _glob_one(
        result_dir,
        [
            f"{p}curvature*valid*mask*.csv",
            f"{p}curv*valid*mask*.csv",
        ],
        label="curvature valid mask CSV",
        required=False,
    )

    derivative_valid_mask = _glob_one(
        result_dir,
        [
            f"{p}derivative*valid*mask*.csv",
            f"{p}deriv*valid*mask*.csv",
        ],
        label="derivative valid mask CSV",
        required=False,
    )

    derivative_failure_reason = _glob_one(
        result_dir,
        [
            f"{p}derivative*failure*.csv",
            f"{p}deriv*failure*.csv",
            f"{p}failure*reason*.csv",
        ],
        label="derivative failure reason CSV",
        required=False,
    )

    metric_domain_mask = _glob_one(
        result_dir,
        [
            f"{p}metric*domain*mask*.csv",
        ],
        label="metric domain mask CSV",
        required=False,
    )

    smooth_valid_mask = _glob_one(
        result_dir,
        [
            f"{p}smooth*valid*mask*.csv",
        ],
        label="smooth valid mask CSV",
        required=False,
    )

    return GeometryCsvFiles(
        result_dir=result_dir,
        prefix=prefix,
        outer_curvature=outer_curvature,
        inner_curvature=inner_curvature,
        outer_derivatives=outer_derivatives,
        inner_derivatives=inner_derivatives,
        curvature_valid_mask=curvature_valid_mask,
        derivative_valid_mask=derivative_valid_mask,
        derivative_failure_reason=derivative_failure_reason,
        metric_domain_mask=metric_domain_mask,
        smooth_valid_mask=smooth_valid_mask,
    )


def read_csv_checked(path: Optional[Path], *, label: str) -> Optional[pd.DataFrame]:
    if path is None:
        return None

    try:
        df = pd.read_csv(path)
    except Exception as exc:
        raise RuntimeError(f"Failed to read {label}: {path}") from exc

    if df.empty:
        print(f"[WARN] CSV is empty: {label}: {path}")
    else:
        print(f"[INFO] Loaded {label}: {len(df)} rows, {len(df.columns)} columns")
        validate_required_columns(df, label=label)

    return df


def load_geometry_csvs(files: GeometryCsvFiles) -> LoadedGeometryCsvs:
    return LoadedGeometryCsvs(
        files=files,
        outer_curvature=read_csv_checked(files.outer_curvature, label="outer curvature"),
        inner_curvature=read_csv_checked(files.inner_curvature, label="inner curvature"),
        outer_derivatives=read_csv_checked(files.outer_derivatives, label="outer derivatives"),
        inner_derivatives=read_csv_checked(files.inner_derivatives, label="inner derivatives"),
        curvature_valid_mask=read_csv_checked(files.curvature_valid_mask, label="curvature valid mask"),
        derivative_valid_mask=read_csv_checked(files.derivative_valid_mask, label="derivative valid mask"),
        derivative_failure_reason=read_csv_checked(files.derivative_failure_reason, label="derivative failure reason"),
        metric_domain_mask=read_csv_checked(files.metric_domain_mask, label="metric domain mask"),
        smooth_valid_mask=read_csv_checked(files.smooth_valid_mask, label="smooth valid mask"),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="CapDAT geometry diagnostic CSV discovery and loading."
    )

    parser.add_argument(
        "--result-dir",
        required=True,
        type=Path,
        help="Directory containing CapDAT geometry CSV outputs.",
    )

    parser.add_argument(
        "--prefix",
        default="",
        help="Artifact prefix used by CapDAT geometry outputs.",
    )

    parser.add_argument(
        "--surface",
        choices=["outer", "inner", "both"],
        default="both",
        help="Surface to inspect.",
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if required surface curvature/derivative CSVs are missing.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    files = discover_geometry_csv_files(
        result_dir=args.result_dir,
        prefix=args.prefix,
        surface=args.surface,
        strict=args.strict,
    )

    _ = load_geometry_csvs(files)

    print("[SUCCESS] File discovery and CSV loading completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
