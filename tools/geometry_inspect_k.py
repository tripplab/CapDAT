#!/usr/bin/env python3
"""
CapDAT geometry K diagnostics.

Includes CSV discovery/loading + merged-surface validation and basic
Gaussian-curvature heatmap generation for --surface outer|inner|both.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scipy.stats import spearmanr as scipy_spearmanr
except Exception:  # pragma: no cover - optional dependency
    scipy_spearmanr = None


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
    curvature_summary: Optional[Path] = None


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
    curvature_summary: Optional[pd.DataFrame] = None


COLUMN_ALIAS_MAP: dict[str, tuple[str, ...]] = {
    "point_index": ("point_index", "point_idx", "idx", "index", "vertex_index", "vertex_idx"),
    "i": ("i", "ix", "grid_i"),
    "j": ("j", "iy", "grid_j"),
    "x": ("x", "coord_x", "px", "pos_x"),
    "y": ("y", "coord_y", "py", "pos_y"),
    "z": ("z", "coord_z", "pz", "pos_z"),
    "k": ("k", "curvature_k", "principal_k", "k_value", "k_raw"),
    "valid": ("valid", "is_valid", "mask", "valid_mask"),
    "reason": ("reason", "failure_reason", "error_reason", "status_reason"),
    "derivative_valid": ("derivative_valid",),
    "curvature_valid": ("curvature_valid",),
    "metric_domain": ("metric_domain",),
    "fit_attempted": ("fit_attempted", "derivative_fit_attempted"),
    "thickness_valid": ("thickness_valid",),
}

REQUIRED_COLUMNS_BY_LABEL: dict[str, tuple[str, ...]] = {
    "outer curvature": ("k",),
    "inner curvature": ("k",),
    "outer derivatives": (),
    "inner derivatives": (),
    "curvature valid mask": (),
    "derivative valid mask": (),
    "derivative failure reason": (),
    "metric domain mask": (),
    "smooth valid mask": (),
    "curvature summary": (),
}

ALTERNATE_REQUIRED_COLUMNS_BY_LABEL: dict[str, tuple[tuple[str, ...], ...]] = {
    "outer curvature": (("point_index",), ("i", "j")),
    "inner curvature": (("point_index",), ("i", "j")),
    "outer derivatives": (("point_index",), ("i", "j")),
    "inner derivatives": (("point_index",), ("i", "j")),
    "curvature valid mask": (("point_index", "valid"), ("i", "j", "curvature_valid")),
    "derivative valid mask": (("point_index", "valid"), ("i", "j", "derivative_valid")),
    "derivative failure reason": (("point_index", "reason"), ("i", "j", "fit_attempted")),
    "metric domain mask": (("point_index", "valid"), ("i", "j", "metric_domain")),
    "smooth valid mask": (("point_index", "valid"), ("i", "j", "thickness_valid")),
    "curvature summary": (),
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
    alternate_groups = ALTERNATE_REQUIRED_COLUMNS_BY_LABEL.get(label, ())
    resolved_aliases = resolve_column_aliases(df)
    missing = [name for name in required_columns if name not in resolved_aliases]

    if missing:
        raise ValueError(
            f"Missing required columns for {label}: {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    if alternate_groups:
        for group in alternate_groups:
            missing_group = [name for name in group if name not in resolved_aliases]
            if not missing_group:
                break
        else:
            expected = " OR ".join(str(list(group)) for group in alternate_groups)
            raise ValueError(
                f"Missing required alternative columns for {label}. "
                f"Expected one of: {expected}. "
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

    curvature_summary = _glob_one(
        result_dir,
        [
            f"{p}curvature*summary*.csv",
        ],
        label="curvature summary CSV",
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
        curvature_summary=curvature_summary,
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
        curvature_summary=read_csv_checked(files.curvature_summary, label="curvature summary"),
    )


def _build_rename_map(df: pd.DataFrame) -> dict[str, str]:
    aliases = resolve_column_aliases(df)
    rename_map: dict[str, str] = {}
    for canonical, actual in aliases.items():
        if actual != canonical and canonical not in df.columns:
            rename_map[actual] = canonical
    return rename_map


def _standardize_columns(df: pd.DataFrame, *, label: str) -> pd.DataFrame:
    rename_map = _build_rename_map(df)
    if rename_map:
        return df.rename(columns=rename_map).copy()
    return df.copy()


def merge_surface_geometry(
    curvature_df: Optional[pd.DataFrame],
    derivative_df: Optional[pd.DataFrame],
    *,
    surface: str,
) -> Optional[pd.DataFrame]:
    if curvature_df is None or derivative_df is None:
        print(f"[WARN] Skipping {surface} merge; missing curvature or derivatives table")
        return None

    cdf = _standardize_columns(curvature_df, label=f"{surface} curvature")
    ddf = _standardize_columns(derivative_df, label=f"{surface} derivatives")

    merge_keys = ["i", "j"]
    missing_curv = [k for k in merge_keys if k not in cdf.columns]
    missing_deriv = [k for k in merge_keys if k not in ddf.columns]
    if missing_curv or missing_deriv:
        raise ValueError(
            f"Cannot merge {surface} tables by i,j. "
            f"Missing in curvature: {missing_curv}; missing in derivatives: {missing_deriv}"
        )

    merged = cdf.merge(ddf, on=merge_keys, how="inner", suffixes=("_curv", "_deriv"), validate="one_to_one")

    # Merge inputs share several coordinate/domain columns. Normalize those
    # duplicated columns back to canonical names so downstream checks can use a
    # stable schema.
    for base_col in ("x", "y", "metric_domain", "derivative_valid", "curvature_valid"):
        curv_col = f"{base_col}_curv"
        deriv_col = f"{base_col}_deriv"
        if base_col in merged.columns:
            continue
        if curv_col in merged.columns and deriv_col in merged.columns:
            if not merged[curv_col].equals(merged[deriv_col]):
                raise ValueError(
                    f"{surface} merge mismatch for duplicated column '{base_col}' "
                    f"between curvature and derivatives tables"
                )
            merged[base_col] = merged[curv_col]
        elif curv_col in merged.columns:
            merged[base_col] = merged[curv_col]
        elif deriv_col in merged.columns:
            merged[base_col] = merged[deriv_col]

    expected_rows = min(len(cdf), len(ddf))
    actual_rows = len(merged)
    print(f"[INFO] {surface} merge row counts: curvature={len(cdf)}, derivatives={len(ddf)}, merged={actual_rows}")
    if actual_rows != expected_rows:
        raise ValueError(
            f"Row-count validation failed for {surface} merge by i,j: "
            f"expected {expected_rows} rows (min of inputs), got {actual_rows}"
        )

    for required_col in ("dz_dx", "dz_dy", "x", "y", "k"):
        if required_col not in merged.columns:
            raise ValueError(f"Missing required column '{required_col}' after {surface} merge")

    merged["K_raw"] = merged["k"]
    merged["slope_mag"] = np.sqrt(merged["dz_dx"] ** 2 + merged["dz_dy"] ** 2)
    merged["J"] = np.sqrt(1.0 + merged["dz_dx"] ** 2 + merged["dz_dy"] ** 2)
    merged["abs_K"] = merged["K_raw"].abs()
    merged["r"] = np.sqrt(merged["x"] ** 2 + merged["y"] ** 2)

    if "curvature_valid" not in merged.columns:
        merged["curvature_valid"] = True
    if "K_qc_warn_flag" not in merged.columns:
        merged["K_qc_warn_flag"] = 0

    merged["K_qc_clean"] = merged["curvature_valid"].astype(bool) & (merged["K_qc_warn_flag"] == 0)
    print(f"[INFO] {surface} derived columns computed: slope_mag, J, abs_K, r, K_qc_clean")
    return merged


def dataframe_to_grid(df: pd.DataFrame, value_col: str) -> tuple[np.ndarray, list[float]]:
    if value_col not in df.columns:
        raise ValueError(f"Missing value column for grid conversion: {value_col}")

    nx = int(df["i"].max()) + 1
    ny = int(df["j"].max()) + 1
    grid = np.full((ny, nx), np.nan, dtype=float)

    i_idx = df["i"].to_numpy(dtype=int)
    j_idx = df["j"].to_numpy(dtype=int)
    values = df[value_col].to_numpy(dtype=float)
    grid[j_idx, i_idx] = values

    extent = [float(df["x"].min()), float(df["x"].max()), float(df["y"].min()), float(df["y"].max())]
    return grid, extent


def _prepare_k_plot_df(
    surface_df: pd.DataFrame,
    surface_name: str,
) -> pd.DataFrame:
    required = ["i", "j", "x", "y", "k", "curvature_valid"]
    missing = [col for col in required if col not in surface_df.columns]
    if missing:
        raise ValueError(f"Missing required columns for {surface_name} plot: {missing}")

    plot_df = surface_df.copy()
    plot_df["k_plot"] = plot_df["k"].astype(float)
    invalid_mask = (plot_df["curvature_valid"].astype(int) == 0) | (~np.isfinite(plot_df["k_plot"]))
    plot_df.loc[invalid_mask, "k_plot"] = np.nan

    return plot_df


def _compute_vlim(plot_df: pd.DataFrame, k_percentile: float) -> Optional[float]:
    k_valid = plot_df["k_plot"].to_numpy(dtype=float)
    k_valid = k_valid[np.isfinite(k_valid)]
    if k_valid.size == 0:
        return None

    vlim = float(np.nanpercentile(np.abs(k_valid), k_percentile))
    if not np.isfinite(vlim) or vlim <= 0:
        vlim = float(np.nanmax(np.abs(k_valid)))
    if not np.isfinite(vlim) or vlim <= 0:
        return None
    return vlim


def _compute_shared_vlim(
    surface_dfs: list[pd.DataFrame],
    surface_names: list[str],
    k_percentile: float,
) -> Optional[float]:
    vlms: list[float] = []
    for df, name in zip(surface_dfs, surface_names):
        plot_df = _prepare_k_plot_df(df, name)
        vlim = _compute_vlim(plot_df, k_percentile)
        if vlim is not None:
            vlms.append(vlim)
    if not vlms:
        return None
    return float(max(vlms))


def plot_k_heatmap_base(
    ax: plt.Axes,
    surface_df: pd.DataFrame,
    surface_name: str,
    k_percentile: float,
    *,
    vlim: Optional[float] = None,
) -> tuple[plt.AxesImage, pd.DataFrame]:
    plot_df = _prepare_k_plot_df(surface_df, surface_name)
    resolved_vlim = vlim if vlim is not None else _compute_vlim(plot_df, k_percentile)
    if resolved_vlim is None:
        raise ValueError(f"No valid finite K values for {surface_name}; cannot render base heatmap")

    grid, extent = dataframe_to_grid(plot_df, "k_plot")
    im = ax.imshow(
        grid,
        origin="lower",
        interpolation="nearest",
        aspect="equal",
        cmap="coolwarm_r",
        vmin=-resolved_vlim,
        vmax=resolved_vlim,
        extent=extent,
    )
    return im, plot_df


def plot_k_heatmap(
    surface_df: pd.DataFrame,
    surface_name: str,
    out_dir: Path,
    k_percentile: float,
) -> None:
    vlim = _compute_shared_vlim([surface_df], [surface_name], k_percentile)
    if vlim is None:
        print(f"[WARN] No valid finite K values for {surface_name}; skipping plot")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    im, _ = plot_k_heatmap_base(ax, surface_df, surface_name, k_percentile, vlim=vlim)

    title_surface = "Outer" if surface_name == "outer" else "Inner"
    ax.set_title(f"{title_surface} surface: Gaussian curvature K, all valid nodes")
    ax.set_xlabel("x [Å]")
    ax.set_ylabel("y [Å]")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("K [Å^-2]")

    out_path = out_dir / f"{surface_name}_K_heatmap_all_valid.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def plot_k_heatmap_both(
    outer_df: pd.DataFrame,
    inner_df: pd.DataFrame,
    out_dir: Path,
    k_percentile: float,
) -> None:
    vlim = _compute_shared_vlim([outer_df, inner_df], ["outer", "inner"], k_percentile)
    if vlim is None:
        print("[WARN] Missing valid K values for outer/inner; skipping combined plot")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    panels = [
        (axes[0], outer_df, "outer", "Outer surface: Gaussian curvature K, all valid nodes"),
        (axes[1], inner_df, "inner", "Inner surface: Gaussian curvature K, all valid nodes"),
    ]
    for ax, surface_df, surface_name, title in panels:
        im, _ = plot_k_heatmap_base(ax, surface_df, surface_name, k_percentile, vlim=vlim)
        ax.set_title(title)
        ax.set_xlabel("x [Å]")
        ax.set_ylabel("y [Å]")

    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.92)
    cbar.set_label("K [Å^-2]")

    out_path = out_dir / "K_heatmap_all_valid.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def _plot_qc_overlay_points(ax: plt.Axes, plot_df: pd.DataFrame) -> tuple[int, int]:
    qc_warn = (plot_df["curvature_valid"].astype(int) == 1) & (plot_df["K_qc_warn_flag"].astype(int) != 0)
    qc_df = plot_df.loc[qc_warn]
    valid_nodes = int((plot_df["curvature_valid"].astype(int) == 1).sum())
    warn_nodes = int(len(qc_df))
    if warn_nodes > 0:
        ax.scatter(
            qc_df["x"],
            qc_df["y"],
            marker="s",
            color="black",
            s=8,
            alpha=0.65,
            linewidths=0,
            label="K QC warn",
        )
        ax.legend(loc="upper right", fontsize=8)
    return warn_nodes, valid_nodes


def plot_k_heatmap_with_qc_overlay(
    surface_df: pd.DataFrame,
    surface_name: str,
    out_dir: Path,
    k_percentile: float,
    *,
    strict: bool,
    vlim: Optional[float] = None,
) -> None:
    required = ["i", "j", "x", "y", "k", "curvature_valid"]
    missing = [col for col in required if col not in surface_df.columns]
    if missing:
        raise ValueError(f"Missing required columns for {surface_name} QC overlay plot: {missing}")
    if "K_qc_warn_flag" not in surface_df.columns:
        msg = f"[WARN] Missing K_qc_warn_flag for {surface_name}; skipping QC overlay plot"
        if strict:
            raise ValueError(msg)
        print(msg)
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    im, plot_df = plot_k_heatmap_base(ax, surface_df, surface_name, k_percentile, vlim=vlim)
    warn_nodes, valid_nodes = _plot_qc_overlay_points(ax, plot_df)

    title_surface = "Outer" if surface_name == "outer" else "Inner"
    ax.set_title(f"{title_surface} surface: K heatmap with QC warnings")
    ax.set_xlabel("x [Å]")
    ax.set_ylabel("y [Å]")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("K [Å^-2]")

    out_path = out_dir / f"{surface_name}_K_heatmap_qc_overlay.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    fraction = (warn_nodes / valid_nodes) if valid_nodes > 0 else 0.0
    print(f"[INFO] Saved: {out_path}")
    print(
        f"[INFO] {surface_name} QC overlay: "
        f"qc_warn_nodes={warn_nodes}, curvature_valid_nodes={valid_nodes}, fraction={fraction:.2e}"
    )


def plot_k_heatmap_with_qc_overlay_both(
    outer_df: pd.DataFrame,
    inner_df: pd.DataFrame,
    out_dir: Path,
    k_percentile: float,
    *,
    strict: bool,
) -> None:
    if "K_qc_warn_flag" not in outer_df.columns or "K_qc_warn_flag" not in inner_df.columns:
        missing_surfaces = [
            name for name, df in (("outer", outer_df), ("inner", inner_df)) if "K_qc_warn_flag" not in df.columns
        ]
        msg = (
            "[WARN] Missing K_qc_warn_flag for "
            + ",".join(missing_surfaces)
            + "; skipping combined QC overlay plot"
        )
        if strict:
            raise ValueError(msg)
        print(msg)
        return

    vlim = _compute_shared_vlim([outer_df, inner_df], ["outer", "inner"], k_percentile)
    if vlim is None:
        print("[WARN] Missing valid K values for outer/inner; skipping combined QC overlay plot")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, surface_df, surface_name, title in (
        (axes[0], outer_df, "outer", "Outer surface: K heatmap with QC warnings"),
        (axes[1], inner_df, "inner", "Inner surface: K heatmap with QC warnings"),
    ):
        im, plot_df = plot_k_heatmap_base(ax, surface_df, surface_name, k_percentile, vlim=vlim)
        warn_nodes, valid_nodes = _plot_qc_overlay_points(ax, plot_df)
        fraction = (warn_nodes / valid_nodes) if valid_nodes > 0 else 0.0
        print(
            f"[INFO] {surface_name} QC overlay: "
            f"qc_warn_nodes={warn_nodes}, curvature_valid_nodes={valid_nodes}, fraction={fraction:.2e}"
        )
        ax.set_title(title)
        ax.set_xlabel("x [Å]")
        ax.set_ylabel("y [Å]")

    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.92)
    cbar.set_label("K [Å^-2]")
    out_path = out_dir / "K_heatmap_qc_overlay.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def plot_condition_indicator_heatmap(
    surface_df: pd.DataFrame,
    surface_name: str,
    out_dir: Path,
) -> None:
    required = ["i", "j", "x", "y", "fit_condition_indicator", "derivative_valid"]
    missing = [col for col in required if col not in surface_df.columns]
    if missing:
        raise ValueError(f"Missing required columns for {surface_name} condition plot: {missing}")

    plot_df = surface_df.copy()
    plot_df["fit_condition_indicator"] = plot_df["fit_condition_indicator"].astype(float)
    plot_df["derivative_valid"] = plot_df["derivative_valid"].astype(int)

    invalid_mask = (
        (plot_df["derivative_valid"] == 0)
        | (~np.isfinite(plot_df["fit_condition_indicator"]))
        | (plot_df["fit_condition_indicator"] <= 0.0)
    )
    plot_df.loc[invalid_mask, "fit_condition_indicator"] = np.nan

    valid_condition = plot_df["fit_condition_indicator"].to_numpy(dtype=float)
    finite_condition = valid_condition[np.isfinite(valid_condition)]
    derivative_valid_count = int((plot_df["derivative_valid"] == 1).sum())
    finite_count = int(finite_condition.size)
    if finite_count == 0:
        print(f"[WARN] No finite positive condition indicator values for {surface_name}; skipping plot")
        return

    median_val = float(np.nanmedian(finite_condition))
    p95_val = float(np.nanpercentile(finite_condition, 95))
    p99_val = float(np.nanpercentile(finite_condition, 99))
    max_val = float(np.nanmax(finite_condition))
    print(
        f"[INFO] {surface_name} condition indicator: "
        f"derivative_valid={derivative_valid_count}, finite={finite_count}, "
        f"median={median_val:.6g}, p95={p95_val:.6g}, p99={p99_val:.6g}, max={max_val:.6g}"
    )

    plot_df["log10_condition_indicator"] = np.nan
    finite_mask = np.isfinite(plot_df["fit_condition_indicator"].to_numpy(dtype=float))
    plot_df.loc[finite_mask, "log10_condition_indicator"] = np.log10(
        plot_df.loc[finite_mask, "fit_condition_indicator"]
    )

    log_values = plot_df["log10_condition_indicator"].to_numpy(dtype=float)
    log_values = log_values[np.isfinite(log_values)]
    vmin = float(np.nanpercentile(log_values, 1))
    vmax = float(np.nanpercentile(log_values, 99))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin = float(np.nanmin(log_values))
        vmax = float(np.nanmax(log_values))

    grid, extent = dataframe_to_grid(plot_df, "log10_condition_indicator")
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        grid,
        origin="lower",
        interpolation="nearest",
        aspect="equal",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
        extent=extent,
    )
    title_surface = "Outer" if surface_name == "outer" else "Inner"
    ax.set_title(f"{title_surface} surface: Stage 8 fit condition indicator")
    ax.set_xlabel("x [Å]")
    ax.set_ylabel("y [Å]")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log10(condition indicator)")

    out_path = out_dir / f"{surface_name}_S8_condition_indicator_heatmap.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def _prepare_residual_plot_df(surface_df: pd.DataFrame, surface_name: str, residual_col: str) -> pd.DataFrame:
    required = ["i", "j", "x", "y", "derivative_valid", residual_col]
    missing = [col for col in required if col not in surface_df.columns]
    if missing:
        raise ValueError(f"Missing required columns for {surface_name} residual plot: {missing}")
    plot_df = surface_df.copy()
    plot_df["derivative_valid"] = plot_df["derivative_valid"].astype(int)
    plot_df[residual_col] = plot_df[residual_col].astype(float)
    invalid_mask = (
        (plot_df["derivative_valid"] == 0)
        | (~np.isfinite(plot_df[residual_col]))
        | (plot_df[residual_col] < 0.0)
    )
    plot_df.loc[invalid_mask, residual_col] = np.nan
    return plot_df


def _compute_residual_vmax(values: np.ndarray, residual_percentile: float) -> Optional[float]:
    if values.size == 0:
        return None
    vmax = float(np.nanpercentile(values, residual_percentile))
    if (not np.isfinite(vmax)) or (vmax <= 0.0):
        vmax = float(np.nanmax(values))
    if (not np.isfinite(vmax)) or (vmax <= 0.0):
        return None
    return vmax


def plot_residual_heatmap(
    surface_df: pd.DataFrame,
    surface_name: str,
    residual_col: str,
    residual_label: str,
    out_dir: Path,
    residual_percentile: float,
    overlay_qc: bool = False,
) -> None:
    plot_df = _prepare_residual_plot_df(surface_df, surface_name, residual_col)
    values = plot_df[residual_col].to_numpy(dtype=float)
    finite_values = values[np.isfinite(values)]
    derivative_valid_count = int((plot_df["derivative_valid"] == 1).sum())
    finite_count = int(finite_values.size)
    if finite_count == 0:
        print(f"[WARN] No finite residual values for {surface_name} {residual_col}; skipping plot")
        return
    mean_val = float(np.nanmean(finite_values))
    median_val = float(np.nanmedian(finite_values))
    p95_val = float(np.nanpercentile(finite_values, 95))
    p99_val = float(np.nanpercentile(finite_values, 99))
    max_val = float(np.nanmax(finite_values))
    print(
        f"[INFO] {surface_name} {residual_col}: derivative_valid={derivative_valid_count}, "
        f"finite={finite_count}, mean={mean_val:.6g}, median={median_val:.6g}, "
        f"p95={p95_val:.6g}, p99={p99_val:.6g}, max={max_val:.6g}"
    )
    vmax = _compute_residual_vmax(finite_values, residual_percentile)
    if vmax is None:
        print(f"[WARN] Invalid residual scaling for {surface_name} {residual_col}; skipping plot")
        return
    grid, extent = dataframe_to_grid(plot_df, residual_col)
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        grid, origin="lower", interpolation="nearest", aspect="equal",
        cmap="magma", vmin=0.0, vmax=vmax, extent=extent
    )
    title_surface = "Outer" if surface_name == "outer" else "Inner"
    title_suffix = " + K QC warnings" if overlay_qc else ""
    ax.set_title(f"{title_surface} surface: Stage 8 {residual_col.replace('_', ' ')}{title_suffix}")
    ax.set_xlabel("x [Å]")
    ax.set_ylabel("y [Å]")
    if overlay_qc and "K_qc_warn_flag" in plot_df.columns and "curvature_valid" in plot_df.columns:
        qc_df = plot_df.loc[(plot_df["curvature_valid"].astype(int) == 1) & (plot_df["K_qc_warn_flag"].astype(int) != 0)]
        if len(qc_df) > 0:
            ax.scatter(qc_df["x"], qc_df["y"], marker="s", color="black", s=8, alpha=0.65, linewidths=0)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(residual_label)
    suffix = "_qc_overlay" if overlay_qc else ""
    out_path = out_dir / f"{surface_name}_{residual_col}_heatmap{suffix}.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def plot_residual_heatmaps_combined(
    surface_frames: list[tuple[str, pd.DataFrame]],
    out_dir: Path,
    residual_percentile: float,
    overlay_qc: bool,
) -> None:
    panels: list[tuple[str, str, str, pd.DataFrame, float]] = []
    for surface_name, surface_df in surface_frames:
        for residual_col, residual_label in (
            ("fit_rms_residual", "fit RMS residual [Å]"),
            ("fit_max_abs_residual", "fit max abs residual [Å]"),
        ):
            plot_df = _prepare_residual_plot_df(surface_df, surface_name, residual_col)
            values = plot_df[residual_col].to_numpy(dtype=float)
            finite_values = values[np.isfinite(values)]
            if finite_values.size == 0:
                print(f"[WARN] No finite residual values for {surface_name} {residual_col}; skipping panel")
                continue
            vmax = _compute_residual_vmax(finite_values, residual_percentile)
            if vmax is None:
                print(f"[WARN] Invalid residual scaling for {surface_name} {residual_col}; skipping panel")
                continue
            panels.append((surface_name, residual_col, residual_label, plot_df, vmax))
    if not panels:
        return
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    axes_list = axes.ravel()
    for idx, (surface_name, residual_col, residual_label, plot_df, vmax) in enumerate(panels[:4]):
        ax = axes_list[idx]
        grid, extent = dataframe_to_grid(plot_df, residual_col)
        im = ax.imshow(grid, origin="lower", interpolation="nearest", aspect="equal", cmap="magma", vmin=0.0, vmax=vmax, extent=extent)
        title_surface = "Outer" if surface_name == "outer" else "Inner"
        title_suffix = " + QC warnings" if overlay_qc else ""
        ax.set_title(f"{title_surface} surface: Stage 8 {residual_col.replace('_', ' ')}{title_suffix}")
        ax.set_xlabel("x [Å]")
        ax.set_ylabel("y [Å]")
        if overlay_qc and "K_qc_warn_flag" in plot_df.columns and "curvature_valid" in plot_df.columns:
            qc_df = plot_df.loc[(plot_df["curvature_valid"].astype(int) == 1) & (plot_df["K_qc_warn_flag"].astype(int) != 0)]
            if len(qc_df) > 0:
                ax.scatter(qc_df["x"], qc_df["y"], marker="s", color="black", s=8, alpha=0.65, linewidths=0)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(residual_label)
    for idx in range(len(panels), 4):
        axes_list[idx].axis("off")
    out_name = "fit_rms_max_abs_residual_qc_overlay.png" if overlay_qc else "fit_rms_max_abs_residual_heatmap.png"
    out_path = out_dir / out_name
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def print_residual_summary_stats(surface_df: pd.DataFrame, surface_name: str) -> None:
    for residual_col in ("fit_rms_residual", "fit_max_abs_residual"):
        plot_df = _prepare_residual_plot_df(surface_df, surface_name, residual_col)
        values = plot_df[residual_col].to_numpy(dtype=float)
        finite_values = values[np.isfinite(values)]
        derivative_valid_count = int((plot_df["derivative_valid"] == 1).sum())
        finite_count = int(finite_values.size)
        if finite_count == 0:
            print(f"[WARN] {surface_name} {residual_col}: no finite derivative-valid residual values")
            continue
        print(
            f"[INFO] {surface_name} {residual_col}: derivative_valid={derivative_valid_count}, "
            f"finite={finite_count}, mean={np.nanmean(finite_values):.6g}, median={np.nanmedian(finite_values):.6g}, "
            f"p95={np.nanpercentile(finite_values,95):.6g}, p99={np.nanpercentile(finite_values,99):.6g}, "
            f"max={np.nanmax(finite_values):.6g}"
        )
        if "curvature_valid" in plot_df.columns and "K_qc_warn_flag" in plot_df.columns:
            qc_mask = (plot_df["curvature_valid"].astype(int) == 1) & (plot_df["K_qc_warn_flag"].astype(int) != 0)
            qc_values = plot_df.loc[qc_mask, residual_col].to_numpy(dtype=float)
            qc_values = qc_values[np.isfinite(qc_values)]
            median_all = float(np.nanmedian(finite_values))
            median_qc = float(np.nanmedian(qc_values)) if qc_values.size > 0 else np.nan
            ratio = (median_qc / median_all) if (np.isfinite(median_qc) and median_all > 0.0) else np.nan
            short = "RMS" if residual_col == "fit_rms_residual" else "max-abs"
            print(
                f"[INFO] {surface_name} {short} residual QC enrichment: "
                f"median_all={median_all:.6g}, median_qc={median_qc:.6g}, ratio={ratio:.6g}"
            )


def plot_neighbor_heatmap(
    surface_df: pd.DataFrame,
    surface_name: str,
    value_col: str,
    label: str,
    out_dir: Path,
    percentile: float,
    overlay_qc: bool = False,
) -> None:
    required = ["i", "j", "x", "y", "derivative_valid", value_col]
    missing = [col for col in required if col not in surface_df.columns]
    if missing:
        print(f"[WARN] Missing required columns for {surface_name} {value_col}: {missing}; skipping plot")
        return
    plot_df = surface_df.copy()
    plot_df["derivative_valid"] = plot_df["derivative_valid"].astype(int)
    plot_df[value_col] = pd.to_numeric(plot_df[value_col], errors="coerce")
    invalid_mask = (
        (plot_df["derivative_valid"] != 1)
        | (~np.isfinite(plot_df[value_col]))
        | (plot_df[value_col] <= 0.0)
    )
    plot_df.loc[invalid_mask, value_col] = np.nan
    values = plot_df[value_col].to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        print(f"[WARN] No finite derivative-valid {value_col} values for {surface_name}; skipping plot")
        return
    vmin = float(np.nanpercentile(finite, 1.0))
    vmax = float(np.nanpercentile(finite, percentile))
    if (not np.isfinite(vmin)) or (not np.isfinite(vmax)) or (vmax <= vmin):
        vmin = float(np.nanmin(finite))
        vmax = float(np.nanmax(finite))
    if (not np.isfinite(vmin)) or (not np.isfinite(vmax)) or (vmax <= vmin):
        const_value = float(np.nanmedian(finite))
        if value_col == "neighbor_max_radius":
            print(f"[INFO] {surface_name} {value_col} is constant (={const_value:.6g} Å); heatmap skipped.")
        else:
            print(f"[INFO] {surface_name} {value_col} is constant (={const_value:.6g}); heatmap skipped.")
        return
    grid, extent = dataframe_to_grid(plot_df, value_col)
    fig, ax = plt.subplots(figsize=(6, 5))
    cmap = "viridis" if value_col == "neighbor_count" else "magma"
    im = ax.imshow(grid, origin="lower", interpolation="nearest", aspect="equal", cmap=cmap, vmin=vmin, vmax=vmax, extent=extent)
    title_surface = "Outer" if surface_name == "outer" else "Inner"
    if value_col == "neighbor_count":
        base_title = f"{title_surface} surface: neighbor count (Stage 8)"
    else:
        base_title = f"{title_surface} surface: neighbor max radius (Å)"
    ax.set_title(base_title + (" + K QC warnings" if overlay_qc else ""))
    ax.set_xlabel("x [Å]")
    ax.set_ylabel("y [Å]")
    if overlay_qc and "K_qc_warn_flag" in plot_df.columns and "curvature_valid" in plot_df.columns:
        qc_df = plot_df.loc[(plot_df["curvature_valid"].astype(int) == 1) & (plot_df["K_qc_warn_flag"].astype(int) != 0)]
        if len(qc_df) > 0:
            ax.scatter(qc_df["x"], qc_df["y"], marker="s", color="black", s=8, alpha=0.6, linewidths=0)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)
    suffix = "_qc_overlay" if overlay_qc else ""
    stem = "neighbor_count" if value_col == "neighbor_count" else "neighbor_radius"
    out_path = out_dir / f"{surface_name}_{stem}{suffix}.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def plot_scale_radial_profiles(surface_df: pd.DataFrame, surface_name: str, out_dir: Path) -> None:
    required = ["x", "y", "neighbor_count", "neighbor_max_radius", "derivative_valid", "curvature_valid", "K_qc_warn_flag"]
    missing = [col for col in required if col not in surface_df.columns]
    if missing:
        print(f"[WARN] Missing required columns for {surface_name} radial profiles: {missing}; skipping")
        return
    df = surface_df.copy()
    df["derivative_valid"] = df["derivative_valid"].astype(int)
    df["curvature_valid"] = df["curvature_valid"].astype(int)
    df["K_qc_warn_flag"] = df["K_qc_warn_flag"].astype(int)
    df["neighbor_count"] = pd.to_numeric(df["neighbor_count"], errors="coerce")
    df["neighbor_max_radius"] = pd.to_numeric(df["neighbor_max_radius"], errors="coerce")
    df.loc[(df["derivative_valid"] != 1) | (~np.isfinite(df["neighbor_count"])) | (df["neighbor_count"] <= 0.0), "neighbor_count"] = np.nan
    df.loc[(df["derivative_valid"] != 1) | (~np.isfinite(df["neighbor_max_radius"])) | (df["neighbor_max_radius"] <= 0.0), "neighbor_max_radius"] = np.nan
    df["r"] = np.sqrt(np.square(pd.to_numeric(df["x"], errors="coerce")) + np.square(pd.to_numeric(df["y"], errors="coerce")))
    df = df[np.isfinite(df["r"])].copy()
    if df.empty:
        print(f"[WARN] No finite radial coordinates for {surface_name}; skipping radial profiles")
        return
    r_max = float(df["r"].max())
    bin_count = int(min(25, max(15, np.sqrt(max(len(df), 1)).astype(int))))
    bins = np.linspace(0.0, r_max, bin_count + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    df["r_bin"] = pd.cut(df["r"], bins=bins, include_lowest=True, labels=False)
    grouped = df.groupby("r_bin", observed=True)
    count_mean = grouped["neighbor_count"].mean().reindex(range(bin_count)).to_numpy(dtype=float)
    count_median = grouped["neighbor_count"].median().reindex(range(bin_count)).to_numpy(dtype=float)
    rad_mean = grouped["neighbor_max_radius"].mean().reindex(range(bin_count)).to_numpy(dtype=float)
    rad_median = grouped["neighbor_max_radius"].median().reindex(range(bin_count)).to_numpy(dtype=float)
    qc_num = grouped.apply(lambda g: int(((g["curvature_valid"] == 1) & (g["K_qc_warn_flag"] != 0)).sum())).reindex(range(bin_count), fill_value=0)
    qc_den = grouped.apply(lambda g: int((g["curvature_valid"] == 1).sum())).reindex(range(bin_count), fill_value=0)
    qc_fraction = np.divide(qc_num.to_numpy(dtype=float), qc_den.to_numpy(dtype=float), out=np.full(bin_count, np.nan), where=qc_den.to_numpy(dtype=float) > 0.0)
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True, constrained_layout=True)
    axes[0].plot(centers, count_mean, label="mean", lw=1.8)
    axes[0].plot(centers, count_median, label="median", lw=1.8, ls="--")
    axes[0].set_ylabel("count")
    axes[0].set_title(("Outer" if surface_name == "outer" else "Inner") + " surface: scale radial profiles")
    axes[0].legend()
    axes[1].plot(centers, rad_mean, label="mean", lw=1.8)
    axes[1].plot(centers, rad_median, label="median", lw=1.8, ls="--")
    axes[1].set_ylabel("radius [Å]")
    axes[1].legend()
    axes[2].plot(centers, qc_fraction, lw=1.8, color="black")
    axes[2].set_ylabel("fraction")
    axes[2].set_xlabel("r [Å]")
    out_path = out_dir / f"{surface_name}_scale_radial_profiles.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def compute_scale_enrichment(surface_df: pd.DataFrame, surface_name: str) -> None:
    required = ["neighbor_count", "neighbor_max_radius", "derivative_valid", "curvature_valid", "K_qc_warn_flag"]
    missing = [col for col in required if col not in surface_df.columns]
    if missing:
        print(f"[WARN] Missing required columns for {surface_name} scale enrichment: {missing}")
        return
    df = surface_df.copy()
    df["derivative_valid"] = df["derivative_valid"].astype(int)
    df["curvature_valid"] = df["curvature_valid"].astype(int)
    df["K_qc_warn_flag"] = df["K_qc_warn_flag"].astype(int)
    for col in ("neighbor_count", "neighbor_max_radius"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
        df.loc[(df["derivative_valid"] != 1) | (~np.isfinite(df[col])) | (df[col] <= 0.0), col] = np.nan
    qc_mask = (df["curvature_valid"] == 1) & (df["K_qc_warn_flag"] != 0)
    nonqc_mask = (df["curvature_valid"] == 1) & (df["K_qc_warn_flag"] == 0)
    for col in ("neighbor_count", "neighbor_max_radius"):
        qc_vals = df.loc[qc_mask, col].to_numpy(dtype=float)
        nonqc_vals = df.loc[nonqc_mask, col].to_numpy(dtype=float)
        qc_vals = qc_vals[np.isfinite(qc_vals)]
        nonqc_vals = nonqc_vals[np.isfinite(nonqc_vals)]
        med_qc = float(np.nanmedian(qc_vals)) if qc_vals.size > 0 else np.nan
        med_nonqc = float(np.nanmedian(nonqc_vals)) if nonqc_vals.size > 0 else np.nan
        ratio = med_qc / med_nonqc if (np.isfinite(med_qc) and np.isfinite(med_nonqc) and med_nonqc != 0.0) else np.nan
        label = "neighbor_count" if col == "neighbor_count" else "neighbor_radius"
        print(f"[INFO] {surface_name} {label} enrichment: median_qc={med_qc:.6g}, median_nonqc={med_nonqc:.6g}, ratio={ratio:.6g}")


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return np.nan
    if np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def safe_rank_corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return np.nan
    if scipy_spearmanr is not None:
        rho, _ = scipy_spearmanr(x, y)
        return float(rho)
    x_rank = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    y_rank = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    return safe_corr(x_rank, y_rank)


def describe_group(values: np.ndarray) -> dict[str, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {"median": np.nan, "mean": np.nan, "p95": np.nan, "max": np.nan}
    return {
        "median": float(np.nanmedian(finite)),
        "mean": float(np.nanmean(finite)),
        "p95": float(np.nanpercentile(finite, 95)),
        "max": float(np.nanmax(finite)),
    }


def safe_ratio(a: float, b: float) -> float:
    if not np.isfinite(a) or not np.isfinite(b) or b == 0.0:
        return np.nan
    return float(a / b)


def _build_residual_k_analysis_subset(surface_df: pd.DataFrame) -> pd.DataFrame:
    analysis_mask = (
        (surface_df["curvature_valid"].astype(int) == 1)
        & (surface_df["derivative_valid"].astype(int) == 1)
        & np.isfinite(pd.to_numeric(surface_df["abs_K"], errors="coerce"))
        & np.isfinite(pd.to_numeric(surface_df["fit_rms_residual"], errors="coerce"))
        & (pd.to_numeric(surface_df["fit_rms_residual"], errors="coerce") >= 0.0)
    )
    analysis_df = surface_df.loc[analysis_mask].copy()
    analysis_df["abs_K"] = pd.to_numeric(analysis_df["abs_K"], errors="coerce")
    analysis_df["fit_rms_residual"] = pd.to_numeric(analysis_df["fit_rms_residual"], errors="coerce")
    return analysis_df


def _build_second_derivative_subset(surface_df: pd.DataFrame, deriv_col: str) -> pd.DataFrame:
    analysis_mask = (
        (surface_df["curvature_valid"].astype(int) == 1)
        & (surface_df["derivative_valid"].astype(int) == 1)
        & np.isfinite(pd.to_numeric(surface_df["k"], errors="coerce"))
        & np.isfinite(pd.to_numeric(surface_df["abs_K"], errors="coerce"))
        & np.isfinite(pd.to_numeric(surface_df[deriv_col], errors="coerce"))
    )
    analysis_df = surface_df.loc[analysis_mask].copy()
    analysis_df["abs_K"] = pd.to_numeric(analysis_df["abs_K"], errors="coerce")
    analysis_df[deriv_col] = pd.to_numeric(analysis_df[deriv_col], errors="coerce")
    abs_col = f"abs_{deriv_col}"
    analysis_df[abs_col] = np.abs(analysis_df[deriv_col].to_numpy(dtype=float))
    return analysis_df


def analyze_absK_vs_second_derivative(surface_df: pd.DataFrame, surface_name: str, deriv_col: str, out_dir: Path) -> tuple[pd.DataFrame, float, float]:
    required = ["curvature_valid", "derivative_valid", "k", "abs_K", deriv_col, "K_qc_warn_flag"]
    missing = [col for col in required if col not in surface_df.columns]
    if missing:
        print(f"[WARN] Missing required columns for {surface_name} abs(K)-vs-{deriv_col} analysis: {missing}; skipping")
        return pd.DataFrame(), np.nan, np.nan
    analysis_df = _build_second_derivative_subset(surface_df, deriv_col)
    if len(analysis_df) < 2:
        print(f"[WARN] {surface_name} abs(K)-vs-{deriv_col} subset has <2 rows; skipping")
        return pd.DataFrame(), np.nan, np.nan
    abs_col = f"abs_{deriv_col}"
    x = analysis_df[abs_col].to_numpy(dtype=float)
    y = analysis_df["abs_K"].to_numpy(dtype=float)
    pearson_r = safe_corr(x, y)
    spearman_rho = safe_rank_corr(x, y)
    print(f"[INFO] {surface_name} abs(K) vs abs({deriv_col}): Pearson r={pearson_r:.6g}, Spearman rho={spearman_rho:.6g}, n={len(analysis_df)}")
    qc_mask = analysis_df["K_qc_warn_flag"].astype(int) != 0
    med_qc = describe_group(analysis_df.loc[qc_mask, abs_col].to_numpy(dtype=float))["median"]
    med_nonqc = describe_group(analysis_df.loc[~qc_mask, abs_col].to_numpy(dtype=float))["median"]
    print(f"[INFO] {surface_name} abs({deriv_col}) QC enrichment: median_qc={med_qc:.6g}, median_nonqc={med_nonqc:.6g}, ratio={safe_ratio(med_qc, med_nonqc):.6g}")
    return analysis_df, pearson_r, spearman_rho


def compute_k_num_like(surface_df: pd.DataFrame) -> pd.DataFrame:
    needed = ["d2z_dx2", "d2z_dy2", "d2z_dxdy"]
    missing = [c for c in needed if c not in surface_df.columns]
    if missing:
        raise ValueError(f"Missing second-derivative columns for K_num_like: {missing}")
    out = surface_df.copy()
    dxx = pd.to_numeric(out["d2z_dx2"], errors="coerce")
    dyy = pd.to_numeric(out["d2z_dy2"], errors="coerce")
    dxy = pd.to_numeric(out["d2z_dxdy"], errors="coerce")
    out["K_num_like"] = dxx * dyy - (dxy ** 2)
    out["abs_K_num_like"] = np.abs(out["K_num_like"].to_numpy(dtype=float))
    return out


def plot_k_num_like_heatmaps_both(outer_df: pd.DataFrame, inner_df: pd.DataFrame, out_dir: Path, overlay_qc: bool = False) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    any_plotted = False
    for ax, (surface_name, df) in zip(axes, (("outer", outer_df), ("inner", inner_df))):
        req = ["i", "j", "x", "y", "derivative_valid", "K_num_like"]
        if overlay_qc:
            req += ["curvature_valid", "K_qc_warn_flag"]
        if any(c not in df.columns for c in req):
            ax.set_title(f"{surface_name}: missing cols")
            ax.axis("off")
            continue
        plot_df = df[req].copy()
        plot_df["K_num_like"] = pd.to_numeric(plot_df["K_num_like"], errors="coerce")
        plot_df.loc[(plot_df["derivative_valid"].astype(int) != 1) | (~np.isfinite(plot_df["K_num_like"])), "K_num_like"] = np.nan
        finite = plot_df["K_num_like"].to_numpy(dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            ax.set_title(f"{surface_name}: no valid values")
            ax.axis("off")
            continue
        any_plotted = True
        vlim = float(np.nanpercentile(np.abs(finite), 99))
        if not np.isfinite(vlim) or vlim <= 0:
            ax.set_title(f"{surface_name}: invalid vlim")
            ax.axis("off")
            continue
        grid, extent = dataframe_to_grid(plot_df, "K_num_like")
        im = ax.imshow(grid, origin="lower", extent=extent, cmap="coolwarm", vmin=-vlim, vmax=vlim, interpolation="nearest")
        if overlay_qc:
            qc_df = plot_df.loc[(plot_df["curvature_valid"].astype(int) == 1) & (plot_df["K_qc_warn_flag"].astype(int) != 0)]
            if len(qc_df) > 0:
                ax.scatter(qc_df["x"], qc_df["y"], marker="s", color="black", s=8, alpha=0.6, linewidths=0)
        ax.set_title(("Outer" if surface_name == "outer" else "Inner") + " surface: K numerator-like term" + (" + K QC warnings" if overlay_qc else ""))
        ax.set_xlabel("x [Å]")
        ax.set_ylabel("y [Å]")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("K numerator-like term [Å^-2]")
    if not any_plotted:
        plt.close(fig)
        print("[WARN] No valid data for K_num_like heatmaps")
        return
    out_path = out_dir / ("K_num_like_qc_overlay.png" if overlay_qc else "K_num_like_heatmap.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def plot_k_num_like_heatmap(surface_df: pd.DataFrame, surface_name: str, out_dir: Path, overlay_qc: bool = False) -> None:
    if surface_name == "outer":
        plot_k_num_like_heatmaps_both(surface_df, pd.DataFrame(), out_dir, overlay_qc=overlay_qc)
    else:
        plot_k_num_like_heatmaps_both(pd.DataFrame(), surface_df, out_dir, overlay_qc=overlay_qc)


def analyze_absK_vs_k_num_like(surface_df: pd.DataFrame, surface_name: str) -> pd.DataFrame:
    required = ["curvature_valid", "derivative_valid", "abs_K", "abs_K_num_like", "K_qc_warn_flag"]
    if any(c not in surface_df.columns for c in required):
        return pd.DataFrame()
    mask = (
        (surface_df["curvature_valid"].astype(int) == 1)
        & (surface_df["derivative_valid"].astype(int) == 1)
        & np.isfinite(pd.to_numeric(surface_df["abs_K"], errors="coerce"))
        & np.isfinite(pd.to_numeric(surface_df["abs_K_num_like"], errors="coerce"))
    )
    analysis_df = surface_df.loc[mask].copy()
    if len(analysis_df) < 2:
        print(f"[WARN] {surface_name} abs(K)-vs-abs(K_num_like) subset has <2 rows")
        return pd.DataFrame()
    x = pd.to_numeric(analysis_df["abs_K_num_like"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(analysis_df["abs_K"], errors="coerce").to_numpy(dtype=float)
    print(f"[INFO] {surface_name} abs(K) vs abs(K_num_like): Pearson r={safe_corr(x, y):.6g}, Spearman rho={safe_rank_corr(x, y):.6g}, n={len(analysis_df)}")
    qc_mask = (analysis_df["curvature_valid"].astype(int) == 1) & (analysis_df["K_qc_warn_flag"].astype(int) != 0)
    nonqc_mask = (analysis_df["curvature_valid"].astype(int) == 1) & (analysis_df["K_qc_warn_flag"].astype(int) == 0)
    qc_vals = pd.to_numeric(analysis_df.loc[qc_mask, "abs_K_num_like"], errors="coerce").to_numpy(dtype=float)
    nonqc_vals = pd.to_numeric(analysis_df.loc[nonqc_mask, "abs_K_num_like"], errors="coerce").to_numpy(dtype=float)
    med_qc = float(np.nanmedian(qc_vals)) if np.isfinite(qc_vals).any() else np.nan
    med_nonqc = float(np.nanmedian(nonqc_vals)) if np.isfinite(nonqc_vals).any() else np.nan
    p95_qc = float(np.nanpercentile(qc_vals[np.isfinite(qc_vals)], 95)) if np.isfinite(qc_vals).any() else np.nan
    p95_nonqc = float(np.nanpercentile(nonqc_vals[np.isfinite(nonqc_vals)], 95)) if np.isfinite(nonqc_vals).any() else np.nan
    print(f"[INFO] {surface_name} abs(K_num_like) QC enrichment: median_qc={med_qc:.6g}, median_nonqc={med_nonqc:.6g}, ratio={safe_ratio(med_qc, med_nonqc):.6g}, p95_qc={p95_qc:.6g}, p95_nonqc={p95_nonqc:.6g}")
    return analysis_df

def plot_second_derivative_heatmap(surface_df: pd.DataFrame, surface_name: str, field_col: str, label: str, out_dir: Path, overlay_qc: bool = False) -> None:
    required = ["i", "j", "x", "y", "derivative_valid", field_col]
    if overlay_qc:
        required += ["curvature_valid", "K_qc_warn_flag"]
    missing = [c for c in required if c not in surface_df.columns]
    if missing:
        print(f"[WARN] Missing columns for {surface_name} {field_col} heatmap: {missing}; skipping")
        return
    plot_df = surface_df[required].copy()
    plot_df[field_col] = pd.to_numeric(plot_df[field_col], errors="coerce")
    invalid = (plot_df["derivative_valid"].astype(int) != 1) | (~np.isfinite(plot_df[field_col]))
    plot_df.loc[invalid, field_col] = np.nan
    valid = plot_df[field_col].to_numpy(dtype=float)
    finite = valid[np.isfinite(valid)]
    if finite.size == 0:
        print(f"[WARN] No finite derivative-valid values for {surface_name} {field_col}; skipping")
        return
    vlim = float(np.nanpercentile(np.abs(finite), 99))
    if not np.isfinite(vlim) or vlim <= 0.0:
        print(f"[WARN] Invalid vlim for {surface_name} {field_col}; skipping")
        return
    grid, extent = dataframe_to_grid(plot_df, field_col)
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    im = ax.imshow(grid, origin="lower", extent=extent, cmap="coolwarm", vmin=-vlim, vmax=vlim, interpolation="nearest")
    if overlay_qc:
        qc_df = plot_df.loc[(plot_df["curvature_valid"].astype(int) == 1) & (plot_df["K_qc_warn_flag"].astype(int) != 0)]
        if len(qc_df) > 0:
            ax.scatter(qc_df["x"], qc_df["y"], marker="s", color="black", s=8, alpha=0.6, linewidths=0)
    ax.set_xlabel("x [Å]")
    ax.set_ylabel("y [Å]")
    title_surface = "Outer" if surface_name == "outer" else "Inner"
    ax.set_title(f"{title_surface} surface: {field_col}" + (" + K QC warnings" if overlay_qc else ""))
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)
    suffix = "_qc_overlay" if overlay_qc else ""
    out_path = out_dir / f"{surface_name}_{field_col}_heatmap{suffix}.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")

def analyze_absK_vs_residual(surface_df: pd.DataFrame, surface_name: str, out_dir: Path, *, save_plot: bool = True) -> None:
    required = ["curvature_valid", "derivative_valid", "abs_K", "fit_rms_residual", "K_qc_warn_flag"]
    missing = [col for col in required if col not in surface_df.columns]
    if missing:
        print(f"[WARN] Missing required columns for {surface_name} residual-vs-K analysis: {missing}; skipping")
        return

    total_rows = len(surface_df)
    analysis_df = _build_residual_k_analysis_subset(surface_df)
    qc_mask = analysis_df["K_qc_warn_flag"].astype(int) != 0
    qc_count = int(qc_mask.sum())
    non_qc_count = int((~qc_mask).sum())
    print(f"[INFO] {surface_name} residual-vs-K subset: total={total_rows}, n={len(analysis_df)}, qc={qc_count}, non_qc={non_qc_count}")

    if len(analysis_df) < 2:
        print(f"[WARN] {surface_name} analysis subset has <2 rows; skipping residual-vs-K correlation/plot")
        return

    x = analysis_df["fit_rms_residual"].to_numpy(dtype=float)
    y = analysis_df["abs_K"].to_numpy(dtype=float)
    pearson_r = safe_corr(x, y)
    spearman_rho = safe_rank_corr(x, y)
    print(
        f"[INFO] {surface_name} abs(K) vs RMS residual: "
        f"Pearson r={pearson_r:.6g}, Spearman rho={spearman_rho:.6g}, n={len(analysis_df)}"
    )

    qc_rms = describe_group(analysis_df.loc[qc_mask, "fit_rms_residual"].to_numpy(dtype=float))
    nonqc_rms = describe_group(analysis_df.loc[~qc_mask, "fit_rms_residual"].to_numpy(dtype=float))
    qc_absk = describe_group(analysis_df.loc[qc_mask, "abs_K"].to_numpy(dtype=float))
    nonqc_absk = describe_group(analysis_df.loc[~qc_mask, "abs_K"].to_numpy(dtype=float))
    print(
        f"[INFO] {surface_name} RMS residual median: "
        f"qc={qc_rms['median']:.6g}, nonqc={nonqc_rms['median']:.6g}, "
        f"ratio={safe_ratio(qc_rms['median'], nonqc_rms['median']):.6g}"
    )
    print(
        f"[INFO] {surface_name} abs(K) median: "
        f"qc={qc_absk['median']:.6g}, nonqc={nonqc_absk['median']:.6g}, "
        f"ratio={safe_ratio(qc_absk['median'], nonqc_absk['median']):.6g}"
    )

    p90 = float(np.nanpercentile(x, 90))
    high_residual = x >= p90
    qc_fraction = float(np.mean(high_residual[qc_mask.to_numpy()])) if qc_count > 0 else np.nan
    nonqc_fraction = float(np.mean(high_residual[(~qc_mask).to_numpy()])) if non_qc_count > 0 else np.nan
    print(
        f"[INFO] {surface_name} high-residual enrichment: "
        f"qc_fraction={qc_fraction:.6g}, nonqc_fraction={nonqc_fraction:.6g}, "
        f"ratio={safe_ratio(qc_fraction, nonqc_fraction):.6g}"
    )

    if not save_plot:
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    nonqc_df = analysis_df.loc[~qc_mask]
    qc_df = analysis_df.loc[qc_mask]
    ax.scatter(nonqc_df["fit_rms_residual"], nonqc_df["abs_K"], s=10, alpha=0.35, c="steelblue", linewidths=0, label="non-QC")
    if len(qc_df) > 0:
        ax.scatter(qc_df["fit_rms_residual"], qc_df["abs_K"], s=18, alpha=0.75, c="crimson", linewidths=0, label="QC warning")
    ax.set_xlabel("fit_rms_residual [Å]")
    ax.set_ylabel("abs(K) [Å^-2]")
    ax.set_title(("Outer" if surface_name == "outer" else "Inner") + " surface: abs(K) vs fit RMS residual")
    ax.legend(loc="upper left", fontsize=8)
    annotation = "\n".join(
        [
            f"n={len(analysis_df)}",
            f"Pearson r={pearson_r:.3g}",
            f"Spearman rho={spearman_rho:.3g}",
            f"QC frac={qc_count / len(analysis_df):.3g}",
        ]
    )
    ax.text(0.98, 0.02, annotation, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.6"})
    out_path = out_dir / f"{surface_name}_absK_vs_fit_rms_residual.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def plot_absK_vs_residual_both(outer_df: pd.DataFrame, inner_df: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    panels = [("outer", outer_df, axes[0]), ("inner", inner_df, axes[1])]
    has_points = False
    for surface_name, surface_df, ax in panels:
        analysis_df = _build_residual_k_analysis_subset(surface_df)
        if len(analysis_df) < 2:
            ax.set_title(f"{surface_name}: insufficient valid data")
            ax.axis("off")
            continue
        has_points = True
        qc_mask = analysis_df["K_qc_warn_flag"].astype(int) != 0
        x = analysis_df["fit_rms_residual"].to_numpy(dtype=float)
        y = analysis_df["abs_K"].to_numpy(dtype=float)
        pearson_r = safe_corr(x, y)
        spearman_rho = safe_rank_corr(x, y)
        nonqc_df = analysis_df.loc[~qc_mask]
        qc_df = analysis_df.loc[qc_mask]
        ax.scatter(nonqc_df["fit_rms_residual"], nonqc_df["abs_K"], s=10, alpha=0.35, c="steelblue", linewidths=0, label="non-QC")
        if len(qc_df) > 0:
            ax.scatter(qc_df["fit_rms_residual"], qc_df["abs_K"], s=18, alpha=0.75, c="crimson", linewidths=0, label="QC warning")
        ax.set_xlabel("fit_rms_residual [Å]")
        ax.set_ylabel("abs(K) [Å^-2]")
        ax.set_title(("Outer" if surface_name == "outer" else "Inner") + " surface")
        ax.legend(loc="upper left", fontsize=8)
        qc_fraction = float(qc_mask.mean()) if len(analysis_df) > 0 else np.nan
        annotation = f"n={len(analysis_df)}\nPearson r={pearson_r:.3g}\nSpearman rho={spearman_rho:.3g}\nQC frac={qc_fraction:.3g}"
        ax.text(0.98, 0.02, annotation, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.6"})
    if not has_points:
        plt.close(fig)
        print("[WARN] No valid outer/inner data for combined abs(K)-vs-residual plot")
        return
    out_path = out_dir / "absK_vs_fit_rms_residual.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def plot_absK_vs_second_derivative_both(outer_df: pd.DataFrame, inner_df: pd.DataFrame, out_dir: Path) -> None:
    deriv_cols = ["d2z_dx2", "d2z_dy2", "d2z_dxdy"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for row, (surface_name, surface_df) in enumerate((("outer", outer_df), ("inner", inner_df))):
        for col, deriv_col in enumerate(deriv_cols):
            ax = axes[row, col]
            analysis_df, pearson_r, spearman_rho = analyze_absK_vs_second_derivative(surface_df, surface_name, deriv_col, out_dir)
            if len(analysis_df) < 2:
                ax.set_title(f"{surface_name} {deriv_col}: insufficient data")
                ax.axis("off")
                continue
            qc_mask = analysis_df["K_qc_warn_flag"].astype(int) != 0
            nonqc_df = analysis_df.loc[~qc_mask]
            qc_df = analysis_df.loc[qc_mask]
            ax.scatter(nonqc_df[f"abs_{deriv_col}"], nonqc_df["abs_K"], s=10, alpha=0.35, c="steelblue", linewidths=0, label="non-QC")
            if len(qc_df) > 0:
                ax.scatter(qc_df[f"abs_{deriv_col}"], qc_df["abs_K"], s=18, alpha=0.75, c="crimson", linewidths=0, label="QC warning")
            ax.set_xlabel(f"abs({deriv_col}) [Å^-1]")
            ax.set_ylabel("abs(K) [Å^-2]")
            title_surface = "Outer" if surface_name == "outer" else "Inner"
            ax.set_title(f"{title_surface}: abs(K) vs abs({deriv_col})")
            qc_fraction = float(qc_mask.mean()) if len(analysis_df) > 0 else np.nan
            annotation = f"n={len(analysis_df)}\nPearson r={pearson_r:.3g}\nSpearman rho={spearman_rho:.3g}\nQC frac={qc_fraction:.3g}"
            ax.text(0.98, 0.02, annotation, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                    bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.6"})
            ax.legend(loc="upper left", fontsize=8)
    out_path = out_dir / "absK_vs_abs_d2z_dxdy.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def plot_absK_vs_k_num_like_both(outer_df: pd.DataFrame, inner_df: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    any_points = False
    for ax, (surface_name, df) in zip(axes, (("outer", outer_df), ("inner", inner_df))):
        analysis_df = analyze_absK_vs_k_num_like(df, surface_name)
        if len(analysis_df) < 2:
            ax.set_title(f"{surface_name}: insufficient data")
            ax.axis("off")
            continue
        any_points = True
        qc_mask = analysis_df["K_qc_warn_flag"].astype(int) != 0
        nonqc_df = analysis_df.loc[~qc_mask]
        qc_df = analysis_df.loc[qc_mask]
        ax.scatter(nonqc_df["abs_K_num_like"], nonqc_df["abs_K"], s=10, alpha=0.3, c="steelblue", linewidths=0, label="non-QC")
        if len(qc_df) > 0:
            ax.scatter(qc_df["abs_K_num_like"], qc_df["abs_K"], s=18, alpha=0.75, c="crimson", linewidths=0, label="QC warning")
        ax.set_xlabel("abs(K_num_like) [Å^-2]")
        ax.set_ylabel("abs(K) [Å^-2]")
        ax.set_title(("Outer" if surface_name == "outer" else "Inner") + " surface: abs(K) vs abs(K_num_like)")
        ax.legend(loc="upper left", fontsize=8)
    if not any_points:
        plt.close(fig)
        print("[WARN] No valid data for abs(K) vs abs(K_num_like) plot")
        return
    out_path = out_dir / "absK_vs_abs_K_num_like.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")


def print_k_reconstruction_check(surface_df: pd.DataFrame, surface_name: str) -> None:
    needed = ["dz_dx", "dz_dy", "K_num_like", "k", "curvature_valid", "derivative_valid"]
    if any(c not in surface_df.columns for c in needed):
        return
    df = surface_df.copy()
    mask = (
        (df["curvature_valid"].astype(int) == 1)
        & (df["derivative_valid"].astype(int) == 1)
        & np.isfinite(pd.to_numeric(df["dz_dx"], errors="coerce"))
        & np.isfinite(pd.to_numeric(df["dz_dy"], errors="coerce"))
        & np.isfinite(pd.to_numeric(df["K_num_like"], errors="coerce"))
        & np.isfinite(pd.to_numeric(df["k"], errors="coerce"))
    )
    df = df.loc[mask].copy()
    if len(df) < 2:
        return
    denom = (1.0 + pd.to_numeric(df["dz_dx"]) ** 2 + pd.to_numeric(df["dz_dy"]) ** 2) ** 2
    k_from = pd.to_numeric(df["K_num_like"]) / denom
    k = pd.to_numeric(df["k"])
    diff = np.abs((k_from - k).to_numpy(dtype=float))
    print(f"[INFO] {surface_name} K reconstruction check: median_abs_diff={np.nanmedian(diff):.6g}, max_abs_diff={np.nanmax(diff):.6g}, corr={safe_corr(k_from.to_numpy(dtype=float), k.to_numpy(dtype=float)):.6g}")


def plot_second_derivative_heatmaps_both(outer_df: pd.DataFrame, inner_df: pd.DataFrame, out_dir: Path, overlay_qc: bool = False) -> None:
    fields = [("d2z_dx2", "d2z/dx2 [Å^-1]"), ("d2z_dy2", "d2z/dy2 [Å^-1]"), ("d2z_dxdy", "d2z/dxdy [Å^-1]")]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for row, (surface_name, surface_df) in enumerate((("outer", outer_df), ("inner", inner_df))):
        for col, (field_col, label) in enumerate(fields):
            ax = axes[row, col]
            req = ["i", "j", "x", "y", "derivative_valid", field_col]
            if overlay_qc:
                req += ["curvature_valid", "K_qc_warn_flag"]
            if any(c not in surface_df.columns for c in req):
                ax.set_title(f"{surface_name} {field_col}: missing cols")
                ax.axis("off")
                continue
            df = surface_df[req].copy()
            df[field_col] = pd.to_numeric(df[field_col], errors="coerce")
            df.loc[(df["derivative_valid"].astype(int) != 1) | (~np.isfinite(df[field_col])), field_col] = np.nan
            finite = df[field_col].to_numpy(dtype=float)
            finite = finite[np.isfinite(finite)]
            if finite.size == 0:
                ax.set_title(f"{surface_name} {field_col}: no data")
                ax.axis("off")
                continue
            vlim = float(np.nanpercentile(np.abs(finite), 99))
            if not np.isfinite(vlim) or vlim <= 0.0:
                print(f"[WARN] Invalid vlim for {surface_name} {field_col}; panel skipped")
                ax.set_title(f"{surface_name} {field_col}: invalid vlim")
                ax.axis("off")
                continue
            grid, extent = dataframe_to_grid(df, field_col)
            im = ax.imshow(grid, origin="lower", extent=extent, cmap="coolwarm", vmin=-vlim, vmax=vlim, interpolation="nearest")
            if overlay_qc:
                qc_df = df.loc[(df["curvature_valid"].astype(int) == 1) & (df["K_qc_warn_flag"].astype(int) != 0)]
                if len(qc_df) > 0:
                    ax.scatter(qc_df["x"], qc_df["y"], marker="s", color="black", s=8, alpha=0.6, linewidths=0)
            ax.set_title((("Outer" if surface_name == "outer" else "Inner") + f": {field_col}"))
            ax.set_xlabel("x [Å]")
            ax.set_ylabel("y [Å]")
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(label)
    out_name = "d2z_dxdy_qc_overlay.png" if overlay_qc else "d2z_dxdy_heatmap.png"
    out_path = out_dir / out_name
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved: {out_path}")

def run_surface_merges(loaded: LoadedGeometryCsvs) -> dict[str, Optional[pd.DataFrame]]:
    return {
        "outer": merge_surface_geometry(loaded.outer_curvature, loaded.outer_derivatives, surface="outer"),
        "inner": merge_surface_geometry(loaded.inner_curvature, loaded.inner_derivatives, surface="inner"),
    }


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
        help=(
            "Surface selection for heatmap generation: "
            "outer (outer only), inner (inner only), both (default)."
        ),
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if required surface curvature/derivative CSVs are missing.",
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory for diagnostics output PNG files.",
    )

    parser.add_argument(
        "--k-percentile",
        type=float,
        default=99.0,
        help="Percentile for symmetric K color clipping per surface.",
    )
    parser.add_argument(
        "--residual-percentile",
        type=float,
        default=99.0,
        help="Percentile for residual color clipping with vmin fixed at 0.",
    )
    parser.add_argument(
        "--qc-overlay",
        dest="qc_overlay",
        action="store_true",
        default=True,
        help="Enable K heatmap QC-warning overlay plots (default: enabled).",
    )
    parser.add_argument(
        "--no-qc-overlay",
        dest="qc_overlay",
        action="store_false",
        help="Disable K heatmap QC-warning overlay plots.",
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

    loaded = load_geometry_csvs(files)
    merged_by_surface = run_surface_merges(loaded)

    out_dir = args.out_dir if args.out_dir is not None else args.result_dir / "diagnostics_k"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output directory: {out_dir}")

    if args.surface == "both":
        outer_df = merged_by_surface.get("outer")
        inner_df = merged_by_surface.get("inner")
        if outer_df is None or inner_df is None:
            print("[WARN] Missing merged table for outer/inner; skipping combined plot")
        else:
            plot_k_heatmap_both(outer_df, inner_df, out_dir, args.k_percentile)
            if args.qc_overlay:
                plot_k_heatmap_with_qc_overlay_both(
                    outer_df,
                    inner_df,
                    out_dir,
                    args.k_percentile,
                    strict=args.strict,
                )
            plot_condition_indicator_heatmap(outer_df, "outer", out_dir)
            plot_condition_indicator_heatmap(inner_df, "inner", out_dir)
            print_residual_summary_stats(outer_df, "outer")
            print_residual_summary_stats(inner_df, "inner")
            plot_residual_heatmaps_combined(
                [("outer", outer_df), ("inner", inner_df)],
                out_dir,
                args.residual_percentile,
                overlay_qc=False,
            )
            if args.qc_overlay:
                plot_residual_heatmaps_combined(
                    [("outer", outer_df), ("inner", inner_df)],
                    out_dir,
                    args.residual_percentile,
                    overlay_qc=True,
                )
            for surf_name, surf_df in (("outer", outer_df), ("inner", inner_df)):
                plot_neighbor_heatmap(surf_df, surf_name, "neighbor_count", "neighbor count", out_dir, 99.0, overlay_qc=False)
                plot_neighbor_heatmap(surf_df, surf_name, "neighbor_max_radius", "neighbor max radius [Å]", out_dir, 99.0, overlay_qc=False)
                if args.qc_overlay:
                    plot_neighbor_heatmap(surf_df, surf_name, "neighbor_count", "neighbor count", out_dir, 99.0, overlay_qc=True)
                    plot_neighbor_heatmap(surf_df, surf_name, "neighbor_max_radius", "neighbor max radius [Å]", out_dir, 99.0, overlay_qc=True)
                plot_scale_radial_profiles(surf_df, surf_name, out_dir)
                compute_scale_enrichment(surf_df, surf_name)
                analyze_absK_vs_residual(surf_df, surf_name, out_dir, save_plot=False)
            plot_absK_vs_residual_both(outer_df, inner_df, out_dir)
            plot_second_derivative_heatmaps_both(outer_df, inner_df, out_dir, overlay_qc=False)
            if args.qc_overlay:
                plot_second_derivative_heatmaps_both(outer_df, inner_df, out_dir, overlay_qc=True)
            plot_absK_vs_second_derivative_both(outer_df, inner_df, out_dir)
            outer_df = compute_k_num_like(outer_df)
            inner_df = compute_k_num_like(inner_df)
            plot_k_num_like_heatmaps_both(outer_df, inner_df, out_dir, overlay_qc=False)
            if args.qc_overlay:
                plot_k_num_like_heatmaps_both(outer_df, inner_df, out_dir, overlay_qc=True)
            plot_absK_vs_k_num_like_both(outer_df, inner_df, out_dir)
            print_k_reconstruction_check(outer_df, "outer")
            print_k_reconstruction_check(inner_df, "inner")
    else:
        surface_df = merged_by_surface.get(args.surface)
        if surface_df is None:
            print(f"[WARN] Missing merged table for {args.surface}; skipping plot")
        else:
            plot_k_heatmap(surface_df, args.surface, out_dir, args.k_percentile)
            if args.qc_overlay:
                plot_k_heatmap_with_qc_overlay(
                    surface_df,
                    args.surface,
                    out_dir,
                    args.k_percentile,
                    strict=args.strict,
                )
            plot_condition_indicator_heatmap(surface_df, args.surface, out_dir)
            print_residual_summary_stats(surface_df, args.surface)
            for residual_col, residual_label in (
                ("fit_rms_residual", "fit RMS residual [Å]"),
                ("fit_max_abs_residual", "fit max abs residual [Å]"),
            ):
                plot_residual_heatmap(
                    surface_df, args.surface, residual_col, residual_label, out_dir, args.residual_percentile, False
                )
                if args.qc_overlay:
                    plot_residual_heatmap(
                        surface_df, args.surface, residual_col, residual_label, out_dir, args.residual_percentile, True
                    )
            plot_neighbor_heatmap(surface_df, args.surface, "neighbor_count", "neighbor count", out_dir, 99.0, overlay_qc=False)
            plot_neighbor_heatmap(surface_df, args.surface, "neighbor_max_radius", "neighbor max radius [Å]", out_dir, 99.0, overlay_qc=False)
            if args.qc_overlay:
                plot_neighbor_heatmap(surface_df, args.surface, "neighbor_count", "neighbor count", out_dir, 99.0, overlay_qc=True)
                plot_neighbor_heatmap(surface_df, args.surface, "neighbor_max_radius", "neighbor max radius [Å]", out_dir, 99.0, overlay_qc=True)
            plot_scale_radial_profiles(surface_df, args.surface, out_dir)
            compute_scale_enrichment(surface_df, args.surface)
            analyze_absK_vs_residual(surface_df, args.surface, out_dir)
            for deriv_col, label in (
                ("d2z_dx2", "d2z/dx2 [Å^-1]"),
                ("d2z_dy2", "d2z/dy2 [Å^-1]"),
                ("d2z_dxdy", "d2z/dxdy [Å^-1]"),
            ):
                plot_second_derivative_heatmap(surface_df, args.surface, deriv_col, label, out_dir, overlay_qc=False)
                if args.qc_overlay:
                    plot_second_derivative_heatmap(surface_df, args.surface, deriv_col, label, out_dir, overlay_qc=True)
                analyze_absK_vs_second_derivative(surface_df, args.surface, deriv_col, out_dir)
            surface_df = compute_k_num_like(surface_df)
            plot_k_num_like_heatmap(surface_df, args.surface, out_dir, overlay_qc=False)
            if args.qc_overlay:
                plot_k_num_like_heatmap(surface_df, args.surface, out_dir, overlay_qc=True)
            analyze_absK_vs_k_num_like(surface_df, args.surface)
            print_k_reconstruction_check(surface_df, args.surface)

    success_msg = "[SUCCESS] K heatmaps"
    success_msg += ", QC overlays," if args.qc_overlay else ""
    success_msg += " condition indicators, residual diagnostics, and scale mismatch diagnostics generated."
    print(success_msg)
    print("[SUCCESS] Scale mismatch diagnostics generated.")
    print("[SUCCESS] Residual-vs-K quantitative diagnostics generated.")
    print("[SUCCESS] Second-derivative K-driver diagnostics generated.")
    print("[SUCCESS] K numerator-like diagnostics generated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
