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
        print(f"[INFO] {label} renaming columns: {rename_map}")
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

    success_msg = "[SUCCESS] Basic K heatmaps"
    success_msg += " and QC overlays generated." if args.qc_overlay else " generated."
    print(success_msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
