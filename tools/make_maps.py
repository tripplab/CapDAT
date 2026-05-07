#!/usr/bin/env python3
import argparse
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate composed 2D field maps for radial thickness and outer/inner curvature "
            "for one capsid at a specified cylinder radius."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Example:\n"
            "  python3 tools/make_maps.py --results results --capsid 1cwp "
            "--fold 2_0,2_1,3_0,5_0 --cyl_radius 35\n\n"
            "Outputs:\n"
            "  results/{capsid}/R{cyl_radius}_thickness_map.png\n"
            "  results/{capsid}/R{cyl_radius}_outer_H_map.png\n"
            "  results/{capsid}/R{cyl_radius}_inner_H_map.png\n"
            "  results/{capsid}/R{cyl_radius}_outer_K_map.png\n"
            "  results/{capsid}/R{cyl_radius}_inner_K_map.png"
        ),
    )
    p.add_argument("--results", required=True, help="Results root directory")
    p.add_argument("--capsid", required=True, help="Capsid id (e.g., 1cwp)")
    p.add_argument("--fold", required=True, help="Comma-separated folds in desired panel order")
    p.add_argument("--cyl_radius", required=True, type=int, help="Cylinder radius R")
    p.add_argument("--res", default=500, type=int, help="Grid resolution per dimension (default: 500)")
    return p.parse_args()


def _grid_field(df, value_col, valid_col, res, cyl_radius):
    d = df[["x", "y", value_col, valid_col]].copy()
    d = d.dropna(subset=["x", "y", value_col, valid_col])

    gx = np.linspace(-cyl_radius, cyl_radius, res)
    gy = np.linspace(-cyl_radius, cyl_radius, res)
    X, Y = np.meshgrid(gx, gy)

    valid_mask = d[valid_col].astype(bool).values
    pts_all = d[["x", "y"]].values
    pts_valid = d.loc[valid_mask, ["x", "y"]].values
    vals_valid = d.loc[valid_mask, value_col].values

    if pts_valid.shape[0] < 3:
        raise ValueError(f"Not enough valid points for interpolation of {value_col}.")

    grid_linear = griddata(pts_valid, vals_valid, (X, Y), method="linear")
    grid_nearest = griddata(pts_valid, vals_valid, (X, Y), method="nearest")
    grid_val = np.where(np.isnan(grid_linear), grid_nearest, grid_linear)

    radial_mask = (X ** 2 + Y ** 2) <= float(cyl_radius) ** 2

    # Define support strictly from the convex hull of valid samples only.
    # This avoids nearest-neighbor extrapolation into unsupported peripheral
    # regions, which can otherwise create ring-like artifacts at the boundary.
    support_seed = np.ones(pts_valid.shape[0], dtype=float)
    support_linear = griddata(pts_valid, support_seed, (X, Y), method="linear")
    support_region = np.isfinite(support_linear)

    grid = np.full_like(grid_val, np.nan, dtype=float)
    keep_mask = radial_mask & support_region
    grid[keep_mask] = grid_val[keep_mask]
    extent = [-cyl_radius, cyl_radius, -cyl_radius, cyl_radius]
    return X, Y, grid, extent


def _read_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def _build_panel_data(results_dir, capsid, folds, r, kind, res):
    out = []
    pooled_valid_values = []

    for fold in folds:
        base = Path(results_dir) / capsid / fold
        if not base.exists():
            raise FileNotFoundError(f"Missing required fold directory: {base}")

        if kind == "thickness":
            fp = base / f"R{r}_thickness_radial.csv"
            df = _read_csv(fp)
            value_col = "thickness_radial"
            valid_col = "thickness_valid"
        elif kind == "outer_H":
            fp = base / f"R{r}_outer_curvature.csv"
            df = _read_csv(fp)
            value_col = "H_oriented"
            valid_col = "curvature_valid"
        elif kind == "inner_H":
            fp = base / f"R{r}_inner_curvature.csv"
            df = _read_csv(fp)
            value_col = "H_oriented"
            valid_col = "curvature_valid"
        elif kind == "outer_K":
            fp = base / f"R{r}_outer_curvature.csv"
            df = _read_csv(fp)
            value_col = "K_raw"
            valid_col = "curvature_valid"
        elif kind == "inner_K":
            fp = base / f"R{r}_inner_curvature.csv"
            df = _read_csv(fp)
            value_col = "K_raw"
            valid_col = "curvature_valid"
        else:
            raise ValueError(kind)

        _, _, grid, extent = _grid_field(df, value_col=value_col, valid_col=valid_col, res=res, cyl_radius=r)
        out.append((fold, grid, extent))

        valid_vals = df.loc[df[valid_col].astype(bool), value_col].dropna().values
        if valid_vals.size > 0:
            pooled_valid_values.append(valid_vals)

    if not pooled_valid_values:
        raise ValueError(f"No valid values found for {kind}")

    pooled = np.concatenate(pooled_valid_values)
    vmin = np.percentile(pooled, 1)
    vmax = np.percentile(pooled, 99)

    return out, vmin, vmax


def _render_composed(panels, vmin, vmax, title, outfile, cmap, center_zero=False):
    n = len(panels)
    ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))

    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    })

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.8 * ncols, 4.4 * nrows), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax) if center_zero else None

    mappable = None
    for idx, (fold, grid, extent) in enumerate(panels):
        ax = axes[idx]
        if norm is not None:
            mappable = ax.imshow(grid, origin="lower", cmap=cmap, norm=norm, extent=extent)
        else:
            mappable = ax.imshow(grid, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, extent=extent)

        ax.set_title(f"Fold {fold}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal")
        ax.set_facecolor("black")

    for ax in axes[n:]:
        ax.set_visible(False)

    cbar = fig.colorbar(mappable, ax=axes[:n], shrink=0.92, pad=0.02)
    cbar.set_label(title)

    fig.suptitle(title, fontsize=15)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=300)
    plt.close(fig)


def main():
    args = parse_args()
    global np, pd, plt, LinearSegmentedColormap, TwoSlopeNorm, griddata
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    from scipy.interpolate import griddata

    folds = [x.strip() for x in args.fold.split(",") if x.strip()]
    if not folds:
        raise ValueError("--fold must be a non-empty comma-separated list")

    results_dir = Path(args.results)
    capsid = args.capsid
    r = args.cyl_radius

    rb_cmap = LinearSegmentedColormap.from_list("rbw", ["red", "white", "blue"], N=256)
    rb_cmap.set_bad("black")
    viridis = plt.get_cmap("viridis").copy()
    viridis.set_bad("black")

    jobs = [
        ("thickness", f"R{r}_thickness_map.png", "Thickness (radial) [Ang]", viridis, False),
        ("outer_H", f"R{r}_outer_H_map.png", "Outer H (H_oriented) [Ang^-1]", rb_cmap, True),
        ("inner_H", f"R{r}_inner_H_map.png", "Inner H (H_oriented) [Ang^-1]", rb_cmap, True),
        ("outer_K", f"R{r}_outer_K_map.png", "Outer K (K_raw) [Ang^-2]", rb_cmap, True),
        ("inner_K", f"R{r}_inner_K_map.png", "Inner K (K_raw) [Ang^-2]", rb_cmap, True),
    ]

    for kind, out_name, title, cmap, center_zero in jobs:
        panels, vmin, vmax = _build_panel_data(results_dir, capsid, folds, r, kind, args.res)
        out = results_dir / capsid / out_name
        _render_composed(panels, vmin, vmax, title, out, cmap, center_zero=center_zero)


if __name__ == "__main__":
    main()
