#!/usr/bin/env python3
"""
make_mate_geography_figure.py

Generate the current two-panel MATE geography figure from a state-level DataFrame:

    A. Raw MATE prevalence choropleth
    B. Hierarchical beta-binomial state estimates with 95% empirical-Bayes
       posterior intervals

For the current manuscript, MATE is the former "MATE-EHR v1+" definition:
    retransplantation
    OR index transplant hospitalization >90 days
    OR dialysis/hemodialysis during days 91-365
    OR >=4 unplanned inpatient readmissions through day 365
    OR oxygen dependence / oxygen DME evidence on or after day 90

Expected input columns
----------------------
Required:
    state
    n_patients

And ONE of:
    MATE
    MATE-EHR v1+
    MATE_EHR_v1_plus

Example Python use
------------------
    import pandas as pd
    from make_mate_geography_figure import make_mate_geography_figure

    df = pd.read_csv("state_mate_counts.csv")
    make_mate_geography_figure(df, outdir="PLOTS")

Command-line use
----------------
    python make_mate_geography_figure.py state_mate_counts.csv --outdir PLOTS

Dependencies
------------
    pip install pandas numpy scipy matplotlib geopandas

The U.S. Census state boundary ZIP is downloaded once and cached at:
    ~/.cache/mate/cb_2024_us_state_20m.zip
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import betaln, gammaln
from scipy.stats import beta as beta_dist


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

CENSUS_STATE_URL = (
    "https://www2.census.gov/geo/tiger/GENZ2024/shp/"
    "cb_2024_us_state_20m.zip"
)

US_50_DC = {
    "AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","HI","ID","IL","IN","IA",
    "KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ",
    "NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN","TX","UT","VT",
    "VA","WA","WV","WI","WY","DC"
}

HIGHLIGHT_STATES = ("PA", "TX", "KY", "SC")

# Same manuscript map range used previously.
MAP_VMIN = 0.0
MAP_VMAX = 65.0

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 14,
    "axes.titlesize": 19,
    "axes.labelsize": 16,
    "xtick.labelsize": 10,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
})


# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

def prepare_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize column names, identify the MATE column, filter to 50 states + DC,
    and calculate raw MATE prevalence.
    """
    d = df.copy()

    # State column.
    if "state" not in d.columns:
        if "STATE" in d.columns:
            d = d.rename(columns={"STATE": "state"})
        else:
            raise ValueError("Input DataFrame must contain 'state' or 'STATE'.")

    # Current MATE = former MATE-EHR v1+.
    mate_candidates = ["MATE", "MATE-EHR v1+", "MATE_EHR_v1_plus"]
    mate_col = next((c for c in mate_candidates if c in d.columns), None)
    if mate_col is None:
        raise ValueError(
            "Could not find the MATE count column. Expected one of:\n"
            "  MATE\n"
            "  MATE-EHR v1+\n"
            "  MATE_EHR_v1_plus"
        )

    if "n_patients" not in d.columns:
        raise ValueError("Input DataFrame must contain 'n_patients'.")

    d = d.rename(columns={mate_col: "MATE"})
    d["state"] = d["state"].astype(str).str.strip().str.upper()

    # Excludes *Masked, *Unspecified, GU, PR, etc.
    d = d[d["state"].isin(US_50_DC)].copy()

    if d["state"].duplicated().any():
        dup = sorted(d.loc[d["state"].duplicated(keep=False), "state"].unique())
        raise ValueError(f"Duplicate state rows found: {dup}")

    d["n_patients"] = pd.to_numeric(d["n_patients"], errors="raise")
    d["MATE"] = pd.to_numeric(d["MATE"], errors="raise")

    if (d["n_patients"] <= 0).any():
        raise ValueError("All n_patients values must be positive.")

    if ((d["MATE"] < 0) | (d["MATE"] > d["n_patients"])).any():
        raise ValueError("MATE counts must satisfy 0 <= MATE <= n_patients.")

    d["mate_prevalence"] = 100.0 * d["MATE"] / d["n_patients"]

    return d.sort_values("state").reset_index(drop=True)


# -----------------------------------------------------------------------------
# Census state geometry
# -----------------------------------------------------------------------------

def get_state_boundary_file(cache_dir: str | Path | None = None) -> Path:
    """
    Download/cache the official U.S. Census 1:20m state boundary ZIP.
    """
    if cache_dir is None:
        cache_dir = Path.home() / ".cache" / "mate"
    else:
        cache_dir = Path(cache_dir)

    cache_dir.mkdir(parents=True, exist_ok=True)
    zip_path = cache_dir / "cb_2024_us_state_20m.zip"

    if zip_path.exists() and zip_path.stat().st_size > 0:
        return zip_path

    print("Downloading U.S. Census state boundaries...")
    print(f"  {CENSUS_STATE_URL}")
    print(f"  -> {zip_path}")

    try:
        urllib.request.urlretrieve(CENSUS_STATE_URL, zip_path)
    except Exception as exc:
        zip_path.unlink(missing_ok=True)
        raise RuntimeError(
            "Could not download Census state boundaries.\n"
            "If this machine has no internet access, manually download:\n"
            f"  {CENSUS_STATE_URL}\n"
            "and provide it with --state-boundaries."
        ) from exc

    return zip_path


def prepare_state_geometry(
    df: pd.DataFrame,
    state_boundary_file: str | Path | None = None,
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Read Census state geometry and return:
        lower48 + DC,
        Alaska,
        Hawaii

    Geometry is projected to EPSG:2163 (US National Atlas Equal Area), matching
    the preferred manuscript map appearance.
    """
    boundary_file = (
        Path(state_boundary_file)
        if state_boundary_file is not None
        else get_state_boundary_file()
    )

    states = gpd.read_file(boundary_file)

    if "STUSPS" not in states.columns:
        raise ValueError(
            f"{boundary_file} does not contain the expected Census STUSPS column."
        )

    states = (
        states[states["STUSPS"].isin(df["state"])]
        .rename(columns={"STUSPS": "state"})
        [["state", "geometry"]]
        .copy()
    )

    states = states.merge(df, on="state", how="left")
    states = states.to_crs(epsg=2163)

    lower48 = states[~states["state"].isin(["AK", "HI"])].copy()
    alaska = states[states["state"] == "AK"].copy()
    hawaii = states[states["state"] == "HI"].copy()

    return lower48, alaska, hawaii


# -----------------------------------------------------------------------------
# Hierarchical beta-binomial model
# -----------------------------------------------------------------------------

def fit_beta_binomial(df: pd.DataFrame) -> Tuple[Dict[str, float], pd.DataFrame]:
    """
    Empirical-Bayes hierarchical beta-binomial model:

        k_s | p_s ~ Binomial(n_s, p_s)
        p_s ~ Beta(alpha, beta)

    alpha and beta are estimated by maximum marginal likelihood.

    Returns
    -------
    summary : dict
        Fitted hyperparameters and heterogeneity summaries.
    state_results : DataFrame
        Raw and shrinkage-adjusted estimates with 95% posterior intervals.
    """
    n = df["n_patients"].to_numpy(dtype=float)
    k = df["MATE"].to_numpy(dtype=float)

    def log_choose(nn, kk):
        return gammaln(nn + 1) - gammaln(kk + 1) - gammaln(nn - kk + 1)

    def negative_log_likelihood(log_ab):
        alpha, beta = np.exp(log_ab)
        ll = np.sum(
            log_choose(n, k)
            + betaln(k + alpha, n - k + beta)
            - betaln(alpha, beta)
        )
        return -ll

    # Method-of-moments initialization from raw state proportions.
    p = k / n
    m = float(np.mean(p))
    v = float(np.var(p, ddof=1))
    concentration = max(m * (1.0 - m) / max(v, 1e-8) - 1.0, 2.0)

    init = np.log([
        max(m * concentration, 1e-3),
        max((1.0 - m) * concentration, 1e-3),
    ])

    opt = minimize(
        negative_log_likelihood,
        init,
        method="BFGS",
    )

    if not opt.success:
        # BFGS can occasionally report precision loss despite a stable optimum.
        # Retry with Nelder-Mead.
        opt = minimize(
            negative_log_likelihood,
            init,
            method="Nelder-Mead",
            options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-10},
        )

    alpha, beta = np.exp(opt.x)

    population_mean = alpha / (alpha + beta)
    between_state_sd = np.sqrt(
        alpha * beta /
        ((alpha + beta) ** 2 * (alpha + beta + 1.0))
    )

    # Beta-binomial intraclass correlation / overdispersion.
    rho = 1.0 / (alpha + beta + 1.0)

    # Empirical-Bayes posterior state estimates.
    out = df.copy()
    out["raw_prevalence"] = out["MATE"] / out["n_patients"]
    out["posterior_mean"] = (
        (out["MATE"] + alpha) /
        (out["n_patients"] + alpha + beta)
    )
    out["posterior_lo"] = beta_dist.ppf(
        0.025,
        out["MATE"] + alpha,
        out["n_patients"] - out["MATE"] + beta,
    )
    out["posterior_hi"] = beta_dist.ppf(
        0.975,
        out["MATE"] + alpha,
        out["n_patients"] - out["MATE"] + beta,
    )

    # Model AIC.
    bb_loglik = -negative_log_likelihood(np.log([alpha, beta]))
    aic_bb = 2 * 2 - 2 * bb_loglik

    # Common-rate binomial model for comparison.
    p_common = k.sum() / n.sum()
    binom_loglik = np.sum(
        log_choose(n, k)
        + k * np.log(p_common)
        + (n - k) * np.log(1.0 - p_common)
    )
    aic_binomial = 2 * 1 - 2 * binom_loglik

    summary = {
        "alpha": float(alpha),
        "beta": float(beta),
        "population_mean": float(population_mean),
        "between_state_sd": float(between_state_sd),
        "rho": float(rho),
        "aic_beta_binomial": float(aic_bb),
        "aic_common_binomial": float(aic_binomial),
        "delta_aic": float(aic_binomial - aic_bb),
        "common_rate": float(p_common),
    }

    return summary, out


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def add_state_inset(
    fig,
    gdf: gpd.GeoDataFrame,
    rect,
    column: str,
    cmap,
    norm,
    label: str,
):
    """
    Add compact Alaska/Hawaii inset.
    """
    if gdf.empty:
        return

    ax = fig.add_axes(rect)
    gdf.plot(
        column=column,
        ax=ax,
        cmap=cmap,
        norm=norm,
        edgecolor="black",
        linewidth=0.6,
    )

    xmin, ymin, xmax, ymax = gdf.total_bounds
    dx = xmax - xmin
    dy = ymax - ymin

    ax.set_xlim(xmin - 0.05 * dx, xmax + 0.05 * dx)
    ax.set_ylim(ymin - 0.05 * dy, ymax + 0.05 * dy)
    ax.set_axis_off()
    ax.text(
        0.02, 0.02, label,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def make_mate_geography_figure(
    df: pd.DataFrame,
    outdir: str | Path = "PLOTS",
    state_boundary_file: str | Path | None = None,
    stem: str = "MATE_geography_hierarchical",
) -> Dict[str, Path]:
    """
    Generate the manuscript-ready two-panel MATE figure.

    Panel A:
        raw MATE prevalence choropleth

    Panel B:
        empirical-Bayes beta-binomial shrinkage estimates with 95% intervals

    Saves PNG, PDF, and SVG.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    d = prepare_dataframe(df)
    lower48, alaska, hawaii = prepare_state_geometry(
        d,
        state_boundary_file=state_boundary_file,
    )
    summary, state_results = fit_beta_binomial(d)

    # Merge posterior results with map geometry only if later desired.
    # Current map intentionally shows raw state prevalence.
    cmap = mpl.colormaps["viridis"]
    map_norm = mpl.colors.Normalize(vmin=MAP_VMIN, vmax=MAP_VMAX)

    fig = plt.figure(figsize=(15.5, 6.8))

    # -------------------------------------------------------------------------
    # Panel A: Choropleth
    # -------------------------------------------------------------------------
    ax_map = fig.add_axes([0.035, 0.13, 0.43, 0.77])

    lower48.plot(
        column="mate_prevalence",
        ax=ax_map,
        cmap=cmap,
        norm=map_norm,
        edgecolor="black",
        linewidth=0.55,
    )

    ax_map.set_axis_off()
    ax_map.set_title(
        "MATE prevalence by state",
        fontsize=19,
        fontweight="bold",
        pad=8,
    )

    # Insets are positioned relative to the entire figure.
    add_state_inset(
        fig, alaska,
        rect=[0.055, 0.115, 0.105, 0.16],
        column="mate_prevalence",
        cmap=cmap,
        norm=map_norm,
        label="AK",
    )
    add_state_inset(
        fig, hawaii,
        rect=[0.165, 0.12, 0.065, 0.095],
        column="mate_prevalence",
        cmap=cmap,
        norm=map_norm,
        label="HI",
    )

    # Horizontal colorbar directly beneath Panel A.
    cax = fig.add_axes([0.105, 0.075, 0.29, 0.028])
    sm = mpl.cm.ScalarMappable(norm=map_norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("Prevalence (%)", fontsize=15, labelpad=7)
    cb.ax.tick_params(labelsize=13, length=4)
    cb.outline.set_linewidth(0.8)

    # Panel label.
    fig.text(0.027, 0.925, "A", fontsize=22, fontweight="bold")

    # -------------------------------------------------------------------------
    # Panel B: Hierarchical estimates
    # -------------------------------------------------------------------------
    ax = fig.add_axes([0.515, 0.14, 0.46, 0.73])

    plot_df = (
        state_results
        .sort_values("posterior_mean")
        .reset_index(drop=True)
    )

    x = np.arange(len(plot_df))
    y = 100.0 * plot_df["posterior_mean"].to_numpy()
    lo = 100.0 * plot_df["posterior_lo"].to_numpy()
    hi = 100.0 * plot_df["posterior_hi"].to_numpy()

    ax.errorbar(
        x,
        y,
        yerr=np.vstack([y - lo, hi - y]),
        fmt="o",
        markersize=4.2,
        linewidth=1.0,
        capsize=0,
    )

    pop_mean_pct = 100.0 * summary["population_mean"]
    ax.axhline(
        pop_mean_pct,
        linestyle="--",
        linewidth=1.8,
        label=f"Hierarchical mean = {pop_mean_pct:.1f}%",
    )

    # Highlight manuscript states.
    y_range = max(hi) - min(lo)
    for st in HIGHLIGHT_STATES:
        rows = plot_df.index[plot_df["state"] == st]
        if len(rows) == 0:
            continue
        i = int(rows[0])
        offset = 1.8
        if st == "SC":
            offset = 2.0
        ax.text(
            i,
            y[i] + offset,
            st,
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(
        plot_df["state"],
        rotation=90,
        fontsize=9,
    )

    ax.set_ylabel("Shrinkage-adjusted prevalence (%)", fontsize=15)
    ax.set_xlabel("State / DC", fontsize=15)
    ax.set_title(
        "Hierarchical beta-binomial estimates: MATE",
        fontsize=19,
        fontweight="bold",
        pad=8,
    )
    # Legend, followed immediately by model-summary details beneath it.
    ax.legend(frameon=False, loc="upper left", fontsize=12)

    ax.text(
        0.02,
        0.885,
        (
            f"Between-state SD = {100*summary['between_state_sd']:.1f} pp\n"
            f"ρ = {summary['rho']:.3f}\n"
            f"β-binomial AIC = {summary['aic_beta_binomial']:.1f}\n"
            f"Common-rate binomial AIC = {summary['aic_common_binomial']:.1f}\n"
            f"ΔAIC = {summary['delta_aic']:.1f}"
        ),
        transform=ax.transAxes,
        fontsize=10.5,
        ha="left",
        va="top",
        linespacing=1.22,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=12)

    fig.text(0.498, 0.925, "B", fontsize=22, fontweight="bold")

    # -------------------------------------------------------------------------
    # Save
    # -------------------------------------------------------------------------
    paths = {
        "png": outdir / f"{stem}.png",
        "pdf": outdir / f"{stem}.pdf",
        "svg": outdir / f"{stem}.svg",
    }

    fig.savefig(paths["png"], dpi=500, bbox_inches="tight", facecolor="white")
    fig.savefig(paths["pdf"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["svg"], bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Print manuscript-relevant results.
    total_n = int(d["n_patients"].sum())
    total_mate = int(d["MATE"].sum())
    print(f"State-attributed cohort: {total_mate:,}/{total_n:,} "
          f"({100*total_mate/total_n:.2f}%)")
    print(f"Hierarchical mean: {100*summary['population_mean']:.2f}%")
    print(f"Between-state SD: {100*summary['between_state_sd']:.2f} percentage points")
    print(f"rho: {summary['rho']:.4f}")
    print(f"Beta-binomial AIC: {summary['aic_beta_binomial']:.2f}")
    print(f"Common-rate binomial AIC: {summary['aic_common_binomial']:.2f}")
    print(f"Delta AIC: {summary['delta_aic']:.2f}")

    for st in ("PA", "TX", "KY", "SC"):
        r = state_results.loc[state_results["state"] == st]
        if r.empty:
            continue
        r = r.iloc[0]
        print(
            f"{st}: raw={100*r['raw_prevalence']:.2f}%; "
            f"adjusted={100*r['posterior_mean']:.2f}% "
            f"(95% interval {100*r['posterior_lo']:.2f}%-"
            f"{100*r['posterior_hi']:.2f}%)"
        )

    print("\nSaved:")
    for p in paths.values():
        print(f"  {p}")

    return paths


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate side-by-side MATE prevalence choropleth and "
            "hierarchical beta-binomial state analysis."
        )
    )

    parser.add_argument(
        "csv",
        type=Path,
        help="State-level CSV, e.g. state_mate_counts.csv",
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("PLOTS"),
        help="Output directory (default: PLOTS)",
    )

    parser.add_argument(
        "--state-boundaries",
        type=Path,
        default=None,
        help=(
            "Optional local Census state-boundary ZIP/shapefile. "
            "If omitted, the official Census file is downloaded and cached."
        ),
    )

    parser.add_argument(
        "--stem",
        default="MATE_geography_hierarchical",
        help="Output filename stem.",
    )

    args = parser.parse_args()

    df = pd.read_csv(args.csv)

    make_mate_geography_figure(
        df,
        outdir=args.outdir,
        state_boundary_file=args.state_boundaries,
        stem=args.stem,
    )


if __name__ == "__main__":
    main()

