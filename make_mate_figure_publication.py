#!/usr/bin/env python3
"""Publication-ready MATE geographic figure, reproducible from the repository CSV.

Outputs a US state choropleth and empirical-Bayes beta-binomial state estimates
with larger print-readable fonts. The reported MATE counts use the five-component
oxygen-inclusive MATE column (historically called 'MATE-EHR v1+' in the CSV).
The statistical method is intentionally identical to the original analysis.

Example:
    python make_mate_figure_publication.py state_mate_counts.csv \
        --outdir PLOTS_publication --font-scale 1.0

Dependencies: numpy pandas scipy matplotlib geopandas
Boundaries: 2024 US Census 1:20m cartographic boundaries, cached locally;
            pass --state-boundaries /path/to/cb_2024_us_state_20m.zip
            when running without internet.
"""
from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import betaln, gammaln
from scipy.stats import beta as beta_distribution

CENSUS_URL = (
    "https://www2.census.gov/geo/tiger/GENZ2024/shp/"
    "cb_2024_us_state_20m.zip"
)
US_50_DC = set((
    "AL AK AZ AR CA CO CT DE FL GA HI ID IL IN IA KS KY LA ME MD MA MI MN MS "
    "MO MT NE NV NH NJ NM NY NC ND OH OK OR PA RI SC SD TN TX UT VT VA WA "
    "WV WI WY DC"
).split())
HIGHLIGHT = ("PA", "TX", "KY", "SC")


def load_counts(csv_file: Path) -> pd.DataFrame:
    """Check that the requested outcome is the one used in the paper."""
    df = pd.read_csv(csv_file)
    if "state" not in df and "STATE" in df:
        df = df.rename(columns={"STATE": "state"})
    if "MATE" in df:
        mate_key = "MATE"
    elif "MATE-EHR v1+" in df:
        mate_key = "MATE-EHR v1+"
    elif "MATE_EHR_v1_plus" in df:
        mate_key = "MATE_EHR_v1_plus"
    else:
        raise ValueError("Expected MATE, MATE-EHR v1+, or MATE_EHR_v1_plus column")
    required = {"state", "n_patients", mate_key}
    if not required.issubset(df):
        raise ValueError(f"Missing required columns: {sorted(required - set(df))}")
    df = df.rename(columns={mate_key: "MATE"}).copy()
    df["state"] = df["state"].str.upper().str.strip()
    df = df[df["state"].isin(US_50_DC)].copy()
    if len(df) != 51 or df["state"].duplicated().any():
        raise ValueError("Expected exactly one row for each of the 50 states and DC")
    for key in ("n_patients", "MATE"):
        df[key] = pd.to_numeric(df[key], errors="raise")
        if not np.allclose(df[key].values, np.round(df[key].values)):
            raise ValueError(f"Noninteger {key} values")
    if (df["n_patients"] <= 0).any() or (df["MATE"] < 0).any():
        raise ValueError("Nonpositive denominators or negative case counts")
    if (df["MATE"] > df["n_patients"]).any():
        raise ValueError("MATE case counts exceed denominators")
    df["prevalence"] = 100.0 * df["MATE"] / df["n_patients"]
    return df.sort_values("state").reset_index(drop=True)


def beta_binomial_fit(df: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    """Match the original fit: maximize beta-binomial marginal likelihood."""
    n = df["n_patients"].to_numpy(dtype=float)
    k = df["MATE"].to_numpy(dtype=float)
    log_choose = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)

    def objective(log_ab: np.ndarray) -> float:
        a, b = np.exp(log_ab)
        return -float(np.sum(log_choose + betaln(k + a, n - k + b) - betaln(a, b)))

    p = k / n
    m = float(p.mean())
    v = float(p.var(ddof=1))
    concentration = max(m * (1 - m) / max(v, 1e-8) - 1, 2.0)
    start = np.log([max(m * concentration, 1e-3),
                    max((1 - m) * concentration, 1e-3)])
    opt = minimize(objective, start, method="BFGS")
    if not opt.success:
        opt = minimize(objective, start, method="Nelder-Mead",
                       options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-10})
    if not opt.success and not np.isfinite(opt.fun):
        raise RuntimeError(f"Beta-binomial estimation failed: {opt.message}")
    a, b = np.exp(opt.x)
    mu = float(a / (a + b))
    sd = float(np.sqrt(a * b / ((a + b) ** 2 * (a + b + 1))))
    rho = float(1 / (a + b + 1))
    beta_aic = float(4 + 2 * objective(np.log([a, b])))
    common_p = float(k.sum() / n.sum())
    binomial_ll = float(np.sum(log_choose + k * np.log(common_p) +
                               (n - k) * np.log(1 - common_p)))
    binomial_aic = float(2 - 2 * binomial_ll)
    out = df.copy()
    out["posterior_mean"] = (k + a) / (n + a + b)
    out["posterior_lo"] = beta_distribution.ppf(.025, k + a, n - k + b)
    out["posterior_hi"] = beta_distribution.ppf(.975, k + a, n - k + b)
    return {
        "mean": mu, "sd": sd, "rho": rho,
        "aic_beta": beta_aic, "aic_binomial": binomial_aic,
        "delta_aic": binomial_aic - beta_aic,
    }, out


def resolve_boundaries(user_file: Path | None) -> Path:
    if user_file is not None:
        if not user_file.exists():
            raise FileNotFoundError(f"Boundary file does not exist: {user_file}")
        return user_file
    cache = Path.home() / ".cache" / "mate" / "cb_2024_us_state_20m.zip"
    if not cache.is_file() or cache.stat().st_size == 0:
        cache.parent.mkdir(parents=True, exist_ok=True)
        try:
            print(f"Downloading state boundaries from {CENSUS_URL}")
            urllib.request.urlretrieve(CENSUS_URL, cache)
        except Exception as exc:
            cache.unlink(missing_ok=True)
            raise RuntimeError(
                "Download failed. Manually download the official Census ZIP and "
                "pass its path with --state-boundaries."
            ) from exc
    return cache


def load_boundaries(counts: pd.DataFrame, zip_or_shp: Path) -> tuple:
    geo = gpd.read_file(zip_or_shp)
    if "STUSPS" not in geo.columns:
        raise ValueError("The Census boundary file needs an STUSPS field")
    geo = geo.loc[geo["STUSPS"].isin(US_50_DC), ["STUSPS", "geometry"]].copy()
    geo = geo.rename(columns={"STUSPS": "state"})
    if set(geo["state"]) != US_50_DC:
        raise ValueError("Boundary file is missing one or more states/DC")
    geo = geo.merge(counts[["state", "prevalence"]], on="state", how="left")
    geo = geo.to_crs("EPSG:2163")
    return (geo.loc[~geo["state"].isin(["AK", "HI"])].copy(),
            geo.loc[geo["state"] == "AK"].copy(),
            geo.loc[geo["state"] == "HI"].copy())


def plot_inset(fig, geo, rect, title, cmap, norm, fs):
    ax = fig.add_axes(rect)
    geo.plot(column="prevalence", ax=ax, cmap=cmap, norm=norm,
             edgecolor="#303642", linewidth=0.65)
    xmin, ymin, xmax, ymax = geo.total_bounds
    xpad, ypad = max((xmax - xmin) * .05, 1), max((ymax - ymin) * .04, 1)
    ax.set_xlim(xmin-xpad, xmax+xpad)
    ax.set_ylim(ymin-ypad, ymax+ypad)
    ax.set_axis_off()
    ax.text(.02, -.02, title, ha="left", va="top", fontsize=fs,
            fontweight="bold", transform=ax.transAxes)


def make_figure(counts: pd.DataFrame, outdir: Path, boundary_file: Path,
                *, stem: str, font_scale: float, state_labels: str) -> dict:
    result, estimates = beta_binomial_fit(counts)
    lower48, ak, hi = load_boundaries(counts, boundary_file)
    outdir.mkdir(parents=True, exist_ok=True)
    F = font_scale
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "pdf.fonttype": 42,
        "svg.fonttype": "none", "axes.linewidth": 1.2,
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })
    color = "#17638a"
    cm = mpl.colormaps["viridis"]
    norm = mpl.colors.Normalize(vmin=0, vmax=65)
    fig = plt.figure(figsize=(18.8, 8.5), facecolor="white")
    try:
        # Panel A: original state boundaries, raw 5-component MATE proportion.
        ax_map = fig.add_axes([.025, .275, .355, .59])
        lower48.plot(column="prevalence", ax=ax_map, cmap=cm, norm=norm,
                     edgecolor="#343a45", linewidth=.75)
        ax_map.set_axis_off()
        fig.text(.028, .942, "A", fontsize=25*F, fontweight="bold", va="top")
        fig.text(.205, .918, "MATE prevalence by state", fontsize=19*F,
                 fontweight="bold", ha="center", va="top")
        plot_inset(fig, ak, [.04, .128, .118, .155], "AK", cm, norm, 15*F)
        plot_inset(fig, hi, [.178, .146, .070, .093], "HI", cm, norm, 15*F)
        cb_axis = fig.add_axes([.257, .158, .113, .030])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cm),
                          cax=cb_axis, orientation="horizontal")
        cb.set_ticks([0, 15, 30, 45, 60])
        cb.ax.tick_params(labelsize=15*F, length=5, pad=4)
        cb.set_label("MATE prevalence (%)", fontsize=17*F, labelpad=8)

        # Panel B: wider plot, including all 51 vertically rotated state labels.
        ax = fig.add_axes([.450, .270, .525, .595])
        p = estimates.sort_values("posterior_mean").reset_index(drop=True)
        x = np.arange(len(p))
        y = 100 * p["posterior_mean"].to_numpy()
        lo = 100 * p["posterior_lo"].to_numpy()
        hi = 100 * p["posterior_hi"].to_numpy()
        ax.errorbar(x, y, yerr=np.vstack((y-lo, hi-y)), fmt="o",
                    color=color, ecolor=color, markersize=5.3*F,
                    elinewidth=1.35*F, capsize=0, zorder=3)
        mean_pct = 100*result["mean"]
        ax.axhline(mean_pct, linestyle="--", linewidth=1.75*F,
                   color="#6c7b89", zorder=1)
        # Vertical leader arrows stop just ABOVE each state's upper interval.
        # Each label is directly aligned with its posterior mean/error bar.
        label_heights = {"PA": 46.0, "TX": 56.0, "KY": 63.0, "SC": 78.0}
        for st in HIGHLIGHT:
            idx = p.index[p["state"] == st]
            if len(idx):
                i = int(idx[0])
                arrow_end_y = hi[i] + 1.0  # do not touch the interval
                label_y = max(label_heights[st], arrow_end_y + 10.0)
                ax.annotate(
                    st,
                    xy=(i, arrow_end_y),
                    xycoords="data",
                    xytext=(i, label_y),
                    textcoords="data",
                    ha="center",
                    va="bottom",
                    fontsize=16*F,
                    fontweight="bold",
                    color="#1f2937",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=.98),
                    arrowprops=dict(
                        arrowstyle="-|>",
                        color="#4b5563",
                        lw=1.25*F,
                        shrinkA=5,
                        shrinkB=0,
                        mutation_scale=9*F,
                        connectionstyle="arc3,rad=0",
                    ),
                    zorder=5,
                )
        ax.set_xlim(-1.5, 51.5)
        # Headroom keeps vertical leaders and their labels inside the axes.
        ax.set_ylim(max(0, min(lo)-4), max(85, max(hi)+17))
        ax.set_ylabel("Adjusted MATE prevalence (%)", fontsize=18*F, labelpad=10)
        ax.set_xlabel("States ranked by estimated prevalence", fontsize=17*F, labelpad=12)
        ax.tick_params(axis="y", labelsize=15*F, length=5, width=1.1)
        if state_labels == "all":
            ax.set_xticks(x)
            ax.set_xticklabels(p["state"], rotation=90, fontsize=10.2*F, ha="center", va="top")
            ax.tick_params(axis="x", pad=3, length=3, width=.8)
        else:
            ax.set_xticks([])
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        fig.text(.413, .942, "B", fontsize=25*F, fontweight="bold", va="top")
        fig.text(.713, .918, "Hierarchical state estimates", fontsize=21*F,
                 fontweight="bold", va="top", ha="center")
        # Place the model summary in the free upper-left part of panel B.
        details = (
            f"Mean {mean_pct:.1f}%  |  SD {100*result['sd']:.1f} pp"
            f"  |  ρ {result['rho']:.3f}\n"
            f"AIC: beta-binomial {result['aic_beta']:.1f}, "
            f"binomial {result['aic_binomial']:.1f}, "
            f"Δ {result['delta_aic']:.1f}"
        )
        ax.text(.02, .98, details, transform=ax.transAxes, ha="left", va="top",
                fontsize=14*F, color="#273444", linespacing=1.5)
        fig.text(.71, .025,
                 "Dots: posterior means; bars: 95% posterior intervals; "
                 "dashed line: fitted state mean.",
                 fontsize=11.5*F, color="#394454", ha="center")

        paths = {ext: outdir / f"{stem}.{ext}" for ext in ("png", "pdf", "svg")}
        fig.savefig(paths["png"], dpi=400)
        fig.savefig(paths["pdf"])
        fig.savefig(paths["svg"])
    finally:
        plt.close(fig)

    n = int(counts["n_patients"].sum())
    k = int(counts["MATE"].sum())
    print(f"State-attributed MATE: {k:,}/{n:,} = {100*k/n:.2f}%")
    print(f"Fitted state mean: {100*result['mean']:.2f}%")
    print(f"Between-state SD: {100*result['sd']:.2f} pp; rho: {result['rho']:.4f}")
    print(f"AIC beta-binomial {result['aic_beta']:.2f}, "
          f"binomial {result['aic_binomial']:.2f}, delta {result['delta_aic']:.2f}")
    print("Outputs:")
    for file in paths.values():
        print(f"  {file}")
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="Input state_mate_counts.csv")
    parser.add_argument("--outdir", type=Path, default=Path("PLOTS_publication"))
    parser.add_argument("--state-boundaries", type=Path, default=None,
                        help="Local Census ZIP or shapefile (if offline)")
    parser.add_argument("--stem", default="MATE_geography_publication")
    parser.add_argument("--font-scale", type=float, default=1.10,
                        help="Multiply all font sizes; try 1.0 or 1.2")
    parser.add_argument("--state-labels", choices=("highlighted", "all"),
                        default="all",
                        help="Default: show all 51 state/DC abbreviations at 90 degrees")
    args = parser.parse_args()
    if not 0.75 <= args.font_scale <= 1.5:
        parser.error("--font-scale must be between 0.75 and 1.5")
    df = load_counts(args.csv)
    boundary_file = resolve_boundaries(args.state_boundaries)
    make_figure(df, args.outdir, boundary_file, stem=args.stem,
                font_scale=args.font_scale, state_labels=args.state_labels)


if __name__ == "__main__":
    main()