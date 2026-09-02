#!/usr/bin/env python3
"""
Generate the five standalone MATE figures from a state-level pandas DataFrame.

Figures
-------
1. Lung transplant recipients represented in COSMOS (n_patients)
2. MATE-EHR v1 prevalence by state
3. MATE-EHR v1+ prevalence by state
4. State-level distribution of MATE-EHR v1 prevalence + logit-normal fit
5. State-level distribution of MATE-EHR v1+ prevalence + logit-normal fit

Expected DataFrame columns
--------------------------
state                    State abbreviation, e.g. "KY", "TX", "PA"
n_patients               Number of transplant recipients
MATE-EHR v1              Number meeting MATE-EHR v1
MATE-EHR v1+             Number meeting MATE-EHR v1+

The function also accepts "STATE" instead of "state".

Primary use from Python
-----------------------
    from make_mate_figures import make_all_figures
    make_all_figures(df, outdir="PLOTS")

Command-line use
----------------
    python make_mate_figures.py state_mate_counts.csv --outdir PLOTS

Dependencies
------------
    pandas
    numpy
    scipy
    matplotlib
    geopandas
    Internet access on first run to download the official U.S. Census state-boundary ZIP

Install if needed:
    pip install pandas numpy scipy matplotlib geopandas
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats


# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

US_50_DC = {
    "AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","HI","ID","IL","IN","IA",
    "KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ",
    "NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN","TX","UT","VT",
    "VA","WA","WV","WI","WY","DC"
}

MARK_STATES = ("PA", "TX", "KY")

# Common scale for the two prevalence choropleths so they are directly comparable.
PREVALENCE_MAP_MIN = 0.0
PREVALENCE_MAP_MAX = 65.0

# Publication typography.
FONT_FAMILY = "DejaVu Sans"
TITLE_SIZE = 20
AXIS_LABEL_SIZE = 17
TICK_SIZE = 15
LEGEND_SIZE = 14
COLORBAR_TICK_SIZE = 16
COLORBAR_LABEL_SIZE = 18


# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize expected state-column naming without changing MATE column names.
    """
    d = df.copy()

    if "state" not in d.columns:
        if "STATE" in d.columns:
            d = d.rename(columns={"STATE": "state"})
        else:
            raise ValueError("DataFrame must contain a 'state' or 'STATE' column.")

    required = ["state", "n_patients", "MATE-EHR v1", "MATE-EHR v1+"]
    missing = [c for c in required if c not in d.columns]
    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing)
            + "\nExpected columns: state, n_patients, MATE-EHR v1, MATE-EHR v1+"
        )

    d["state"] = d["state"].astype(str).str.strip().str.upper()

    # Use only 50 states + DC. This removes *Masked, *Unspecified, GU, PR, etc.
    d = d[d["state"].isin(US_50_DC)].copy()

    # One row per state is expected.
    if d["state"].duplicated().any():
        dup = d.loc[d["state"].duplicated(keep=False), "state"].tolist()
        raise ValueError(f"Duplicate state rows found: {sorted(set(dup))}")

    for col in ["n_patients", "MATE-EHR v1", "MATE-EHR v1+"]:
        d[col] = pd.to_numeric(d[col], errors="raise")

    if (d["n_patients"] <= 0).any():
        bad = d.loc[d["n_patients"] <= 0, "state"].tolist()
        raise ValueError(f"n_patients must be > 0. Invalid states: {bad}")

    d["MATE_v1_prevalence"] = 100.0 * d["MATE-EHR v1"] / d["n_patients"]
    d["MATE_v1plus_prevalence"] = 100.0 * d["MATE-EHR v1+"] / d["n_patients"]

    return d.sort_values("state").reset_index(drop=True)


# -----------------------------------------------------------------------------
# Geography
# -----------------------------------------------------------------------------

CENSUS_STATE_URL = (
    "https://www2.census.gov/geo/tiger/GENZ2024/shp/"
    "cb_2024_us_state_20m.zip"
)


def _get_state_boundary_file(cache_dir: str | Path | None = None) -> Path:
    """
    Return a local Census state-boundary ZIP file.

    The file is downloaded once from the U.S. Census Bureau and cached locally.
    This removes the previous dependency on ``basemap-data``.

    Default cache:
        ~/.cache/mate/cb_2024_us_state_20m.zip
    """
    import urllib.request

    if cache_dir is None:
        cache_dir = Path.home() / ".cache" / "mate"
    else:
        cache_dir = Path(cache_dir)

    cache_dir.mkdir(parents=True, exist_ok=True)
    zip_path = cache_dir / "cb_2024_us_state_20m.zip"

    if zip_path.exists() and zip_path.stat().st_size > 0:
        return zip_path

    print(f"Downloading U.S. state boundaries from Census:\n  {CENSUS_STATE_URL}")
    print(f"Caching at:\n  {zip_path}")

    try:
        urllib.request.urlretrieve(CENSUS_STATE_URL, zip_path)
    except Exception as exc:
        if zip_path.exists():
            zip_path.unlink(missing_ok=True)
        raise RuntimeError(
            "Could not download the U.S. Census state boundary file.\n"
            f"URL: {CENSUS_STATE_URL}\n\n"
            "Either ensure this machine has internet access, or manually download\n"
            "that ZIP file and save it as:\n"
            f"  {zip_path}\n"
        ) from exc

    return zip_path


def _make_state_geometry(
    df: pd.DataFrame,
    state_boundary_file: str | Path | None = None,
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Return lower-48+DC, Alaska, and Hawaii in US National Atlas Equal Area
    projection (EPSG:2163).

    Parameters
    ----------
    df
        Prepared MATE state DataFrame.
    state_boundary_file
        Optional local Census state shapefile or ZIP. If omitted, the script
        automatically downloads/caches the Census 2024 1:20m state boundaries.

    Notes
    -----
    Alaska and Hawaii are subsequently drawn as compact lower-left insets.
    """
    if state_boundary_file is None:
        state_boundary_file = _get_state_boundary_file()
    else:
        state_boundary_file = Path(state_boundary_file)

    # GeoPandas can read the Census ZIP directly.
    states = gpd.read_file(state_boundary_file)

    # Census cartographic boundary files use STUSPS for postal abbreviations.
    if "STUSPS" not in states.columns:
        raise ValueError(
            f"{state_boundary_file} does not contain the expected Census "
            "'STUSPS' state-abbreviation column."
        )

    states = (
        states[states["STUSPS"].isin(df["state"])]
        .rename(columns={"STUSPS": "state"})
        [["state", "geometry"]]
        .copy()
    )

    states = states.merge(df, on="state", how="left")

    # US National Atlas Equal Area: same projection used in the preferred plots.
    states = states.to_crs(epsg=2163)

    lower48 = states[~states["state"].isin(["AK", "HI"])].copy()
    alaska = states[states["state"] == "AK"].copy()
    hawaii = states[states["state"] == "HI"].copy()

    return lower48, alaska, hawaii


# -----------------------------------------------------------------------------
# Output helper
# -----------------------------------------------------------------------------

def _save_figure(fig: mpl.figure.Figure, outdir: Path, stem: str) -> Dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)

    paths = {
        "png": outdir / f"{stem}.png",
        "pdf": outdir / f"{stem}.pdf",
        "svg": outdir / f"{stem}.svg",
    }

    fig.savefig(paths["png"], dpi=500, bbox_inches="tight", facecolor="white")
    fig.savefig(paths["pdf"], bbox_inches="tight", facecolor="white")
    fig.savefig(paths["svg"], bbox_inches="tight", facecolor="white")

    return paths


# -----------------------------------------------------------------------------
# Choropleths
# -----------------------------------------------------------------------------

def _add_state_inset(
    fig: mpl.figure.Figure,
    gdf: gpd.GeoDataFrame,
    rect: Iterable[float],
    column: str,
    cmap,
    norm: mpl.colors.Normalize,
    label: str,
) -> None:
    """
    Add AK or HI as a compact inset. Each inset is independently scaled to fit
    its box, which is the conventional visual treatment for US maps.
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
        linewidth=0.65,
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
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def make_choropleth(
    lower48: gpd.GeoDataFrame,
    alaska: gpd.GeoDataFrame,
    hawaii: gpd.GeoDataFrame,
    column: str,
    title: str,
    colorbar_label: str,
    norm: mpl.colors.Normalize,
    outdir: Path,
    stem: str,
) -> Dict[str, Path]:
    """
    Draw one standalone US choropleth with compact AK/HI insets.
    """
    cmap = mpl.colormaps["viridis"]

    fig = plt.figure(figsize=(11.6, 7.0))

    # Main lower-48 + DC map.
    ax = fig.add_axes([0.055, 0.11, 0.79, 0.80])
    lower48.plot(
        column=column,
        ax=ax,
        cmap=cmap,
        norm=norm,
        edgecolor="black",
        linewidth=0.65,
    )
    ax.set_axis_off()
    ax.set_title(title, pad=10, fontweight="bold", fontsize=TITLE_SIZE)

    # Conventional compact insets.
    _add_state_inset(
        fig, alaska,
        rect=[0.065, 0.115, 0.18, 0.18],
        column=column, cmap=cmap, norm=norm, label="AK",
    )
    _add_state_inset(
        fig, hawaii,
        rect=[0.255, 0.12, 0.10, 0.11],
        column=column, cmap=cmap, norm=norm, label="HI",
    )

    # Compact, readable vertical colorbar.
    cax = fig.add_axes([0.875, 0.23, 0.026, 0.54])
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.ax.tick_params(labelsize=COLORBAR_TICK_SIZE, width=0.9, length=4)
    cb.outline.set_linewidth(0.9)
    cb.set_label(
        colorbar_label,
        fontsize=COLORBAR_LABEL_SIZE,
        labelpad=12,
    )

    paths = _save_figure(fig, outdir, stem)
    plt.close(fig)
    return paths


# -----------------------------------------------------------------------------
# Distribution fitting
# -----------------------------------------------------------------------------

def _fit_state_prevalence(prevalence_percent: np.ndarray) -> pd.DataFrame:
    """
    Compare logit-normal, beta, and Gaussian fits using AIC.

    Input is prevalence in percentage points (e.g. 25.7), converted internally
    to proportions.
    """
    p = np.asarray(prevalence_percent, dtype=float) / 100.0

    if np.any((p <= 0.0) | (p >= 1.0)):
        raise ValueError(
            "State prevalences must lie strictly between 0 and 100 "
            "for beta/logit-normal fitting."
        )

    # Logit-normal.
    z = np.log(p / (1.0 - p))
    mu_z, sigma_z = stats.norm.fit(z)
    ll_logn = np.sum(
        stats.norm.logpdf(z, mu_z, sigma_z)
        - np.log(p)
        - np.log(1.0 - p)
    )
    aic_logn = 2 * 2 - 2 * ll_logn

    # Beta.
    alpha, beta, _, _ = stats.beta.fit(p, floc=0, fscale=1)
    ll_beta = np.sum(stats.beta.logpdf(p, alpha, beta))
    aic_beta = 2 * 2 - 2 * ll_beta

    # Gaussian on the raw proportion.
    mu, sigma = stats.norm.fit(p)
    ll_norm = np.sum(stats.norm.logpdf(p, mu, sigma))
    aic_norm = 2 * 2 - 2 * ll_norm

    result = pd.DataFrame(
        {
            "model": ["Logit-normal", "Beta", "Gaussian"],
            "AIC": [aic_logn, aic_beta, aic_norm],
        }
    ).sort_values("AIC", ignore_index=True)

    result["delta_AIC"] = result["AIC"] - result["AIC"].min()

    # Store parameters in attrs for plotting.
    result.attrs["logit_mu"] = float(mu_z)
    result.attrs["logit_sigma"] = float(sigma_z)
    result.attrs["beta_alpha"] = float(alpha)
    result.attrs["beta_beta"] = float(beta)
    result.attrs["normal_mu"] = float(mu)
    result.attrs["normal_sigma"] = float(sigma)

    return result


def _logitnormal_pdf(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """
    Logit-normal PDF on x in (0, 1).
    """
    z = np.log(x / (1.0 - x))
    return stats.norm.pdf(z, mu, sigma) / (x * (1.0 - x))


def make_distribution_figure(
    df: pd.DataFrame,
    prevalence_column: str,
    title: str,
    x_label: str,
    outdir: Path,
    stem: str,
    bins: np.ndarray,
    mark_states: Tuple[str, ...] = MARK_STATES,
) -> Tuple[Dict[str, Path], pd.DataFrame]:
    """
    Draw one standalone histogram with the best-fitting logit-normal density and
    dashed markers for PA, TX, and KY.

    The candidate model comparison is printed and returned.
    """
    values = df[prevalence_column].to_numpy(dtype=float)
    fits = _fit_state_prevalence(values)

    mu_z = fits.attrs["logit_mu"]
    sigma_z = fits.attrs["logit_sigma"]

    fig = plt.figure(figsize=(7.4, 4.9))
    ax = fig.add_axes([0.12, 0.15, 0.83, 0.77])

    # Histogram in density units.
    ax.hist(
        values,
        bins=bins,
        density=True,
        alpha=0.24,
        edgecolor="black",
        linewidth=1.1,
        label="State prevalence",
    )

    # Logit-normal density in density-per-percentage-point units.
    xmin = max(0.001, values.min() / 100.0 - 0.04)
    xmax = min(0.999, values.max() / 100.0 + 0.06)
    x = np.linspace(xmin, xmax, 900)
    y = _logitnormal_pdf(x, mu_z, sigma_z) / 100.0

    ax.plot(
        x * 100.0,
        y,
        linewidth=3.1,
        label="Logit-normal fit",
    )

    # State markers.
    state_values = {}
    for st in mark_states:
        row = df.loc[df["state"] == st]
        if row.empty:
            continue
        state_values[st] = float(row.iloc[0][prevalence_column])

    for st, val in state_values.items():
        ax.axvline(val, linestyle="--", linewidth=2.0)

    # Label placement:
    # PA is left of the central cluster; TX/KY may be nearly coincident for v1+.
    ymax = max(ax.get_ylim()[1], float(np.nanmax(y)) * 1.10)
    ax.set_ylim(0, ymax)

    if "PA" in state_values:
        ax.text(
            state_values["PA"] + 0.5,
            0.88 * ymax,
            f"PA  {state_values['PA']:.1f}%",
            rotation=90,
            va="top",
            ha="left",
            fontsize=15,
            fontweight="bold",
        )

    if "TX" in state_values:
        tx_val = state_values["TX"]
        tx_offset = 0.5

        # If KY is very close, separate the TX and KY text horizontally.
        if "KY" in state_values and abs(state_values["KY"] - tx_val) < 3.0:
            tx_offset = 0.9

        ax.text(
            tx_val + tx_offset,
            0.73 * ymax,
            f"TX  {tx_val:.1f}%",
            rotation=90,
            va="top",
            ha="left",
            fontsize=15,
            fontweight="bold",
        )

    if "KY" in state_values:
        ky_val = state_values["KY"]
        ky_offset = 0.5

        if "TX" in state_values and abs(ky_val - state_values["TX"]) < 3.0:
            ky_offset = 2.1

        ax.text(
            ky_val + ky_offset,
            0.56 * ymax,
            f"KY  {ky_val:.1f}%",
            rotation=90,
            va="top",
            ha="left",
            fontsize=15,
            fontweight="bold",
        )

    ax.set_xlabel(x_label, fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("Density", fontsize=AXIS_LABEL_SIZE)
    ax.set_title(title, fontsize=TITLE_SIZE)
    ax.tick_params(labelsize=TICK_SIZE, width=1.2, length=5)
    ax.legend(frameon=False, loc="upper right", fontsize=LEGEND_SIZE)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    paths = _save_figure(fig, outdir, stem)
    plt.close(fig)

    return paths, fits


# -----------------------------------------------------------------------------
# Main public function
# -----------------------------------------------------------------------------

def make_all_figures(
    df: pd.DataFrame,
    outdir: str | Path = "PLOTS",
    state_boundary_file: str | Path | None = None,
) -> Dict[str, Dict[str, Path]]:
    """
    Generate all five standalone figures.

    Parameters
    ----------
    df
        State-level DataFrame containing:
            state
            n_patients
            MATE-EHR v1
            MATE-EHR v1+
    outdir
        Output directory.
    state_boundary_file
        Optional local Census state-boundary shapefile/ZIP. If omitted, the
        script downloads and caches the official Census state ZIP automatically.

    Returns
    -------
    dict
        Dictionary mapping figure names to PNG/PDF/SVG output paths.
    """
    outdir = Path(outdir)
    d = _normalize_columns(df)

    lower48, alaska, hawaii = _make_state_geometry(d, state_boundary_file)

    # Global publication settings.
    plt.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.size": 16,
            "axes.titlesize": TITLE_SIZE,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "xtick.labelsize": TICK_SIZE,
            "ytick.labelsize": TICK_SIZE,
            "legend.fontsize": LEGEND_SIZE,
        }
    )

    outputs: Dict[str, Dict[str, Path]] = {}

    # 1. Cohort size map.
    n_norm = mpl.colors.Normalize(
        vmin=0,
        vmax=float(d["n_patients"].max()),
    )

    outputs["n_patients"] = make_choropleth(
        lower48,
        alaska,
        hawaii,
        column="n_patients",
        title="Lung transplant recipients represented in COSMOS",
        colorbar_label="Recipients, n",
        norm=n_norm,
        outdir=outdir,
        stem="MATE_choropleth_n_patients",
    )

    # Same prevalence scale for v1 and v1+.
    prevalence_norm = mpl.colors.Normalize(
        vmin=PREVALENCE_MAP_MIN,
        vmax=PREVALENCE_MAP_MAX,
    )

    # 2. MATE-EHR v1 prevalence map.
    outputs["v1_prevalence"] = make_choropleth(
        lower48,
        alaska,
        hawaii,
        column="MATE_v1_prevalence",
        title="MATE-EHR v1 prevalence by state",
        colorbar_label="Prevalence (%)",
        norm=prevalence_norm,
        outdir=outdir,
        stem="MATE_choropleth_v1_prevalence",
    )

    # 3. MATE-EHR v1+ prevalence map.
    outputs["v1plus_prevalence"] = make_choropleth(
        lower48,
        alaska,
        hawaii,
        column="MATE_v1plus_prevalence",
        title="MATE-EHR v1+ prevalence by state",
        colorbar_label="Prevalence (%)",
        norm=prevalence_norm,
        outdir=outdir,
        stem="MATE_choropleth_v1plus_prevalence",
    )

    # 4. MATE-EHR v1 state distribution.
    v1_paths, v1_fit = make_distribution_figure(
        d,
        prevalence_column="MATE_v1_prevalence",
        title="State-level variation in MATE-v1",
        x_label="MATE-v1 prevalence across states (%)",
        outdir=outdir,
        stem="MATE_v1_logitnormal_fit",
        bins=np.arange(10, 66, 5),
    )
    outputs["v1_distribution"] = v1_paths

    # 5. MATE-EHR v1+ state distribution.
    v1p_paths, v1p_fit = make_distribution_figure(
        d,
        prevalence_column="MATE_v1plus_prevalence",
        title="State-level variation in MATE-EHR v1+",
        x_label="MATE-EHR v1+ prevalence across states (%)",
        outdir=outdir,
        stem="MATE_v1plus_logitnormal_fit",
        bins=np.arange(15, 71, 5),
    )
    outputs["v1plus_distribution"] = v1p_paths

    # Console summary.
    print("\nMATE-EHR v1 distribution fit:")
    print(v1_fit.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    print("\nMATE-EHR v1+ distribution fit:")
    print(v1p_fit.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    print("\nState values marked on distribution figures:")
    for st in MARK_STATES:
        row = d.loc[d["state"] == st]
        if row.empty:
            continue
        r = row.iloc[0]
        print(
            f"{st}: "
            f"v1={r['MATE_v1_prevalence']:.2f}% "
            f"({int(r['MATE-EHR v1'])}/{int(r['n_patients'])}); "
            f"v1+={r['MATE_v1plus_prevalence']:.2f}% "
            f"({int(r['MATE-EHR v1+'])}/{int(r['n_patients'])})"
        )

    print("\nWrote figures to:", outdir.resolve())
    return outputs


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate the five standalone MATE state-level figures."
    )
    parser.add_argument(
        "csv",
        type=Path,
        help="CSV containing state, n_patients, MATE-EHR v1, MATE-EHR v1+.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("PLOTS"),
        help="Output directory (default: PLOTS).",
    )
    parser.add_argument(
        "--state-boundaries",
        type=Path,
        default=None,
        help=(
            "Optional local Census state-boundary ZIP/shapefile. "
            "If omitted, the Census file is downloaded and cached automatically."
        ),
    )
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    make_all_figures(
        df,
        outdir=args.outdir,
        state_boundary_file=args.state_boundaries,
    )


if __name__ == "__main__":
    main()
