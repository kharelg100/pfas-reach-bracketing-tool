# -*- coding: utf-8 -*-
"""
PFAS Reach Bracketing Tool — Python Edition
===========================================
Version: 1.1.1-public-release
Author: Dr. Gehendra Kharel, STREAM Lab (Sustainable Tools for Risk Evaluation and Climate-Water Modeling), Texas Christian University
Archived software DOI: https://doi.org/10.5281/zenodo.20464907
Concept DOI recommended for citation: https://doi.org/10.5281/zenodo.18937739

Purpose
-------
Machine-readable reach experiment analysis for PFAS surface-water reconnaissance.
The tool converts long-format PFAS results and a reach_experiment_table.csv file
into censoring-aware site totals, reach contrasts, sensitivity diagnostics, and
figures. It is aligned with the revised Water Research X Research Letter
manuscript, in which:

  * ΣPFAS40 = summed concentration of the 40 target PFAS analytes.
  * Non-detects are handled by multiple imputation (MI), with 5,000 replicates
    by default for manuscript-level analyses.
  * Analytes with >=2 detections are imputed from a fitted lognormal
    distribution truncated at [0, RL].
  * Analytes with <2 detections use Uniform(0, RL), avoiding arbitrary
    single-detect variance assumptions.
  * Reach experiments are read from a machine-readable table rather than hard-
    coded into the script.
  * f_trib is reported only as an apparent concentration-space screening
    diagnostic, not as a flow fraction, source-apportionment fraction, or PFAS
    load fraction.

Inputs
------
Place these in ./inputs/ relative to this script, or edit CONFIG below:

Required:
  1. pfas_long_filtered_standardized.csv
     Required columns:
       site_id, campaign, analyte_std, result_value, detect_flag, rl
     Optional columns retained if present:
       qualifier, campaign_date, client_sample_id, sample_id, result_text, unit

Optional:
  2. sites_metadata.csv OR sites_metadata_gk2v3.csv
     Required/recognized columns:
       Site ID or site_id, Lat/latitude, Lon/longitude, Reach_group, Bracket_role

  3. reach_experiment_table.csv
     Required columns:
       experiment_id, experiment_type, downstream_site, upstream_sites
     Optional but recommended:
       reach_group, tributary_sites, baseline_sites, all_sites, experiment_name, notes

Outputs
-------
Written to ./outputs/:
  * PFAS_Reach_Bracketing_Results.xlsx
  * validation_report.txt
  * SiteTotals.csv
  * ReachExperiments.csv
  * Sensitivity_RL2.csv
  * Convergence.csv
  * ImputationParams.csv
  * Fig_SiteConcentrations.png/svg
  * Fig_Sensitivity_MI_vs_RL2.png/svg
  * Fig_ConvergenceDiag.png/svg
  * Fig_DetectionHeatmap.png/svg
  * Fig_NetworkMap.png/svg (only when coordinates are available)


Methodological references
-------------------------
This tool builds on established methods and guidance for censored environmental
chemistry data, multiple imputation, reproducible workflows, PFAS LC-MS/MS
analysis, and PFAS stream loading/source-screening:

  [1] Helsel, D.R. (2012). Statistics for Censored Environmental Data Using
      Minitab and R, 2nd ed. Wiley.
  [2] Rubin, D.B. (1987). Multiple Imputation for Nonresponse in Surveys.
      Wiley.
  [3] Lubin, J.H. et al. (2004). Epidemiologic evaluation of measurement data
      in the presence of detection limits. Environmental Health Perspectives
      112(17), 1691-1696. doi:10.1289/ehp.7199.
  [4] Helsel, D.R. (2006). Fabricating data: How substituting values for
      nondetects can ruin results. Chemosphere 65(11), 2434-2439.
      doi:10.1016/j.chemosphere.2006.04.051.
  [5] Antweiler, R.C. & Taylor, H.E. (2008). Evaluation of statistical
      treatments of left-censored environmental data using coincident
      uncensored data sets: I. Summary statistics. Environmental Science &
      Technology 42(10), 3732-3738. doi:10.1021/es071301c.
  [6] Wilkinson, M.D. et al. (2016). The FAIR Guiding Principles for scientific
      data management and stewardship. Scientific Data 3, 160018.
      doi:10.1038/sdata.2016.18.
  [7] Goble, C.A. et al. (2025). Applying the FAIR principles to computational
      workflows. Scientific Data 12, 337. doi:10.1038/s41597-025-04451-9.
  [8] U.S. EPA (2024). Method 1633: Analysis of PFAS in aqueous, solid,
      biosolids, and tissue samples by LC-MS/MS, Revision A.
  [9] Spearman, C. (1904). The proof and measurement of association between
      two things. American Journal of Psychology 15(1), 72-101.
      doi:10.2307/1412159.
  [10] Woodward, E.E. et al. (2024). Using a time-of-travel sampling approach
       to quantify PFAS stream loading and source inputs in a mixed-source,
       urban catchment. ACS ES&T Water 4(10), 4356-4367.
       doi:10.1021/acsestwater.4c00288.
  [11] Asare, P.T. et al. (2026). Handling left-censored PFAS data in
       surface-water reconnaissance: A reproducible workflow applied to ten
       sites in the Trinity River headwaters, Texas. PLOS Water, e0000554.
       doi:10.1371/journal.pwat.0000554.

Dependencies
------------
Standard Anaconda stack:
  numpy, pandas, scipy, matplotlib, openpyxl

Public-release notes
--------------------
This script performs local computation only and does not transmit data. It is
intended for transparent screening and prioritization, not regulatory compliance
or source apportionment. Concentration-space f_trib diagnostics require no flow
data and must not be interpreted as discharge or load fractions.

License
-------
License: MIT. See the accompanying LICENSE file for the full license text.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import json
import math
import sys
import warnings
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
from scipy.stats import norm, spearmanr

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────
DEFAULT_N_IMPUTE = 5000
DEFAULT_SEED = 20260227
DEFAULT_DPI = 300

ANALYTE_TOTAL_LABEL = "ΣPFAS40"
UNIT_LABEL = "ng L-1"

# Recognized filename aliases for user convenience
PFAS_FILENAMES = ["pfas_long_filtered_standardized.csv", "pfas_results.csv", "PFAS_results.csv"]
SITES_FILENAMES = ["sites_metadata.csv", "sites_metadata_gk2v3.csv", "Sampling_Sites_GK2v3.csv"]
REACH_FILENAMES = ["reach_experiment_table.csv", "GK2v3_reach_experiments.csv"]

# Colors chosen for high contrast and print readability
COL_UP = "#2563eb"
COL_TRIB = "#6b21a8"
COL_DOWN = "#ea580c"
COL_GRAY = "#6b7280"
COL_DARK = "#111827"
COL_LIGHT = "#f3f4f6"

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
})


@dataclass
class RunConfig:
    base_dir: Path
    input_dir: Path
    output_dir: Path
    n_impute: int = DEFAULT_N_IMPUTE
    seed: int = DEFAULT_SEED
    dpi: int = DEFAULT_DPI
    require_exact_40: bool = False


def _print(msg: str) -> None:
    print(msg, flush=True)


# ──────────────────────────────────────────────────────────────────────────────
# PATHS, NORMALIZATION, AND VALIDATION
# ──────────────────────────────────────────────────────────────────────────────
def find_first_existing(directory: Path, candidates: Iterable[str]) -> Optional[Path]:
    for name in candidates:
        p = directory / name
        if p.exists():
            return p
    return None


def normalize_site_id(x) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    s = s.replace(" ", "")
    # Convert TW005 to TW05 only if user supplied a leading-zero variant
    if s.startswith("TW") and s[2:].isdigit():
        # retain two or more digits when originally present; TW251 remains TW251
        n = int(s[2:])
        if n < 100:
            return f"TW{n:02d}"
    return s


def split_pipe(x) -> List[str]:
    if pd.isna(x) or str(x).strip() == "":
        return []
    return [normalize_site_id(s) for s in str(x).split("|") if str(s).strip()]


def fmt_interval(med: float, lo: float, hi: float, dec: int = 1) -> str:
    if not np.isfinite(med):
        return ""
    return f"{med:.{dec}f} [{lo:.{dec}f}, {hi:.{dec}f}]"


def validate_required_columns(df: pd.DataFrame, required: Iterable[str], name: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")


def read_pfas(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(path)
    # Normalize column names only for spaces/case aliases; preserve analyte_std and user columns
    rename = {}
    for c in df.columns:
        c_clean = c.strip()
        if c != c_clean:
            rename[c] = c_clean
    df = df.rename(columns=rename)

    required = ["site_id", "campaign", "analyte_std", "result_value", "detect_flag", "rl"]
    validate_required_columns(df, required, "PFAS results CSV")

    raw_rows = len(df)
    df = df.copy()
    df["site_id"] = df["site_id"].map(normalize_site_id)
    df["campaign"] = df["campaign"].astype(str).str.strip()
    df["analyte_std"] = df["analyte_std"].astype(str).str.strip()
    df["result_value"] = pd.to_numeric(df["result_value"], errors="coerce")
    df["rl"] = pd.to_numeric(df["rl"], errors="coerce")
    df["detect_flag"] = pd.to_numeric(df["detect_flag"], errors="coerce").fillna(0).astype(int)

    # Results below RL but flagged as detect are allowed; result_value missing with detect=1 is invalid.
    problems = []
    invalid_detect = df[~df["detect_flag"].isin([0, 1])]
    if len(invalid_detect):
        problems.append(f"detect_flag contains values other than 0/1 in {len(invalid_detect)} rows")
    invalid_rl = df[~np.isfinite(df["rl"]) | (df["rl"] <= 0)]
    if len(invalid_rl):
        problems.append(f"rl is missing or <=0 in {len(invalid_rl)} rows")
    invalid_detect_values = df[(df["detect_flag"] == 1) & (~np.isfinite(df["result_value"]))]
    if len(invalid_detect_values):
        problems.append(f"detect_flag=1 but result_value missing in {len(invalid_detect_values)} rows")
    invalid_nd_values = df[(df["detect_flag"] == 0) & (df["result_value"].notna()) & (df["result_value"] < 0)]
    if len(invalid_nd_values):
        problems.append(f"negative result_value in {len(invalid_nd_values)} non-detect rows")

    duplicate_key = ["campaign", "site_id", "analyte_std"]
    dups = df[df.duplicated(duplicate_key, keep=False)].sort_values(duplicate_key)
    if len(dups):
        # Keep first but report; duplicate rows compromise reproducibility if not noted.
        problems.append(f"duplicate campaign-site-analyte rows found: {len(dups)} rows; keeping first occurrence")
        df = df.drop_duplicates(duplicate_key, keep="first")

    if problems:
        msg = "PFAS input validation warnings:\n  - " + "\n  - ".join(problems)
        _print(msg)

    validation = pd.DataFrame({
        "check": ["raw_rows", "deduplicated_rows", "n_campaigns", "n_sites", "n_analytes", "warnings"],
        "value": [raw_rows, len(df), df["campaign"].nunique(), df["site_id"].nunique(), df["analyte_std"].nunique(), "; ".join(problems) if problems else "none"],
    })
    return df, validation


def read_sites(path: Optional[Path], df: pd.DataFrame) -> pd.DataFrame:
    if path is None or not path.exists():
        sites = pd.DataFrame({
            "site_id": sorted(df["site_id"].unique()),
            "Lat": np.nan,
            "Lon": np.nan,
            "Reach_group": "",
            "Bracket_role": "",
        })
        return sites

    sites = pd.read_csv(path)
    # Handle common aliases
    rename_map = {}
    for c in sites.columns:
        cl = c.lower().strip()
        if cl in ["site id", "site_id", "site"]:
            rename_map[c] = "site_id"
        elif cl in ["lat", "latitude"]:
            rename_map[c] = "Lat"
        elif cl in ["lon", "long", "longitude"]:
            rename_map[c] = "Lon"
        elif cl in ["reach_group", "reach group"]:
            rename_map[c] = "Reach_group"
        elif cl in ["bracket_role", "bracket role", "role"]:
            rename_map[c] = "Bracket_role"
    sites = sites.rename(columns=rename_map)
    if "site_id" not in sites.columns:
        raise ValueError("Sites metadata must include 'Site ID' or 'site_id'.")
    sites["site_id"] = sites["site_id"].map(normalize_site_id)
    for c in ["Lat", "Lon", "Reach_group", "Bracket_role"]:
        if c not in sites.columns:
            sites[c] = np.nan if c in ["Lat", "Lon"] else ""
    sites["Lat"] = pd.to_numeric(sites["Lat"], errors="coerce")
    sites["Lon"] = pd.to_numeric(sites["Lon"], errors="coerce")
    sites = sites.drop_duplicates("site_id", keep="first")
    return sites


def read_reach(path: Optional[Path]) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    reach = pd.read_csv(path)
    required = ["experiment_id", "experiment_type", "downstream_site", "upstream_sites"]
    validate_required_columns(reach, required, "Reach experiment CSV")
    for c in ["reach_group", "experiment_name", "tributary_sites", "baseline_sites", "all_sites", "notes"]:
        if c not in reach.columns:
            reach[c] = ""
    for c in ["downstream_site"]:
        reach[c] = reach[c].map(normalize_site_id)
    for c in ["upstream_sites", "tributary_sites", "baseline_sites", "all_sites"]:
        reach[c] = reach[c].fillna("").astype(str).apply(lambda s: "|".join(split_pipe(s)))
    reach["experiment_type"] = reach["experiment_type"].fillna("").astype(str).str.strip()
    reach["experiment_id"] = reach["experiment_id"].fillna("").astype(str).str.strip()
    return reach


# ──────────────────────────────────────────────────────────────────────────────
# MULTIPLE IMPUTATION
# ──────────────────────────────────────────────────────────────────────────────
def build_wide(df: pd.DataFrame):
    analytes = sorted(df["analyte_std"].unique())
    index_cols = ["campaign", "site_id"]
    value_wide = df.pivot_table(index=index_cols, columns="analyte_std", values="result_value", aggfunc="first").reindex(columns=analytes)
    rl_wide = df.pivot_table(index=index_cols, columns="analyte_std", values="rl", aggfunc="first").reindex(columns=analytes)
    detect_wide = df.pivot_table(index=index_cols, columns="analyte_std", values="detect_flag", aggfunc="first").reindex(columns=analytes)
    # Missing rows after pivot are treated as non-detect only if RL available; otherwise invalid for that analyte/sample.
    detect_wide = detect_wide.fillna(0).astype(int)
    return analytes, value_wide, rl_wide, detect_wide


def fit_lognormal_params(value_wide, detect_wide, analytes) -> pd.DataFrame:
    rows = []
    for a in analytes:
        det_vals = value_wide[a][detect_wide[a] == 1].dropna().astype(float).values
        det_vals = det_vals[np.isfinite(det_vals) & (det_vals > 0)]
        if len(det_vals) >= 2:
            logs = np.log(det_vals)
            mu = float(logs.mean())
            sigma = float(logs.std(ddof=1))
            if not np.isfinite(sigma) or sigma < 1e-6:
                sigma = 0.5
            mode = "truncated_lognormal"
        else:
            mu = np.nan
            sigma = np.nan
            mode = "uniform_0_RL"
        rows.append({"analyte_std": a, "n_detect": int(len(det_vals)), "imputation_mode": mode, "mu_log": mu, "sigma_log": sigma})
    return pd.DataFrame(rows)


def run_multiple_imputation(df: pd.DataFrame, config: RunConfig):
    analytes, value_wide, rl_wide, detect_wide = build_wide(df)
    if config.require_exact_40 and len(analytes) != 40:
        raise ValueError(f"Expected exactly 40 analytes; found {len(analytes)}")

    params = fit_lognormal_params(value_wide, detect_wide, analytes)
    param_map = {r["analyte_std"]: r for _, r in params.iterrows()}

    n_samples = len(value_wide.index)
    totals = np.zeros((n_samples, config.n_impute), dtype=np.float64)
    detect_only = np.zeros(n_samples, dtype=np.float64)
    rl2_totals = np.zeros(n_samples, dtype=np.float64)
    rng = np.random.default_rng(config.seed)

    for a in analytes:
        obs = value_wide[a].to_numpy(dtype=float)
        rl = rl_wide[a].to_numpy(dtype=float)
        det = detect_wide[a].to_numpy(dtype=int)
        det_mask = (det == 1) & np.isfinite(obs) & (obs >= 0)
        nd_mask = ~det_mask

        vals = np.zeros((n_samples, config.n_impute), dtype=np.float64)
        if det_mask.any():
            vals[det_mask, :] = obs[det_mask, None]
            detect_only[det_mask] += obs[det_mask]
            rl2_totals[det_mask] += obs[det_mask]

        if nd_mask.any():
            upper = rl[nd_mask]
            invalid = ~np.isfinite(upper) | (upper <= 0)
            if invalid.any():
                upper = upper.copy()
                upper[invalid] = 0.0
            rl2_totals[nd_mask] += upper / 2.0
            pinfo = param_map[a]
            if pinfo["imputation_mode"] == "truncated_lognormal" and np.all(upper >= 0):
                mu = float(pinfo["mu_log"])
                sigma = float(pinfo["sigma_log"])
                positive = upper > 0
                draw = np.zeros((nd_mask.sum(), config.n_impute), dtype=np.float64)
                if positive.any():
                    logu = np.log(np.maximum(upper[positive], 1e-12))
                    p = norm.cdf((logu - mu) / sigma)
                    p = np.clip(p, 1e-12, 1.0)
                    u = rng.random((positive.sum(), config.n_impute)) * p[:, None]
                    z = mu + sigma * norm.ppf(u)
                    d = np.exp(z)
                    d = np.minimum(d, upper[positive, None])
                    d = np.maximum(d, 0.0)
                    draw[positive, :] = d
            else:
                draw = rng.random((nd_mask.sum(), config.n_impute)) * upper[:, None]
            vals[nd_mask, :] = draw
        totals += vals

    sample_index = value_wide.index
    site_stats = pd.DataFrame({
        "campaign": [i[0] for i in sample_index],
        "site_id": [i[1] for i in sample_index],
        "Sigma_detect_only": detect_only,
        "Sigma_RL2": rl2_totals,
        "Sigma_MI_median": np.median(totals, axis=1),
        "Sigma_MI_p2.5": np.quantile(totals, 0.025, axis=1),
        "Sigma_MI_p97.5": np.quantile(totals, 0.975, axis=1),
    })
    det_counts = detect_wide.sum(axis=1).astype(int).to_numpy()
    site_stats["n_detect"] = det_counts
    site_stats["n_analytes"] = len(analytes)
    site_stats["detect_freq"] = [f"{n}/{len(analytes)} ({n/len(analytes)*100:.1f}%)" for n in det_counts]
    site_stats["Sigma_MI_median_[95%II]"] = site_stats.apply(lambda r: fmt_interval(r.Sigma_MI_median, r["Sigma_MI_p2.5"], r["Sigma_MI_p97.5"]), axis=1)
    return site_stats, totals, sample_index, analytes, params, value_wide, rl_wide, detect_wide


# ──────────────────────────────────────────────────────────────────────────────
# DIAGNOSTICS AND REACH ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────
def convergence_diagnostics(totals: np.ndarray, n_impute: int) -> pd.DataFrame:
    ref_med = np.median(totals, axis=1)
    ref_lo = np.quantile(totals, 0.025, axis=1)
    ref_hi = np.quantile(totals, 0.975, axis=1)
    rows = []
    for N in [250, 500, 1000, 2000, 5000, 7500, 10000]:
        if N >= n_impute:
            continue
        medN = np.median(totals[:, :N], axis=1)
        loN = np.quantile(totals[:, :N], 0.025, axis=1)
        hiN = np.quantile(totals[:, :N], 0.975, axis=1)
        rows.append({
            "N": N,
            "max_abs_diff_median_ngL": float(np.max(np.abs(medN - ref_med))),
            "median_abs_diff_median_ngL": float(np.median(np.abs(medN - ref_med))),
            "max_abs_diff_II95_lo_ngL": float(np.max(np.abs(loN - ref_lo))),
            "max_abs_diff_II95_hi_ngL": float(np.max(np.abs(hiN - ref_hi))),
        })
    return pd.DataFrame(rows)


def sensitivity_rl2(site_stats: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sens = site_stats.copy()
    sens["AbsDiff_MI_minus_RL2"] = (sens["Sigma_MI_median"] - sens["Sigma_RL2"]).abs()
    sens["PctDiff_MI_minus_RL2"] = np.where(sens["Sigma_MI_median"] > 0, (sens["Sigma_MI_median"] - sens["Sigma_RL2"]) / sens["Sigma_MI_median"] * 100, np.nan)
    if len(sens) > 1:
        rho, pval = spearmanr(sens["Sigma_MI_median"], sens["Sigma_RL2"])
    else:
        rho, pval = np.nan, np.nan
    summary = pd.DataFrame({
        "metric": ["Spearman_rho_MI_vs_RL2", "Spearman_pvalue", "Max_abs_MI_minus_RL2_ngL", "Median_abs_MI_minus_RL2_ngL"],
        "value": [rho, pval, sens["AbsDiff_MI_minus_RL2"].max(), sens["AbsDiff_MI_minus_RL2"].median()],
    })
    return sens, summary


def build_distribution_lookup(sample_index, totals):
    return {(sample_index[i][0], sample_index[i][1]): totals[i, :] for i in range(len(sample_index))}


def dist_summary(arr: np.ndarray, dec: int = 1) -> Tuple[float, float, float, str]:
    if arr is None or len(arr) == 0:
        return np.nan, np.nan, np.nan, ""
    med = float(np.nanmedian(arr))
    lo = float(np.nanquantile(arr, 0.025))
    hi = float(np.nanquantile(arr, 0.975))
    return med, lo, hi, fmt_interval(med, lo, hi, dec)


def compute_ftrib(up_dist: np.ndarray, trib_dist: np.ndarray, dn_dist: np.ndarray) -> Tuple[float, float, float, str]:
    """Compute concentration-space f_trib when downstream lies between upstream and tributary.

    Returns median, lo, hi, status string. Does not clip invalid values.
    """
    denom = trib_dist - up_dist
    with np.errstate(divide="ignore", invalid="ignore"):
        f = (dn_dist - up_dist) / denom
    valid = np.isfinite(f) & (denom > 0) & (f >= 0) & (f <= 1)
    valid_rate = float(np.mean(valid)) if len(f) else 0.0
    if valid_rate < 0.95:
        # Identify common cause for reviewer-proof table notes.
        med_up, med_trib, med_dn = np.nanmedian(up_dist), np.nanmedian(trib_dist), np.nanmedian(dn_dist)
        if med_dn > med_trib:
            reason = "n/a: C_dn > C_trib"
        elif med_dn < med_up:
            reason = "n/a: C_dn < C_up"
        elif med_trib <= med_up:
            reason = "n/a: C_trib <= C_up"
        else:
            reason = "n/a: outside two-endmember range"
        return np.nan, np.nan, np.nan, f"{reason}; valid={valid_rate:.2f}"
    f_valid = f[valid]
    med, lo, hi, txt = dist_summary(f_valid, dec=3)
    return med, lo, hi, f"f_trib={txt}; valid={valid_rate:.2f}"


def analyze_reaches(reach: pd.DataFrame, dist_lookup: Dict[Tuple[str, str], np.ndarray]) -> pd.DataFrame:
    if reach.empty:
        return pd.DataFrame()
    campaigns = sorted({k[0] for k in dist_lookup.keys()})
    rows = []
    for _, r in reach.iterrows():
        exp_id = str(r["experiment_id"])
        exp_type = str(r["experiment_type"]).strip()
        dn = normalize_site_id(r["downstream_site"])
        ups = split_pipe(r.get("upstream_sites", ""))
        tribs = split_pipe(r.get("tributary_sites", ""))
        baselines = split_pipe(r.get("baseline_sites", ""))
        all_sites = split_pipe(r.get("all_sites", "")) or sorted(set([dn] + ups + tribs + baselines))
        for camp in campaigns:
            present_sites = {s for (c, s) in dist_lookup.keys() if c == camp}
            complete = set(all_sites).issubset(present_sites)
            if not complete:
                continue
            dn_dist = dist_lookup.get((camp, dn))
            if dn_dist is None:
                continue
            primary_up = ups[0] if ups else None
            up_dist = dist_lookup.get((camp, primary_up)) if primary_up else None
            delta = None
            delta_summary = ""
            delta_med = np.nan
            delta_lo = np.nan
            delta_hi = np.nan
            if up_dist is not None:
                delta = dn_dist - up_dist
                delta_med, delta_lo, delta_hi, delta_summary = dist_summary(delta, dec=1)
            diag_type = ""
            f_med = f_lo = f_hi = np.nan
            f_status = ""
            if exp_type == "tributary_increment" and up_dist is not None and tribs:
                trib_dist = dist_lookup.get((camp, tribs[0]))
                if trib_dist is not None:
                    f_med, f_lo, f_hi, f_status = compute_ftrib(up_dist, trib_dist, dn_dist)
                    diag_type = "f_trib_concentration_space"
            elif exp_type == "confluence_mix" and len(ups + tribs) >= 2:
                endmember_sites = ups + tribs
                dists = [dist_lookup.get((camp, u)) for u in endmember_sites]
                if all(d is not None for d in dists):
                    lo_end = np.minimum.reduce(dists)
                    hi_end = np.maximum.reduce(dists)
                    outside = (dn_dist < lo_end) | (dn_dist > hi_end)
                    outside_prob = float(np.mean(outside))
                    with np.errstate(divide="ignore", invalid="ignore"):
                        f_mix = (dn_dist - lo_end) / (hi_end - lo_end)
                    valid = np.isfinite(f_mix) & (f_mix >= 0) & (f_mix <= 1)
                    if valid.any():
                        f_med, f_lo, f_hi, f_txt = dist_summary(f_mix[valid], dec=3)
                    else:
                        f_txt = ""
                    f_status = f"outside_range_prob={outside_prob:.2f}; apparent_mix_fraction={f_txt if f_txt else 'n/a'}; endmembers={'|'.join(endmember_sites)}"
                    diag_type = "confluence_range_check"
            else:
                if exp_type == "source_reach":
                    diag_type = "linear_segment_delta"
                    f_status = "not applicable to source_reach; use ΔΣPFAS40"
                elif exp_type == "downstream_propagation":
                    diag_type = "downstream_propagation_delta"
                    f_status = "not applicable to downstream_propagation; use ΔΣPFAS40"

            base_record = {
                "campaign": camp,
                "experiment_id": exp_id,
                "reach_group": r.get("reach_group", ""),
                "experiment_type": exp_type,
                "experiment_name": r.get("experiment_name", ""),
                "contrast_role": "primary",
                "downstream_site": dn,
                "upstream_sites": "|".join(ups),
                "tributary_sites": "|".join(tribs),
                "baseline_sites": "|".join(baselines),
                "all_sites": "|".join(all_sites),
                "complete_required_sites": complete,
                "delta_SigmaPFAS40_median": delta_med,
                "delta_SigmaPFAS40_p2.5": delta_lo,
                "delta_SigmaPFAS40_p97.5": delta_hi,
                "delta_SigmaPFAS40_[95%II]": delta_summary,
                "diagnostic_type": diag_type,
                "f_trib_or_mix_median": f_med,
                "f_trib_or_mix_p2.5": f_lo,
                "f_trib_or_mix_p97.5": f_hi,
                "diagnostic_note": f_status,
            }
            rows.append(base_record)

            # Optional cumulative/baseline context rows. These are useful for
            # repeated sentinel comparisons (e.g., TW08-TW09) and for making the
            # conventional-vs-reach interpretation in the manuscript auditable.
            for b in baselines:
                bdist = dist_lookup.get((camp, b))
                if bdist is None:
                    continue
                bdelta = dn_dist - bdist
                b_med, b_lo, b_hi, b_txt = dist_summary(bdelta, dec=1)
                rows.append({
                    **base_record,
                    "experiment_id": f"{exp_id}__baseline_{b}",
                    "experiment_type": f"{exp_type}_baseline",
                    "experiment_name": f"{r.get('experiment_name', '')} — cumulative vs baseline {b}",
                    "contrast_role": "cumulative_vs_baseline",
                    "upstream_sites": b,
                    "tributary_sites": "",
                    "delta_SigmaPFAS40_median": b_med,
                    "delta_SigmaPFAS40_p2.5": b_lo,
                    "delta_SigmaPFAS40_p97.5": b_hi,
                    "delta_SigmaPFAS40_[95%II]": b_txt,
                    "diagnostic_type": "baseline_context_delta",
                    "f_trib_or_mix_median": np.nan,
                    "f_trib_or_mix_p2.5": np.nan,
                    "f_trib_or_mix_p97.5": np.nan,
                    "diagnostic_note": "cumulative downstream-minus-baseline contrast; not a mixing fraction",
                })
    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURES
# ──────────────────────────────────────────────────────────────────────────────
def save_fig(fig, out_dir: Path, stem: str, dpi: int) -> None:
    for ext in ["png", "svg"]:
        fig.savefig(out_dir / f"{stem}.{ext}", dpi=dpi, bbox_inches="tight", facecolor="white")


def role_color(role: str) -> str:
    r = str(role).lower() if role is not None else ""
    if "tributary" in r:
        return COL_TRIB
    if "upstream" in r or "baseline" in r or "pre-source" in r:
        return COL_UP
    if "downstream" in r or "integrator" in r or "post" in r or "outlet" in r or "confluence" in r:
        return COL_DOWN
    return COL_GRAY


def make_site_concentration_figure(site_stats: pd.DataFrame, sites: pd.DataFrame, out_dir: Path, dpi: int) -> None:
    plot = site_stats.merge(sites[["site_id", "Reach_group", "Bracket_role"]], on="site_id", how="left", suffixes=("", "_site"))
    plot = plot.sort_values("Sigma_MI_median", ascending=False)
    fig, ax = plt.subplots(figsize=(max(9, 0.42 * len(plot)), 5.6))
    x = np.arange(len(plot))
    y = plot["Sigma_MI_median"].to_numpy()
    yerr = np.vstack([y - plot["Sigma_MI_p2.5"].to_numpy(), plot["Sigma_MI_p97.5"].to_numpy() - y])
    colors = [role_color(r) for r in plot.get("Bracket_role", [""] * len(plot))]
    ax.bar(x, y, color=colors, edgecolor="#222222", linewidth=0.3)
    ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor="#111111", elinewidth=0.8, capsize=2)
    if y.max() / max(y[y > 0].min(), 1e-9) > 50:
        ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n{s}" for c, s in zip(plot["campaign"], plot["site_id"])], rotation=90, fontsize=7)
    ax.set_ylabel(f"{ANALYTE_TOTAL_LABEL} ({UNIT_LABEL})")
    ax.set_title("Censoring-aware site totals")
    patches = [mpatches.Patch(color=COL_UP, label="Upstream/baseline"), mpatches.Patch(color=COL_TRIB, label="Tributary"), mpatches.Patch(color=COL_DOWN, label="Downstream/integrator")]
    ax.legend(handles=patches, loc="upper right", fontsize=8, frameon=False)
    ax.grid(axis="y", alpha=0.25)
    save_fig(fig, out_dir, "Fig_SiteConcentrations", dpi)
    plt.close(fig)


def make_sensitivity_figure(sens: pd.DataFrame, out_dir: Path, dpi: int) -> None:
    plot = sens.sort_values("Sigma_MI_median")
    fig, ax = plt.subplots(figsize=(7.5, max(5.0, 0.28 * len(plot))))
    y = np.arange(len(plot))
    ax.hlines(y, plot["Sigma_RL2"], plot["Sigma_MI_median"], color="#9ca3af", lw=1)
    ax.scatter(plot["Sigma_RL2"], y, label="RL/2", color="#6b7280", s=22)
    ax.scatter(plot["Sigma_MI_median"], y, label="MI median", color="#2563eb", s=22)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c} {s}" for c, s in zip(plot["campaign"], plot["site_id"])], fontsize=7)
    ax.set_xlabel(f"{ANALYTE_TOTAL_LABEL} ({UNIT_LABEL})")
    ax.set_title("MI versus RL/2 substitution")
    ax.legend(frameon=False)
    ax.grid(axis="x", alpha=0.2)
    save_fig(fig, out_dir, "Fig_Sensitivity_MI_vs_RL2", dpi)
    plt.close(fig)


def make_convergence_figure(conv: pd.DataFrame, out_dir: Path, dpi: int) -> None:
    if conv.empty:
        return
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(conv["N"], conv["max_abs_diff_median_ngL"], marker="o", label="Median")
    ax.plot(conv["N"], conv["max_abs_diff_II95_lo_ngL"], marker="o", label="95% II lower")
    ax.plot(conv["N"], conv["max_abs_diff_II95_hi_ngL"], marker="o", label="95% II upper")
    ax.set_xlabel("Number of MI replicates used")
    ax.set_ylabel(f"Maximum absolute difference vs full run ({UNIT_LABEL})")
    ax.set_title("Monte Carlo convergence diagnostic")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    save_fig(fig, out_dir, "Fig_ConvergenceDiag", dpi)
    plt.close(fig)


def make_detection_heatmap(df: pd.DataFrame, out_dir: Path, dpi: int) -> None:
    mat = df.pivot_table(index=["campaign", "site_id"], columns="analyte_std", values="detect_flag", aggfunc="max").fillna(0).astype(int)
    # Keep detected analytes first for readability
    cols = mat.sum(axis=0).sort_values(ascending=False).index.tolist()
    mat = mat[cols]
    labels = [f"{c} {s}" for c, s in mat.index]
    fig, ax = plt.subplots(figsize=(max(9, 0.23 * len(cols)), max(4.5, 0.22 * len(labels))))
    ax.imshow(mat.values, aspect="auto", interpolation="nearest", cmap=mpl.colors.ListedColormap(["#f3f4f6", "#111827"]))
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, rotation=90, fontsize=6)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title("Detected PFAS analytes by site-campaign sample")
    ax.set_xlabel("Analyte")
    ax.set_ylabel("Sample")
    save_fig(fig, out_dir, "Fig_DetectionHeatmap", dpi)
    plt.close(fig)


def make_network_map(site_stats: pd.DataFrame, sites: pd.DataFrame, out_dir: Path, dpi: int) -> None:
    if not {"Lat", "Lon"}.issubset(sites.columns) or sites[["Lat", "Lon"]].isna().all().any():
        return
    # Plot maximum/only concentration by site for overview
    site_max = site_stats.groupby("site_id", as_index=False)["Sigma_MI_median"].max()
    plot = sites.merge(site_max, on="site_id", how="left")
    if plot["Lat"].isna().all() or plot["Lon"].isna().all():
        return
    fig, ax = plt.subplots(figsize=(7, 6))
    vals = plot["Sigma_MI_median"].fillna(0).to_numpy()
    positive = vals[vals > 0]
    norm = LogNorm(vmin=max(positive.min(), 1) if len(positive) else 1, vmax=max(vals.max(), 2))
    sc = ax.scatter(plot["Lon"], plot["Lat"], c=np.maximum(vals, 1), norm=norm, s=90, cmap="viridis", edgecolor="black", linewidth=0.5)
    for _, row in plot.iterrows():
        ax.text(row["Lon"], row["Lat"], f" {row['site_id']}", fontsize=7, va="center")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Site network overview")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(f"Max {ANALYTE_TOTAL_LABEL} ({UNIT_LABEL})")
    ax.grid(alpha=0.25)
    save_fig(fig, out_dir, "Fig_NetworkMap", dpi)
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# EXPORTS
# ──────────────────────────────────────────────────────────────────────────────
def write_excel(out_path: Path, sheets: Dict[str, pd.DataFrame], readme_lines: List[str]) -> None:
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        pd.DataFrame({"README": readme_lines}).to_excel(writer, sheet_name="README", index=False)
        for name, table in sheets.items():
            safe_name = name[:31]
            table.to_excel(writer, sheet_name=safe_name, index=False)
        # Light formatting with openpyxl through writer internals (output only; not part of analysis)
        wb = writer.book
        for ws in wb.worksheets:
            ws.freeze_panes = "A2"
            for col in ws.columns:
                max_len = 0
                col_letter = col[0].column_letter
                for cell in col[:200]:
                    try:
                        max_len = max(max_len, len(str(cell.value)) if cell.value is not None else 0)
                    except Exception:
                        pass
                ws.column_dimensions[col_letter].width = min(max(max_len + 2, 10), 42)


def write_validation_report(path: Path, lines: List[str]) -> None:
    path.write_text("\n".join(lines), encoding="utf-8")


def export_csvs(out_dir: Path, tables: Dict[str, pd.DataFrame]) -> None:
    for name, table in tables.items():
        table.to_csv(out_dir / f"{name}.csv", index=False)


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────
def run(config: RunConfig) -> Dict[str, pd.DataFrame]:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    pfas_path = find_first_existing(config.input_dir, PFAS_FILENAMES)
    sites_path = find_first_existing(config.input_dir, SITES_FILENAMES)
    reach_path = find_first_existing(config.input_dir, REACH_FILENAMES)

    if pfas_path is None:
        raise FileNotFoundError(f"No PFAS results CSV found in {config.input_dir}. Expected one of: {PFAS_FILENAMES}")

    report = []
    report.append("PFAS Reach Bracketing Tool validation report")
    report.append(f"Version: 1.1.1-public-release")
    report.append(f"PFAS input: {pfas_path}")
    report.append(f"Sites input: {sites_path if sites_path else 'not provided'}")
    report.append(f"Reach input: {reach_path if reach_path else 'not provided'}")
    report.append(f"N_IMPUTE: {config.n_impute}")
    report.append(f"Seed: {config.seed}")

    _print("Loading and validating inputs...")
    df, input_validation = read_pfas(pfas_path)
    sites = read_sites(sites_path, df)
    reach = read_reach(reach_path)

    # Filter to declared sites when sites metadata supplied; explicitly report excluded rows.
    excluded = pd.DataFrame()
    if sites_path is not None:
        allowed = set(sites["site_id"])
        excluded = df[~df["site_id"].isin(allowed)].copy()
        if len(excluded):
            report.append(f"Excluded rows not in site metadata: {len(excluded)}")
        df = df[df["site_id"].isin(allowed)].copy()

    n_analytes = df["analyte_std"].nunique()
    report.append(f"Included rows: {len(df)}")
    report.append(f"Included site-campaign samples: {df[['campaign','site_id']].drop_duplicates().shape[0]}")
    report.append(f"Included sites: {df['site_id'].nunique()}")
    report.append(f"Included analytes: {n_analytes}")
    if n_analytes != 40:
        report.append("WARNING: analyte count is not 40; ΣPFAS40 label may not be appropriate for this input.")

    _print(f"Running censoring-aware MI ({config.n_impute:,} replicates)...")
    site_stats, totals, sample_index, analytes, params, value_wide, rl_wide, detect_wide = run_multiple_imputation(df, config)
    site_stats = site_stats.merge(sites[["site_id", "Reach_group", "Bracket_role", "Lat", "Lon"]], on="site_id", how="left")

    _print("Computing diagnostics...")
    conv = convergence_diagnostics(totals, config.n_impute)
    sens, sens_summary = sensitivity_rl2(site_stats)
    report.append(f"Spearman rho MI vs RL/2: {sens_summary.loc[sens_summary.metric=='Spearman_rho_MI_vs_RL2','value'].iloc[0]:.4f}")
    report.append(f"Max |MI - RL/2|: {sens_summary.loc[sens_summary.metric=='Max_abs_MI_minus_RL2_ngL','value'].iloc[0]:.3f} ng/L")
    if len(conv):
        last = conv.iloc[-1]
        report.append(f"Convergence last N={int(last['N'])}: max |Δ median|={last['max_abs_diff_median_ngL']:.3f} ng/L")

    _print("Analyzing reach experiments...")
    dist_lookup = build_distribution_lookup(sample_index, totals)
    reach_results = analyze_reaches(reach, dist_lookup)
    report.append(f"Reach-campaign results: {len(reach_results)}")

    # Known manuscript-alignment checks when GK2v3 example dataset is present.
    def get_site(camp, site):
        rr = site_stats[(site_stats["campaign"] == camp) & (site_stats["site_id"] == site)]
        return float(rr["Sigma_MI_median"].iloc[0]) if len(rr) else np.nan
    tw91 = get_site("Campaign 2", "TW91")
    if np.isfinite(tw91):
        report.append(f"Manuscript check: Campaign 2 TW91 ΣPFAS40 median = {tw91:.1f} ng/L")
    if not reach_results.empty and "experiment_id" in reach_results.columns:
        wf1 = reach_results[reach_results["experiment_id"].eq("WF1_FarmersBranch") & reach_results["campaign"].eq("Campaign 2")]
        if len(wf1):
            report.append(f"Manuscript check: WF1_FarmersBranch C2 ΔΣPFAS40 = {wf1['delta_SigmaPFAS40_median'].iloc[0]:.1f} ng/L; {wf1['diagnostic_note'].iloc[0]}")

    tables = {
        "InputValidation": input_validation,
        "SiteTotals": site_stats,
        "ReachExperiments": reach_results,
        "Sensitivity_RL2": sens,
        "Sensitivity_Summary": sens_summary,
        "Convergence": conv,
        "ImputationParams": params,
    }
    if len(excluded):
        tables["ExcludedRows"] = excluded
    export_csvs(config.output_dir, tables)

    readme = [
        "PFAS Reach Bracketing Tool results workbook.",
        f"Generated by PFAS_Reach_Bracketing_Tool.py v1.1.1-public-release.",
        f"N_IMPUTE = {config.n_impute}; RNG seed = {config.seed}.",
        "ΣPFAS40 is the summed concentration of the target PFAS analytes in the input dataset.",
        "95% II = 2.5th to 97.5th percentile interval across multiple-imputation replicates.",
        "f_trib is an apparent concentration-space screening diagnostic only; it is not a flow, source-apportionment, or load fraction.",
        "Values outside 0-1 indicate that downstream concentration lies outside the two-endmember range and f_trib is not applicable.",
    ]
    _print("Writing workbook...")
    write_excel(config.output_dir / "PFAS_Reach_Bracketing_Results.xlsx", tables, readme)
    write_validation_report(config.output_dir / "validation_report.txt", report)

    _print("Generating figures...")
    make_site_concentration_figure(site_stats, sites, config.output_dir, config.dpi)
    make_sensitivity_figure(sens, config.output_dir, config.dpi)
    make_convergence_figure(conv, config.output_dir, config.dpi)
    make_detection_heatmap(df, config.output_dir, config.dpi)
    make_network_map(site_stats, sites, config.output_dir, config.dpi)

    _print("Done.")
    _print(f"Outputs written to: {config.output_dir.resolve()}")
    return tables


def parse_args(argv: Optional[List[str]] = None) -> RunConfig:
    parser = argparse.ArgumentParser(description="PFAS Reach Bracketing Tool — Python Edition")
    parser.add_argument("--base-dir", default=".", help="Base directory containing inputs/ and outputs/ folders")
    parser.add_argument("--input-dir", default=None, help="Input directory; default: <base-dir>/inputs, or <base-dir>/data if inputs/ is absent")
    parser.add_argument("--output-dir", default=None, help="Output directory; default: <base-dir>/outputs")
    parser.add_argument("--n-impute", type=int, default=DEFAULT_N_IMPUTE, help="Number of multiple-imputation replicates (default 5000)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed for reproducibility")
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI, help="Figure resolution")
    parser.add_argument("--require-exact-40", action="store_true", help="Fail if input does not contain exactly 40 analytes")
    args, unknown = parser.parse_known_args(argv)
    base = Path(args.base_dir).resolve()
    if args.input_dir:
        input_dir = Path(args.input_dir).resolve()
    else:
        input_dir = base / "inputs"
        # Browser bundles store demo inputs in ./data; use that automatically when
        # ./inputs is absent so the Python edition runs out of the box from the
        # same public-release package.
        if not input_dir.exists() and (base / "data").exists():
            input_dir = base / "data"
    output_dir = Path(args.output_dir).resolve() if args.output_dir else base / "outputs"
    return RunConfig(base_dir=base, input_dir=input_dir, output_dir=output_dir, n_impute=args.n_impute, seed=args.seed, dpi=args.dpi, require_exact_40=args.require_exact_40)


if __name__ == "__main__":
    cfg = parse_args()
    try:
        run(cfg)
    except Exception as exc:
        print("\nERROR: PFAS Reach Bracketing Tool failed.", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        raise

