"""
PFAS Reach Bracketing Tool — Python Edition
=============================================
Machine-readable reach experiment analysis with multiple imputation.

Works in:
  • Spyder / Anaconda  — run as a script (F5)
  • Jupyter Notebook   — run cells top-to-bottom
  • Google Colab       — upload CSVs via files.upload() or mount Drive

Inputs (place in ./inputs/ or set paths below):
  1. pfas_long_filtered_standardized.csv  [REQUIRED]
  2. sites_metadata.csv                   [OPTIONAL]
  3. reach_experiment_table.csv           [OPTIONAL]

Outputs (written to ./outputs/):
  • PFAS_Reach_Bracketing_Results.xlsx    (5-sheet workbook)
  • Fig_SiteConcentrations.svg / .png
  • Fig_NetworkMap.svg / .png
  • Fig_Sensitivity_MI_vs_RL2.svg / .png
  • Fig_ConvergenceDiag.svg / .png

Methodological References:
  [1] Helsel, D.R. (2012) Statistics for Censored Environmental Data
      Using Minitab and R, 2nd ed. John Wiley & Sons.
  [2] Rubin, D.B. (1987) Multiple Imputation for Nonresponse in
      Surveys. John Wiley & Sons.
  [3] Lubin, J.H. et al. (2004) Epidemiologic quantification of exposure
      to environmental pollutants with left-censored data. Environ.
      Health Perspect. 112(17), 1691–1696.
  [4] Helsel, D.R. (2006) Fabricating data: How substituting values for
      nondetects can ruin results. Chemosphere 65(11), 2434–2439.
  [5] Antweiler, R.C. & Taylor, H.E. (2008) Evaluation of statistical
      treatments of left-censored environmental data. Environ. Sci.
      Technol. 42(10), 3732–3738.
  [6] US EPA (2024) Method 1633: Analysis of PFAS in aqueous, solid,
      biosolids, and tissue samples by LC-MS/MS. EPA-821-R-24-002.
  [7] Spearman, C. (1904) The proof and measurement of association
      between two things. Am. J. Psychol. 15(1), 72–101.
  [8] Woodward, D.S. et al. (2024) Time-of-travel synoptic survey with
      concurrent streamflow measurement for PFAS stream loading. ACS
      ES&T Water 4(10), 4356–4367.

© 2025–2026 Dr. Gehendra Kharel, Texas Christian University.
Contact: g.kharel@tcu.edu
"""

# %% ═══════════════════════════════════════════
#  IMPORTS
# ═══════════════════════════════════════════════
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from scipy.stats import norm
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# For Colab: uncomment below to upload files interactively
# from google.colab import files
# uploaded = files.upload()


# %% ═══════════════════════════════════════════
#  CONFIGURATION
# ═══════════════════════════════════════════════
N_IMPUTE    = 1000      # MI replicates (increase to 5000–10000 for publication)
RNG_SEED    = 20260227  # reproducibility seed
DPI         = 300       # figure resolution

# ── File paths (edit as needed) ──
BASE_DIR    = Path(".")
IN_DIR      = BASE_DIR / "inputs"
OUT_DIR     = BASE_DIR / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PFAS_PATH   = IN_DIR / "pfas_long_filtered_standardized.csv"
SITES_PATH  = IN_DIR / "sites_metadata.csv"
REACH_PATH  = IN_DIR / "reach_experiment_table.csv"

# TCU color palette
TCU         = "#4d1979"
TCU2        = "#6b2fa0"
TCU3        = "#8b5fbf"
TCUL        = "#f3edf9"
COL_UP      = "#2563eb"   # upstream / baseline
COL_TRIB    = TCU          # tributary
COL_DOWN    = "#ea580c"    # downstream / integrator
COL_GRAY    = "#8c95a5"


# %% ═══════════════════════════════════════════
#  HELPER FUNCTIONS
# ═══════════════════════════════════════════════
def split_pipe(s):
    """Split pipe-delimited string into list."""
    if pd.isna(s) or str(s).strip() == "":
        return []
    return [x.strip() for x in str(s).split("|") if x.strip()]


def fmt_interval(med, lo, hi, dec=1):
    """Format median [lo, hi] string."""
    if np.isnan(med):
        return ""
    return f"{med:.{dec}f} [{lo:.{dec}f}, {hi:.{dec}f}]"


def role_color(role):
    """Return color based on bracket role string."""
    if not role or pd.isna(role):
        return COL_GRAY
    r = str(role).lower()
    if any(k in r for k in ["upstream", "baseline", "pre-source"]):
        return COL_UP
    if "tributary" in r:
        return COL_TRIB
    if any(k in r for k in ["downstream", "integrator", "outlet", "confluence", "post-"]):
        return COL_DOWN
    return COL_GRAY


def role_label(role):
    """Short label for bracket role."""
    if not role or pd.isna(role):
        return ""
    r = str(role).lower()
    if "upstream" in r or "baseline" in r:
        return "Upstream"
    if "tributary" in r:
        return "Tributary"
    if "post-" in r or "post " in r:
        return "Post-confluence"
    if "downstream" in r:
        return "Downstream"
    if "outlet" in r:
        return "Outlet"
    if "integrator" in r:
        return "Integrator"
    if "confluence" in r:
        return "Confluence"
    return ""


# %% ═══════════════════════════════════════════
#  LOAD DATA
# ═══════════════════════════════════════════════
print("Loading data...")

# ── PFAS results (required) ──
if not PFAS_PATH.exists():
    raise FileNotFoundError(
        f"PFAS results CSV not found at {PFAS_PATH}\n"
        "Place your file in ./inputs/ or update PFAS_PATH above."
    )
df = pd.read_csv(PFAS_PATH)
df["result_value"] = pd.to_numeric(df["result_value"], errors="coerce")
df["rl"] = pd.to_numeric(df["rl"], errors="coerce")
df["detect_flag"] = df["detect_flag"].fillna(0).astype(int)
print(f"  PFAS results: {len(df):,} rows, {df['site_id'].nunique()} sites, "
      f"{df['analyte_std'].nunique()} analytes")

# ── Sites metadata (optional) ──
if SITES_PATH.exists():
    sites_meta = pd.read_csv(SITES_PATH)
    print(f"  Sites metadata: {len(sites_meta)} sites loaded")
else:
    sites_meta = pd.DataFrame({
        "Site ID": sorted(df["site_id"].unique()),
        "Lat": np.nan, "Lon": np.nan,
        "Reach_group": "", "Bracket_role": ""
    })
    print(f"  Sites metadata: auto-generated from {len(sites_meta)} unique site IDs")

# ── Reach experiments (optional) ──
if REACH_PATH.exists():
    reach = pd.read_csv(REACH_PATH)
    print(f"  Reach experiments: {len(reach)} experiments loaded")
    HAS_REACH = True
else:
    reach = pd.DataFrame()
    HAS_REACH = False
    print("  Reach experiments: not provided (skipping reach analysis)")

# Filter to declared sites
allowed_sites = sorted(sites_meta["Site ID"].unique())
n_before = len(df)
df = df[df["site_id"].isin(allowed_sites)].copy()
n_after = len(df)
if n_before != n_after:
    print(f"  Filtered: {n_before - n_after} rows from sites not in metadata")


# %% ═══════════════════════════════════════════
#  MULTIPLE IMPUTATION ENGINE
# ═══════════════════════════════════════════════
print(f"\nRunning MI with {N_IMPUTE:,} replicates (seed={RNG_SEED})...")

analytes = sorted(df["analyte_std"].unique())
idx_cols = ["campaign", "site_id"]

value_wide = df.pivot_table(index=idx_cols, columns="analyte_std",
                            values="result_value", aggfunc="first").reindex(columns=analytes)
rl_wide = df.pivot_table(index=idx_cols, columns="analyte_std",
                         values="rl", aggfunc="first").reindex(columns=analytes)
detect_wide = df.pivot_table(index=idx_cols, columns="analyte_std",
                             values="detect_flag", aggfunc="first").reindex(columns=analytes)

index = value_wide.index
n_samples = len(index)
rng = np.random.default_rng(RNG_SEED)

# Fit analyte-specific lognormal parameters (Helsel 2012 [1]; Lubin et al. 2004 [3])
params = {}
for a in analytes:
    det_vals = value_wide[a][detect_wide[a] == 1].dropna().values
    det_vals = det_vals[np.isfinite(det_vals) & (det_vals > 0)]
    if len(det_vals) >= 2:
        logs = np.log(det_vals)
        mu = float(logs.mean())
        sigma = float(logs.std(ddof=1))
        if sigma < 1e-6:
            sigma = 0.5
        params[a] = (mu, sigma, "lognorm", len(det_vals))
    else:
        params[a] = (None, None, "uniform", len(det_vals))

n_lognorm = sum(1 for v in params.values() if v[2] == "lognorm")
n_uniform = sum(1 for v in params.values() if v[2] == "uniform")
print(f"  Distribution fits: {n_lognorm} lognorm, {n_uniform} uniform")

# ── Impute totals (Rubin 1987 [2]; non-detects handled per Helsel 2012 [1]) ──
totals = np.zeros((n_samples, N_IMPUTE), dtype=np.float64)
rl2_totals = np.zeros(n_samples)
detect_only = np.zeros(n_samples)

for a in analytes:
    obs = value_wide[a].values.astype(float)
    rl_vals = rl_wide[a].values.astype(float)
    det_mask = (detect_wide[a].values.astype(float) == 1) & np.isfinite(obs)

    vals = np.zeros((n_samples, N_IMPUTE))

    # Detected: hold fixed
    if det_mask.any():
        vals[det_mask, :] = obs[det_mask, None]
        rl2_totals[det_mask] += obs[det_mask]
        detect_only[det_mask] += obs[det_mask]

    # Non-detect: impute
    nd_mask = ~det_mask
    if nd_mask.any():
        mu, sigma, mode, _ = params[a]
        upper = rl_vals[nd_mask]

        # RL/2 substitution
        valid_rl = np.where(np.isfinite(upper) & (upper > 0), upper, 0)
        rl2_totals[nd_mask] += valid_rl / 2

        if mode == "lognorm":
            logu = np.log(np.maximum(upper, 1e-10))
            p = norm.cdf((logu - mu) / sigma)
            p = np.clip(p, 1e-10, 1.0)
            u = rng.random((nd_mask.sum(), N_IMPUTE)) * p[:, None]
            z = mu + sigma * norm.ppf(u)
            draw = np.exp(z)
            draw = np.minimum(draw, upper[:, None])
            draw = np.maximum(draw, 0)
        else:
            draw = rng.random((nd_mask.sum(), N_IMPUTE)) * upper[:, None]

        vals[nd_mask, :] = draw

    totals += vals

# ── Summarize ──
site_stats = pd.DataFrame({
    "campaign": [i[0] for i in index],
    "site_id": [i[1] for i in index],
    "MI_median": np.median(totals, axis=1),
    "MI_lo": np.quantile(totals, 0.025, axis=1),
    "MI_hi": np.quantile(totals, 0.975, axis=1),
    "RL2": rl2_totals,
    "DetectOnly": detect_only,
})
site_stats["MI_interval"] = site_stats.apply(
    lambda r: fmt_interval(r["MI_median"], r["MI_lo"], r["MI_hi"]), axis=1)

# Add metadata
meta = sites_meta.rename(columns={"Site ID": "site_id"})
site_stats = site_stats.merge(meta[["site_id", "Reach_group", "Bracket_role"]],
                               on="site_id", how="left")

print(f"  Computed ΣPFAS for {n_samples} samples")
print(f"  Concentration range: {site_stats['MI_median'].min():.1f} – "
      f"{site_stats['MI_median'].max():.1f} ng/L")


# %% ═══════════════════════════════════════════
#  CONVERGENCE DIAGNOSTICS
# ═══════════════════════════════════════════════
# Convergence diagnostics (Rubin 1987 [2])
print("\nConvergence diagnostics...")
ref_med = np.median(totals, axis=1)
ref_lo = np.quantile(totals, 0.025, axis=1)
ref_hi = np.quantile(totals, 0.975, axis=1)

conv_Ns = [250, 500, 1000, 2000, 5000]
conv_Ns = [N for N in conv_Ns if N < N_IMPUTE]
conv_rows = []
for N in conv_Ns:
    medN = np.median(totals[:, :N], axis=1)
    loN = np.quantile(totals[:, :N], 0.025, axis=1)
    hiN = np.quantile(totals[:, :N], 0.975, axis=1)
    conv_rows.append({
        "N": N,
        "max_diff_median": float(np.max(np.abs(medN - ref_med))),
        "max_diff_lo": float(np.max(np.abs(loN - ref_lo))),
        "max_diff_hi": float(np.max(np.abs(hiN - ref_hi))),
    })
conv_df = pd.DataFrame(conv_rows)
if len(conv_df):
    print(f"  Max |Δ median| at N={conv_Ns[-1]}: {conv_rows[-1]['max_diff_median']:.3f} ng/L")


# %% ═══════════════════════════════════════════
#  RL/2 SENSITIVITY
# ═══════════════════════════════════════════════
# RL/2 sensitivity (cf. Helsel 2006 [4]; Antweiler & Taylor 2008 [5])
print("\nRL/2 sensitivity...")
from scipy.stats import spearmanr

rho, _ = spearmanr(site_stats["MI_median"], site_stats["RL2"])  # Spearman 1904 [7]
max_diff = np.max(np.abs(site_stats["MI_median"] - site_stats["RL2"]))
print(f"  Spearman ρ = {rho:.4f}")
print(f"  Max |MI − RL/2| = {max_diff:.3f} ng/L")


# %% ═══════════════════════════════════════════
#  REACH EXPERIMENT ANALYSIS
#  Concentration-based mixing fraction (f_trib) as screening diagnostic;
#  does not replace discharge-based load calculations (Woodward et al. 2024 [8])
# ═══════════════════════════════════════════════
row_map = {(index[i][0], index[i][1]): i for i in range(n_samples)}

def get_dist(campaign, site):
    key = (campaign, site)
    if key in row_map:
        return totals[row_map[key], :]
    return None

reach_results = []

if HAS_REACH:
    print("\nAnalyzing reach experiments...")
    for _, r in reach.iterrows():
        ds = r["downstream_site"]
        ups = split_pipe(r.get("upstream_sites"))
        trs = split_pipe(r.get("tributary_sites"))
        bls = split_pipe(r.get("baseline_sites"))
        all_sites = split_pipe(r.get("all_sites"))
        exp_type = r["experiment_type"]

        for camp in ["Campaign 1", "Campaign 2"]:
            have = set(site_stats[site_stats["campaign"] == camp]["site_id"])
            if not set(all_sites).issubset(have):
                continue

            dd = get_dist(camp, ds)
            if dd is None:
                continue

            delta_txt, diag, delta_med, f_med = "", "", np.nan, np.nan

            if exp_type in ["tributary_increment", "downstream_propagation"]:
                up = ups[0] if ups else None
                if up:
                    ud = get_dist(camp, up)
                    if ud is not None:
                        delta = dd - ud
                        delta_med = float(np.median(delta))
                        delta_txt = fmt_interval(
                            np.median(delta), np.quantile(delta, 0.025),
                            np.quantile(delta, 0.975))
                if trs:
                    td = get_dist(camp, trs[0])
                    ud2 = get_dist(camp, up) if up else None
                    if td is not None and ud2 is not None:
                        denom = td - ud2
                        numer = dd - ud2
                        f = np.where(denom != 0, numer / denom, np.nan)
                        f_clip = np.clip(f, 0, 1)
                        valid_frac = float(np.nanmean((f >= 0) & (f <= 1)))
                        f_med = float(np.nanmedian(f_clip))
                        diag = f"f≈{f_med:.3f} (valid {valid_frac:.2f})"

            elif exp_type == "source_reach":
                up = ups[0] if ups else None
                if up:
                    ud = get_dist(camp, up)
                    if ud is not None:
                        delta = dd - ud
                        delta_med = float(np.median(delta))
                        delta_txt = fmt_interval(
                            np.median(delta), np.quantile(delta, 0.025),
                            np.quantile(delta, 0.975))
                        diag = f"segment: {ds}–{up}"

            elif exp_type == "confluence_mix" and len(ups) >= 2:
                u1 = get_dist(camp, ups[0])
                u2 = get_dist(camp, ups[1])
                if u1 is not None and u2 is not None:
                    hi = np.maximum(u1, u2)
                    lo = np.minimum(u1, u2)
                    outside = (dd < lo) | (dd > hi)
                    outside_prob = float(np.mean(outside))
                    denom = hi - lo
                    f = np.where(denom != 0, (dd - lo) / denom, np.nan)
                    f_clip = np.clip(f, 0, 1)
                    f_med = float(np.nanmedian(f_clip))
                    diag = f"outside-range={outside_prob:.2f}; f≈{f_med:.3f}"

            reach_results.append({
                "experiment_id": r["experiment_id"],
                "reach_group": r["reach_group"],
                "experiment_type": exp_type,
                "campaign": camp,
                "delta": delta_txt,
                "delta_median": delta_med,
                "f_trib": f_med,
                "diagnostic": diag,
            })

    reach_df = pd.DataFrame(reach_results)
    print(f"  {len(reach_df)} reach-campaign combinations computed")
else:
    reach_df = pd.DataFrame()


# %% ═══════════════════════════════════════════
#  EXPORT EXCEL WORKBOOK
# ═══════════════════════════════════════════════
print("\nExporting Excel workbook...")
xl_path = OUT_DIR / "PFAS_Reach_Bracketing_Results.xlsx"

params_df = pd.DataFrame([
    {"analyte": a, "n_detect": params[a][3], "mode": params[a][2],
     "mu_log": params[a][0], "sigma_log": params[a][1]}
    for a in analytes
]).sort_values(["mode", "n_detect"], ascending=[True, False])

sens_df = site_stats[["site_id", "campaign", "MI_median", "RL2", "DetectOnly"]].copy()
sens_df["Diff_MI_RL2"] = sens_df["MI_median"] - sens_df["RL2"]
sens_df["Pct_diff"] = np.where(
    sens_df["MI_median"] > 0,
    sens_df["Diff_MI_RL2"] / sens_df["MI_median"] * 100, 0)

with pd.ExcelWriter(xl_path, engine="openpyxl") as writer:
    site_stats.to_excel(writer, sheet_name="SiteTotals", index=False)
    if len(reach_df):
        reach_df.to_excel(writer, sheet_name="ReachExperiments", index=False)
    sens_df.to_excel(writer, sheet_name="Sensitivity_RL2", index=False)
    conv_df.to_excel(writer, sheet_name="Convergence", index=False)
    params_df.to_excel(writer, sheet_name="ImputationParams", index=False)

print(f"  Saved: {xl_path}")


# %% ═══════════════════════════════════════════
#  FIGURE 1: SITE CONCENTRATION BAR CHART
# ═══════════════════════════════════════════════
print("\nGenerating figures...")

sorted_stats = site_stats.sort_values("MI_median", ascending=False).reset_index(drop=True)
n_bars = len(sorted_stats)

fig1, ax1 = plt.subplots(figsize=(max(10, n_bars * 0.6), 5.5))

use_log = (sorted_stats["MI_hi"].max() / sorted_stats["MI_lo"][sorted_stats["MI_lo"] > 0].min()) > 50

for i, (_, row) in enumerate(sorted_stats.iterrows()):
    col = role_color(row.get("Bracket_role"))
    ax1.bar(i, row["MI_median"], width=0.65, color=col, alpha=0.8, zorder=3)
    ax1.errorbar(i, row["MI_median"],
                 yerr=[[row["MI_median"] - row["MI_lo"]],
                       [row["MI_hi"] - row["MI_median"]]],
                 fmt="none", ecolor=COL_GRAY, capsize=3, linewidth=1.2, zorder=4)

labels = [f"{r['site_id']}\n({r['campaign'].replace('Campaign ', 'C')})"
          for _, r in sorted_stats.iterrows()]
ax1.set_xticks(range(n_bars))
ax1.set_xticklabels(labels, fontsize=9, fontweight="bold", rotation=40, ha="right")

if use_log:
    ax1.set_yscale("log")
ax1.set_ylabel("ΣPFAS (ng/L)" + (" — log scale" if use_log else ""),
               fontsize=12, fontweight="bold")
ax1.set_title("ΣPFAS by Site (MI Median + 95% Imputation Interval)",
              fontsize=14, fontweight="bold", color=TCU)
ax1.grid(axis="y", alpha=0.3, linewidth=0.5)
ax1.set_xlim(-0.6, n_bars - 0.4)

# Legend
patches = [
    mpatches.Patch(color=COL_UP, alpha=0.8, label="Upstream / Baseline"),
    mpatches.Patch(color=COL_TRIB, alpha=0.8, label="Tributary"),
    mpatches.Patch(color=COL_DOWN, alpha=0.8, label="Downstream / Integrator"),
]
ax1.legend(handles=patches, fontsize=9, frameon=True, loc="upper right")
plt.tight_layout()

for ext in ["svg", "png"]:
    fig1.savefig(OUT_DIR / f"Fig_SiteConcentrations.{ext}", dpi=DPI, facecolor="white")
plt.show()
print(f"  Saved: Fig_SiteConcentrations.svg/.png")


# %% ═══════════════════════════════════════════
#  FIGURE 2: NETWORK MAP
# ═══════════════════════════════════════════════
has_coords = sites_meta["Lat"].notna().any() and sites_meta["Lon"].notna().any()

if has_coords:
    fig2, ax2 = plt.subplots(figsize=(9, 7))
    ax2.set_facecolor("#fafbfc")

    # Build stats lookup
    stats_map = {}
    for _, r in site_stats.iterrows():
        if r["site_id"] not in stats_map:
            stats_map[r["site_id"]] = {}
        stats_map[r["site_id"]][r["campaign"]] = r

    # Concentration colormap
    cmap = LinearSegmentedColormap.from_list("pfas", [
        (0.0, "#2166ac"), (0.15, "#67a9cf"), (0.35, "#f5a623"),
        (0.6, "#e8491d"), (1.0, "#8b0000")])
    meds = site_stats["MI_median"].values
    cnorm = LogNorm(vmin=max(1, meds[meds > 0].min() * 0.5),
                    vmax=meds.max() * 1.5)

    # Draw edges from reach table
    if HAS_REACH:
        for _, r in reach.iterrows():
            ds_row = sites_meta[sites_meta["Site ID"] == r["downstream_site"]]
            if ds_row.empty or pd.isna(ds_row.iloc[0]["Lat"]):
                continue
            ds_xy = (ds_row.iloc[0]["Lon"], ds_row.iloc[0]["Lat"])

            for up in split_pipe(r.get("upstream_sites")):
                up_row = sites_meta[sites_meta["Site ID"] == up]
                if up_row.empty or pd.isna(up_row.iloc[0]["Lat"]):
                    continue
                ax2.plot([up_row.iloc[0]["Lon"], ds_xy[0]],
                         [up_row.iloc[0]["Lat"], ds_xy[1]],
                         color=TCU3, linewidth=2, alpha=0.4, zorder=1)

            for tr in split_pipe(r.get("tributary_sites")):
                tr_row = sites_meta[sites_meta["Site ID"] == tr]
                if tr_row.empty or pd.isna(tr_row.iloc[0]["Lat"]):
                    continue
                ax2.plot([tr_row.iloc[0]["Lon"], ds_xy[0]],
                         [tr_row.iloc[0]["Lat"], ds_xy[1]],
                         color=TCU, linewidth=1.5, alpha=0.5,
                         linestyle="--", dashes=(5, 3), zorder=1)

    # Plot sites
    for _, s in sites_meta.iterrows():
        if pd.isna(s["Lat"]) or pd.isna(s["Lon"]):
            continue
        sid = s["Site ID"]
        sm = stats_map.get(sid, {})
        best = sm.get("Campaign 2") if "Campaign 2" in sm else sm.get("Campaign 1")
        med = best["MI_median"] if best is not None else None

        r_size = 80
        if med and med > 0:
            t = np.clip((np.log10(med) - np.log10(cnorm.vmin)) /
                        (np.log10(cnorm.vmax) - np.log10(cnorm.vmin)), 0, 1)
            r_size = 60 + t * 300
            color = cmap(cnorm(med))
        else:
            color = "#cccccc"

        ax2.scatter(s["Lon"], s["Lat"], s=r_size, c=[color],
                    edgecolors="white", linewidths=1.5, zorder=5)

        rl = role_label(s.get("Bracket_role"))
        if rl:
            ax2.text(s["Lon"], s["Lat"] + 0.003, rl, fontsize=8,
                     fontstyle="italic", ha="center", va="bottom",
                     color=role_color(s.get("Bracket_role")), zorder=6)

        ax2.text(s["Lon"], s["Lat"] - 0.003, sid, fontsize=10,
                 fontweight="bold", ha="center", va="top", color=TCU, zorder=6)

        if med is not None:
            label = f"{med:.0f}" if med < 1000 else f"{med/1000:.1f}k"
            ax2.text(s["Lon"], s["Lat"] - 0.007, f"{label} ng/L",
                     fontsize=8, ha="center", va="top", color="#3f4a5a", zorder=6)

    ax2.set_xlabel("Longitude", fontsize=11)
    ax2.set_ylabel("Latitude", fontsize=11)
    ax2.set_title("Sampling Network — ΣPFAS Concentration",
                  fontsize=14, fontweight="bold", color=TCU)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=TCU3, linewidth=2, alpha=0.5, label="Mainstem connection"),
        Line2D([0], [0], color=TCU, linewidth=1.5, linestyle="--", alpha=0.6,
               label="Tributary connection"),
    ]
    ax2.legend(handles=legend_elements, fontsize=9, loc="lower right", frameon=True)
    plt.tight_layout()

    for ext in ["svg", "png"]:
        fig2.savefig(OUT_DIR / f"Fig_NetworkMap.{ext}", dpi=DPI, facecolor="white")
    plt.show()
    print(f"  Saved: Fig_NetworkMap.svg/.png")
else:
    print("  Network map: skipped (no coordinates in sites metadata)")


# %% ═══════════════════════════════════════════
#  FIGURE 3: MI vs RL/2 SENSITIVITY (DUMBBELL)
# ═══════════════════════════════════════════════
sorted_sens = site_stats.sort_values("MI_median", ascending=True).reset_index(drop=True)
n_rows = len(sorted_sens)

fig3, ax3 = plt.subplots(figsize=(10, max(4, n_rows * 0.42)))

for i, (_, row) in enumerate(sorted_sens.iterrows()):
    mi_val = row["MI_median"]
    rl_val = row["RL2"]
    diff_pct = abs(mi_val - rl_val) / mi_val * 100 if mi_val > 0 else 0

    # Gap bar color
    if diff_pct < 1:
        gap_col = "#d1fae5"   # green
    elif diff_pct < 5:
        gap_col = "#e8d5f5"   # light purple
    elif diff_pct < 15:
        gap_col = "#c084fc"   # medium purple
    else:
        gap_col = "#9333ea"   # deep purple

    lo_val = min(mi_val, rl_val)
    hi_val = max(mi_val, rl_val)

    # Alternating row background
    if i % 2 == 0:
        ax3.axhspan(i - 0.4, i + 0.4, color=TCUL, alpha=0.3, zorder=0)

    # Gap bar
    if hi_val > 0 and lo_val > 0:
        ax3.barh(i, hi_val - lo_val, left=lo_val, height=0.35,
                 color=gap_col, alpha=0.6, zorder=2)

    # RL/2 dot (open)
    ax3.scatter(rl_val, i, s=70, facecolors="white", edgecolors=COL_GRAY,
                linewidths=2, zorder=4)
    # MI dot (filled)
    ax3.scatter(mi_val, i, s=90, c=TCU, edgecolors="white",
                linewidths=1.5, zorder=5)

    # % annotation
    pct_col = "#16a34a" if diff_pct < 1 else ("#5f6b7a" if diff_pct < 5 else "#9333ea")
    ax3.text(ax3.get_xlim()[1] if i == 0 else sorted_sens["RL2"].max() * 2.8,
             i, f"  {diff_pct:.1f}%", fontsize=10, fontweight="bold",
             color=pct_col, va="center", zorder=6)

ylabels = [f"{r['site_id']} ({r['campaign'].replace('Campaign ', 'C')})"
           for _, r in sorted_sens.iterrows()]
ax3.set_yticks(range(n_rows))
ax3.set_yticklabels(ylabels, fontsize=10, fontweight="bold")
ax3.set_xscale("log")
ax3.set_xlabel("ΣPFAS (ng/L, log scale)", fontsize=12, fontweight="bold")
ax3.set_title(f"Censoring Sensitivity: MI Median vs RL/2 Substitution\n"
              f"Spearman ρ = {rho:.4f}",
              fontsize=13, fontweight="bold", color=TCU)
ax3.grid(axis="x", alpha=0.3, linewidth=0.5)
ax3.set_ylim(-0.6, n_rows - 0.4)

# Legend
legend_elements = [
    plt.scatter([], [], s=90, c=TCU, edgecolors="white", linewidths=1.5,
                label="MI Median"),
    plt.scatter([], [], s=70, facecolors="white", edgecolors=COL_GRAY,
                linewidths=2, label="RL/2"),
    mpatches.Patch(color="#e8d5f5", alpha=0.6, label="Gap between methods"),
]
ax3.legend(handles=legend_elements, fontsize=9, loc="lower right", frameon=True)

# Recalculate xlim to make room for % annotations
xmax = max(sorted_sens["MI_median"].max(), sorted_sens["RL2"].max()) * 4
ax3.set_xlim(right=xmax)

plt.tight_layout()
for ext in ["svg", "png"]:
    fig3.savefig(OUT_DIR / f"Fig_Sensitivity_MI_vs_RL2.{ext}", dpi=DPI, facecolor="white")
plt.show()
print(f"  Saved: Fig_Sensitivity_MI_vs_RL2.svg/.png")


# %% ═══════════════════════════════════════════
#  FIGURE 4: CONVERGENCE DIAGNOSTICS
# ═══════════════════════════════════════════════
if len(conv_df):
    fig4, ax4 = plt.subplots(figsize=(7, 4.5))

    ax4.plot(conv_df["N"], conv_df["max_diff_median"], "o-",
             color=TCU, linewidth=2, markersize=8, label="|Δ median|", zorder=3)
    ax4.plot(conv_df["N"], conv_df["max_diff_lo"], "s--",
             color=COL_UP, linewidth=1.5, markersize=6, label="|Δ 95% II lower|", zorder=3)
    ax4.plot(conv_df["N"], conv_df["max_diff_hi"], "^--",
             color=COL_DOWN, linewidth=1.5, markersize=6, label="|Δ 95% II upper|", zorder=3)

    ax4.axhline(y=1.0, color="#dc2626", linewidth=1, linestyle=":", alpha=0.5, zorder=1)
    ax4.text(conv_df["N"].iloc[0], 1.05, "1 ng/L threshold", fontsize=8,
             color="#dc2626", alpha=0.7)

    ax4.set_xlabel("Number of MI Replicates (N)", fontsize=12, fontweight="bold")
    ax4.set_ylabel("Max |Difference| vs Full Run (ng/L)", fontsize=11, fontweight="bold")
    ax4.set_title(f"Imputation Convergence (reference N = {N_IMPUTE:,})",
                  fontsize=13, fontweight="bold", color=TCU)
    ax4.legend(fontsize=9, frameon=True)
    ax4.grid(alpha=0.3, linewidth=0.5)
    ax4.set_yscale("log")
    plt.tight_layout()

    for ext in ["svg", "png"]:
        fig4.savefig(OUT_DIR / f"Fig_ConvergenceDiag.{ext}", dpi=DPI, facecolor="white")
    plt.show()
    print(f"  Saved: Fig_ConvergenceDiag.svg/.png")


# %% ═══════════════════════════════════════════
#  PRINT SUMMARY
# ═══════════════════════════════════════════════
print("\n" + "=" * 60)
print("  PFAS REACH BRACKETING TOOL — ANALYSIS COMPLETE")
print("=" * 60)
print(f"  Sites analyzed:       {site_stats['site_id'].nunique()}")
print(f"  Total samples:        {n_samples}")
print(f"  Analytes:             {len(analytes)}")
print(f"  MI replicates:        {N_IMPUTE:,}")
print(f"  Lognorm fits:         {n_lognorm}")
print(f"  Uniform fits:         {n_uniform}")
if HAS_REACH:
    print(f"  Reach experiments:    {len(reach_df)}")
print(f"  Spearman ρ (MI/RL2):  {rho:.4f}")
print(f"  Concentration range:  {site_stats['MI_median'].min():.1f} – "
      f"{site_stats['MI_median'].max():.1f} ng/L")
print(f"\n  Outputs → {OUT_DIR.resolve()}")
print("=" * 60)
print("  © 2025–2026 Dr. Gehendra Kharel, TCU | g.kharel@tcu.edu")
print("=" * 60)
