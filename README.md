# PFAS Reach Bracketing Tool

**Machine-readable reach experiment analysis with censoring-aware multiple imputation for PFAS surface water networks.**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10.5281/zenodo.18937740.svg)](https://doi.org/10.5281/zenodo.10.5281/zenodo.18937740)
[![License: MIT](https://img.shields.io/badge/License-MIT-purple.svg)](LICENSE)

> **Live Tool:** [https://gkharel.github.io/pfas-reach-bracketing-tool/](https://gkharel.github.io/pfas-reach-bracketing-tool/)

---

## Overview

An interactive, open-access browser-based tool that implements the machine-readable reach experiment framework for identifying PFAS sources in surface water networks. The tool performs censoring-aware multiple imputation (MI), computes reach-level concentration contrasts and tributary mixing fractions, evaluates sensitivity to censoring treatment (MI vs RL/2), and provides convergence diagnostics — all within the user's browser with no data uploaded to any server.

**Associated manuscript:**  
Kharel, G and Harvey, O. Machine-readable reach bracketing reveals persistent PFAS inputs from headwater tributaries. *ES&T Letters* (submitted 2026).

## Quick Start

### Browser Tool (no installation)
1. Download [`index.html`](index.html) or visit the [live tool](https://gkharel.github.io/pfas-reach-bracketing-tool/)
2. Click **"▶ Run with demo data"** to see all features instantly
3. Upload your own CSVs to analyze your network

### Python Script (Spyder / Jupyter / Colab)
1. Place your CSVs in `./inputs/`
2. Run `python PFAS_Reach_Bracketing_Tool.py`
3. Outputs appear in `./outputs/`

```bash
# Dependencies (all standard Anaconda)
pip install numpy pandas matplotlib scipy openpyxl
```

## Input Data

| File | Required? | Description |
|------|-----------|-------------|
| `pfas_results.csv` | **Yes** | One row per analyte per sample: `site_id, campaign, analyte_std, result_value, detect_flag, rl` |
| `sites_metadata.csv` | Optional | One row per site: `Site ID, Lat, Lon, Reach_group, Bracket_role` |
| `reach_experiment_table.csv` | Optional | One row per experiment: `experiment_id, experiment_type, downstream_site, upstream_sites, tributary_sites, all_sites` |

Download blank templates from the [`templates/`](templates/) folder or from within the browser tool.

## Features

- **Multiple Imputation** — 250 to 10,000 replicates with analyte-specific truncated lognormal distributions
- **Reach Experiments** — Automated ΔΣPFAS contrasts and f_trib mixing fraction diagnostics from a machine-readable experiment table
- **Sensitivity Analysis** — MI vs RL/2 substitution comparison with Spearman ρ and per-site dumbbell chart
- **Convergence Diagnostics** — Monte Carlo stability assessment across replicate subsets
- **Network Visualization** — Concentration-encoded site map with reach connectivity
- **Publication-Quality Exports** — SVG figures and 5-sheet Excel workbook
- **Zero Installation** — Single HTML file, runs in any modern browser
- **Complete Privacy** — All computation client-side; no data transmitted anywhere

## Repository Structure

```
pfas-reach-bracketing-tool/
├── index.html                          # Browser tool (GitHub Pages entry point)
├── python/
│   └── PFAS_Reach_Bracketing_Tool.py   # Python script (Spyder/Jupyter/Colab)
├── templates/
│   ├── pfas_results_template.csv       # PFAS results CSV template
│   ├── sites_metadata_template.csv     # Sites metadata CSV template
│   └── reach_experiment_template.csv   # Reach experiment CSV template
├── docs/
│   └── SETUP_GUIDE.md                  # Hosting setup instructions
├── LICENSE                             # MIT License
├── CITATION.cff                        # Citation metadata
├── .zenodo.json                        # Zenodo DOI metadata
├── netlify.toml                        # Netlify deployment config
└── README.md                           # This file
```

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{kharel2026pfas_tool,
  author    = {Kharel, Gehendra},
  title     = {{PFAS Reach Bracketing Tool: Machine-readable reach experiment 
                analysis with multiple imputation}},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.18937740},
  url       = {https://github.com/gkharel/pfas-reach-bracketing-tool}
}
```

## Methodological References

The tool implements established statistical methods for censored environmental data. Key references:

1. **Helsel, D.R.** (2012) *Statistics for Censored Environmental Data Using Minitab and R*, 2nd ed. John Wiley & Sons. — Foundational text for censored data methods.
2. **Rubin, D.B.** (1987) *Multiple Imputation for Nonresponse in Surveys*. John Wiley & Sons. — Theoretical basis for MI.
3. **Lubin, J.H. et al.** (2004) Epidemiologic quantification of exposure to environmental pollutants with left-censored data. *Environ. Health Perspect.* 112(17), 1691–1696.
4. **Helsel, D.R.** (2006) Fabricating data: How substituting values for nondetects can ruin results. *Chemosphere* 65(11), 2434–2439.
5. **Antweiler, R.C. & Taylor, H.E.** (2008) Evaluation of statistical treatments of left-censored environmental data. *Environ. Sci. Technol.* 42(10), 3732–3738.
6. **US EPA** (2024) Method 1633: Analysis of PFAS in aqueous, solid, biosolids, and tissue samples by LC-MS/MS. EPA-821-R-24-002.
7. **Spearman, C.** (1904) The proof and measurement of association between two things. *Am. J. Psychol.* 15(1), 72–101.
8. **Woodward, D.S. et al.** (2024) Time-of-travel synoptic survey with concurrent streamflow measurement for PFAS stream loading. *ACS ES&T Water* 4(10), 4356–4367.

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contact

**Dr. Gehendra Kharel**  
Department of Environmental Sciences  
Texas Christian University  
Fort Worth, TX 76129, USA  
📧 [g.kharel@tcu.edu](mailto:g.kharel@tcu.edu)

---
© 2025–2026 Dr. Gehendra Kharel, Texas Christian University. All rights reserved.
