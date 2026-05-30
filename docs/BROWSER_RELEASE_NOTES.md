# Browser Release Notes

## v1.1.1 browser + Python bundle (final deployment package)

This release adds the companion Python script to the browser package and exposes it through a download link in the tool interface. The author attribution has been updated to Dr. Gehendra Kharel, STREAM Lab (Sustainable Tools for Risk Evaluation and Climate-Water Modeling), Texas Christian University.


This browser release is intended to replace the earlier GitHub Pages version at:

https://kharelg100.github.io/pfas-reach-bracketing-tool/

Recommended update steps:

1. Replace the existing `index.html` with this release's `index.html`.
2. Add the `assets/`, `data/`, `docs/`, `validation/`, and `code/` folders to the GitHub Pages branch.
3. Confirm that `data/reach_experiment_table.csv` includes the manuscript schema fields:
   `experiment_id, reach_group, experiment_type, downstream_site, upstream_sites, tributary_sites, baseline_sites, all_sites`.
4. Confirm the site runs with the bundled demo data.
5. Confirm that the browser page has a working “Download Python script” link pointing to `code/PFAS_Reach_Bracketing_Tool.py`.
6. Update Zenodo only after GitHub Pages and the Python edition point to the same release.

Manuscript consistency checks:

- Default MI replicates: 5,000.
- `ΣPFAS40`: summed concentration of 40 target PFAS analytes.
- `f_trib`: apparent concentration-space tributary fraction only; not a flow or load fraction.
- Site IDs normalized as `TW08`, `TW09`, etc.
- Campaign 1 TW04 rows excluded because TW04 is not part of the GK2v3 network.


Additional final-deployment updates:

- Added a visible methodological references section in the browser interface and README to acknowledge the established literature underlying censored-data handling, multiple imputation, FAIR/reproducible workflow design, EPA Method 1633A, and PFAS stream loading/source-screening.
- Kept the public version label at **1.1.1** because intermediate packages were not deployed.
- Kept the full visual design in both `index.html` and `PFAS_Reach_Bracketing_Tool.html`; `index.html` is now the full tool page rather than a plain redirect.


GitHub Pages live tool: https://kharelg100.github.io/pfas-reach-bracketing-tool/
GitHub source repository: https://github.com/kharelg100/pfas-reach-bracketing-tool


Validation note: When running the manuscript demo data, the summary metric cards should display TW91 ΣPFAS40, WF1 ΔΣPFAS40, and WF1 ftrib. WF1 ftrib appears in the summary metrics and in the ReachExperiments.csv output.
