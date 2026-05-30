# Changelog

## v1.1.1 browser + Python bundle (final deployment package)

- Bundles the companion Python edition (`code/PFAS_Reach_Bracketing_Tool.py`) with the GitHub Pages browser tool.
- Adds a browser-interface download placeholder/link for users who prefer scripted reuse in Spyder/Anaconda, Jupyter, or Google Colab.
- Updates tool attribution to Dr. Gehendra Kharel, STREAM Lab (Sustainable Tools for Risk Evaluation and Climate-Water Modeling), Texas Christian University.
- Keeps the manuscript-aligned browser defaults and reach schema from the v1.1.1 browser release.
- Adds browser and README methodological references acknowledging established censored-data, multiple-imputation, FAIR workflow, EPA Method 1633A, and PFAS stream-loading/source-screening literature.
- Keeps both `index.html` and `PFAS_Reach_Bracketing_Tool.html` as fully styled tool pages so GitHub Pages and local double-click use retain the visual design.


## v1.1.1 browser release

- Aligns browser defaults with revised manuscript: 5,000 MI replicates by default.
- Adds full reach schema fields in UI text: `reach_group`, `baseline_sites`, and `all_sites`.
- Updates `f_trib` logic: no clipping; reports `n/a` when downstream concentration falls outside the two-endmember range.
- Adds validation checks for required columns, invalid reporting limits, invalid detect flags, duplicate rows, and detected rows with missing concentration values.
- Excludes rows outside declared site metadata and reports exclusions.
- Adds cumulative baseline reach rows where `baseline_sites` are supplied.
- Uses updated `reach_experiment_table.csv`, including `WF0_TW09_to_TW08_Sentinel`.
- Adds downloadable CSV outputs and a browser validation report.


GitHub Pages live tool: https://kharelg100.github.io/pfas-reach-bracketing-tool/
GitHub source repository: https://github.com/kharelg100/pfas-reach-bracketing-tool
