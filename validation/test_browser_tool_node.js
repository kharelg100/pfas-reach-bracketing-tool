const fs = require('fs');
const path = require('path');
const core = require('../assets/pfas_reach_tool.js');

const base = path.join(__dirname, '..');
const pfas = fs.readFileSync(path.join(base, 'data', 'pfas_long_filtered_standardized.csv'), 'utf8');
const sites = fs.readFileSync(path.join(base, 'data', 'sites_metadata.csv'), 'utf8');
const reach = fs.readFileSync(path.join(base, 'data', 'reach_experiment_table.csv'), 'utf8');

const result = core.runAnalysis({
  pfasRows: core.parseCSV(pfas),
  sitesRows: core.parseCSV(sites),
  reachRows: core.parseCSV(reach)
}, {nImpute: 5000, seed: 20260227});

const tw91 = result.siteTotals.find(r => r.campaign === 'Campaign 2' && r.site_id === 'TW91');
const wf1 = result.reachResults.find(r => r.campaign === 'Campaign 2' && r.experiment_id === 'WF1_FarmersBranch');
const rho = result.sensitivitySummary.find(r => r.metric === 'Spearman_rho_MI_vs_RL2').value;

if (!tw91 || Math.abs(tw91.Sigma_MI_median - 15667) > 5) throw new Error('TW91 validation failed');
if (!wf1 || Math.abs(wf1.delta_SigmaPFAS40_median - 728.4) > 5) throw new Error('WF1 validation failed');
if (Math.abs(rho - 1) > 1e-8) throw new Error('Spearman validation failed');

fs.writeFileSync(path.join(base, 'validation', 'browser_validation_report.txt'), result.report);
console.log(result.report);
