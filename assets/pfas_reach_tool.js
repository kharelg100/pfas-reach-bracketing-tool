/* PFAS Reach Bracketing Tool — browser edition v1.1.1
 * Machine-readable reach experiment analysis with multiple imputation.
 * MIT License. No data are uploaded; all computation occurs in the browser.
 */
(function(){
  'use strict';

  const UNIT = 'ng L-1';
  const DEFAULT_N_IMPUTE = 5000;
  const DEFAULT_SEED = 20260227;
  const REQUIRED_PFAS_COLUMNS = ['site_id','campaign','analyte_std','result_value','detect_flag','rl'];
  const REQUIRED_REACH_COLUMNS = ['experiment_id','experiment_type','downstream_site','upstream_sites'];

  function normalizeSiteId(x){
    if (x === null || x === undefined) return '';
    let s = String(x).trim().toUpperCase().replace(/\s+/g,'');
    if (/^TW\d+$/.test(s)) {
      const n = parseInt(s.slice(2), 10);
      if (n < 100) return 'TW' + String(n).padStart(2,'0');
    }
    return s;
  }
  function splitPipe(x){
    if (x === null || x === undefined || String(x).trim()==='') return [];
    return String(x).split('|').map(v => normalizeSiteId(v)).filter(Boolean);
  }
  function num(x){
    if (x === null || x === undefined || String(x).trim()==='') return NaN;
    const v = Number(String(x).replace(/,/g,''));
    return Number.isFinite(v) ? v : NaN;
  }
  function isFiniteNumber(x){ return typeof x === 'number' && Number.isFinite(x); }
  function median(arr){ return quantile(arr, 0.5); }
  function quantile(arr, q){
    const a = Array.from(arr).filter(Number.isFinite).sort((x,y)=>x-y);
    if (!a.length) return NaN;
    const pos = (a.length - 1) * q;
    const lo = Math.floor(pos), hi = Math.ceil(pos), frac = pos - lo;
    return lo === hi ? a[lo] : a[lo]*(1-frac) + a[hi]*frac;
  }
  function fmt(x, d=1){ return Number.isFinite(x) ? x.toFixed(d) : ''; }
  function fmtInterval(med, lo, hi, d=1){ return Number.isFinite(med) ? `${fmt(med,d)} [${fmt(lo,d)}, ${fmt(hi,d)}]` : ''; }
  function escapeHtml(s){ return String(s ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c])); }

  function parseCSV(text){
    if (text && text.charCodeAt(0) === 0xFEFF) text = text.slice(1);
    const rows = [];
    let row = [], field = '', inQuotes = false;
    for (let i=0; i<text.length; i++){
      const c = text[i], next = text[i+1];
      if (inQuotes){
        if (c === '"'){
          if (next === '"'){ field += '"'; i++; }
          else inQuotes = false;
        } else field += c;
      } else {
        if (c === '"') inQuotes = true;
        else if (c === ','){ row.push(field); field = ''; }
        else if (c === '\n') { row.push(field); rows.push(row); row = []; field = ''; }
        else if (c === '\r') { /* ignore */ }
        else field += c;
      }
    }
    row.push(field); rows.push(row);
    while (rows.length && rows[rows.length-1].every(v => String(v).trim()==='')) rows.pop();
    if (!rows.length) return [];
    const headers = rows[0].map(h => String(h).trim());
    return rows.slice(1).filter(r => r.some(v => String(v).trim() !== '')).map(r => {
      const obj = {};
      headers.forEach((h,i)=> obj[h] = r[i] !== undefined ? r[i] : '');
      return obj;
    });
  }
  function toCSV(rows){
    if (!rows || !rows.length) return '';
    const headers = Object.keys(rows[0]);
    const esc = v => {
      const s = String(v ?? '');
      return /[",\n\r]/.test(s) ? '"' + s.replace(/"/g,'""') + '"' : s;
    };
    return [headers.join(','), ...rows.map(r => headers.map(h => esc(r[h])).join(','))].join('\n');
  }
  function downloadText(name, content, type='text/plain'){
    const blob = new Blob([content], {type});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(()=>URL.revokeObjectURL(url), 5000);
  }

  // Deterministic PRNG (Mulberry32)
  function makeRng(seed){
    let t = (Number(seed) >>> 0) || 1;
    return function(){
      t += 0x6D2B79F5;
      let r = Math.imul(t ^ (t >>> 15), 1 | t);
      r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
      return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
    };
  }
  // Abramowitz-Stegun erf approximation
  function erf(x){
    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);
    const a1=0.254829592, a2=-0.284496736, a3=1.421413741, a4=-1.453152027, a5=1.061405429, p=0.3275911;
    const t = 1/(1+p*x);
    const y = 1 - (((((a5*t+a4)*t)+a3)*t+a2)*t+a1)*t*Math.exp(-x*x);
    return sign*y;
  }
  function normCdf(x){ return 0.5 * (1 + erf(x / Math.SQRT2)); }
  // Acklam inverse normal CDF approximation
  function normInv(p){
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    const a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00];
    const b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01];
    const c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00];
    const d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00];
    const plow=0.02425, phigh=1-plow;
    let q, r;
    if (p < plow){
      q = Math.sqrt(-2*Math.log(p));
      return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    } else if (p <= phigh){
      q = p - 0.5; r = q*q;
      return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    } else {
      q = Math.sqrt(-2*Math.log(1-p));
      return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  }

  function rankAverage(values){
    const arr = values.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);
    const ranks = new Array(values.length);
    let i=0;
    while (i<arr.length){
      let j=i+1;
      while (j<arr.length && arr[j].v === arr[i].v) j++;
      const avg = (i + j - 1)/2 + 1;
      for (let k=i;k<j;k++) ranks[arr[k].i] = avg;
      i = j;
    }
    return ranks;
  }
  function pearson(x,y){
    const n = x.length; if (n<2) return NaN;
    const mx = x.reduce((a,b)=>a+b,0)/n, my = y.reduce((a,b)=>a+b,0)/n;
    let num=0, dx=0, dy=0;
    for (let i=0;i<n;i++){ const vx=x[i]-mx, vy=y[i]-my; num+=vx*vy; dx+=vx*vx; dy+=vy*vy; }
    return (dx>0 && dy>0) ? num/Math.sqrt(dx*dy) : NaN;
  }
  function spearman(x,y){ return pearson(rankAverage(x), rankAverage(y)); }

  function validateAndNormalize(pfasRows, sitesRows, reachRows){
    const validation = [];
    const columns = pfasRows.length ? Object.keys(pfasRows[0]) : [];
    const missing = REQUIRED_PFAS_COLUMNS.filter(c => !columns.includes(c));
    if (missing.length) throw new Error('PFAS CSV missing required columns: '+missing.join(', '));
    const seen = new Set();
    let duplicates = 0, invalidRL=0, invalidFlag=0, detectedMissing=0;
    const pfas = pfasRows.map((r, idx) => {
      const site_id = normalizeSiteId(r.site_id);
      const campaign = String(r.campaign ?? '').trim();
      const analyte_std = String(r.analyte_std ?? '').trim();
      const result_value = num(r.result_value);
      const rl = num(r.rl);
      const detect_flag = Number(r.detect_flag);
      const key = [campaign, site_id, analyte_std].join('||');
      if (seen.has(key)) duplicates++; else seen.add(key);
      if (!(detect_flag===0 || detect_flag===1)) invalidFlag++;
      if (!isFiniteNumber(rl) || rl <= 0) invalidRL++;
      if (detect_flag===1 && !isFiniteNumber(result_value)) detectedMissing++;
      return {...r, site_id, campaign, analyte_std, result_value, rl, detect_flag: detect_flag===1 ? 1 : 0, _row_index: idx+2};
    });
    if (duplicates) validation.push({severity:'warning', message:`Duplicate campaign-site-analyte rows: ${duplicates}. First values are used.`});
    if (invalidFlag) validation.push({severity:'error', message:`Invalid detect_flag values: ${invalidFlag}. Required values are 0/1.`});
    if (invalidRL) validation.push({severity:'error', message:`Missing or non-positive reporting limits: ${invalidRL}.`});
    if (detectedMissing) validation.push({severity:'error', message:`Detected rows with missing result_value: ${detectedMissing}.`});

    let sites = [];
    if (sitesRows && sitesRows.length){
      sites = sitesRows.map(r => {
        const site_id = normalizeSiteId(r.site_id ?? r['Site ID'] ?? r.SiteID ?? r.site);
        return {...r, site_id, Lat:num(r.Lat ?? r.latitude ?? r.Latitude), Lon:num(r.Lon ?? r.longitude ?? r.Longitude), Reach_group:r.Reach_group ?? r.reach_group ?? '', Bracket_role:r.Bracket_role ?? r.bracket_role ?? ''};
      }).filter(r=>r.site_id);
    } else {
      const unique = Array.from(new Set(pfas.map(r=>r.site_id))).sort();
      sites = unique.map(s=>({site_id:s, Lat:NaN, Lon:NaN, Reach_group:'', Bracket_role:''}));
      validation.push({severity:'info', message:`No sites metadata supplied; generated ${sites.length} sites from PFAS input.`});
    }
    const allowed = new Set(sites.map(r=>r.site_id));
    const included = pfas.filter(r => allowed.has(r.site_id));
    const excluded = pfas.filter(r => !allowed.has(r.site_id));
    if (excluded.length) validation.push({severity:'info', message:`Excluded ${excluded.length} PFAS rows from sites not present in site metadata.`});

    let reach = [];
    if (reachRows && reachRows.length){
      const rcols = Object.keys(reachRows[0]);
      const rmiss = REQUIRED_REACH_COLUMNS.filter(c => !rcols.includes(c));
      if (rmiss.length) validation.push({severity:'error', message:'Reach CSV missing required columns: '+rmiss.join(', ')});
      reach = reachRows.map(r => {
        const dn = normalizeSiteId(r.downstream_site);
        const ups = splitPipe(r.upstream_sites);
        const tribs = splitPipe(r.tributary_sites);
        const baselines = splitPipe(r.baseline_sites);
        const allSites = splitPipe(r.all_sites).length ? splitPipe(r.all_sites) : Array.from(new Set([dn, ...ups, ...tribs, ...baselines])).filter(Boolean);
        return {...r, experiment_id:String(r.experiment_id ?? '').trim(), reach_group:r.reach_group ?? '', experiment_type:String(r.experiment_type ?? '').trim(), downstream_site:dn, upstream_sites:ups.join('|'), tributary_sites:tribs.join('|'), baseline_sites:baselines.join('|'), all_sites:allSites.join('|')};
      });
    } else validation.push({severity:'info', message:'No reach table supplied; reach-contrast analysis skipped.'});

    if (validation.some(v=>v.severity==='error')) throw new Error(validation.filter(v=>v.severity==='error').map(v=>v.message).join('\n'));
    return {pfas:included, sites, reach, excludedRows:excluded, validation};
  }

  function runAnalysis(input, options={}){
    const nImpute = Number(options.nImpute ?? DEFAULT_N_IMPUTE);
    const seed = Number(options.seed ?? DEFAULT_SEED);
    const rng = makeRng(seed);
    const {pfas, sites, reach, excludedRows, validation} = validateAndNormalize(input.pfasRows, input.sitesRows, input.reachRows);
    const analytes = Array.from(new Set(pfas.map(r=>r.analyte_std))).sort();
    const sampleKeys = Array.from(new Set(pfas.map(r=>r.campaign+'||'+r.site_id))).sort();
    const sampleIndex = new Map(sampleKeys.map((k,i)=>[k,i]));
    const sampleRecords = sampleKeys.map(k=>{ const [campaign,site_id]=k.split('||'); return {campaign, site_id}; });
    const rowMap = new Map();
    pfas.forEach(r => {
      const key = [r.campaign,r.site_id,r.analyte_std].join('||');
      if (!rowMap.has(key)) rowMap.set(key,r);
    });
    const nSamples = sampleKeys.length;
    const totals = Array.from({length:nSamples}, () => new Float64Array(nImpute));
    const rl2Totals = new Float64Array(nSamples);
    const detectOnly = new Float64Array(nSamples);
    const detectCounts = new Int32Array(nSamples);
    const params = [];

    for (const analyte of analytes){
      const detVals = pfas.filter(r => r.analyte_std===analyte && r.detect_flag===1 && isFiniteNumber(r.result_value) && r.result_value>0).map(r=>r.result_value);
      let mode='uniform', mu=NaN, sigma=NaN;
      if (detVals.length >= 2){
        const logs = detVals.map(Math.log);
        mu = logs.reduce((a,b)=>a+b,0)/logs.length;
        const variance = logs.reduce((a,b)=>a+(b-mu)*(b-mu),0)/(logs.length-1);
        sigma = Math.sqrt(Math.max(variance, 0));
        if (sigma < 1e-6) sigma = 0.5;
        mode='lognorm';
      }
      params.push({analyte, n_detect:detVals.length, mode, mu_log:Number.isFinite(mu)?mu:'', sigma_log:Number.isFinite(sigma)?sigma:''});
      for (let si=0; si<nSamples; si++){
        const {campaign, site_id} = sampleRecords[si];
        const row = rowMap.get([campaign,site_id,analyte].join('||'));
        if (!row) continue;
        const detected = row.detect_flag === 1 && isFiniteNumber(row.result_value);
        if (detected){
          const val = row.result_value;
          detectCounts[si] += 1;
          rl2Totals[si] += val;
          detectOnly[si] += val;
          const arr = totals[si];
          for (let k=0;k<nImpute;k++) arr[k] += val;
        } else {
          const rl = isFiniteNumber(row.rl) && row.rl>0 ? row.rl : 0;
          rl2Totals[si] += rl/2;
          const arr = totals[si];
          if (mode === 'lognorm' && rl>0){
            const p = Math.max(1e-12, Math.min(1, normCdf((Math.log(Math.max(rl,1e-10)) - mu) / sigma)));
            for (let k=0;k<nImpute;k++){
              const u = rng() * p;
              let draw = Math.exp(mu + sigma*normInv(Math.max(1e-12, Math.min(1-1e-12, u))));
              if (draw > rl) draw = rl;
              if (draw < 0 || !Number.isFinite(draw)) draw = 0;
              arr[k] += draw;
            }
          } else {
            for (let k=0;k<nImpute;k++) arr[k] += rng()*rl;
          }
        }
      }
    }

    const sitesById = new Map(sites.map(s=>[s.site_id,s]));
    const siteTotals = sampleRecords.map((s, i) => {
      const med=median(totals[i]), lo=quantile(totals[i],0.025), hi=quantile(totals[i],0.975);
      const m = sitesById.get(s.site_id) || {};
      return {
        campaign:s.campaign, site_id:s.site_id, Reach_group:m.Reach_group||m.reach_group||'', Bracket_role:m.Bracket_role||m.bracket_role||'',
        Sigma_MI_median:med, Sigma_MI_p2_5:lo, Sigma_MI_p97_5:hi, Sigma_MI_95II:fmtInterval(med,lo,hi,1),
        Sigma_RL2:rl2Totals[i], Sigma_DetectOnly:detectOnly[i], n_detect:detectCounts[i], n_analytes:analytes.length,
        detection_frequency:`${detectCounts[i]}/${analytes.length}`
      };
    }).sort((a,b)=> a.campaign.localeCompare(b.campaign) || a.site_id.localeCompare(b.site_id));

    const distLookup = new Map(sampleRecords.map((s,i)=>[s.campaign+'||'+s.site_id, totals[i]]));
    const reachResults = analyzeReaches(reach, distLookup, siteTotals);
    const sensitivity = siteTotals.map(r => ({...r, AbsDiff_MI_minus_RL2:Math.abs(r.Sigma_MI_median-r.Sigma_RL2), Diff_MI_minus_RL2:r.Sigma_MI_median-r.Sigma_RL2, PctDiff_MI_minus_RL2:(r.Sigma_MI_median>0?(r.Sigma_MI_median-r.Sigma_RL2)/r.Sigma_MI_median*100:NaN)}));
    const rho = spearman(siteTotals.map(r=>r.Sigma_MI_median), siteTotals.map(r=>r.Sigma_RL2));
    const sensitivitySummary = [{metric:'Spearman_rho_MI_vs_RL2', value:rho}, {metric:'Max_abs_MI_minus_RL2_ngL', value:Math.max(...sensitivity.map(r=>r.AbsDiff_MI_minus_RL2))}, {metric:'Median_abs_MI_minus_RL2_ngL', value:median(sensitivity.map(r=>r.AbsDiff_MI_minus_RL2))}];
    const convergence = convergenceDiagnostics(totals, nImpute);
    const validationRows = validation.concat([{severity:'info', message:`Included PFAS rows: ${pfas.length}`},{severity:'info', message:`Included site-campaign samples: ${nSamples}`},{severity:'info', message:`Included sites: ${new Set(siteTotals.map(r=>r.site_id)).size}`},{severity:'info', message:`Included analytes: ${analytes.length}`}]);
    const report = buildValidationReport(validationRows, siteTotals, reachResults, sensitivitySummary, convergence, excludedRows, nImpute, seed);
    return {siteTotals, reachResults, sensitivity, sensitivitySummary, convergence, params, excludedRows, validationRows, report, totals, sampleRecords, analytes, nImpute, seed};
  }

  function summarizeDistribution(arr, d=1){
    const med=median(arr), lo=quantile(arr,0.025), hi=quantile(arr,0.975);
    return {med, lo, hi, text:fmtInterval(med,lo,hi,d)};
  }
  function ftrib(up, trib, dn){
    const vals=[]; let valid=0;
    for (let i=0;i<dn.length;i++){
      const denom = trib[i] - up[i];
      const f = (dn[i] - up[i]) / denom;
      if (Number.isFinite(f) && denom > 0 && f >= 0 && f <= 1){ vals.push(f); valid++; }
    }
    const validRate = dn.length ? valid/dn.length : 0;
    if (validRate < 0.95){
      const mu=median(up), mt=median(trib), md=median(dn);
      let reason='n/a: outside two-endmember range';
      if (md > mt) reason = 'n/a: C_dn > C_trib';
      else if (md < mu) reason = 'n/a: C_dn < C_up';
      else if (mt <= mu) reason = 'n/a: C_trib <= C_up';
      return {med:NaN, lo:NaN, hi:NaN, note:`${reason}; valid=${validRate.toFixed(2)}`};
    }
    const s = summarizeDistribution(vals, 3);
    return {med:s.med, lo:s.lo, hi:s.hi, note:`f_trib=${s.text}; valid=${validRate.toFixed(2)}`};
  }
  function analyzeReaches(reach, distLookup, siteTotals){
    if (!reach || !reach.length) return [];
    const campaigns = Array.from(new Set(siteTotals.map(r=>r.campaign))).sort();
    const out=[];
    for (const r of reach){
      const expType = String(r.experiment_type || '').trim();
      const dn = normalizeSiteId(r.downstream_site);
      const ups = splitPipe(r.upstream_sites), tribs = splitPipe(r.tributary_sites), baselines = splitPipe(r.baseline_sites);
      const allSites = splitPipe(r.all_sites).length ? splitPipe(r.all_sites) : Array.from(new Set([dn,...ups,...tribs,...baselines])).filter(Boolean);
      for (const campaign of campaigns){
        const present = new Set(siteTotals.filter(s=>s.campaign===campaign).map(s=>s.site_id));
        if (!allSites.every(s=>present.has(s))) continue;
        const dnDist = distLookup.get(campaign+'||'+dn); if (!dnDist) continue;
        const primaryUp = ups.length ? ups[0] : '';
        const upDist = primaryUp ? distLookup.get(campaign+'||'+primaryUp) : null;
        let delta = null, ds={med:NaN,lo:NaN,hi:NaN,text:''};
        if (upDist){
          delta = new Float64Array(dnDist.length);
          for (let i=0;i<dnDist.length;i++) delta[i] = dnDist[i] - upDist[i];
          ds = summarizeDistribution(delta,1);
        }
        let diagType='', f={med:NaN,lo:NaN,hi:NaN,note:''};
        if (expType === 'tributary_increment' && upDist && tribs.length){
          const tribDist = distLookup.get(campaign+'||'+tribs[0]);
          if (tribDist && ds.med >= 0){ f = ftrib(upDist, tribDist, dnDist); diagType='f_trib_concentration_space'; }
          else { diagType='segment_delta'; f.note='not reported for negative segment delta; use ΔΣPFAS40'; }
        } else if (expType === 'confluence_mix'){
          const endmembers = [...ups, ...tribs];
          const dists = endmembers.map(s=>distLookup.get(campaign+'||'+s));
          if (dists.length >= 2 && dists.every(Boolean)){
            let outside=0, fvals=[];
            for (let i=0;i<dnDist.length;i++){
              let lo=Infinity, hi=-Infinity;
              for (const d of dists){ if (d[i]<lo) lo=d[i]; if (d[i]>hi) hi=d[i]; }
              if (dnDist[i] < lo || dnDist[i] > hi) outside++;
              const f0 = (dnDist[i] - lo)/(hi-lo);
              if (Number.isFinite(f0) && f0>=0 && f0<=1) fvals.push(f0);
            }
            const fs = fvals.length ? summarizeDistribution(fvals,3) : {med:NaN,lo:NaN,hi:NaN,text:'n/a'};
            f={med:fs.med,lo:fs.lo,hi:fs.hi,note:`outside_range_prob=${(outside/dnDist.length).toFixed(2)}; apparent_mix_fraction=${fs.text}; endmembers=${endmembers.join('|')}`};
            diagType='confluence_range_check';
          }
        } else if (expType === 'source_reach') { diagType='linear_segment_delta'; f.note='not applicable to source_reach; use ΔΣPFAS40'; }
        else if (expType === 'downstream_propagation') { diagType='downstream_propagation_delta'; f.note='not applicable to downstream_propagation; use ΔΣPFAS40'; }
        const base = {campaign, experiment_id:r.experiment_id, reach_group:r.reach_group||'', experiment_type:expType, experiment_name:r.experiment_name||'', contrast_role:'primary', downstream_site:dn, upstream_sites:ups.join('|'), tributary_sites:tribs.join('|'), baseline_sites:baselines.join('|'), all_sites:allSites.join('|'), complete_required_sites:true, delta_SigmaPFAS40_median:ds.med, delta_SigmaPFAS40_p2_5:ds.lo, delta_SigmaPFAS40_p97_5:ds.hi, delta_SigmaPFAS40_95II:ds.text, diagnostic_type:diagType, f_trib_or_mix_median:f.med, f_trib_or_mix_p2_5:f.lo, f_trib_or_mix_p97_5:f.hi, diagnostic_note:f.note};
        out.push(base);
        for (const b of baselines){
          const bdist = distLookup.get(campaign+'||'+b); if (!bdist) continue;
          const bd = new Float64Array(dnDist.length); for (let i=0;i<bd.length;i++) bd[i]=dnDist[i]-bdist[i];
          const bs = summarizeDistribution(bd,1);
          out.push({...base, experiment_id:`${r.experiment_id}__baseline_${b}`, experiment_type:`${expType}_baseline`, experiment_name:`${r.experiment_name||''} - cumulative vs baseline ${b}`, contrast_role:'cumulative_vs_baseline', upstream_sites:b, tributary_sites:'', delta_SigmaPFAS40_median:bs.med, delta_SigmaPFAS40_p2_5:bs.lo, delta_SigmaPFAS40_p97_5:bs.hi, delta_SigmaPFAS40_95II:bs.text, diagnostic_type:'baseline_context_delta', f_trib_or_mix_median:NaN, f_trib_or_mix_p2_5:NaN, f_trib_or_mix_p97_5:NaN, diagnostic_note:'cumulative downstream-minus-baseline contrast; not a mixing fraction'});
        }
      }
    }
    return out;
  }
  function convergenceDiagnostics(totals, nImpute){
    const Ns=[250,500,1000,2000,5000,7500,10000].filter(N=>N<nImpute);
    const ref = totals.map(a=>({med:median(a), lo:quantile(a,0.025), hi:quantile(a,0.975)}));
    return Ns.map(N=>{
      let maxMed=0, maxLo=0, maxHi=0, medDiffs=[];
      totals.forEach((a,i)=>{
        const slice = Array.from(a.slice(0,N));
        const medN=median(slice), loN=quantile(slice,0.025), hiN=quantile(slice,0.975);
        const dm=Math.abs(medN-ref[i].med), dl=Math.abs(loN-ref[i].lo), dh=Math.abs(hiN-ref[i].hi);
        if (dm>maxMed) maxMed=dm; if (dl>maxLo) maxLo=dl; if (dh>maxHi) maxHi=dh; medDiffs.push(dm);
      });
      return {N, max_abs_diff_median_ngL:maxMed, median_abs_diff_median_ngL:median(medDiffs), max_abs_diff_II95_lo_ngL:maxLo, max_abs_diff_II95_hi_ngL:maxHi};
    });
  }
  function buildValidationReport(rows, siteTotals, reachResults, sensSummary, convergence, excludedRows, nImpute, seed){
    const maxDiff = sensSummary.find(r=>r.metric==='Max_abs_MI_minus_RL2_ngL')?.value;
    const rho = sensSummary.find(r=>r.metric==='Spearman_rho_MI_vs_RL2')?.value;
    const lines = [];
    lines.push('PFAS Reach Bracketing Tool validation report');
    lines.push(`Generated: ${new Date().toISOString()}`);
    lines.push(`Version: 1.1.1-browser-python-bundle`);
    lines.push(`MI replicates: ${nImpute}`);
    lines.push(`Random seed: ${seed}`);
    lines.push(`Site-campaign samples: ${siteTotals.length}`);
    lines.push(`Unique sites: ${new Set(siteTotals.map(r=>r.site_id)).size}`);
    lines.push(`Reach-campaign output rows: ${reachResults.length}`);
    lines.push(`Excluded rows outside site metadata: ${excludedRows.length}`);
    lines.push(`Spearman rho MI vs RL/2: ${Number.isFinite(rho)?rho.toFixed(4):'n/a'}`);
    lines.push(`Max |MI - RL/2|: ${Number.isFinite(maxDiff)?maxDiff.toFixed(3):'n/a'} ng L-1`);
    if (convergence.length){
      const last = convergence[convergence.length-1];
      lines.push(`Convergence last N=${last.N}: max |Δ median|=${last.max_abs_diff_median_ngL.toFixed(3)} ng L-1`);
    }
    const tw91 = siteTotals.find(r=>r.campaign==='Campaign 2' && r.site_id==='TW91');
    if (tw91) lines.push(`Campaign 2 TW91 ΣPFAS40 median=${tw91.Sigma_MI_median.toFixed(1)} ng L-1`);
    const wf1 = reachResults.find(r=>r.campaign==='Campaign 2' && r.experiment_id==='WF1_FarmersBranch');
    if (wf1) lines.push(`WF1_FarmersBranch C2 ΔΣPFAS40=${wf1.delta_SigmaPFAS40_median.toFixed(1)} ng L-1; ${wf1.diagnostic_note}`);
    lines.push(''); lines.push('Input validation messages:');
    rows.forEach(r=>lines.push(`- ${r.severity}: ${r.message}`));
    lines.push(''); lines.push('f_trib definition: apparent concentration-space tributary fraction under conservative two-endmember mixing; not a flow, source-apportionment, or load fraction.');
    return lines.join('\n');
  }

  // ────────────────────────────────────────────────────────────────────────
  // Rendering helpers
  function renderTable(el, rows, maxRows=200){
    if (!rows || !rows.length){ el.innerHTML = '<p class="small">No rows to display.</p>'; return; }
    const headers = Object.keys(rows[0]);
    const display = rows.slice(0,maxRows);
    let html = '<table><thead><tr>'+headers.map(h=>`<th>${escapeHtml(h)}</th>`).join('')+'</tr></thead><tbody>';
    html += display.map(r=>'<tr>'+headers.map(h=>`<td>${escapeHtml(formatCell(r[h]))}</td>`).join('')+'</tr>').join('');
    html += '</tbody></table>';
    if (rows.length > maxRows) html += `<p class="small">Showing first ${maxRows} of ${rows.length} rows. Download CSV for full table.</p>`;
    el.innerHTML = html;
  }
  function formatCell(v){
    if (typeof v === 'number') return Number.isFinite(v) ? (Math.abs(v)>=1000 ? v.toFixed(1) : v.toFixed(3).replace(/\.000$/,'').replace(/(\.\d*?)0+$/,'$1')) : '';
    return v ?? '';
  }
  function metric(label, value){ return `<div class="metric"><div class="label">${escapeHtml(label)}</div><div class="value">${escapeHtml(value)}</div></div>`; }
  function renderSummary(result){
    const tw91 = result.siteTotals.find(r=>r.campaign==='Campaign 2' && r.site_id==='TW91');
    const wf1 = result.reachResults.find(r=>r.campaign==='Campaign 2' && r.experiment_id==='WF1_FarmersBranch');
    const rho = result.sensitivitySummary.find(r=>r.metric==='Spearman_rho_MI_vs_RL2')?.value;
    const maxDiff = result.sensitivitySummary.find(r=>r.metric==='Max_abs_MI_minus_RL2_ngL')?.value;
    const wf1F = wf1 && Number.isFinite(wf1.f_trib_or_mix_median) ? `${wf1.f_trib_or_mix_median.toFixed(3)} [${wf1.f_trib_or_mix_p2_5.toFixed(3)}, ${wf1.f_trib_or_mix_p97_5.toFixed(3)}]` : 'n/a';
    return [
      metric('Site-campaign samples', String(result.siteTotals.length)),
      metric('Analytes', String(result.analytes.length)),
      metric('TW91 C2 ΣPFAS40', tw91 ? `${tw91.Sigma_MI_median.toFixed(1)} ng L-1` : 'n/a'),
      metric('WF1 ΔΣPFAS40', wf1 ? `${wf1.delta_SigmaPFAS40_median.toFixed(1)} ng L-1` : 'n/a'),
      metric('WF1 ftrib', wf1F),
      metric('MI vs RL/2 Spearman ρ', Number.isFinite(rho) ? rho.toFixed(4) : 'n/a'),
      metric('Max |MI − RL/2|', Number.isFinite(maxDiff) ? `${maxDiff.toFixed(3)} ng L-1` : 'n/a'),
      metric('MI replicates', result.nImpute.toLocaleString()),
      metric('Seed', String(result.seed))
    ].join('');
  }
  function plotSiteTotals(rows){
    const sorted = [...rows].sort((a,b)=>b.Sigma_MI_median-a.Sigma_MI_median);
    const w=900, h=420, ml=80, mr=20, mt=20, mb=100;
    const max = Math.max(...sorted.map(r=>r.Sigma_MI_p97_5));
    const minPositive = Math.min(...sorted.map(r=>Math.max(0.1,r.Sigma_MI_p2_5)).filter(v=>v>0));
    const log = max/minPositive > 50;
    const y = val => {
      if (log){ const a=Math.log10(minPositive), b=Math.log10(max*1.15); return mt+(h-mt-mb)*(1-(Math.log10(Math.max(val,minPositive))-a)/(b-a)); }
      return mt+(h-mt-mb)*(1-val/(max*1.15));
    };
    const bw=(w-ml-mr)/sorted.length*0.72;
    let svg=`<svg viewBox="0 0 ${w} ${h}" width="100%" height="${h}" role="img" aria-label="Site totals">`;
    svg += `<line x1="${ml}" y1="${h-mb}" x2="${w-mr}" y2="${h-mb}" stroke="#111827"/><line x1="${ml}" y1="${mt}" x2="${ml}" y2="${h-mb}" stroke="#111827"/>`;
    const ticks = log ? [10,100,1000,10000] : [0,max/4,max/2,3*max/4,max];
    ticks.filter(t=>t>=minPositive && t<=max*1.15).forEach(t=>{ const yy=y(t); svg += `<line x1="${ml-4}" y1="${yy}" x2="${w-mr}" y2="${yy}" stroke="#e5e7eb"/><text x="${ml-8}" y="${yy+4}" text-anchor="end" font-size="11">${t>=1000?t.toFixed(0):t.toFixed(1).replace(/\.0$/,'')}</text>`; });
    sorted.forEach((r,i)=>{ const x=ml+i*(w-ml-mr)/sorted.length+(w-ml-mr)/sorted.length/2; const y0=y(r.Sigma_MI_median), yb=y(0.1); const color = r.Bracket_role && r.Bracket_role.toLowerCase().includes('tributary') ? '#6b21a8' : (r.Bracket_role && r.Bracket_role.toLowerCase().includes('downstream') ? '#ea580c' : '#2563eb'); svg += `<rect x="${x-bw/2}" y="${y0}" width="${bw}" height="${Math.max(1,yb-y0)}" fill="${color}" opacity="0.82"/>`; const lo=y(r.Sigma_MI_p2_5), hi=y(r.Sigma_MI_p97_5); svg += `<line x1="${x}" y1="${hi}" x2="${x}" y2="${lo}" stroke="#111827"/><line x1="${x-4}" y1="${hi}" x2="${x+4}" y2="${hi}" stroke="#111827"/><line x1="${x-4}" y1="${lo}" x2="${x+4}" y2="${lo}" stroke="#111827"/>`; svg += `<text transform="translate(${x},${h-mb+14}) rotate(60)" font-size="10" text-anchor="start">${escapeHtml(r.site_id)} ${escapeHtml(r.campaign.replace('Campaign ', 'C'))}</text>`; });
    svg += `<text x="${ml}" y="14" font-size="13" font-weight="700">ΣPFAS40 medians with 95% imputation intervals (${UNIT})${log?' — log scale':''}</text></svg>`;
    return svg;
  }
  
function plotSensitivity(rows){
    const sorted = [...rows].sort((a,b)=>b.Sigma_MI_median-a.Sigma_MI_median);
    const n = sorted.length;
    const rowH = 34, w = 900, mt = 78, mb = 68, ml = 128, mr = 125;
    const ph = Math.max(1, n * rowH);
    const h = mt + ph + mb;
    const vals = [];
    rows.forEach(r=>{ if (r.Sigma_MI_median>0) vals.push(r.Sigma_MI_median); if (r.Sigma_RL2>0) vals.push(r.Sigma_RL2); });
    const minVal = Math.max(0.1, Math.min(...vals) * 0.45);
    const maxVal = Math.max(...vals) * 1.8;
    const logMin = Math.log10(minVal), logMax = Math.log10(maxVal);
    const x = v => ml + (Math.log10(Math.max(v, minVal)) - logMin) / (logMax - logMin) * (w - ml - mr);
    const rho = spearman(rows.map(r=>r.Sigma_MI_median), rows.map(r=>r.Sigma_RL2));
    let svg = `<svg id="rsv" viewBox="0 0 ${w} ${h}" width="100%" height="${h}" role="img" aria-label="Censoring sensitivity comparing MI median and RL/2 substitution">`;
    svg += `<rect x="0" y="0" width="${w}" height="${h}" fill="white"/>`;
    svg += `<text x="${w/2}" y="25" text-anchor="middle" font-size="17" font-weight="800" fill="#4d1979">Censoring sensitivity: MI median vs RL/2 substitution</text>`;
    svg += `<text x="${w/2}" y="45" text-anchor="middle" font-size="12" fill="#5b6472">Spearman ρ = ${Number.isFinite(rho)?rho.toFixed(4):'n/a'} | longer bars indicate larger method disagreement</text>`;
    const legendY = mt - 14;
    svg += `<circle cx="${ml+6}" cy="${legendY}" r="6.2" fill="#4d1979" stroke="white" stroke-width="1.4"/>`;
    svg += `<text x="${ml+20}" y="${legendY+4}" font-size="11" font-weight="700" fill="#111827">MI median</text>`;
    svg += `<circle cx="${ml+117}" cy="${legendY}" r="5.2" fill="white" stroke="#8c95a5" stroke-width="2"/>`;
    svg += `<text x="${ml+131}" y="${legendY+4}" font-size="11" font-weight="700" fill="#111827">RL/2</text>`;
    svg += `<rect x="${ml+176}" y="${legendY-5}" width="24" height="10" rx="4" fill="#e8d5f5"/>`;
    svg += `<text x="${ml+207}" y="${legendY+4}" font-size="11" font-weight="700" fill="#111827">Gap</text>`;
    svg += `<text x="${w-mr+18}" y="${legendY-1}" font-size="10" font-weight="800" fill="#4d1979">% diff</text>`;
    svg += `<text x="${w-mr+66}" y="${legendY-1}" font-size="10" font-weight="800" fill="#4d1979">Δ ng L-1</text>`;
    const ticks = [1,10,100,1000,10000,100000].filter(t=>t>=minVal && t<=maxVal);
    ticks.forEach(t=>{
      const xx = x(t);
      svg += `<line x1="${xx}" y1="${mt}" x2="${xx}" y2="${mt+ph}" stroke="#e5e7eb" stroke-width="0.8"/>`;
      svg += `<text x="${xx}" y="${mt+ph+22}" text-anchor="middle" font-size="11" font-weight="700" fill="#5b6472">${t>=1000?(t/1000)+'k':t}</text>`;
    });
    svg += `<text x="${(ml+w-mr)/2}" y="${h-12}" text-anchor="middle" font-size="13" font-weight="800" fill="#374151">ΣPFAS40 (ng L-1, log scale)</text>`;
    sorted.forEach((r,i)=>{
      const y = mt + i*rowH + rowH/2;
      if (i%2===0) svg += `<rect x="${ml}" y="${y-rowH/2}" width="${w-ml-mr}" height="${rowH}" fill="#f3edf9" opacity="0.30"/>`;
      const role = r.Bracket_role || '';
      const roleCol = roleColor(role);
      svg += `<text x="${ml-10}" y="${y-2}" text-anchor="end" font-size="11.5" font-weight="800" fill="${roleCol}">${escapeHtml(r.site_id)}</text>`;
      svg += `<text x="${ml-10}" y="${y+11}" text-anchor="end" font-size="8.8" fill="#8c95a5">${escapeHtml(String(r.campaign).replace('Campaign ','C'))}</text>`;
      const xMI = x(r.Sigma_MI_median), xRL = x(r.Sigma_RL2);
      const lo = Math.min(xMI, xRL), hi = Math.max(xMI, xRL);
      const diff = r.Sigma_MI_median - r.Sigma_RL2;
      const absDiff = Math.abs(diff);
      const absPct = r.Sigma_MI_median>0 ? absDiff / r.Sigma_MI_median * 100 : 0;
      const gapColor = absPct < 0.5 ? '#d1fae5' : absPct < 2 ? '#e8d5f5' : absPct < 5 ? '#c084fc' : '#9333ea';
      if ((hi-lo) > 2) svg += `<rect x="${lo}" y="${y-6}" width="${hi-lo}" height="12" rx="6" fill="${gapColor}" opacity="0.60"/>`;
      svg += `<circle cx="${xRL}" cy="${y}" r="5.4" fill="white" stroke="#8c95a5" stroke-width="2.1"><title>${escapeHtml(r.site_id)} ${escapeHtml(r.campaign)} RL/2 ${r.Sigma_RL2.toFixed(2)}</title></circle>`;
      svg += `<circle cx="${xMI}" cy="${y}" r="6.4" fill="#4d1979" stroke="white" stroke-width="1.5"><title>${escapeHtml(r.site_id)} ${escapeHtml(r.campaign)} MI ${r.Sigma_MI_median.toFixed(2)}</title></circle>`;
      const pctCol = absPct < 1 ? '#16a34a' : absPct < 3 ? '#5b6472' : '#9333ea';
      svg += `<text x="${w-mr+18}" y="${y+4}" font-size="11.2" font-weight="800" fill="${pctCol}">${absPct.toFixed(1)}%</text>`;
      svg += `<text x="${w-mr+68}" y="${y+4}" font-size="10.5" fill="#64748b">${absDiff.toFixed(1)}</text>`;
    });
    svg += `</svg>`;
    return svg;
  }

  function roleColor(role){
    const r = String(role||'').toLowerCase();
    if (r.includes('tributary')) return '#6b21a8';
    if (r.includes('downstream') || r.includes('integrator') || r.includes('outlet') || r.includes('post-')) return '#ea580c';
    if (r.includes('upstream') || r.includes('baseline') || r.includes('pre-source')) return '#2563eb';
    return '#5b6472';
  }

  function howToReadSiteTotals(){
    return `<div class="chart-note"><strong>How to read this chart:</strong> Each bar is the censoring-aware MI median ΣPFAS40 for one site-campaign sample, and the whisker shows the 95% imputation interval. Colors indicate site role: blue for upstream/baseline, purple for tributary end-member, and orange for downstream/integrator positions. Taller bars identify priority locations, while wider intervals indicate greater uncertainty from non-detects.</div>`;
  }

  function sensitivityStatsHtml(rows){
    const rho = spearman(rows.map(r=>r.Sigma_MI_median), rows.map(r=>r.Sigma_RL2));
    const diffs = rows.map(r=>Math.abs(r.Sigma_MI_median - r.Sigma_RL2));
    const pcts = rows.map(r=>r.Sigma_MI_median>0 ? Math.abs(r.Sigma_MI_median-r.Sigma_RL2)/r.Sigma_MI_median*100 : 0);
    return `<div class="summary-grid compact-summary">${metric('Spearman ρ', Number.isFinite(rho)?rho.toFixed(4):'n/a')}${metric('Max |MI − RL/2|', `${Math.max(...diffs).toFixed(3)} ng L-1`)}${metric('Median |MI − RL/2|', `${median(diffs).toFixed(3)} ng L-1`)}${metric('Max % difference', `${Math.max(...pcts).toFixed(2)}%`)}</div>`;
  }

  function howToReadSensitivity(){
    return `<div class="chart-note"><strong>How to read this chart:</strong> Each row compares two ways of handling non-detected PFAS for the same sample. The filled purple dot is the MI median; the open gray circle is the simpler RL/2 substitution. The colored bar between them shows the gap: green indicates negligible disagreement (&lt;1%), while purple indicates larger disagreement. Samples with many non-detects tend to show wider gaps because more analyte values are imputed. Samples dominated by detections show near-zero gaps because both methods use the same measured concentrations.</div>`;
  }

  // UI
  function $(id){ return typeof document !== 'undefined' ? document.getElementById(id) : null; }
  async function readFileInput(input){
    const f = input.files && input.files[0];
    return f ? await f.text() : '';
  }
  function setStatus(msg, type='info'){ const el=$('status'); if (el){ el.className='status '+type; el.textContent=msg; } }
  async function runFromInputs(useDemo=false){
    const runBtn = $('runBtn'), demoBtn = $('demoBtn');
    try{
      if (runBtn) runBtn.disabled = true;
      if (demoBtn) demoBtn.disabled = true;
      setStatus('Running analysis. This may take a few seconds for 5,000 imputations...', 'info');
      await new Promise(r=>setTimeout(r,25));
      let pfasText='', sitesText='', reachText='';
      if (useDemo){
        if (!window.PFAS_DEMO_DATA) throw new Error('Demo data are not available. Confirm that data/demo_data.js is present or use the standalone HTML file.');
        pfasText=window.PFAS_DEMO_DATA.pfas; sitesText=window.PFAS_DEMO_DATA.sites; reachText=window.PFAS_DEMO_DATA.reach;
      } else {
        pfasText = await readFileInput($('pfasFile'));
        sitesText = await readFileInput($('sitesFile'));
        reachText = await readFileInput($('reachFile'));
        if (!pfasText) throw new Error('Please upload a PFAS Results CSV or use demo data.');
      }
      const input = {pfasRows:parseCSV(pfasText), sitesRows:sitesText?parseCSV(sitesText):[], reachRows:reachText?parseCSV(reachText):[]};
      const nImputeEl = $('nImpute'), seedEl = $('seed');
      const result = runAnalysis(input, {nImpute:Number(nImputeEl ? nImputeEl.value : DEFAULT_N_IMPUTE), seed:Number(seedEl ? seedEl.value : DEFAULT_SEED)});
      window.PFAS_LAST_RESULT = result;
      renderResult(result);
      setStatus('Analysis complete. Results are shown below and can be downloaded as CSV files.', 'ok');
      const card = $('resultsCard');
      if (card && card.scrollIntoView) card.scrollIntoView({behavior:'smooth', block:'start'});
    } catch(e){
      console.error(e);
      setStatus(e.message || String(e), 'err');
      const card = $('resultsCard');
      const placeholder = $('resultsPlaceholder');
      if (card) card.classList.remove('hidden');
      if (placeholder) { placeholder.className='status err'; placeholder.textContent = e.message || String(e); }
    } finally {
      if (runBtn) runBtn.disabled = false;
      if (demoBtn) demoBtn.disabled = false;
    }
  }
  function renderResult(result){
    const card = $('resultsCard');
    if (card) card.classList.remove('hidden');
    const placeholder = $('resultsPlaceholder');
    if (placeholder) placeholder.style.display = 'none';
    $('summary').innerHTML = renderSummary(result);
    renderTable($('siteTable'), result.siteTotals.map(cleanRow));
    renderTable($('reachTable'), result.reachResults.map(cleanRow));
    renderTable($('sensitivityTable'), result.sensitivity.map(cleanRow));
    renderTable($('convergenceTable'), result.convergence.map(cleanRow));
    renderTable($('paramsTable'), result.params.map(cleanRow));
    renderTable($('validationTable'), result.validationRows);
    $('validationReport').textContent = result.report;
    $('figSiteTotals').innerHTML = plotSiteTotals(result.siteTotals) + howToReadSiteTotals();
    $('figSensitivity').innerHTML = plotSensitivity(result.sensitivity) + sensitivityStatsHtml(result.sensitivity) + howToReadSensitivity();
    const downloads = [
      ['SiteTotals.csv', toCSV(result.siteTotals.map(cleanRow))],
      ['ReachExperiments.csv', toCSV(result.reachResults.map(cleanRow))],
      ['Sensitivity_RL2.csv', toCSV(result.sensitivity.map(cleanRow))],
      ['Convergence.csv', toCSV(result.convergence.map(cleanRow))],
      ['ImputationParams.csv', toCSV(result.params.map(cleanRow))],
      ['ExcludedRows.csv', toCSV(result.excludedRows.map(cleanRow))],
      ['validation_report.txt', result.report],
      ['results.json', JSON.stringify(result, (k,v)=> (v instanceof Float64Array ? Array.from(v) : v), 2)]
    ];
    $('downloadRow').innerHTML = '';
    downloads.forEach(([name, content])=>{
      const b=document.createElement('button'); b.className='light'; b.textContent='⬇ '+name; b.addEventListener('click',()=>downloadText(name, content, name.endsWith('.json')?'application/json':(name.endsWith('.csv')?'text/csv':'text/plain'))); $('downloadRow').appendChild(b);
    });
  }
  function cleanRow(row){
    const out={};
    for (const [k,v] of Object.entries(row)){
      if (v instanceof Float64Array) continue;
      out[k] = typeof v === 'number' ? (Number.isFinite(v) ? v : '') : v;
    }
    return out;
  }
  function initUi(){
    const runBtn = $('runBtn'), demoBtn = $('demoBtn'), resetBtn = $('resetBtn');
    if (!runBtn || !demoBtn) return;
    runBtn.addEventListener('click',()=>runFromInputs(false));
    demoBtn.addEventListener('click',()=>runFromInputs(true));
    if (resetBtn) resetBtn.addEventListener('click',()=>{
      ['pfasFile','sitesFile','reachFile'].forEach(id=>{ const el=$(id); if (el) el.value=''; });
      const placeholder = $('resultsPlaceholder');
      if (placeholder) { placeholder.style.display='block'; placeholder.className='status info'; placeholder.textContent='Results will appear here after running the analysis or manuscript demo data.'; }
      ['summary','downloadRow','siteTable','reachTable','sensitivityTable','convergenceTable','paramsTable','validationTable','figSiteTotals','figSensitivity'].forEach(id=>{ const el=$(id); if (el) el.innerHTML=''; });
      const vr = $('validationReport'); if (vr) vr.textContent='';
      setStatus('Upload PFAS CSV or run the demo data. All computation is local to your browser.','info');
    });
    document.querySelectorAll('.tab').forEach(tab=>tab.addEventListener('click',()=>{
      document.querySelectorAll('.tab').forEach(t=>t.classList.remove('active'));
      document.querySelectorAll('.tab-panel').forEach(p=>p.classList.remove('active'));
      tab.classList.add('active');
      const panel = $('panel-'+tab.dataset.tab);
      if (panel) panel.classList.add('active');
    }));
    try {
      const url = new URL(window.location.href);
      if (url.searchParams.get('demo') === '1' || url.hash === '#demo') {
        setTimeout(()=>runFromInputs(true), 100);
      }
    } catch(e) { /* ignore URL parsing errors for local file paths */ }
  }

  const core = {parseCSV,toCSV,runAnalysis,normalizeSiteId,splitPipe,ftrib,quantile,median,spearman};
  if (typeof window !== 'undefined') {
    window.PFASReachToolCore = core;
    window.PFASReachToolRunDemo = ()=>runFromInputs(true);
    window.PFASReachToolRun = ()=>runFromInputs(false);
    if (document.readyState === 'loading') window.addEventListener('DOMContentLoaded', initUi);
    else initUi();
  }
  if (typeof module !== 'undefined' && module.exports) module.exports = core;
})();
