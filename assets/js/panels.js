/* StarNest panels (Phase 2): front ends for the Phase-1 Python endpoints. Loaded after app.js (uses its globals: probeBackend, apiCall, apiOpen,
   getTaskFilename, showToast, getPlotlyLayout, PLOTLY_CONFIG, Plotly).
   TERMINOLOGY  panel = one card mounted into a <div id="panel-..."> placeholder of index.html; backend-only = needs the Python server
                (the in-browser engine has no equivalent), so without it the card shows a notice and its buttons do nothing.
   WORKFLOW     mount(): inject the card HTML -> refreshNotices(): probe the backend, show / hide the notice -> button click ->
                call(): probe, POST/GET, render tiles / table / plot -> errors are shown inside the card (and as a toast). */
(function () {
  'use strict';

  // ===================== GLOBAL CONFIGURATION =====================
  const CFG = {
    BORDERLINE_HALF_WIDTH: 0.10,   // |PMemb - cutoff| below this counts as "borderline"
    UNSTABLE_FRAC: 0.9,            // mirrors main.py: member in < 90 % of the resampling draws
    TABLE_ROWS: 15,                // rows shown in the on-page tables (full tables are downloadable / in the API reply)
    DEFAULT_N_RESAMPLE: 100, DEFAULT_MIN_CLUSTER: 10,
    DEFAULT_F_BIN: 0.5, DEFAULT_N_MC: 50, DEFAULT_N_SYNTH: 10000,
    DEFAULT_ETA: 10, DEFAULT_N_BOOT: 500, DEFAULT_N_RANDOM: 100,
    PLOT_HEIGHT: 340,
    TIMEOUT_MS: 600000             // MC mass assignment / bootstrap can take minutes
  };
  // =================================================================

  const $ = (id) => document.getElementById(id);
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const num = (v, d = 3) => (v === null || v === undefined || !Number.isFinite(Number(v))) ? '—' : Number(v).toFixed(d);
  const pm = (v, e, d = 3) => `${num(v, d)} <span class="text-gray-500">&plusmn; ${num(e, d)}</span>`;
  const INPUT = 'bg-white/[0.04] border border-white/[0.08] rounded-lg px-2 py-1.5 text-sm text-white w-full';
  const BTN = 'py-2 px-4 bg-gradient-to-r from-accent-cyan to-accent-violet rounded-xl text-sm font-semibold hover:shadow-[0_0_20px_rgba(0,212,255,0.3)] transition-all duration-300 flex items-center justify-center gap-2';
  const BTN2 = 'py-2 px-3 bg-white/[0.04] hover:bg-white/[0.08] border border-white/[0.1] rounded-lg text-xs font-medium text-gray-200 transition-all';

  const tile = (label, value, sub = '') => `<div class="rounded-xl p-3 border border-white/[0.06] bg-white/[0.03]"><div class="text-[10px] uppercase tracking-wider text-gray-500">${label}</div><div class="text-lg font-semibold text-white mt-0.5">${value}</div>${sub ? `<div class="text-[11px] text-gray-400 mt-0.5">${sub}</div>` : ''}</div>`;
  const tiles = (arr) => `<div class="grid grid-cols-2 md:grid-cols-4 gap-3">${arr.join('')}</div>`;
  const warn = (msgs, tone = 'amber') => msgs.length ? `<div class="rounded-lg p-3 border ${tone === 'red' ? 'border-red-500/40 bg-red-500/[0.07] text-red-200' : 'border-amber-500/40 bg-amber-500/[0.07] text-amber-200'} text-xs space-y-1">${msgs.map((m) => `<div><i class="fas fa-triangle-exclamation mr-1"></i>${esc(m)}</div>`).join('')}</div>` : '';
  const table = (cols, rows) => `<div class="overflow-x-auto"><table class="w-full text-[11px] text-gray-300"><thead><tr class="text-gray-500 text-left">${cols.map((c) => `<th class="py-1 pr-3 font-medium">${esc(c.h)}</th>`).join('')}</tr></thead><tbody>${rows.map((r) => `<tr class="border-t border-white/[0.05]">${cols.map((c) => `<td class="py-1 pr-3">${c.f(r)}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`;

  function card(id, icon, title, subtitle, body) {
    return `<div id="${id}" class="rounded-2xl p-4 border border-white/[0.08] bg-white/[0.03] backdrop-blur-sm space-y-3" data-panel="1">
      <div class="flex flex-wrap items-center justify-between gap-2">
        <h4 class="text-sm font-semibold text-white flex items-center gap-2"><i class="fas ${icon} text-accent-cyan"></i>${title}</h4>
        <span class="text-[10px] text-accent-violet bg-accent-violet/10 border border-accent-violet/20 px-2 py-0.5 rounded-full">Python backend only</span>
      </div>
      <p class="text-[11px] text-gray-400 leading-relaxed">${subtitle}</p>
      <div class="pn-notice hidden rounded-lg p-2.5 border border-amber-500/30 bg-amber-500/[0.06] text-[11px] text-amber-200"><i class="fas fa-plug mr-1"></i>This panel needs the Python backend. The in-browser engine has no equivalent, so nothing is computed here. Start the backend (<code>uvicorn main:app</code>) or use the Streamlit app in bridge mode.</div>
      ${body}
      <div class="pn-err hidden rounded-lg p-2.5 border border-red-500/40 bg-red-500/[0.07] text-xs text-red-200"></div>
      <div class="pn-out space-y-3"></div>
    </div>`;
  }

  // ---------- shared plumbing ----------
  const isoProvenance = () => {                       // what isochrone is behind the numbers: shown on every isochrone-dependent panel
    const m = window._lastCmdMeta, p = m && m.workflow && m.workflow.parsec;
    if (!p || p.state !== 'done') return warn(['No isochrone has been fitted yet: press "Apply Extinction & Fetch PARSEC" on the Isochrone page first. Without it nothing here is meaningful.']);
    const msgs = [];
    if (String(p.grid_source || '').toUpperCase() === 'MOCK') msgs.push('The isochrone grid in use is a MOCK test grid, NOT PARSEC. Every mass, A_V and DM below is invalid.');
    return (msgs.length ? warn(msgs, 'red') : '') + `<div class="text-[11px] text-gray-500">Isochrone: ${esc(p.service || 'CMD 3.8')}${p.interpolated ? ' (bilinear blend of grid nodes)' : ''}, logAge ${num(p.logAge, 2)}, [M/H] ${num(p.MH, 2)}, A<sub>V</sub> sent ${num(p.A_V_sent, 3)}.</div>`;
  };

  async function call(root, endpoint, opts) {
    const err = root.querySelector('.pn-err'), out = root.querySelector('.pn-out');
    err.classList.add('hidden');
    if (!(await probeBackend())) { refreshNotices(); showToast('This panel needs the Python backend.', 'warning'); return null; }
    try {
      return await apiCall(endpoint, { ...opts, timeoutMs: CFG.TIMEOUT_MS });
    } catch (e) {
      err.textContent = String(e.message || e); err.classList.remove('hidden'); out.innerHTML = '';
      showToast(String(e.message || e), 'error');
      return null;
    }
  }
  const post = (root, endpoint, body) => call(root, endpoint, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
  const busy = (btn, on) => { btn.disabled = on; btn.classList.toggle('opacity-60', on); };
  const optNum = (id) => { const v = parseFloat($(id).value); return Number.isFinite(v) ? v : undefined; };
  const needFile = (task) => { const f = getTaskFilename(task); if (!f) showToast('Run the GMM first (or upload a custom filtered CSV).', 'warning'); return f; };
  const needStage2 = () => { if (!(window._lastCmdMeta && window._lastCmdMeta.workflow && window._lastCmdMeta.workflow.parsec && window._lastCmdMeta.workflow.parsec.state === 'done')) { showToast('Press "Apply Extinction & Fetch PARSEC" on the Isochrone page first.', 'warning'); return false; } return true; };
  const plot = (div, traces, layout) => Plotly.newPlot(div, traces, getPlotlyLayout({ height: CFG.PLOT_HEIGHT, ...layout }), PLOTLY_CONFIG);

  async function refreshNotices() {
    const alive = await probeBackend();
    document.querySelectorAll('[data-panel="1"] .pn-notice').forEach((n) => n.classList.toggle('hidden', !!alive));
    document.querySelectorAll('[data-panel="1"] button[data-act]').forEach((b) => { b.disabled = !alive; b.classList.toggle('opacity-50', !alive); b.classList.toggle('cursor-not-allowed', !alive); });
  }

  // ---------- 1. GMM reliability (Membership page) ----------
  function mountReliability(host) {
    host.innerHTML = card('pn-rel', 'fa-shield-halved', 'Membership reliability',
      'Re-draws the stars within their Gaia errors <em>n</em> times and re-computes PMemb (scatter of the probability, how often a member stays a member), and cross-checks the member list with an independent HDBSCAN clustering. Needs a finished GMM run in this session.',
      `<div class="grid grid-cols-2 md:grid-cols-4 gap-3 text-xs text-gray-400 items-end">
        <label class="flex flex-col gap-1">Resamples<input id="rel-n" type="number" min="10" max="2000" step="10" value="${CFG.DEFAULT_N_RESAMPLE}" class="${INPUT}"></label>
        <label class="flex flex-col gap-1">HDBSCAN min cluster size<input id="rel-mcs" type="number" min="3" step="1" value="${CFG.DEFAULT_MIN_CLUSTER}" class="${INPUT}"></label>
        <label class="flex flex-col gap-1">Seed (blank = default)<input id="rel-seed" type="number" step="1" class="${INPUT}"></label>
        <button data-act="rel" class="${BTN}"><i class="fas fa-vial"></i> Check reliability</button></div>`);
    host.querySelector('[data-act="rel"]').addEventListener('click', async (ev) => {
      const btn = ev.currentTarget;                 // currentTarget is null after the first await

      const root = $('pn-rel'); busy(btn, true);
      const body = { n_resample: optNum('rel-n'), min_cluster_size: optNum('rel-mcs') }; const s = optNum('rel-seed'); if (s !== undefined) body.seed = s;
      const r = await post(root, '/gmm-diagnostics', body); busy(btn, false); if (!r) return;
      const cut = parseFloat(($('input-cutoff') || {}).value) || 0.5;
      const rows = r.table || [], mem = rows.filter((x) => x.PMemb >= cut);
      const border = rows.filter((x) => Math.abs(x.PMemb - cut) < CFG.BORDERLINE_HALF_WIDTH).length;
      const h = r.hdbscan || {}, m = r.members || {};
      const worst = mem.slice().sort((a, c) => a.PMemb_frac_member - c.PMemb_frac_member).slice(0, CFG.TABLE_ROWS);
      root.querySelector('.pn-out').innerHTML =
        tiles([tile('Members', m.n ?? '—', `PMemb &ge; ${cut}`), tile('Mean &sigma;(PMemb)', num(m.mean_pmemb_std, 3), `max ${num(m.max_pmemb_std, 3)}`),
               tile('Unstable members', m.n_unstable ?? '—', esc(m.n_unstable_def || '')), tile(`Borderline (|PMemb&minus;${cut}| &lt; ${CFG.BORDERLINE_HALF_WIDTH})`, border, `of ${rows.length} fitted stars`)]) +
        tiles([tile('HDBSCAN agreement (Jaccard)', num(h.jaccard, 2), 'GMM &cap; HDBSCAN / GMM &cup; HDBSCAN'), tile('GMM members', h.n_gmm ?? '—'), tile('HDBSCAN members', h.n_hdbscan ?? '—', `${h.n_groups ?? '—'} groups`), tile('In both', h.n_both ?? '—')]) +
        warn([...(m.n_unstable > 0 ? [`${m.n_unstable} member(s) are members in < ${CFG.UNSTABLE_FRAC * 100}% of the error draws: treat them as uncertain.`] : []), ...(h.jaccard < 0.6 ? ['Low agreement with HDBSCAN: the Gaussian-mixture membership may be shaped by the model assumption, not only the data.'] : [])]) +
        `<div id="rel-plot"></div><div class="text-[11px] text-gray-400">Least stable members (row = position in the fitted sample, not the CSV row):</div>` +
        table([{ h: 'row', f: (x) => x._row }, { h: 'PMemb', f: (x) => num(x.PMemb, 3) }, { h: 'mean', f: (x) => num(x.PMemb_mean, 3) }, { h: '&sigma;', f: (x) => num(x.PMemb_std, 3) },
               { h: 'frac. member', f: (x) => num(x.PMemb_frac_member, 2) }, { h: 'HDBSCAN', f: (x) => (x.HDBSCAN_in_cluster ? 'yes' : 'no') }], worst.map((x, i) => ({ ...x, _row: rows.indexOf(x) })));
      plot('rel-plot', [{ x: rows.map((x) => x.PMemb), y: rows.map((x) => x.PMemb_std), mode: 'markers', type: 'scatter', marker: { size: 5, color: rows.map((x) => (x.HDBSCAN_in_cluster ? '#00d4ff' : '#6b7280')), opacity: 0.8 }, text: rows.map((x) => `frac member ${num(x.PMemb_frac_member, 2)}`), hoverinfo: 'x+y+text' }],
        { title: { text: 'PMemb scatter vs PMemb (cyan = also in HDBSCAN cluster)', font: { size: 12, color: '#9ca3af' } }, xaxis: { ...getPlotlyLayout().xaxis, title: 'PMemb' }, yaxis: { ...getPlotlyLayout().yaxis, title: '&sigma;(PMemb) over draws' }, shapes: [{ type: 'line', x0: cut, x1: cut, y0: 0, y1: 1, yref: 'paper', line: { color: '#f59e0b', dash: 'dot' } }] });
    });
  }

  // ---------- 2. Cluster-component rule (Membership page; read by app.js when it posts /run-gmm) ----------
  function mountRule(host) {
    host.innerHTML = card('pn-rule', 'fa-bullseye', 'Cluster-component rule',
      'Which Gaussian of the mixture is "the cluster"? <b>Tightest</b> = smallest covariance determinant (default; can pick a dense field blob). <b>Closest to literature</b> = component whose mean is nearest your literature proper motion and parallax. These fields are sent with the next <i>Run GMM</i>.',
      `<div class="grid grid-cols-2 md:grid-cols-5 gap-3 text-xs text-gray-400 items-end">
        <label class="flex flex-col gap-1 col-span-2 md:col-span-1">Rule<select id="rule-sel" class="bg-space-800 border border-white/[0.08] rounded-lg px-2 py-1.5 text-sm text-white"><option value="smallest_det" selected>Tightest (smallest determinant)</option><option value="closest_to_literature">Closest to literature</option></select></label>
        <label class="flex flex-col gap-1">Lit. pmRA (mas/yr)<input id="rule-pmra" type="number" step="0.001" class="${INPUT}" disabled></label>
        <label class="flex flex-col gap-1">Lit. pmDE (mas/yr)<input id="rule-pmde" type="number" step="0.001" class="${INPUT}" disabled></label>
        <label class="flex flex-col gap-1">Lit. parallax (mas)<input id="rule-plx" type="number" step="0.001" class="${INPUT}" disabled></label>
        <div class="flex flex-col gap-1"><input id="rule-name" type="text" placeholder="Cluster name, e.g. UPK 545" class="${INPUT}"><button data-act="hr" type="button" class="${BTN2}"><i class="fas fa-book-open mr-1"></i>Prefill (Hunt &amp; Reffert 2023)</button></div></div>
        <div class="text-[11px] text-gray-500">The prefill only fills the boxes above for you to inspect; nothing is applied until you press <i>Run GMM</i>.</div>`);
    const sel = $('rule-sel');
    const sync = () => ['rule-pmra', 'rule-pmde', 'rule-plx'].forEach((i) => { $(i).disabled = sel.value !== 'closest_to_literature'; });
    sel.addEventListener('change', sync);
    host.querySelector('[data-act="hr"]').addEventListener('click', async (ev) => {
      const btn = ev.currentTarget;                 // currentTarget is null after the first await

      const root = $('pn-rule'), name = $('rule-name').value.trim(); if (!name) { showToast('Type a cluster name first.', 'warning'); return; }
      busy(btn, true);
      const r = await call(root, '/hunt-reffert?name=' + encodeURIComponent(name), { method: 'GET' }); busy(btn, false); if (!r) return;
      if (r.pmRA != null) $('rule-pmra').value = r.pmRA; if (r.pmDE != null) $('rule-pmde').value = r.pmDE; if (r.plx != null) $('rule-plx').value = r.plx;
      if (r.pmRA != null && r.plx != null) { sel.value = 'closest_to_literature'; sync(); }
      root.querySelector('.pn-out').innerHTML = `<div class="text-[11px] text-gray-300">Hunt &amp; Reffert (2023) for ${esc(name)}: pmRA ${num(r.pmRA)}, pmDE ${num(r.pmDE)}, parallax ${num(r.plx)}, distance ${num(r.distance_pc, 1)} pc.${r.missing && r.missing.length ? ` Missing: ${esc(r.missing.join(', '))}.` : ''} Distance is shown for information only (set the window centre yourself).</div>`;
    });
    // app.js merges this into the /run-gmm request body
    window.starnestClusterRuleParams = () => {
      const rule = sel.value, out = { cluster_rule: rule };
      if (rule === 'closest_to_literature') {
        const a = parseFloat($('rule-pmra').value), d = parseFloat($('rule-pmde').value), p = parseFloat($('rule-plx').value);
        if (Number.isFinite(a) && Number.isFinite(d)) out.lit_pm = [a, d];
        if (Number.isFinite(p)) out.lit_plx = p;
      }
      return out;
    };
  }

  // ---------- 3. Likelihood isochrone fit (Isochrone page) ----------
  function mountIsoFit(host) {
    host.innerHTML = card('pn-isofit', 'fa-chart-line', 'Likelihood fit (A<sub>V</sub>, DM, binaries, outliers)',
      'Fits the isochrone to the member CMD with a mixture likelihood: single-star sequence + equal-mass-binary sequence + uniform outliers, with Gaussian priors on A<sub>V</sub> (dust map) and DM (parallax). The nearest-point fit quality shown on this page is a separate, older diagnostic and is unchanged.',
      `<div class="grid grid-cols-2 md:grid-cols-4 gap-3 text-xs text-gray-400 items-end">
        <label class="flex flex-col gap-1">Dust source<select id="ifit-dust" class="bg-space-800 border border-white/[0.08] rounded-lg px-2 py-1.5 text-sm text-white"><option value="sfd" selected>SFD (IRSA value in use)</option><option value="bayestar">Bayestar19 (needs dustmaps)</option><option value="free">Free (no dust prior)</option></select></label>
        <label class="flex items-center gap-2 mt-4"><input id="ifit-bin" type="checkbox" checked> Fit binary sequence</label>
        <div></div><button data-act="ifit" class="${BTN}"><i class="fas fa-bullseye"></i> Fit</button></div>`);
    host.querySelector('[data-act="ifit"]').addEventListener('click', async (ev) => {
      const btn = ev.currentTarget;                 // currentTarget is null after the first await

      const root = $('pn-isofit'); const f = needFile('isochrone'); if (!f || !needStage2()) return;
      busy(btn, true);
      const r = await post(root, '/fit-isochrone-likelihood', { filename: f, dust_source: $('ifit-dust').value, fit_binaries: $('ifit-bin').checked }); busy(btn, false); if (!r) return;
      const msgs = [];
      if (!r.converged) msgs.push('The optimiser did not converge: do not quote these numbers.');
      if (r.f_out > 0.45 || r.f_out < 1e-3) msgs.push(`Outlier fraction ${num(r.f_out, 3)} sits at a boundary of its allowed range: it is poorly constrained by this sample.`);
      if (r.dust_source === 'free') msgs.push('Dust is free: A_V and DM are degenerate with age and metallicity; only the DM prior holds the fit.');
      root.querySelector('.pn-out').innerHTML = isoProvenance() + warn(msgs) +
        tiles([tile('A<sub>V</sub>', pm(r.A_V, r.A_V_err, 3), `E(B&minus;V) = ${num(r.E_BV, 3)}`), tile('Distance modulus', pm(r.DM, r.DM_err, 3), `d = ${num(r.d_pc, 1)} pc`),
               tile('Binary-sequence fraction', num(r.f_bin, 3), `${r.n_binary_seq ?? '—'} stars with p &gt; 0.5`), tile('Outlier fraction', num(r.f_out, 3), `${r.n_outlier ?? '—'} stars with p &gt; 0.5`)]) +
        `<div class="text-[11px] text-gray-400">N = ${r.n}, parameters ${r.n_par}, BIC ${num(r.bic, 1)}, dust source <b>${esc(r.dust_source)}</b>${r.av_prior != null ? `, A<sub>V</sub> prior ${num(r.av_prior, 3)}` : ''}, DM prior ${num(r.dm_prior, 3)}. Errors on A<sub>V</sub> and DM are conditional on the best-fit mixture fractions. ${esc(r.legacy_note || '')}</div>`;
    });
  }

  // ---------- 4. Velocities (Kinematics page) ----------
  function mountVelocities(host) {
    host.innerHTML = card('pn-vel', 'fa-gauge-high', 'Velocities',
      'Tangential velocities v = 4.74 &mu;/&varpi;, projected into v<sub>out</sub> (radial expansion) and v<sub>az</sub> (rotation) about the cluster centre, with the perspective-expansion correction, and a maximum-likelihood 1-D velocity dispersion that deconvolves each star&rsquo;s own measurement covariance.',
      `<div class="grid grid-cols-2 md:grid-cols-5 gap-3 text-xs text-gray-400 items-end">
        <label class="flex flex-col gap-1">Parallax used<select id="vel-plx" class="bg-space-800 border border-white/[0.08] rounded-lg px-2 py-1.5 text-sm text-white"><option value="star" selected>Each star&rsquo;s own</option><option value="cluster">Cluster median</option></select></label>
        <label class="flex flex-col gap-1">Virial &eta;<input id="vel-eta" type="number" step="0.5" min="1" value="${CFG.DEFAULT_ETA}" class="${INPUT}"></label>
        <label class="flex flex-col gap-1">Bootstrap draws<input id="vel-boot" type="number" step="50" min="50" value="${CFG.DEFAULT_N_BOOT}" class="${INPUT}"></label>
        <label class="flex items-center gap-2 mt-4"><input id="vel-persp" type="checkbox" checked> Perspective correction</label>
        <button data-act="vel" class="${BTN}"><i class="fas fa-play"></i> Compute</button></div>`);
    host.querySelector('[data-act="vel"]').addEventListener('click', async (ev) => {
      const btn = ev.currentTarget;                 // currentTarget is null after the first await

      const root = $('pn-vel'); const f = needFile('kinematics'); if (!f) return; busy(btn, true);
      const r = await post(root, '/compute-velocities', { filename: f, plx_mode: $('vel-plx').value, eta: optNum('vel-eta'), n_boot: optNum('vel-boot'), perspective: $('vel-persp').checked }); busy(btn, false); if (!r) return;
      const notes = (r.notes || []).slice();
      if (r.mass_msun && r.sigma_vir_kms) notes.push(`Virial check: sigma_1D / sigma_vir = ${num(r.sigma_1d_kms / r.sigma_vir_kms, 2)} (mass ${num(r.mass_msun, 0)} Msun, eta ${r.eta}). A ratio well above 1 suggests expansion / unbound or binary-inflated dispersion.`);
      root.querySelector('.pn-out').innerHTML =
        tiles([tile('&sigma;<sub>1D</sub> (km/s)', pm(r.sigma_1d_kms, r.sigma_1d_err, 3), r.mixture_used ? 'outlier-mixture model' : 'single Gaussian'), tile('v<sub>out</sub> (km/s)', pm(r.v_out_kms, r.v_out_err, 3), 'expansion &gt; 0'),
               tile('v<sub>az</sub> (km/s)', pm(r.v_az_kms, r.v_az_err, 3), 'rotation'), tile('Cluster RV (km/s)', num(r.rv_cluster_kms, 2), `${r.n_rv ?? 0} Gaia RVs`)]) +
        `<div class="text-[11px] text-gray-400">N = ${r.n}; parallax mode <b>${esc(r.plx_mode)}</b>; cluster parallax ${num(r.plx_cluster_mas, 3)} mas. Mixture vs single Gaussian: &Delta;BIC = ${num(r.delta_bic, 1)} (${r.mixture_used ? 'mixture preferred, &Delta;BIC &gt; 10' : 'single Gaussian kept'}). Median measurement error ${num(r.median_meas_err_kms, 3)} km/s.${r.sigma_vir_kms ? ` &sigma;<sub>vir</sub> = ${num(r.sigma_vir_kms, 3)} km/s.` : ''}</div>` + warn(notes);
    });
  }

  // ---------- 5. Binaries & Mass (IMF page) ----------
  function mountBinaries(host) {
    host.innerHTML = card('pn-bin', 'fa-star-half-stroke', 'Binaries &amp; total mass',
      'Monte-Carlo masses that allow for unresolved binaries: each star is matched to synthetic single / binary systems built from the isochrone (flux addition, companion mass uniform in 0&ndash;m<sub>1</sub>) and gets a binary probability and mass ratio q. Total mass M<sub>T</sub> = observed + unseen low-mass stars (with companions) + white dwarfs.',
      `<div class="grid grid-cols-2 md:grid-cols-4 gap-3 text-xs text-gray-400 items-end">
        <label class="flex flex-col gap-1">Assumed binary fraction<input id="bin-f" type="number" step="0.05" min="0" max="1" value="${CFG.DEFAULT_F_BIN}" class="${INPUT}"></label>
        <label class="flex flex-col gap-1">Monte-Carlo repeats<input id="bin-mc" type="number" step="10" min="1" value="${CFG.DEFAULT_N_MC}" class="${INPUT}"></label>
        <label class="flex flex-col gap-1">Synthetic systems<input id="bin-syn" type="number" step="1000" min="1000" value="${CFG.DEFAULT_N_SYNTH}" class="${INPUT}"></label>
        <button data-act="bin" class="${BTN}"><i class="fas fa-play"></i> Compute</button></div>`);
    host.querySelector('[data-act="bin"]').addEventListener('click', async (ev) => {
      const btn = ev.currentTarget;                 // currentTarget is null after the first await

      const root = $('pn-bin'); const f = needFile('imf'); if (!f || !needStage2()) return; busy(btn, true);
      const r = await post(root, '/compute-binaries', { filename: f, f_bin: optNum('bin-f'), n_mc: optNum('bin-mc'), n_synth: optNum('bin-syn') }); busy(btn, false); if (!r) return;
      const m = r.mass, msgs = [...(r.warnings || []), ...(r.note ? [r.note] : []), r.caveat].filter(Boolean);
      const rows = (r.table || []).slice(0, CFG.TABLE_ROWS);
      root.querySelector('.pn-out').innerHTML = isoProvenance() + warn(msgs) +
        tiles([tile('Stars', r.n), tile('Binary candidates', r.n_binary, `observed fraction ${num(r.observed_binary_fraction, 2)} (prior-driven)`), tile('Unmatched', r.n_unmatched, 'no synthetic neighbour'),
               tile('M<sub>T</sub> (M&#9737;)', m ? num(m.M_T, 1) : '—', m ? `slope &Gamma; = ${num(m.gamma, 2)}` : 'not computed')]) +
        (m ? `<div class="text-[11px] text-gray-300">M<sub>T</sub> = M<sub>obs</sub> ${num(m.M_obs, 1)} + M<sub>unseen</sub> ${num(m.M_unseen, 1)} &times; (1 + f<sub>bin</sub> q&#772;) [f<sub>bin</sub> = ${num(m.f_bin, 2)}, q&#772; = ${num(m.q_mean, 2)}] + M<sub>WD</sub> ${num(m.M_WD, 1)} M&#9737;. Slope fitted down to m<sub>min</sub> = ${num(m.m_min, 2)} M&#9737;; unseen / WD terms are extrapolations of that slope.</div>` : '') +
        `<div class="text-[11px] text-gray-400">First ${rows.length} rows of ${r.n} (full table: button below).</div>` +
        table([{ h: 'G', f: (x) => num(x.Gmag, 2) }, { h: 'MC_Mass', f: (x) => num(x.MC_Mass, 3) }, { h: 'MC_CompMass', f: (x) => num(x.MC_CompMass, 3) }, { h: 'is_binary', f: (x) => (x.is_binary ? 'yes' : 'no') },
               { h: 'q', f: (x) => num(x.q, 2) }, { h: 'binary_prob', f: (x) => num(x.binary_prob, 2) }, { h: 'match_dist', f: (x) => num(x.match_dist, 3) }], rows) +
        `<button type="button" class="${BTN2}" id="bin-dl"><i class="fas fa-download mr-1"></i>Download full table (CSV)</button>`;
      const dl = $('bin-dl'); if (dl) dl.addEventListener('click', () => apiOpen('/download-binaries'));
    });
  }

  // ---------- 6. Dynamics (Kinematics page) ----------
  function mountDynamics(host) {
    host.innerHTML = card('pn-dyn', 'fa-sitemap', 'Dynamics',
      'Mass segregation ratio &Lambda;<sub>MSR</sub> (Allison et al. 2009: &Lambda; &gt; 1 = massive stars more concentrated than random stars), tidal radius, relaxation time versus age, and &mdash; if <i>galpy</i> is installed &mdash; an orbit integration checking whether the members are stretched along the orbit (a tidal tail).',
      `<div class="grid grid-cols-2 md:grid-cols-5 gap-3 text-xs text-gray-400 items-end">
        <label class="flex flex-col gap-1">Random draws<input id="dyn-nr" type="number" step="50" min="20" value="${CFG.DEFAULT_N_RANDOM}" class="${INPUT}"></label>
        <label class="flex flex-col gap-1">Cluster mass, M&#9737; (blank = auto)<input id="dyn-mass" type="number" step="1" min="1" class="${INPUT}"></label>
        <label class="flex flex-col gap-1">sigma 1-D, km/s (blank = auto)<input id="dyn-sig" type="number" step="0.01" min="0.01" class="${INPUT}"></label>
        <label class="flex items-center gap-2 mt-4"><input id="dyn-orbit" type="checkbox" checked> Orbit / tail check</label>
        <button data-act="dyn" class="${BTN}"><i class="fas fa-play"></i> Compute</button></div>`);
    host.querySelector('[data-act="dyn"]').addEventListener('click', async (ev) => {
      const btn = ev.currentTarget;                 // currentTarget is null after the first await

      const root = $('pn-dyn'); const f = needFile('kinematics'); if (!f || !needStage2()) return; busy(btn, true);
      const body = { filename: f, n_random: optNum('dyn-nr'), orbit: $('dyn-orbit').checked }; const m0 = optNum('dyn-mass'), s0 = optNum('dyn-sig'); if (m0) body.mass = m0; if (s0) body.sigma_1d = s0;
      const r = await post(root, '/compute-dynamics', body); busy(btn, false); if (!r) return;
      const rel = r.relaxation, msgs = (r.notes || []).slice();
      if (/lower limit/.test(r.mass_source || '')) msgs.push('Mass = sum of the observed stars only (a lower limit): the tidal radius and relaxation time are underestimated. Run Binaries & total mass first, or type a mass.');
      const msr = (r.msr || []).filter((x) => Number.isFinite(x.Lambda));
      root.querySelector('.pn-out').innerHTML = isoProvenance() + warn(msgs) +
        tiles([tile('Tidal radius', `${num(r.r_t_pc, 1)} pc`, `mass ${num(r.mass_msun, 0)} M&#9737; (${esc(r.mass_source)})`), tile('Half-number radius', `${num(r.r_half_number_pc, 2)} pc`, `d = ${num(r.d_pc, 1)} pc`),
               tile('Relaxation time', rel ? `${num(rel.t_relax_myr, 0)} Myr` : '—', rel ? `crossing ${num(rel.t_cross_myr, 1)} Myr, N = ${rel.n_stars}` : 'needs &sigma;<sub>1D</sub>'),
               tile('Age / t<sub>relax</sub>', rel && rel.age_over_trelax != null ? num(rel.age_over_trelax, 2) : '—', rel && rel.age_over_trelax > 1 ? 'dynamically evolved' : 'not yet relaxed')]) +
        `<div id="dyn-msr"></div>${r.orbit ? `<div class="text-[11px] text-gray-300">Orbit check: elongation ${num(r.orbit.elongation, 2)}, angle to orbit ${num(r.orbit.angle_deg, 1)}&deg; &rarr; <b>${r.orbit.tail_like ? 'tail-like' : 'no tail signature'}</b>.</div><div id="dyn-orbit-plot"></div>` : ''}`;
      if (msr.length) plot('dyn-msr', [{ x: msr.map((x) => x.N), y: msr.map((x) => x.Lambda), mode: 'lines+markers', type: 'scatter', line: { color: '#00d4ff' }, error_y: { type: 'data', symmetric: false, array: msr.map((x) => x.Lambda_hi - x.Lambda), arrayminus: msr.map((x) => x.Lambda - x.Lambda_lo), color: '#00d4ff' } }],
        { title: { text: 'Mass segregation ratio vs N of heaviest stars', font: { size: 12, color: '#9ca3af' } }, xaxis: { ...getPlotlyLayout().xaxis, title: 'N heaviest stars', type: 'log' }, yaxis: { ...getPlotlyLayout().yaxis, title: '&Lambda;<sub>MSR</sub>' }, shapes: [{ type: 'line', x0: 0, x1: 1, xref: 'paper', y0: 1, y1: 1, line: { color: '#f59e0b', dash: 'dot' } }] });
      if (r.orbit && r.orbit.orbit_path_pc) plot('dyn-orbit-plot', [{ x: r.orbit.orbit_path_pc.map((p) => p[0]), y: r.orbit.orbit_path_pc.map((p) => p[1]), mode: 'lines', type: 'scatter', name: 'orbit', line: { color: '#a78bfa' } }, { x: [0], y: [0], mode: 'markers', name: 'Sun', marker: { color: '#f59e0b', size: 9 } }],
        { title: { text: `Orbit &plusmn;50 Myr, heliocentric X-Y (pc)`, font: { size: 12, color: '#9ca3af' } }, xaxis: { ...getPlotlyLayout().xaxis, title: 'X (pc)' }, yaxis: { ...getPlotlyLayout().yaxis, title: 'Y (pc)', scaleanchor: 'x' } });
    });
  }

  // ---------- mount ----------
  const MOUNTS = { 'panel-gmm-reliability': mountReliability, 'panel-cluster-rule': mountRule, 'panel-iso-likelihood': mountIsoFit, 'panel-velocities': mountVelocities, 'panel-binaries': mountBinaries, 'panel-dynamics': mountDynamics };
  Object.entries(MOUNTS).forEach(([id, fn]) => { const h = $(id); if (h) fn(h); });
  refreshNotices();
  window.starnestPanelsRefresh = refreshNotices;       // app.js may call it after the engine badge changes
})();
