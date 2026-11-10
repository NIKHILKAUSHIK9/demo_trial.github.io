/* StarNest UI logic: in-browser GMM / IMF engine, IRSA+CMD client, plotting, page navigation. Extracted from index.html. */
    'use strict';

    window.globalFileName = null;
    window.gmmFilteredFileName = null;
    window.rawFileName = null;   // the uploaded RAW catalogue; GMM always re-runs on this, never on members_latest.csv
    const PAGE_ORDER = ['home', 'membership', 'isochrone', 'kinematics', 'imf'];
    let currentPage = 'home';
    let isTransitioning = false;
    // Backend address: set window.STARNEST_API_BASE to override. Served over http(s) by the backend -> 'api' next to the page
    // ('/api' at the site root, '/starnest/api' when the backend is mounted under /starnest); opened from a file -> the local backend.
    const API_BASE = window.STARNEST_API_BASE || ((location.protocol === 'http:' || location.protocol === 'https:')
      ? location.origin + location.pathname.replace(/[^/]*$/, '') + 'api' : 'http://localhost:8000/api');
    // Visitor id for a backend on another origin (page hosted elsewhere, e.g. Streamlit): sent as X-Session-Id / ?sid=
    const SESSION_ID = (() => {
      let v = ''; try { v = sessionStorage.getItem('sn_sid') || ''; } catch (e) {}
      if (!v) { v = (crypto.randomUUID ? crypto.randomUUID() : String(Math.random()).slice(2) + Date.now()).replace(/-/g, ''); try { sessionStorage.setItem('sn_sid', v); } catch (e) {} }
      return v;
    })();
    const apiOpen = async (path) => {
      if (window.STARNEST_BRIDGE) {       // no URL a new tab could open: fetch through the bridge and save the Blob
        try {
          const r = await fetch(API_BASE + path); if (!r.ok) throw new Error('download failed (' + r.status + ')');
          const m = /filename="?([^";]+)"?/i.exec(r.headers.get('content-disposition') || '');
          const a = document.createElement('a'); a.href = URL.createObjectURL(await r.blob()); a.download = m ? m[1] : (path.split('/').pop() || 'download');
          document.body.appendChild(a); a.click(); a.remove();
        } catch (e) { if (typeof showToast === 'function') showToast(String(e.message || e), 'error'); }
        return;
      }
      window.open(API_BASE + path + (path.includes('?') ? '&' : '?') + 'sid=' + SESSION_ID, '_blank');
    };
    { const _f = window.fetch.bind(window);
      window.fetch = (u, o) => { if (typeof u === 'string' && u.startsWith(API_BASE)) { o = { ...(o || {}), headers: { ...((o || {}).headers || {}), 'X-Session-Id': SESSION_ID } }; } return _f(u, o); }; }

    const taskSources = {
      isochrone: { mode: 'gmm', customFile: null, customLabel: null },
      kinematics: { mode: 'gmm', customFile: null, customLabel: null },
      imf: { mode: 'gmm', customFile: null, customLabel: null }
    };

    function getTaskFilename(task) {
      const cfg = taskSources[task];
      if (cfg && cfg.mode === 'custom' && cfg.customFile) {
        return cfg.customFile;
      }
      return window.gmmFilteredFileName || null;   // never fall back to the RAW file (field stars)
    }

    function getTaskDisplayLabel(task) {
      const cfg = taskSources[task];
      if (cfg && cfg.mode === 'custom' && cfg.customLabel) {
        return `Custom: ${cfg.customLabel}`;
      }
      return 'GMM Filtered File';
    }

    const STAR_TYPE_COLORS = {
      "Main Sequence": "#00FFFF",
      "Pre-Main Sequence": "#BA8E23",
      "White Dwarf": "#000000",
      "Red Giant / Supergiant": "#FF0000"
    };

    const CLUSTER_CENTERS = {
      "melotte111": [186.64917, 26.107],
      "ngc2910": [142.625, -52.9167],
      "haffner14": [116.182, -28.389],
      "rp147": [289.175, -16.300],
      "ubc52": [64.741, 52.373],
      "ubc7": [108.300, -37.350],
      "upk545": [130.838, -58.725]
    };

    // In-memory store for uploaded or sample cluster rows
    const localStore = {
      files: {},
      latestMembers: [],
      latestExtinction: [],
      latestIsochrone: [],
      latestImfBins: [],
      latestImfChi2: [],
      customIsochrone: null,
      activeClusterName: '',
      cache: {},
      ver: {},            // content version per stored file name (cache invalidation)
      verCounter: 0,
      isoVer: 0,          // bumps when a custom isochrone is (re)loaded
      parsec: null,       // {text, name}: last PARSEC file fetched by the browser (for download)
      analysis: null,     // {srcFile, srcVer, rows, id}: extinction-corrected table shared by Kinematics + IMF
      engine: 'unknown'   // 'python' | 'browser'
    };

    // Draw the analytic stand-in curve when no real isochrone is available? It is NOT PARSEC. Default: off.
    const ALLOW_PLACEHOLDER_ISOCHRONE = false;
    const API_TIMEOUT_MS = 120000;       // real backend: PARSEC / IRSA calls can be slow
    const API_PROBE_TIMEOUT_MS = window.STARNEST_PROBE_TIMEOUT_MS || 1200;   // decides whether the backend exists (hosts that wake up slowly raise it)

    // StarNest module constants (exact match to gmm.py, extinction_isochrone.py, imf.py)
    const STARNEST_REQUIRED_COLS = ['Plx', 'e_Plx', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE', 'Gmag', 'BP-RP', 'RUWE', 'gofAL', 'Nper'];
    const STARNEST_AV_PER_EBV = 2.742;
    const STARNEST_AG_AV = 0.789, STARNEST_ABP_AV = 1.002, STARNEST_ARP_AV = 0.589;
    const STARNEST_BPRP0_RANGE = [-0.8, 5.0];
    const STARNEST_LOGAGE_FLOOR = 6.6;
    const STARNEST_MIN_MASSES = 20;
    const STARNEST_M_MIN_FLOOR = 0.08, STARNEST_M_MAX_CAP = 50.0;

    function showModuleError(boxId, textId, msg) {
      const box = document.getElementById(boxId);
      const txt = document.getElementById(textId);
      if (!box || !txt) return;
      if (msg) {
        txt.textContent = msg;
        box.classList.remove('hidden');
      } else {
        txt.textContent = '';
        box.classList.add('hidden');
      }
    }

    function validateColumnsJS(rows) {
      if (!rows || !rows.length) return STARNEST_REQUIRED_COLS;
      const cols = Object.keys(rows[0]);
      return STARNEST_REQUIRED_COLS.filter(req => !cols.includes(req));
    }

    // ===== TEAM PHOTO LOADER =====
    // Put each photo in  images/team/  (or images/ or next to index.html). The file name is the member's name, e.g.
    //   neelam_panwar.jpg  richa_kundu.jpg  nikhil_kaushik.jpg  pryanshu_basur.jpg  monu.jpg  shruti_prasad.jpg
    // Accepted: .jpg .jpeg .png .webp ; the alternative spellings listed in each <img data-photo="..."> are tried too.
    // Use lower-case names: GitHub Pages / Linux servers are case-sensitive. Not found -> initials badge stays visible.
    const TEAM_PHOTO_CONFIG = {
      dirs: ['images/team/', 'images/', ''],        // searched in this order, relative to index.html
      exts: ['jpg', 'jpeg', 'png', 'webp']
    };
    function loadTeamPhotos() {
      document.querySelectorAll('img[data-photo]').forEach(img => {
        const embedded = window.STARNEST_PHOTOS && img.dataset.photo.split(',').map(st => window.STARNEST_PHOTOS[st]).find(Boolean);   // photos inlined by a host page (Streamlit)
        if (embedded) { img.src = embedded; img.style.display = ''; if (img.nextElementSibling) img.nextElementSibling.style.display = 'none'; return; }
        const cands = [];
        TEAM_PHOTO_CONFIG.dirs.forEach(d => img.dataset.photo.split(',').forEach(stem => TEAM_PHOTO_CONFIG.exts.forEach(x => cands.push(`${d}${stem}.${x}`))));
        const tryNext = i => {
          if (i >= cands.length) return;                                   // nothing found: keep the initials badge
          const probe = new Image();
          probe.onload = () => {
            img.src = cands[i]; img.style.display = '';
            if (img.nextElementSibling) img.nextElementSibling.style.display = 'none';
          };
          probe.onerror = () => tryNext(i + 1);
          probe.src = cands[i];
        };
        tryNext(0);
      });
    }

    function escapeHTML(v) {
      return String(v).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
    }

    function medianOf(arr) {
      const a = arr.filter(Number.isFinite).sort((x, y) => x - y);
      if (!a.length) return NaN;
      const m = Math.floor(a.length / 2);
      return a.length % 2 ? a[m] : 0.5 * (a[m - 1] + a[m]);
    }

    function minMax(arr) {          // safe for very large arrays (no spread -> no stack overflow)
      let lo = Infinity, hi = -Infinity;
      for (let i = 0; i < arr.length; i++) { const v = arr[i]; if (v < lo) lo = v; if (v > hi) hi = v; }
      return [lo, hi];
    }

    function uiEbv() {              // E(B-V) currently typed in the Isochrone page (0 if blank)
      const el = document.getElementById('input-ebv');
      const v = el ? parseFloat(el.value) : NaN;
      return Number.isFinite(v) ? v : 0;
    }

    function pruneCache(max = 150) {            // slider drags create many cache keys: keep the newest `max`
      const keys = Object.keys(localStore.cache);
      if (keys.length > max) keys.slice(0, keys.length - max).forEach(k => delete localStore.cache[k]);
    }

    function commitMembers(members) {
      localStore.latestMembers = members;
      localStore.files['members_latest.csv'] = members;
      localStore.ver['members_latest.csv'] = ++localStore.verCounter;
      localStore.analysis = null;
      return localStore.ver['members_latest.csv'];
    }

    function restoreMembers(entry) {
      localStore.latestMembers = entry.members;
      localStore.files['members_latest.csv'] = entry.members;
      localStore.ver['members_latest.csv'] = entry.verId;
      localStore.analysis = null;
    }

    function mlMass(lum) {
      // piecewise main-sequence mass-luminosity relation (solar units):
      //   L = 0.23 M^2.3 (M < 0.43) | L = M^4 (0.43-2) | L = 1.4 M^3.5 (M > 2)
      if (!(lum > 0)) return NaN;
      if (lum < 0.033) return Math.pow(lum / 0.23, 1 / 2.3);
      if (lum < 16) return Math.pow(lum, 1 / 4.0);
      return Math.pow(lum / 1.4, 1 / 3.5);
    }

    function approxMass(mg, color) {
      if (!Number.isFinite(mg) || !Number.isFinite(color)) return NaN;
      const mBol = mg + (-0.012 - 0.25 * color - 0.15 * color * color);
      return mlMass(Math.pow(10, 0.4 * (4.74 - mBol)));
    }

    // star type + mass on one row. User-supplied values are kept; values this app derived are always refreshed.
    function assignDerived(copy, mg, bprpCorr) {
      const existingType = copy.Star_Type || copy['type of stars'] || copy.star_type;
      if (existingType && STAR_TYPE_COLORS[existingType] && copy.Star_Type_src !== 'auto') {
        copy.Star_Type = existingType; copy.Star_Type_src = 'user';
      } else {
        copy.Star_Type = classifyStarRow(mg, bprpCorr); copy.Star_Type_src = 'auto';
      }
      copy['type of stars'] = copy.Star_Type;
      const given = getVal(copy, ['MC_Mass', 'Mass', 'mass']);
      if (Number.isFinite(given) && copy.MC_Mass_src !== 'approx_ML') {
        if (!Number.isFinite(Number(copy.MC_Mass))) copy.MC_Mass = given;
        copy.MC_Mass_src = copy.MC_Mass_src || 'user';
      } else {
        const m = approxMass(mg, bprpCorr);
        if (Number.isFinite(m)) { copy.MC_Mass = Number(m.toFixed(4)); copy.MC_Mass_src = 'approx_ML'; }
      }
    }

    // Rows for Kinematics / IMF: the extinction-corrected table of the CURRENT members, else a fresh enrichment.
    function getAnalysisRows(filename) {
      const a = localStore.analysis;
      if (a && a.srcFile === filename && a.srcVer === (localStore.ver[filename] || 0)) return a.rows;
      const rows = localStore.files[filename] || localStore.latestMembers || [];
      const plx = rows.map(r => getVal(r, ['Plx', 'parallax', 'plx'])).filter(v => Number.isFinite(v) && v > 0);
      const med = medianOf(plx);
      const dm = Number.isFinite(med) ? 5 * Math.log10(1000 / med) - 5 : 7.66;
      return enrichPhotometryAndClassify(rows, dm, STARNEST_AV_PER_EBV * uiEbv());
    }

    // ---- CSV parsing: RFC-4180 quotes, auto delimiter (comma / tab / semicolon / whitespace) ----
    function splitCSVLine(line, delim) {
      if (delim === 'ws') return line.trim().split(/\s+/);
      const out = []; let cur = '', q = false;
      for (let i = 0; i < line.length; i++) {
        const ch = line[i];
        if (q) {
          if (ch === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; }
          else cur += ch;
        } else if (ch === '"') q = true;
        else if (ch === delim) { out.push(cur); cur = ''; }
        else cur += ch;
      }
      out.push(cur);
      return out;
    }

    const NUM_RE = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;

    function parseCSVText(text) {
      const cleaned = text.replace(/^﻿/, '');
      const all = cleaned.split(/\r?\n/);
      const lines = all.filter(l => l.trim() && !l.trim().startsWith('#'));
      if (lines.length < 2 && !all.some(l => l.trim().startsWith('#'))) return [];
      let headerLine = lines[0];
      let dataLines = lines.slice(1);
      const first = lines[0] || '';
      const cnt = c => (first.match(new RegExp(c === '\t' ? '\\t' : '\\' + c, 'g')) || []).length;
      let delim = ',';
      if (cnt(',') === 0) delim = cnt('\t') > 0 ? '\t' : (cnt(';') > 0 ? ';' : 'ws');
      // PARSEC-style whitespace tables keep the header in a leading "# col col ..." comment line
      const firstTokens = first.trim().split(/\s+/);
      if (delim === 'ws' && firstTokens.every(t => NUM_RE.test(t))) {
        const commentHdr = [...all].reverse().find(l => l.trim().startsWith('#') && /[A-Za-z]/.test(l));
        if (commentHdr) { headerLine = commentHdr.replace(/^\s*#\s*/, ''); dataLines = lines; }
      }
      const headers = splitCSVLine(headerLine, delim).map(h => h.trim());
      const rows = [];
      for (const line of dataLines) {
        const vals = splitCSVLine(line, delim);
        if (vals.length < headers.length * 0.5) continue;
        const obj = {};
        headers.forEach((h, idx) => {
          const raw = vals[idx] !== undefined ? vals[idx].trim() : '';
          obj[h] = NUM_RE.test(raw) ? Number(raw) : raw;
        });
        rows.push(obj);
      }
      return rows;
    }

    // Isochrone upload (PARSEC download or processed CSV): colour + absolute G, optional age/label selection
    function parseIsochroneRows(rows, logAge) {
      if (!rows.length) return null;
      let sel = rows;
      const ages = [...new Set(rows.map(r => getVal(r, ['logAge', 'logage', 'log_age'])).filter(Number.isFinite))];
      let usedAge = null;
      if (ages.length > 1 && Number.isFinite(logAge)) {
        usedAge = ages.reduce((b, a) => Math.abs(a - logAge) < Math.abs(b - logAge) ? a : b);
        sel = rows.filter(r => Math.abs(getVal(r, ['logAge', 'logage', 'log_age']) - usedAge) < 1e-6);
      } else if (ages.length === 1) usedAge = ages[0];
      if (sel.some(r => Number.isFinite(getVal(r, ['label'])))) {
        sel = sel.filter(r => { const l = getVal(r, ['label']); return !Number.isFinite(l) || l <= 3; });
      }
      const color = [], mag = [];
      sel.forEach(r => {
        let c = getVal(r, ['bp-rp_corr', 'BP-RP', 'color', 'bprp']);
        if (!Number.isFinite(c)) c = getVal(r, ['G_BPmag', 'G_BP']) - getVal(r, ['G_RPmag', 'G_RP']);
        const m = getVal(r, ['M_G_corr', 'M_G', 'Gmag', 'mag']);
        if (Number.isFinite(c) && Number.isFinite(m)) { color.push(c); mag.push(m); }
      });
      return color.length ? { color, mag, source: 'custom', logAge: usedAge } : null;
    }

    function getVal(row, candidates) {
      for (const key of candidates) {
        if (row[key] !== undefined && row[key] !== '' && !Number.isNaN(Number(row[key]))) {
          return Number(row[key]);
        }
      }
      const lowerMap = {};
      Object.keys(row).forEach(k => { lowerMap[k.toLowerCase()] = row[k]; });
      for (const key of candidates) {
        const val = lowerMap[key.toLowerCase()];
        if (val !== undefined && val !== '' && !Number.isNaN(Number(val))) {
          return Number(val);
        }
      }
      return NaN;
    }

    function classifyStarRow(mg, bprpCorr) {
      if (Number.isNaN(mg) || Number.isNaN(bprpCorr)) return "Main Sequence";
      if (mg > 10.0 && bprpCorr < 0.6) return "White Dwarf";
      if (mg < 3.0 && bprpCorr > 0.8) return "Red Giant / Supergiant";
      if (mg > 4.0 && mg < (3.0 * bprpCorr + 2.0)) return "Pre-Main Sequence";
      return "Main Sequence";
    }

    function enrichPhotometryAndClassify(rows, distMod = 7.66, av = STARNEST_AV_PER_EBV * uiEbv()) {
      return rows.map((r, idx) => {
        const copy = { ...r, _idx: idx };
        let mg = getVal(copy, ['M_G_corr', 'M_G', 'mg_corr']);
        let bprpCorr = getVal(copy, ['bp_rp_corr', 'bp-rp_corr', 'bprp_corr']);
        const gmag = getVal(copy, ['Gmag', 'phot_g_mean_mag', 'G']);
        const bpmag = getVal(copy, ['BPmag', 'phot_bp_mean_mag']);
        const rpmag = getVal(copy, ['RPmag', 'phot_rp_mean_mag']);
        const ag = getVal(copy, ['A_G', 'ag']);
        const bprp = getVal(copy, ['BP-RP', 'bp_rp', 'bprp']);

        const extAv = Number.isFinite(av) ? av : 0;
        const extG = Number.isFinite(ag) ? ag : extAv * STARNEST_AG_AV;
        const extBP = extAv * STARNEST_ABP_AV;
        const extRP = extAv * STARNEST_ARP_AV;

        if (!Number.isFinite(mg) && Number.isFinite(gmag)) mg = gmag - distMod - extG;
        if (!Number.isFinite(bprpCorr)) {
          if (Number.isFinite(bpmag) && Number.isFinite(rpmag)) bprpCorr = (bpmag - extBP) - (rpmag - extRP);
          else if (Number.isFinite(bprp)) bprpCorr = bprp - (extBP - extRP);
        }
        copy.M_G_corr = mg;
        copy.bp_rp_corr = bprpCorr;
        copy['bp-rp_corr'] = bprpCorr;
        copy['BP-RP'] = Number.isFinite(bprp) ? bprp : bprpCorr;
        copy.bprp_corr_ok = Number.isFinite(bprpCorr) && bprpCorr >= STARNEST_BPRP0_RANGE[0] && bprpCorr <= STARNEST_BPRP0_RANGE[1];
        assignDerived(copy, mg, bprpCorr);
        return copy;
      });
    }

    function wlinregressJS(x, y, sig) {
      // weighted least squares y = a + b x with weights 1/sig^2; error scaled by sqrt(chi2_red) when > 1
      const n = x.length;
      const w = sig.map(s => (Number.isFinite(s) && s > 0) ? 1 / (s * s) : 0);
      const sw = w.reduce((a, b) => a + b, 0);
      if (n < 3 || sw <= 0) return linregressJS(x, y);
      const mx = x.reduce((a, xi, i) => a + w[i] * xi, 0) / sw;
      const my = y.reduce((a, yi, i) => a + w[i] * yi, 0) / sw;
      let sxx = 0, sxy = 0, syy = 0;
      for (let i = 0; i < n; i++) { sxx += w[i] * (x[i] - mx) ** 2; sxy += w[i] * (x[i] - mx) * (y[i] - my); syy += w[i] * (y[i] - my) ** 2; }
      const slope = sxx > 0 ? sxy / sxx : 0;
      const intercept = my - slope * mx;
      let chi2 = 0;
      for (let i = 0; i < n; i++) chi2 += w[i] * (y[i] - (intercept + slope * x[i])) ** 2;
      const chi2r = n > 2 ? chi2 / (n - 2) : 1;
      const se = sxx > 0 ? Math.sqrt(Math.max(1, chi2r) / sxx) : 0;
      const r2 = (sxx > 0 && syy > 0) ? (sxy * sxy) / (sxx * syy) : 0;
      return { slope, intercept, r2, se };
    }

    function linregressJS(x, y) {
      const n = x.length;
      if (n < 2) return { slope: -1.35, intercept: 1.0, r2: 0, se: 0 };
      const meanX = x.reduce((a, b) => a + b, 0) / n;
      const meanY = y.reduce((a, b) => a + b, 0) / n;
      let sxx = 0, syy = 0, sxy = 0;
      for (let i = 0; i < n; i++) {
        const dx = x[i] - meanX;
        const dy = y[i] - meanY;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
      }
      const slope = sxx > 0 ? sxy / sxx : 0;
      const intercept = meanY - slope * meanX;
      const r2 = (sxx > 0 && syy > 0) ? (sxy * sxy) / (sxx * syy) : 0;
      let ssr = 0;
      for (let i = 0; i < n; i++) {
        const pred = slope * x[i] + intercept;
        ssr += (y[i] - pred) ** 2;
      }
      const se = (n > 2 && sxx > 0) ? Math.sqrt((ssr / (n - 2)) / sxx) : 0;
      return { slope, intercept, r2, se };
    }

    function trapzJS(y, x) {
      let sum = 0;
      for (let i = 0; i < x.length - 1; i++) {
        sum += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i]);
      }
      return sum;
    }

    function det3x3(cov) {
      const a = cov[0][0], b = cov[0][1], c = cov[0][2];
      const d = cov[1][0], e = cov[1][1], f = cov[1][2];
      const g = cov[2][0], h = cov[2][1], i = cov[2][2];
      return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    }

    function inv3x3(cov) {
      const a = cov[0][0], b = cov[0][1], c = cov[0][2];
      const d = cov[1][0], e = cov[1][1], f = cov[1][2];
      const g = cov[2][0], h = cov[2][1], i = cov[2][2];
      const det = det3x3(cov);
      const safeDet = Math.abs(det) < 1e-300 ? 1e-300 : det;
      return {
        det: Math.abs(safeDet),
        inv: [
          [(e * i - f * h) / safeDet, (c * h - b * i) / safeDet, (b * f - c * e) / safeDet],
          [(f * g - d * i) / safeDet, (a * i - c * g) / safeDet, (c * d - a * f) / safeDet],
          [(d * h - e * g) / safeDet, (b * g - a * h) / safeDet, (a * e - b * d) / safeDet]
        ]
      };
    }

    function mulberry32(seed) {
      let a = seed >>> 0;
      return function () {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
      };
    }

    // M-step from responsibilities (weights nk, means, full covariances + reg_covar on the diagonal)
    function gmmMStep(X, resp, K, reg) {
      const N = X.length;
      const pis = [], mus = [], covs = [];
      for (let k = 0; k < K; k++) {
        let nk = 0;
        for (let n = 0; n < N; n++) nk += resp[n][k];
        nk = Math.max(nk, 1e-10);
        const mu = [0, 0, 0];
        for (let n = 0; n < N; n++) { const r = resp[n][k]; mu[0] += r * X[n][0]; mu[1] += r * X[n][1]; mu[2] += r * X[n][2]; }
        mu[0] /= nk; mu[1] /= nk; mu[2] /= nk;
        const c = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        for (let n = 0; n < N; n++) {
          const r = resp[n][k];
          const d0 = X[n][0] - mu[0], d1 = X[n][1] - mu[1], d2 = X[n][2] - mu[2];
          c[0][0] += r * d0 * d0; c[0][1] += r * d0 * d1; c[0][2] += r * d0 * d2;
          c[1][1] += r * d1 * d1; c[1][2] += r * d1 * d2; c[2][2] += r * d2 * d2;
        }
        c[0][0] = c[0][0] / nk + reg; c[1][1] = c[1][1] / nk + reg; c[2][2] = c[2][2] / nk + reg;
        c[0][1] = c[1][0] = c[0][1] / nk; c[0][2] = c[2][0] = c[0][2] / nk; c[1][2] = c[2][1] = c[1][2] / nk;
        pis.push(nk / N); mus.push(mu); covs.push(c);
      }
      return { pis, mus, covs };
    }

    // E-step: responsibilities + total log-likelihood (log-sum-exp, no underflow)
    function gmmEStep(X, p, K) {
      const N = X.length;
      const inv = p.covs.map(inv3x3);
      const logConst = inv.map((iv, k) => Math.log(p.pis[k] + 1e-300) - 0.5 * (3 * Math.log(2 * Math.PI) + Math.log(iv.det)));
      const resp = new Array(N);
      let ll = 0;
      const lp = new Array(K);
      for (let n = 0; n < N; n++) {
        const x = X[n];
        let mx = -Infinity;
        for (let k = 0; k < K; k++) {
          const dx0 = x[0] - p.mus[k][0], dx1 = x[1] - p.mus[k][1], dx2 = x[2] - p.mus[k][2];
          const m = inv[k].inv;
          const quad = dx0 * (m[0][0] * dx0 + m[0][1] * dx1 + m[0][2] * dx2) +
                       dx1 * (m[1][0] * dx0 + m[1][1] * dx1 + m[1][2] * dx2) +
                       dx2 * (m[2][0] * dx0 + m[2][1] * dx1 + m[2][2] * dx2);
          lp[k] = logConst[k] - 0.5 * Math.max(0, quad);
          if (lp[k] > mx) mx = lp[k];
        }
        let sum = 0;
        for (let k = 0; k < K; k++) sum += Math.exp(lp[k] - mx);
        const lse = mx + Math.log(sum);
        ll += lse;
        const row = new Array(K);
        for (let k = 0; k < K; k++) row[k] = Math.exp(lp[k] - lse);
        resp[n] = row;
      }
      return { resp, ll };
    }

    function gmmSingleRun(X, K, rng, maxIter, tol, reg) {
      const N = X.length;
      // k-means++ seeding, then a few Lloyd steps -> hard responsibilities (as scikit-learn's default init)
      const cent = [[...X[Math.floor(rng() * N)]]];
      while (cent.length < K) {
        const d2 = X.map(x => Math.min(...cent.map(c => (x[0] - c[0]) ** 2 + (x[1] - c[1]) ** 2 + (x[2] - c[2]) ** 2)));
        const tot = d2.reduce((a, b) => a + b, 0);
        let t = rng() * tot, idx = 0;
        for (; idx < N - 1; idx++) { t -= d2[idx]; if (t <= 0) break; }
        cent.push([...X[idx]]);
      }
      let lab = new Array(N).fill(0);
      for (let it = 0; it < 8; it++) {
        for (let n = 0; n < N; n++) {
          let best = 0, bd = Infinity;
          for (let k = 0; k < K; k++) { const d = (X[n][0] - cent[k][0]) ** 2 + (X[n][1] - cent[k][1]) ** 2 + (X[n][2] - cent[k][2]) ** 2; if (d < bd) { bd = d; best = k; } }
          lab[n] = best;
        }
        const sums = Array.from({ length: K }, () => [0, 0, 0, 0]);
        for (let n = 0; n < N; n++) { const s = sums[lab[n]]; s[0] += X[n][0]; s[1] += X[n][1]; s[2] += X[n][2]; s[3]++; }
        for (let k = 0; k < K; k++) if (sums[k][3] > 0) cent[k] = [sums[k][0] / sums[k][3], sums[k][1] / sums[k][3], sums[k][2] / sums[k][3]];
      }
      let resp = Array.from({ length: N }, (_, n) => { const r = new Array(K).fill(0); r[lab[n]] = 1; return r; });
      let params = gmmMStep(X, resp, K, reg);
      let prev = -Infinity, iters = 0, converged = false;
      for (; iters < maxIter; iters++) {
        const e = gmmEStep(X, params, K);
        const meanLL = e.ll / N;
        if (Math.abs(meanLL - prev) < tol) { converged = true; break; }
        prev = meanLL;
        params = gmmMStep(X, e.resp, K, reg);
      }
      const fin = gmmEStep(X, params, K);        // final E-step: responsibilities & LL consistent with final parameters
      return { params, resp: fin.resp, ll: fin.ll, converged, iters };
    }

    // GMM_N_INIT / GMM_TOL / GMM_MAX_ITER mirror N_INIT / GMM_TOL / GMM_MAX_ITER in gmm.py (tight tol: EM must reach the true optimum or members/non-members shift)
    const GMM_N_INIT = 5, GMM_TOL = 1e-6, GMM_MAX_ITER = 1000, GMM_REG = 1e-6;
    const GMM_COLAB_TOL = 1e-3, GMM_COLAB_MAX_ITER = 100;   // sklearn GaussianMixture defaults = what the Colab script uses; 'converged' option uses GMM_TOL / GMM_MAX_ITER
    const GMM_MIN_STARS_PER_COMPONENT = 5;    // = MIN_COMPONENT_STARS in gmm.py. A 3-D covariance needs > 3 points to be full-rank; fewer = collapsed likelihood spike
    function fitFullGMM3D(scaled, K, nInit = GMM_N_INIT, maxIter = GMM_MAX_ITER, tol = GMM_TOL, reg = GMM_REG) {
      const N = scaled.length;
      let best = null;
      for (let r = 0; r < nInit; r++) {
        const run = gmmSingleRun(scaled, K, mulberry32(42 + 1000 * K + r), maxIter, tol, reg);
        // a component sitting on a handful of stars has a near-singular covariance -> huge density -> BIC "wins" -> tiny fake cluster
        if (Math.min(...run.params.pis) * N < GMM_MIN_STARS_PER_COMPONENT) continue;
        if (!best || run.ll > best.ll) best = run;
      }
      if (!best) return { K, bic: Infinity, cIdx: -1, probs: null, converged: false, degenerate: true };
      const nParams = K * 9 + (K - 1);              // K*(3 means + 6 covariance) + (K-1) weights
      const bic = -2.0 * best.ll + nParams * Math.log(N);
      const dets = best.params.covs.map(c => Math.abs(det3x3(c)));
      let cIdx = 0;
      for (let k = 1; k < K; k++) if (dets[k] < dets[cIdx]) cIdx = k;    // tightest blob = cluster
      return { K, bic, cIdx, probs: best.resp.map(r => r[cIdx]), converged: best.converged, params: best.params, resp: best.resp };
    }

    // ---- 3-D cluster map (mirrors kin3d() in main.py) --------------------------------------------------------------------
    // TERMS: sigma_d = d * e_Plx/Plx is the line-of-sight distance error of one star; sky position is known to ~mas, so the
    // error is purely radial -> a round cluster is smeared into a needle pointing at the Sun. Modes:
    //   precision : keep only stars with e_Plx/Plx <= maxFrac, distance 1000/Plx
    //   shrink    : empirical-Bayes parallax  plx* = (plx/s^2 + plx_c/t^2)/(1/s^2 + 1/t^2), t = plx_c^2 * depth/1000 (assumes cluster depth)
    //   raw       : 1000/Plx for every star (tiny Plx -> millions of pc)
    const ICRS_TO_GAL = [[-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
                         [ 0.4941094278755837, -0.4448296299600112,  0.7469822444972189],
                         [-0.8676661490190047, -0.1980763734312015,  0.4559837761750669]];
    function kin3dCompute(ra, dec, plx, eplx, o = {}) {
      const mode = o.mode || 'precision';
      const maxFrac = Number.isFinite(o.maxFrac) && o.maxFrac > 0 ? o.maxFrac : 0.05;
      const depth = Number.isFinite(o.depthPc) && o.depthPc > 0 ? o.depthPc : 5;
      const frame = o.frame || 'galactic', centre = o.centre !== false;
      const n = ra.length, hasErr = eplx.some(e => Number.isFinite(e) && e > 0);
      const all = ra.map((_, i) => i);
      let keep = all, plxUse = plx.slice(), note = '';
      if (mode === 'precision') {
        if (!hasErr) note = 'No e_Plx column in this file: precision cut skipped (raw 1/Plx).';
        else keep = all.filter(i => Number.isFinite(eplx[i]) && eplx[i] > 0 && eplx[i] / plx[i] <= maxFrac);
      } else if (mode === 'shrink') {
        if (!hasErr) note = 'No e_Plx column in this file: shrinkage skipped (raw 1/Plx).';
        else {
          let sw = 0, swp = 0;
          all.forEach(i => { if (Number.isFinite(eplx[i]) && eplx[i] > 0) { const w = 1 / (eplx[i] * eplx[i]); sw += w; swp += w * plx[i]; } });
          const pc = swp / sw, wt = 1 / Math.pow(pc * pc * depth / 1000, 2);
          all.forEach(i => { if (Number.isFinite(eplx[i]) && eplx[i] > 0) { const w = 1 / (eplx[i] * eplx[i]); plxUse[i] = (plx[i] * w + pc * wt) / (w + wt); } });
        }
      }
      const P = keep.map(i => {
        const d = 1000 / plxUse[i], a = ra[i] * Math.PI / 180, b = dec[i] * Math.PI / 180;
        const v = [d * Math.cos(b) * Math.cos(a), d * Math.cos(b) * Math.sin(a), d * Math.sin(b)];
        return frame === 'galactic' ? ICRS_TO_GAL.map(r => r[0] * v[0] + r[1] * v[1] + r[2] * v[2]) : v;
      });
      if (P.length < 3) throw new Error(`3D map: only ${P.length} stars pass the distance filter (mode ${mode}${mode === 'precision' ? ', sigma_plx/plx <= ' + maxFrac : ''}). Loosen it or use another mode.`);
      const med = [0, 1, 2].map(k => medianOf(P.map(p => p[k])));
      let diag = null;
      if (P.length >= 4) {                       // principal axis of the cloud vs line of sight (needle test)
        const mu = [0, 1, 2].map(k => P.reduce((t, p) => t + p[k], 0) / P.length);     // covariance about the MEAN (same as np.cov in main.py)
        const Q = P.map(p => [p[0] - mu[0], p[1] - mu[1], p[2] - mu[2]]);
        const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        Q.forEach(q => { for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) C[a][b] += q[a] * q[b] / Q.length; });
        const mv = v => [0, 1, 2].map(a => C[a][0] * v[0] + C[a][1] * v[1] + C[a][2] * v[2]);
        let v = [1, 0.7, 0.4];
        for (let it = 0; it < 300; it++) { const w = mv(v), nn = Math.hypot(w[0], w[1], w[2]) || 1; v = w.map(x => x / nn); }
        const mvv = mv(v), l1 = v[0] * mvv[0] + v[1] * mvv[1] + v[2] * mvv[2];
        const rest = Math.max((C[0][0] + C[1][1] + C[2][2] - l1) / 2, 1e-12), mn = Math.hypot(med[0], med[1], med[2]);
        diag = { major_sd_pc: Math.sqrt(l1), minor_sd_pc: Math.sqrt(rest), elongation: Math.sqrt(l1 / rest),
                 angle_to_los_deg: Math.acos(Math.min(1, Math.abs((v[0] * med[0] + v[1] * med[1] + v[2] * med[2]) / mn))) * 180 / Math.PI };
        // Observed radial spread vs the spread parallax errors alone would produce (sigma_d = d*e_Plx/Plx, per star, RMS)
        const los = med.map(x => x / mn), mlos = mv(los);
        diag.radial_sd_pc = Math.sqrt(Math.max(los[0] * mlos[0] + los[1] * mlos[1] + los[2] * mlos[2], 0));
        diag.transverse_sd_pc = Math.sqrt(Math.max((C[0][0] + C[1][1] + C[2][2] - diag.radial_sd_pc ** 2) / 2, 0));
        if (mode !== 'shrink' && hasErr) {
          const fin = keep.filter(i => Number.isFinite(eplx[i]) && eplx[i] > 0);
          if (fin.length) {
            diag.expected_err_sd_pc = Math.sqrt(fin.reduce((t, i) => t + Math.pow(1000 / plxUse[i] * eplx[i] / plxUse[i], 2), 0) / fin.length);
            diag.implied_depth_pc = Math.sqrt(Math.max(diag.radial_sd_pc ** 2 - diag.expected_err_sd_pc ** 2, 0));
          }
        }
      }
      const sh = centre ? med : [0, 0, 0];
      return {
        keep, x: P.map(p => p[0] - sh[0]), y: P.map(p => p[1] - sh[1]), z: P.map(p => p[2] - sh[2]),
        d: keep.map(i => 1000 / plxUse[i]), frac: keep.map(i => eplx[i] / plx[i]),
        meta: { mode, frame, centred: centre, max_frac: maxFrac, depth_pc: depth, n_used: P.length, n_dropped: n - P.length,
                median_d_pc: medianOf(keep.map(i => 1000 / plxUse[i])), note, diag }
      };
    }

    // =====================================================================
    // BROWSER-SIDE IRSA DUST + CMD 3.8 CLIENT (no Python backend; direct fetch from this page)
    // Works ONLY if the two servers answer with CORS permission (Access-Control-Allow-Origin). If they do not,
    // the browser blocks the reply and every call below ends in NetBlockedError -> shown in the workflow tracker.
    // Mirrors extinction_isochrone.py (same form fields, same parsing, same checks).
    // =====================================================================
    const BNET_CFG = {
      IRSA_URL: 'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust',
      IRSA_MIN_INTERVAL_MS: 300,
      IRSA_ROUND_DEG: 0.05,            // per_star: stars in the same 0.05 deg cell share ONE IRSA query
      IRSA_MAX_QUERIES: 150,
      CMD_URL: 'https://stev.oapd.inaf.it/cgi-bin/cmd',      // submit CGI of the CMD form page /cgi-bin/cmd_3.8
      CMD_VERSION: '3.8',
      PHOTSYS: 'gaiaEDR3',             // Gaia DR3 photometry = EDR3 photometry (Vega)
      PHOTSYS_VERSION: 'YBCnewVega',
      CMD_MIN_INTERVAL_MS: 2000,       // be gentle with the shared server
      TIMEOUT_MS: 180000,
      LOGAGE_RANGE: [6.6, 10.13], MH_RANGE: [-2.5, 0.8], AV_RANGE: [0, 20],
      ISO_MAX_LABEL: 7,                // 0 PMS, 1 MS, 2 SGB, 3 RGB, 4 core-He, 5-7 E-AGB (TP-AGB / post-AGB dropped)
      CHECK_DG_RANGE: [0.60, 1.00],    // expected dG / A_V
      CHECK_DC_RANGE: [0.25, 0.65],    // expected d(BP-RP) / A_V
      FIT_NEAR_MAG: 0.2,
      RA_KEYS: ['RA_ICRS', 'RAJ2000', '_RAJ2000', 'ra', 'RA'], DE_KEYS: ['DE_ICRS', 'DEJ2000', '_DEJ2000', 'dec', 'DEC'],
      ISO_COLS: ['Gmag', 'G_BPmag', 'G_RPmag', 'Mass', 'logL', 'logTe', 'logg'],
      ISO_EXTRA: ['logAge', 'MH', 'Mini', 'label'],
      ISO_ALIASES: { Gmag: ['Gmag'], G_BPmag: ['G_BPmag', 'G_BPmag_bright', 'BPmag'], G_RPmag: ['G_RPmag', 'RPmag'], Mass: ['Mass'], logL: ['logL'], logTe: ['logTe'], logg: ['logg'] },
      CMD_BASE_PARAMS: {
        cmd_version: '3.8', track_omegai: '0.00', track_parsec: 'parsec_CAF09_v1.2S', track_colibri: 'parsec_CAF09_v1.2S_S_LMC_08_web',
        track_postagb: 'no', n_inTPC: '10', eta_reimers: '0.2', kind_interp: '1', kind_postagb: '-1', dust_sourceM: 'dpmod60alox40',
        dust_sourceC: 'AMCSIC15', kind_mag: '2', kind_dust: '0', extinction_coeff: 'constant', extinction_curve: 'cardelli', kind_LPV: '3',
        imf_file: 'tab_imf/imf_kroupa_orig.dat', output_kind: '0', output_evstage: '1', lf_maginf: '-15', lf_magsup: '20',
        lf_deltamag: '0.5', sim_mtot: '1.0e4', submit_form: 'Submit'
      }
    };
    class NetBlockedError extends Error {}        // browser refused / could not read the reply (CORS, offline, DNS)
    class CmdRequestError extends Error {}        // server answered but rejected the request or returned no isochrone
    class NetworkNotAllowed extends Error {}      // cache-only call and the value is not cached
    const BNET = { irsa: new Map(), parsec: new Map(), irsaLast: 0, cmdLast: 0, cmdChain: Promise.resolve() };
    const bnetSleep = ms => new Promise(r => setTimeout(r, ms));

    async function bnetFetch(url, opts, what) {
      const ctl = new AbortController();
      const timer = setTimeout(() => ctl.abort(), BNET_CFG.TIMEOUT_MS);
      try {
        return await fetch(url, { ...(opts || {}), signal: ctl.signal });
      } catch (e) {
        if (e && e.name === 'AbortError') throw new Error(`${what} timed out after ${BNET_CFG.TIMEOUT_MS / 1000} s`);
        let host = url; try { host = new URL(url).host; } catch (_) { /* keep url */ }
        throw new NetBlockedError(`${what}: the browser could not read the reply from ${host} (the server does not allow cross-origin reads, or you are offline / the host is down)`);
      } finally { clearTimeout(timer); }
    }

    // ---------------- IRSA ----------------
    function bnetParseIrsa(xml) {
      if (/<results[^>]*status\s*=\s*"error"/i.test(xml) || /<message>/i.test(xml) && !/<meanValueSFD>/i.test(xml)) {
        const m = xml.match(/<message>([\s\S]*?)<\/message>/i);
        throw new CmdRequestError('IRSA Dust answered with an error: ' + (m ? m[1].trim() : xml.replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').slice(0, 200)));
      }
      const m = xml.match(/<meanValueSFD>\s*([-+0-9.eE]+)/i);          // = astroquery column 'ext SFD mean' (mean over the IRSA region, >= 2 deg)
      if (!m || !Number.isFinite(parseFloat(m[1]))) throw new CmdRequestError('IRSA reply has no <meanValueSFD> value: ' + xml.replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').slice(0, 200));
      return parseFloat(m[1]);
    }
    function bnetCell(ra, dec) { const g = BNET_CFG.IRSA_ROUND_DEG; return [Number((Math.round(ra / g) * g).toFixed(4)), Number((Math.round(dec / g) * g).toFixed(4))]; }
    async function bnetIrsaCached(cell, cacheOnly) {
      const key = cell.join(',');
      if (BNET.irsa.has(key)) return BNET.irsa.get(key);
      if (cacheOnly) throw new NetworkNotAllowed('IRSA value not cached');
      const wait = BNET_CFG.IRSA_MIN_INTERVAL_MS - (Date.now() - BNET.irsaLast);
      if (wait > 0) await bnetSleep(wait);
      let xml;
      try {
        const r = await bnetFetch(`${BNET_CFG.IRSA_URL}?locstr=${cell[0]}+${cell[1]}+equ+j2000`, { method: 'GET' }, 'IRSA Dust query');
        if (!r.ok) throw new Error(`IRSA Dust answered HTTP ${r.status}`);
        xml = await r.text();
      } finally { BNET.irsaLast = Date.now(); }
      const v = bnetParseIrsa(xml);
      BNET.irsa.set(key, v);
      return v;
    }
    // returns {values: E(B-V) per row (NaN if failed), queries (network calls made), failed}
    async function bnetEbv(rows, mode, cacheOnly) {
      const ra = rows.map(r => getVal(r, BNET_CFG.RA_KEYS)), de = rows.map(r => getVal(r, BNET_CFG.DE_KEYS));
      const fin = a => a.filter(Number.isFinite);
      if (!fin(ra).length || !fin(de).length) throw new Error('No RA/Dec columns in the members file: IRSA needs coordinates.');
      if (mode === 'center') {
        const cell = bnetCell(medianOf(fin(ra)), medianOf(fin(de)));
        const had = BNET.irsa.has(cell.join(','));
        const v = await bnetIrsaCached(cell, cacheOnly);
        return { values: rows.map(() => v), queries: had ? 0 : 1, failed: 0 };
      }
      if (mode === 'per_star') {
        const cells = rows.map((r, i) => (Number.isFinite(ra[i]) && Number.isFinite(de[i])) ? bnetCell(ra[i], de[i]) : null);
        const uniq = [...new Set(cells.filter(Boolean).map(c => c.join(',')))].sort();
        const todo = uniq.filter(k => !BNET.irsa.has(k));
        if (cacheOnly && todo.length) throw new NetworkNotAllowed('IRSA values not cached');
        if (todo.length > BNET_CFG.IRSA_MAX_QUERIES) throw new Error(`per_star needs ${todo.length} IRSA queries (limit ${BNET_CFG.IRSA_MAX_QUERIES}, cell ${BNET_CFG.IRSA_ROUND_DEG} deg). Use E(B-V) mode center or fixed.`);
        const vals = {}; let failed = 0, firstErr = null;
        for (const k of uniq) {
          try { vals[k] = await bnetIrsaCached(k.split(',').map(Number), false); }
          catch (e) { if (e instanceof NetBlockedError && !Object.keys(vals).length) throw e; vals[k] = NaN; failed++; firstErr = firstErr || e; }
        }
        if (failed === uniq.length && firstErr) throw firstErr;
        return { values: cells.map(c => c ? vals[c.join(',')] : NaN), queries: todo.length, failed };
      }
      throw new Error('Unknown E(B-V) mode: ' + mode);
    }

    // ---------------- CMD 3.8 ----------------
    function bnetCmdParams(logage, mh, av) {
      const C = BNET_CFG;
      if (!(logage >= C.LOGAGE_RANGE[0] && logage <= C.LOGAGE_RANGE[1])) throw new Error(`logAge ${logage} outside the PARSEC range ${C.LOGAGE_RANGE}`);
      if (!(mh >= C.MH_RANGE[0] && mh <= C.MH_RANGE[1])) throw new Error(`[M/H] ${mh} outside the accepted range ${C.MH_RANGE}`);
      if (!(av >= C.AV_RANGE[0] && av <= C.AV_RANGE[1])) throw new Error(`A_V ${av} outside ${C.AV_RANGE} mag`);
      return {
        ...C.CMD_BASE_PARAMS,
        photsys_file: `YBC_tab_mag_odfnew/tab_mag_${C.PHOTSYS}.dat`, photsys_version: C.PHOTSYS_VERSION,
        isoc_isagelog: '1', isoc_lagelow: logage.toFixed(4), isoc_lageupp: logage.toFixed(4), isoc_dlage: '0',
        isoc_ismetlog: '1', isoc_metlow: mh.toFixed(4), isoc_metupp: mh.toFixed(4), isoc_dmet: '0',
        extinction_av: av.toFixed(4),
        isoc_agelow: '1.0e9', isoc_ageupp: '1.0e10', isoc_dage: '0.0', isoc_zlow: '0.0152', isoc_zupp: '0.03', isoc_dz: '0.0'
      };
    }
    async function bnetInflate(buf) {
      const u8 = new Uint8Array(buf);
      const run = async fmt => new Uint8Array(await new Response(new Blob([u8]).stream().pipeThrough(new DecompressionStream(fmt))).arrayBuffer());
      if (u8[0] === 0x1f && u8[1] === 0x8b) return run('gzip');
      if (u8[0] === 0x78 && [0x01, 0x5e, 0x9c, 0xda].includes(u8[1])) { try { return await run('deflate'); } catch (_) { /* plain text that starts with x */ } }
      return u8;
    }
    async function bnetFetchParsecRaw(logage, mh, av) {
      const params = bnetCmdParams(logage, mh, av);
      const prev = BNET.cmdChain;                                      // one request at a time + minimum spacing
      let release; BNET.cmdChain = new Promise(r => { release = r; });
      await prev;
      try {
        const wait = BNET_CFG.CMD_MIN_INTERVAL_MS - (Date.now() - BNET.cmdLast);
        if (wait > 0) await bnetSleep(wait);
        const qs = new URLSearchParams(params).toString();
        const r = await bnetFetch(`${BNET_CFG.CMD_URL}?${qs}`, { method: 'POST' }, 'CMD 3.8 request');   // same as ezpadova: POST with the form fields in the query string
        if (!r.ok) throw new Error(`CMD server answered HTTP ${r.status}`);
        const page = await r.text();
        const m = page.match(/output\d+/);
        if (!m) throw new CmdRequestError('CMD did not return an output file. Server said: ' + page.replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').trim().slice(0, 300));
        const src = `${new URL(BNET_CFG.CMD_URL).origin}/tmp/${m[0]}.dat`;
        const r2 = await bnetFetch(src, { method: 'GET' }, 'CMD output download');
        if (!r2.ok) throw new Error(`CMD output file ${src} answered HTTP ${r2.status}`);
        const data = await bnetInflate(await r2.arrayBuffer());
        const text = new TextDecoder('utf-8').decode(data);
        if (!text.trim()) throw new Error('CMD output file is empty');
        return { text, params, info: { source_url: src } };
      } finally { BNET.cmdLast = Date.now(); release(); }
    }
    // header = the last '#' line before the first data row; returns {names, rows: number[][]}
    function bnetParseParsec(text) {
      let header = null; const rows = [];
      for (const line of text.split(/\r?\n/)) {
        const s = line.trim(); if (!s) continue;
        if (s.startsWith('#')) { if (!rows.length && s.includes('logAge') && (s.includes('Gmag') || s.includes('mbolmag'))) header = s.replace(/^#+/, '').trim().split(/\s+/); continue; }
        rows.push(s.split(/\s+/));
      }
      if (!header || !rows.length) throw new CmdRequestError('Could not find the column header / data rows in the PARSEC file (is it an isochrone table in a Gaia system?).');
      const w = header.length;
      return { names: header, rows: rows.filter(r => r.length >= w).map(r => r.slice(0, w).map(Number)) };
    }
    function bnetSelectIso(tbl) {
      const out = {}, rnd = v => Number.isFinite(v) ? Math.round(v * 1e6) / 1e6 : NaN;
      for (const std of BNET_CFG.ISO_COLS) {
        const found = BNET_CFG.ISO_ALIASES[std].find(a => tbl.names.includes(a));
        if (!found) throw new Error(`PARSEC output lacks '${std}'. Available: ${tbl.names.join(', ')}`);
        const j = tbl.names.indexOf(found); out[std] = tbl.rows.map(r => rnd(r[j]));
      }
      for (const c of BNET_CFG.ISO_EXTRA) if (tbl.names.includes(c)) { const j = tbl.names.indexOf(c); out[c] = tbl.rows.map(r => rnd(r[j])); }
      out.n = tbl.rows.length;
      return out;
    }
    function bnetValidate(iso, logage, mh) {
      if (iso.n < 20) throw new CmdRequestError(`PARSEC file has only ${iso.n} rows`);
      let ok = 0; for (let i = 0; i < iso.n; i++) if (Number.isFinite(iso.Gmag[i]) && Number.isFinite(iso.G_BPmag[i]) && Number.isFinite(iso.G_RPmag[i])) ok++;
      if (!(ok / iso.n > 0.9)) throw new CmdRequestError('PARSEC file has too many non-finite magnitudes');
      if (iso.logAge && Math.abs(medianOf(iso.logAge) - logage) > 0.01) throw new CmdRequestError(`PARSEC returned logAge ${medianOf(iso.logAge).toFixed(3)}, requested ${logage.toFixed(3)}`);
      if (iso.MH && Math.abs(medianOf(iso.MH) - mh) > 0.05) throw new CmdRequestError(`PARSEC returned [M/H] ${medianOf(iso.MH).toFixed(3)}, requested ${mh.toFixed(3)}`);
    }
    function bnetProvenance(raw, params, info) {
      const head = ['# ---- StarNest web app: PARSEC file provenance ----',
        `# retrieved_utc   : ${new Date().toISOString().replace('T', ' ').slice(0, 19)}`,
        `# service         : CMD ${BNET_CFG.CMD_VERSION}  (https://stev.oapd.inaf.it/cgi-bin/cmd_${BNET_CFG.CMD_VERSION})`,
        `# photometric sys : Gaia EDR3 = Gaia DR3 photometry, Vega  (${params.photsys_file}, ${params.photsys_version})`,
        `# logAge / [M/H]  : ${params.isoc_lagelow} / ${params.isoc_metlow}`,
        `# interstellar A_V: ${params.extinction_av} mag (${params.extinction_curve}, coefficients: ${params.extinction_coeff})`,
        `# tracks          : ${params.track_parsec} + ${params.track_colibri}`,
        '# fetched by      : the browser directly (no backend)',
        `# file url        : ${info.source_url}`,
        '# request params  : ' + JSON.stringify(params, Object.keys(params).sort()),
        '# ---- original CMD output follows unchanged ----'];
      return head.join('\n') + '\n' + raw;
    }
    async function bnetGetParsec(logage, mh, av, cacheOnly) {
      const key = [logage.toFixed(4), mh.toFixed(4), av.toFixed(3)].join('|');
      if (BNET.parsec.has(key)) return BNET.parsec.get(key);
      if (cacheOnly) throw new NetworkNotAllowed('PARSEC file not cached');
      const { text, params, info } = await bnetFetchParsecRaw(Number(logage.toFixed(4)), Number(mh.toFixed(4)), Number(av.toFixed(3)));
      const iso = bnetSelectIso(bnetParseParsec(text));
      bnetValidate(iso, logage, mh);
      const out = { iso, raw: bnetProvenance(text, params, info), params, info, key };
      BNET.parsec.set(key, out);
      if (BNET.parsec.size > 24) BNET.parsec.delete(BNET.parsec.keys().next().value);
      return out;
    }
    function bnetInterp(x, xp, fp) {                                    // numpy.interp (xp ascending, clamped at the ends)
      const n = xp.length;
      if (x <= xp[0]) return fp[0];
      if (x >= xp[n - 1]) return fp[n - 1];
      let lo = 0, hi = n - 1;
      while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xp[mid] <= x) lo = mid; else hi = mid; }
      const dx = xp[hi] - xp[lo];
      return dx === 0 ? fp[lo] : fp[lo] + (fp[hi] - fp[lo]) * (x - xp[lo]) / dx;
    }
    // A_V>0 isochrone vs A_V=0 one at equal initial mass: the shift must look like Gaia extinction (guards against CMD ignoring / doubling it)
    function bnetCheckShift(iso0, isoAv, av) {
      const out = { av, dG: null, dBPRP: null, ok: null, message: '' };
      if (av <= 0.0005) { out.ok = true; out.message = 'A_V = 0: nothing to check'; return out; }
      const key = (iso0.Mini && isoAv.Mini) ? 'Mini' : 'Mass';
      const ms = iso => {
        const idx = []; for (let i = 0; i < iso.n; i++) if (!iso.label || iso.label[i] === 1) idx.push(i);
        idx.sort((a, b) => iso[key][a] - iso[key][b]);
        return { k: idx.map(i => iso[key][i]), g: idx.map(i => iso.Gmag[i]), c: idx.map(i => iso.G_BPmag[i] - iso.G_RPmag[i]) };
      };
      const a = ms(iso0), b = ms(isoAv);
      if (a.k.length < 5 || b.k.length < 5) { out.message = 'too few main-sequence points to verify the reddening shift'; return out; }
      const dGs = b.k.map((x, i) => b.g[i] - bnetInterp(x, a.k, a.g)), dCs = b.k.map((x, i) => b.c[i] - bnetInterp(x, a.k, a.c));
      const dG = medianOf(dGs), dC = medianOf(dCs), C = BNET_CFG;
      const ok = dG >= C.CHECK_DG_RANGE[0] * av && dG <= C.CHECK_DG_RANGE[1] * av && dC >= C.CHECK_DC_RANGE[0] * av && dC <= C.CHECK_DC_RANGE[1] * av;
      Object.assign(out, { dG, dBPRP: dC, A_G_over_A_V: dG / av, E_BPRP_over_A_V: dC / av, ok,
        message: ok ? 'reddening applied by CMD as expected' : `UNEXPECTED reddening shift (dG=${dG.toFixed(3)}, d(BP-RP)=${dC.toFixed(3)} for A_V=${av.toFixed(3)}): do not trust this fit` });
      return out;
    }
    // distance (mag, CMD plane) from each star to the nearest point of the polyline isochrone
    function bnetFitQuality(color, mag, ic, im, near) {
      near = near === undefined ? BNET_CFG.FIT_NEAR_MAG : near;
      const pts = []; for (let i = 0; i < ic.length; i++) if (Number.isFinite(ic[i]) && Number.isFinite(im[i])) pts.push([ic[i], im[i]]);
      const st = []; for (let i = 0; i < color.length; i++) if (Number.isFinite(color[i]) && Number.isFinite(mag[i])) st.push([color[i], mag[i]]);
      if (!st.length || pts.length < 3) return { n: st.length, median: null, rms: null, frac_near: null, near };
      const s = [0]; for (let i = 1; i < pts.length; i++) s.push(s[i - 1] + Math.hypot(pts[i][0] - pts[i - 1][0], pts[i][1] - pts[i - 1][1]));
      const N = 3000, dx = [], dy = [], xs = pts.map(p => p[0]), ys = pts.map(p => p[1]);
      if (s[s.length - 1] > 0) for (let k = 0; k < N; k++) { const t = s[s.length - 1] * k / (N - 1); dx.push(bnetInterp(t, s, xs)); dy.push(bnetInterp(t, s, ys)); }
      else { dx.push(xs[0]); dy.push(ys[0]); }
      const d = st.map(([c, m]) => { let best = Infinity; for (let k = 0; k < dx.length; k++) { const v = (c - dx[k]) ** 2 + (m - dy[k]) ** 2; if (v < best) best = v; } return Math.sqrt(best); });
      return { n: d.length, median: medianOf(d), rms: Math.sqrt(d.reduce((a, v) => a + v * v, 0) / d.length), frac_near: d.filter(v => v <= near).length / d.length, near };
    }

    async function localBackendEngine(endpoint, options = {}) {
      if (endpoint === '/run-gmm') {
        const params = JSON.parse(options.body || '{}');
        const cacheKey = `gmm_${ options.body }|v${localStore.ver[params.filename] || 0}|e${uiEbv()}`;
        const hitGmm = localStore.cache[cacheKey];
        if (hitGmm) { restoreMembers(hitGmm); return hitGmm.result; }   // cache hit must restore downstream state too

        const rawRows = localStore.files[params.filename] || [];
        if (!rawRows.length) throw new Error('Uploaded file has no valid rows.');

        const missingCols = validateColumnsJS(rawRows);
        if (missingCols.length > 0) {
          throw new Error(`Missing required columns: [${missingCols.map(c => `'${c}'`).join(', ')}]`);
        }

        const ruweMax = params.ruwe_max !== undefined ? params.ruwe_max : 1.6;
        const epmraMax = params.epmra_max !== undefined ? params.epmra_max : 0.1;
        const epmdeMax = params.epmde_max !== undefined ? params.epmde_max : 0.1;
        const plxErrMin = params.plx_err_min !== undefined ? params.plx_err_min : 5.0;
        const gofalMax = params.gofal_max !== undefined ? params.gofal_max : 3.0;
        const nperMin = params.nper_min !== undefined ? params.nper_min : 10.0;
        const distanceLit = params.distance_lit !== undefined ? params.distance_lit : 341.0;
        const halfWidth = params.half_width !== undefined ? params.half_width : 195.5;
        const cutoff = params.cutoff !== undefined ? params.cutoff : 0.5;
        const plxZp = params.plx_zp !== undefined ? params.plx_zp : 0.0;
        const kParam = params.k;

        // Exact apply_quality_cuts from gmm.py (NaN rows dropped)
        const cleaned = rawRows.filter(r => {
          const ruwe = Number(r.RUWE);
          const epmra = Number(r.e_pmRA);
          const epmde = Number(r.e_pmDE);
          const plx = Number(r.Plx);
          const eplx = Number(r.e_Plx);
          const gofal = Number(r.gofAL);
          const nper = Number(r.Nper);
          return (
            Number.isFinite(ruwe) && ruwe < ruweMax &&
            Number.isFinite(epmra) && epmra < epmraMax &&
            Number.isFinite(plx) && Number.isFinite(eplx) && (plx / eplx) > plxErrMin &&
            Number.isFinite(epmde) && epmde < epmdeMax &&
            Number.isFinite(gofal) && gofal < gofalMax &&
            Number.isFinite(nper) && nper > nperMin
          );
        });

        const plxPositive = cleaned.filter(r => (Number(r.Plx) - plxZp) > 0).map(r => {
          const copy = { ...r };
          copy.distance_pc = 1000.0 / (Number(copy.Plx) - plxZp);
          return copy;
        });

        const workingData = plxPositive.filter(r => Math.abs(r.distance_pc - distanceLit) < halfWidth);
        const validFeatRows = workingData.filter(r =>
          Number.isFinite(Number(r.pmRA)) &&
          Number.isFinite(Number(r.pmDE)) &&
          Number.isFinite(Number(r.Plx))
        );

        const kLo = kParam ? Number(kParam) : 2;
        const kHi = kParam ? Number(kParam) : 6;
        if (validFeatRows.length < 10 * kHi) {
          throw new Error(`Only ${validFeatRows.length} stars in distance window; need >= ${10 * kHi}. Check distance / half-width / quality cuts.`);
        }

        const calcMedian = arr => {
          const s = [...arr].sort((a, b) => a - b);
          const mid = Math.floor(s.length / 2);
          return s.length % 2 !== 0 ? s[mid] : 0.5 * (s[mid - 1] + s[mid]);
        };
        const calcStd = arr => {
          if (arr.length < 2) return 0;
          const m = arr.reduce((a, b) => a + b, 0) / arr.length;
          return Math.sqrt(arr.reduce((acc, x) => acc + (x - m) ** 2, 0) / arr.length);
        };

        const pmraArr = validFeatRows.map(r => Number(r.pmRA));
        const pmdeArr = validFeatRows.map(r => Number(r.pmDE));
        const plxArr = validFeatRows.map(r => Number(r.Plx));

        const med0 = calcMedian(pmraArr), std0 = calcStd(pmraArr);
        const med1 = calcMedian(pmdeArr), std1 = calcStd(pmdeArr);
        const med2 = calcMedian(plxArr), std2 = calcStd(plxArr);
        if (!std0 || !std1 || !std2) {
          throw new Error('Zero variance in a GMM feature column.');
        }

        const scaled = validFeatRows.map((_, i) => [
          (pmraArr[i] - med0) / std0,
          (pmdeArr[i] - med1) / std1,
          (plxArr[i] - med2) / std2
        ]);

        const bicMap = {}, kTable = [];
        const useConverged = params.gmm_opt === 'converged';
        const optTol = useConverged ? GMM_TOL : GMM_COLAB_TOL, optIter = useConverged ? GMM_MAX_ITER : GMM_COLAB_MAX_ITER;
        let bestModel = null;
        for (let kk = kLo; kk <= kHi; kk++) {
          await new Promise(res => setTimeout(res, 0));     // let the loading spinner repaint between K fits (tight-tol EM takes seconds)
          const m = fitFullGMM3D(scaled, kk, GMM_N_INIT, optIter, optTol);   // 'colab' = sklearn stock stop rule (tol 1e-3, 100 it); 'converged' = tol 1e-6, 1000 it
          bicMap[kk] = Number.isFinite(m.bic) ? m.bic : null;
          if (m.probs) {                                   // robustness table: the member count this K would give (cluster = smallest det)
            const nm = m.probs.filter(q => q >= cutoff).length, nn = m.probs.filter(q => q <= 1 - cutoff).length;
            kTable.push({ k: kk, bic: m.bic, members: nm, non_members: nn, min_comp_stars: Math.round(Math.min(...m.params.pis) * m.probs.length), converged: m.converged });
          } else kTable.push({ k: kk, bic: null, members: null, non_members: null, min_comp_stars: null, converged: false, degenerate: true });
          if (Number.isFinite(m.bic) && (!bestModel || m.bic < bestModel.bic)) bestModel = m;
        }
        if (!bestModel) throw new Error(`GMM: every K in ${kLo}..${kHi} produced a component with < ${GMM_MIN_STARS_PER_COMPONENT} stars (only ${validFeatRows.length} stars in the window). Widen the distance half-width or relax the quality cuts.`);

        validFeatRows.forEach((r, i) => {
          r.PMemb = bestModel.probs[i];          // unrounded: rounding to 6 dp flipped stars sitting exactly on the cutoff
          r.gmm_label = r.PMemb >= cutoff ? 1 : 0;
        });

        let members = validFeatRows.filter(r => r.PMemb >= cutoff);
        let nonMembers = validFeatRows.filter(r => r.PMemb <= (1.0 - cutoff));

        let mssMetric = 0;
        if (members.length >= 2 && nonMembers.length >= 2) {
          const fGetters = [r => Number(r.pmRA), r => Number(r.pmDE), r => Number(r.Plx)];
          const metrics = fGetters.map(fn => {
            const sm = calcStd(members.map(fn));
            const sn = calcStd(nonMembers.map(fn));
            return (sn - sm) / Math.max(sm, sn, 1e-7);
          });
          mssMetric = (metrics[0] + metrics[1] + metrics[2]) / 3.0;
        }

        const medMemPlx = members.length ? calcMedian(members.map(r => Number(r.Plx))) : NaN;
        const distFromMemPc = Number.isFinite(medMemPlx) ? 1000.0 / (medMemPlx - plxZp) : distanceLit;
        const dmCluster = 5.0 * Math.log10(distFromMemPc) - 5.0;

        // E(B-V) comes from the Isochrone page input (was hard-coded 0.15 for every cluster)
        const avUI = STARNEST_AV_PER_EBV * uiEbv();
        members = enrichPhotometryAndClassify(members, dmCluster, avUI);
        nonMembers = enrichPhotometryAndClassify(nonMembers, dmCluster, avUI);

        const verId = commitMembers(members);

        const counts = {};
        members.forEach(m => { counts[m.Star_Type] = (counts[m.Star_Type] || 0) + 1; });

        const result = {
          ra_dec: {
            members: {
              ra: members.map(r => getVal(r, ['RA_ICRS', 'RAJ2000', '_RAJ2000', 'ra', 'RA'])),
              dec: members.map(r => getVal(r, ['DE_ICRS', 'DEJ2000', '_DEJ2000', 'dec', 'DEC'])),
              color: members.map(r => STAR_TYPE_COLORS[r.Star_Type] || '#ef4444'),
              star_type: members.map(r => r.Star_Type),
              pmemb: members.map(r => Number(r.PMemb.toFixed(6)))
            },
            field: {
              ra: nonMembers.map(r => getVal(r, ['RA_ICRS', 'RAJ2000', '_RAJ2000', 'ra', 'RA'])),
              dec: nonMembers.map(r => getVal(r, ['DE_ICRS', 'DEJ2000', '_DEJ2000', 'dec', 'DEC']))
            }
          },
          vpd: {
            members: {
              pmra: members.map(r => Number(r.pmRA)),
              pmdec: members.map(r => Number(r.pmDE)),
              color: members.map(r => STAR_TYPE_COLORS[r.Star_Type] || '#ef4444')
            },
            field: {
              pmra: nonMembers.map(r => Number(r.pmRA)),
              pmdec: nonMembers.map(r => Number(r.pmDE))
            }
          },
          plx_Cmd: {
            members_plx: members.map(r => Number(r.Plx)),
            field_plx: nonMembers.map(r => Number(r.Plx)),
            members_color: members.map(r => Number(r['BP-RP'])),
            members_gmag: members.map(r => Number(r.Gmag)),
            field_color: nonMembers.map(r => Number(r['BP-RP'])),
            field_gmag: nonMembers.map(r => Number(r.Gmag))
          },
          stats: {
            n_input: cleaned.length,
            total: workingData.length,
            n_fit: validFeatRows.length,
            k: bestModel.K,
            bic: bicMap,
            cluster_component: bestModel.cIdx,
            converged: bestModel.converged,
            k_table: kTable,
            methods: { engine: 'in-browser JS EM (approximate sklearn equivalent)', optimizer: useConverged ? 'fully converged' : 'colab-compatible', tol: optTol, max_iter: optIter, n_init: GMM_N_INIT, reg_covar: GMM_REG, k_mode: kParam ? 'fixed' : 'BIC over 2..6', cutoff, distance_lit: distanceLit, half_width: halfWidth, plx_zp: plxZp },
            components: bestModel.params.pis.map((w, k) => ({
              k, weight: w, n_hard: bestModel.resp.reduce((a, r) => a + (r.indexOf(Math.max(...r)) === k ? 1 : 0), 0),
              det: det3x3(bestModel.params.covs[k]),
              pmRA: bestModel.params.mus[k][0] * std0 + med0, pmDE: bestModel.params.mus[k][1] * std1 + med1, Plx: bestModel.params.mus[k][2] * std2 + med2,
              sd_pmRA: Math.sqrt(bestModel.params.covs[k][0][0]) * std0, sd_pmDE: Math.sqrt(bestModel.params.covs[k][1][1]) * std1, sd_Plx: Math.sqrt(bestModel.params.covs[k][2][2]) * std2,
              is_cluster: k === bestModel.cIdx
            })),
            members: members.length,
            field: nonMembers.length,
            threshold: cutoff,
            mss_metric: Number(mssMetric.toFixed(4)),
            median_plx: Number.isFinite(medMemPlx) ? Number(medMemPlx.toFixed(4)) : 0,
            distance_from_members_pc: Number(distFromMemPc.toFixed(1))
          },
          classification_counts: counts
        };
        localStore.cache[cacheKey] = { result, members, verId };
        return result;
      }

      if (endpoint === '/compute-cmd') {
        const params = JSON.parse(options.body || '{}');
        const srcName = params.members_filename;
        const cacheKey = `cmd_${options.body}|v${localStore.ver[srcName] || 0}|i${localStore.isoVer}`;
        const hitCmd = localStore.cache[cacheKey];
        if (hitCmd) {
          localStore.analysis = hitCmd.analysis; localStore.latestExtinction = hitCmd.analysis.rows;
          localStore.latestIsochrone = hitCmd.isoRows; localStore.parsec = hitCmd.parsec || null;
          return hitCmd.res;
        }

        if (!(params.age >= STARNEST_LOGAGE_FLOOR)) {
          throw new Error(`logAge ${params.age} below PARSEC floor ${STARNEST_LOGAGE_FLOOR}`);
        }
        const mode = params.ebv_mode || 'center';
        const handling = params.iso_extinction || 'reddened';
        if (handling !== 'reddened' && handling !== 'intrinsic') throw new Error(`iso_extinction must be 'reddened' or 'intrinsic', got '${handling}'`);
        const allowNet = params.allow_network !== false;        // false = automatic refresh: cached IRSA / PARSEC values only, nothing is requested
        const typedGiven = params.ebv !== null && params.ebv !== undefined && params.ebv !== '' && Number.isFinite(Number(params.ebv));
        if (mode === 'fixed' && !typedGiven) throw new Error("mode='fixed' requires ebv_fixed");

        const rows = localStore.files[srcName] || localStore.latestMembers || [];
        if (!rows.length) {
          throw new Error('No member stars loaded. Run GMM Analysis or upload a custom CSV first.');
        }
        const age = Number(params.age), mh = Number(params.metallicity) || 0;
        const wf = { n_members: rows.length, irsa: {}, extinction: {}, parsec: {}, fit: {} };
        const notes = [];

        // ---- step 2: E(B-V) from IRSA (direct browser fetch), else the typed value ----
        const ebvTyped = typedGiven ? Number(params.ebv) : NaN;
        let irsaVals = null;
        if (mode === 'fixed') {
          wf.irsa = { state: 'manual', mode };
        } else {
          try {
            const r = await bnetEbv(rows, mode, !allowNet);
            irsaVals = r.values;
            wf.irsa = { state: 'done', mode, queries: r.queries, failed: r.failed, dropped_stars: 0 };
          } catch (e) {
            if (e instanceof NetworkNotAllowed) {
              if (!typedGiven) throw new Error('Press "Apply Extinction & Fetch PARSEC" to query IRSA (or choose mode fixed and type an E(B-V)).');
              wf.irsa = { state: 'pending', mode };
              notes.push('IRSA / PARSEC have not been queried yet: showing the typed E(B-V). Press "Apply Extinction & Fetch PARSEC".');
            } else {
              wf.irsa = { state: 'failed', mode, error: e.message || String(e) };
              if (!typedGiven) throw new Error(`${e.message || e}. Type an E(B-V) and choose mode fixed.`);
              notes.push(`IRSA query failed (${e.message || e}); the typed E(B-V) = ${ebvTyped} is used instead.`);
            }
          }
        }

        const plxZp = Number(params.plx_zp) || 0.0;
        const useClusterDm = params.use_cluster_dm !== false;
        const medPlx = medianOf(rows.map(r => getVal(r, ['Plx', 'parallax', 'plx']) - plxZp).filter(v => Number.isFinite(v) && v > 0));
        const dPc = Number.isFinite(medPlx) ? (1000.0 / medPlx) : Math.pow(10, ((params.dist_mod || 7.66) + 5) / 5);
        const clusterDm = 5.0 * Math.log10(dPc) - 5.0;

        const enriched = rows.map((r, idx) => {
          const copy = { ...r, _idx: idx };
          // per-star E(B-V): IRSA cell value > a column carried by the FILE (per_star only) > the typed value
          const fileEbv = getVal(r, ['E_BV', 'Extinction']);
          const rowEbv = irsaVals ? irsaVals[idx] : ((mode === 'per_star' && Number.isFinite(fileEbv)) ? fileEbv : ebvTyped);
          const rowAv = STARNEST_AV_PER_EBV * rowEbv;
          copy.E_BV = Number(rowEbv.toFixed(4));
          copy.A_G = Number((STARNEST_AG_AV * rowAv).toFixed(4));
          copy.A_BP = Number((STARNEST_ABP_AV * rowAv).toFixed(4));
          copy.A_RP = Number((STARNEST_ARP_AV * rowAv).toFixed(4));

          const rPlx = getVal(copy, ['Plx', 'parallax', 'plx']) - plxZp;
          const dm = useClusterDm
            ? clusterDm
            : (Number.isFinite(rPlx) && rPlx > 0 ? 5.0 * Math.log10(1000.0 / rPlx) - 5.0 : clusterDm);
          copy.distance_modulus = Number(dm.toFixed(4));

          const g = getVal(copy, ['Gmag', 'phot_g_mean_mag']);
          const bp = getVal(copy, ['BPmag', 'phot_bp_mean_mag']);
          const rp = getVal(copy, ['RPmag', 'phot_rp_mean_mag']);
          const bprp = getVal(copy, ['BP-RP', 'bp_rp']);

          copy.M_G_corr = Number.isFinite(g) ? Number((g - dm - copy.A_G).toFixed(4)) : getVal(r, ['M_G_corr']);
          if (Number.isFinite(bp) && Number.isFinite(rp)) {
            copy.M_BP_corr = Number((bp - dm - copy.A_BP).toFixed(4));
            copy.M_RP_corr = Number((rp - dm - copy.A_RP).toFixed(4));
            copy.bp_rp_corr = Number(((bp - copy.A_BP) - (rp - copy.A_RP)).toFixed(4));
          } else if (Number.isFinite(bprp)) {
            copy.bp_rp_corr = Number((bprp - (copy.A_BP - copy.A_RP)).toFixed(4));
          } else {
            copy.bp_rp_corr = getVal(r, ['bp_rp_corr', 'bp-rp_corr']);
          }
          copy['bp-rp_corr'] = copy.bp_rp_corr;
          // observed plane = where a REDDENED isochrone lives (no extinction removed from the stars)
          copy.M_G_obs = Number.isFinite(g) ? Number((g - dm).toFixed(4)) : NaN;
          copy.bp_rp_obs = (Number.isFinite(bp) && Number.isFinite(rp)) ? Number((bp - rp).toFixed(4)) : (Number.isFinite(bprp) ? bprp : NaN);
          copy.bprp_corr_ok = copy.bp_rp_corr >= STARNEST_BPRP0_RANGE[0] && copy.bp_rp_corr <= STARNEST_BPRP0_RANGE[1];
          assignDerived(copy, copy.M_G_corr, copy.bp_rp_corr);     // star type + mass re-derived with THIS extinction/DM
          return copy;
        }).filter(r => Number.isFinite(r.M_G_corr) && Number.isFinite(r.bp_rp_corr) && Number.isFinite(r.M_G_obs) && Number.isFinite(r.bp_rp_obs));

        if (!enriched.length) {
          throw new Error('No stars left after E(B-V)/photometry NaN removal.');
        }
        if (wf.irsa.state === 'done') wf.irsa.dropped_stars = rows.length - enriched.length;

        const counts = {};
        enriched.forEach(r => { counts[r.Star_Type] = (counts[r.Star_Type] || 0) + 1; });
        const fracOk = enriched.filter(r => r.bprp_corr_ok).length / enriched.length;
        const eList = enriched.map(r => r.E_BV);
        const eMean = eList.reduce((a, v) => a + v, 0) / eList.length;
        const eStd = Math.sqrt(eList.reduce((a, v) => a + (v - eMean) ** 2, 0) / eList.length);
        const avMean = STARNEST_AV_PER_EBV * eMean;
        wf.extinction = { E_BV_mean: eMean, E_BV_std: eStd, A_V_mean: avMean, A_G: STARNEST_AG_AV * avMean, A_BP: STARNEST_ABP_AV * avMean,
                          A_RP: STARNEST_ARP_AV * avMean, AV_PER_EBV: STARNEST_AV_PER_EBV, DM: clusterDm, d_pc: dPc };

        // ---- step 4: PARSEC from CMD 3.8 (direct browser fetch); the MEAN A_V is forwarded ----
        const avSend = handling === 'reddened' ? Number(avMean.toFixed(3)) : 0.0;
        let iso0 = null, isoPlot = null, check = null, used = null;
        try {
          const r0 = await bnetGetParsec(age, mh, 0.0, !allowNet);                    // intrinsic: shift check reference
          iso0 = r0.iso; used = r0;
          if (handling === 'reddened' && avSend > 0.0005) {
            used = await bnetGetParsec(age, mh, avSend, !allowNet);
            check = bnetCheckShift(iso0, used.iso, avSend);
            if (check.ok === false) notes.push('WARNING: ' + check.message);
          }
          isoPlot = used.iso;
          wf.parsec = { state: 'done', rows: isoPlot.n, kb: Number((used.raw.length / 1024).toFixed(1)), logAge: age, MH: mh, A_V_sent: avSend,
                        system: 'Gaia EDR3 = DR3 photometry (Vega)', service: 'CMD 3.8', tls_verified: true, file_url: used.info.source_url, check };
        } catch (e) {
          if (e instanceof NetworkNotAllowed) {
            wf.parsec = { state: 'pending' };
            notes.push('PARSEC not fetched yet: press "Apply Extinction & Fetch PARSEC" (one request per logAge / [M/H] / A_V).');
          } else {
            wf.parsec = { state: 'failed', error: e.message || String(e) };
            notes.push(`PARSEC download failed (${e.message || e}); isochrone not drawn.`);
          }
        }
        const parsecName = `PARSEC_gaiaEDR3_logAge${age.toFixed(2)}_MH${(mh >= 0 ? '+' : '-') + Math.abs(mh).toFixed(2)}_AV${avSend.toFixed(3)}.dat`;
        localStore.parsec = used && isoPlot ? { text: used.raw, name: parsecName } : null;

        // ---- step 5: plot + fit quality ----
        const reddened = handling === 'reddened';
        const cx = reddened ? 'bp_rp_obs' : 'bp_rp_corr', my = reddened ? 'M_G_obs' : 'M_G_corr';
        let isochrone = { color: [], mag: [], rows: [], source: 'none' };
        if (isoPlot) {
          const keep = []; for (let i = 0; i < isoPlot.n; i++) if (!isoPlot.label || isoPlot.label[i] <= BNET_CFG.ISO_MAX_LABEL) keep.push(i);
          const icol = keep.map(i => isoPlot.G_BPmag[i] - isoPlot.G_RPmag[i]), imag = keep.map(i => isoPlot.Gmag[i]);
          isochrone = { color: icol, mag: imag, rows: [], source: 'parsec' };
          const okm = enriched.filter(r => r.bprp_corr_ok && r.Star_Type !== 'White Dwarf');
          wf.fit = { state: 'done', plane: reddened ? 'reddened (observed colour, G-DM)' : 'de-reddened', ...bnetFitQuality(okm.map(r => r[cx]), okm.map(r => r[my]), icol, imag) };
        } else {
          wf.fit = { state: 'pending' };
        }

        const res = {
          empirical: {
            color: enriched.map(r => r[cx]),
            mag: enriched.map(r => r[my]),
            star_type: enriched.map(r => r.Star_Type),
            bprp_corr_ok: enriched.map(r => r.bprp_corr_ok)
          },
          isochrone,
          meta: {
            mode: mode === 'fixed' ? 'fixed' : (wf.irsa.state === 'done' ? mode : `${mode} (IRSA not available: typed E(B-V))`),
            engine: 'browser',
            note: notes.join(' '),
            workflow: wf,
            plot_space: reddened ? 'reddened' : 'intrinsic',
            A_V_sent: avSend,
            n_in: rows.length,
            n_out: enriched.length,
            E_BV_median: medianOf(eList),
            E_BV_mean: eMean,
            E_BV_std: eStd,
            A_V: Number((STARNEST_AV_PER_EBV * medianOf(eList)).toFixed(4)),
            A_V_mean: avMean,
            d_pc: Number(dPc.toFixed(1)),
            DM: Number(clusterDm.toFixed(3)),
            frac_colour_ok: Number(fracOk.toFixed(4))
          },
          classification_counts: counts,
          total_classified: enriched.length
        };
        // members_latest.csv is NOT overwritten: it stays the pure GMM output. Derived table lives in localStore.analysis.
        const analysis = { srcFile: srcName, srcVer: localStore.ver[srcName] || 0, rows: enriched, id: ++localStore.verCounter };
        localStore.analysis = analysis;
        localStore.latestExtinction = enriched;
        localStore.latestIsochrone = [];
        if (wf.irsa.state !== 'failed' && wf.parsec.state !== 'failed') localStore.cache[cacheKey] = { res, analysis, isoRows: [], parsec: localStore.parsec };   // failures are retried on the next Apply
        return res;
      }

      if (endpoint === '/compute-kinematics') {
        const params = JSON.parse(options.body || '{}');
        const base = getAnalysisRows(params.filename);
        const RA_KEYS = ['RA_ICRS', 'RAJ2000', '_RAJ2000', 'ra', 'RA'], DE_KEYS = ['DE_ICRS', 'DEJ2000', '_DEJ2000', 'dec', 'DEC'];
        // Rows lacking position / parallax / proper motion are DROPPED (previously silently given RA=116.18, Dec=-28.39, plx=2.9, pm=0)
        const enriched = base.filter(r =>
          Number.isFinite(getVal(r, RA_KEYS)) && Number.isFinite(getVal(r, DE_KEYS)) &&
          getVal(r, ['Plx', 'parallax', 'plx']) > 0 &&
          Number.isFinite(getVal(r, ['pmRA', 'pmra'])) && Number.isFinite(getVal(r, ['pmDE', 'pmdec', 'pmDEC'])));
        if (!enriched.length) throw new Error('No stars with valid RA/Dec, parallax and proper motion for the kinematics plots.');

        let raList = enriched.map(r => getVal(r, RA_KEYS));
        const [raLo, raHi] = minMax(raList);
        if (raHi - raLo > 180) raList = raList.map(v => v < 180 ? v + 360 : v);     // cluster straddling RA = 0/360
        const decList = enriched.map(r => getVal(r, DE_KEYS));
        const plxList = enriched.map(r => getVal(r, ['Plx', 'parallax', 'plx']));
        const pmraList = enriched.map(r => getVal(r, ['pmRA', 'pmra']));
        const pmdecList = enriched.map(r => getVal(r, ['pmDE', 'pmdec', 'pmDEC']));

        const medPmra = medianOf(pmraList), medPmdec = medianOf(pmdecList);
        const centerRa = medianOf(raList), centerDec = medianOf(decList);

        const colors = [], starTypes = [], uArr = [], vArr = [];
        for (let i = 0; i < enriched.length; i++) {
          const pmraRes = pmraList[i] - medPmra;
          const pmdecRes = pmdecList[i] - medPmdec;
          const vnorm = Math.sqrt(pmraRes * pmraRes + pmdecRes * pmdecRes);
          const denom = vnorm < 1e-8 ? 1.0 : vnorm;
          uArr.push(pmraRes / denom);
          vArr.push(pmdecRes / denom);

          const stype = enriched[i].Star_Type || 'Main Sequence';
          starTypes.push(stype);
          colors.push(STAR_TYPE_COLORS[stype] || 'gray');
        }

        const eplxList = enriched.map(r => getVal(r, ['e_Plx', 'e_plx', 'parallax_error', 'plx_error']));
        const k3 = kin3dCompute(raList, decList, plxList, eplxList, { mode: params.kin3d_mode, maxFrac: Number(params.kin3d_max_frac), depthPc: Number(params.kin3d_depth), frame: params.kin3d_frame, centre: params.kin3d_centre });
        const p3 = { x: k3.x, y: k3.y, z: k3.z, d: k3.d, frac: k3.frac, color: k3.keep.map(i => colors[i]), star_type: k3.keep.map(i => starTypes[i]), meta: k3.meta };

        const rawClusterName = (localStore.activeClusterName || params.filename || 'Cluster')
          .replace(/^Extinctioncorr_/i, '')
          .replace(/newcorrected.*$/i, '')
          .split('_')[0];

        return {
          cluster_name: rawClusterName,
          meta: { n_in: base.length, n_used: enriched.length, n_dropped: base.length - enriched.length },
          positions_3d: p3, meta3d: p3.meta,
          vectors_2d: {
            ra: raList, dec: decList, u: uArr, v: vArr, color: colors, star_type: starTypes,
            center_ra: centerRa, center_dec: centerDec
          }
        };
      }

      if (endpoint === '/compute-imf') {
        const params = JSON.parse(options.body || '{}');
        const aId = (localStore.analysis && localStore.analysis.srcFile === params.filename) ? localStore.analysis.id : 0;
        const cacheKey = `imf_${ options.body }|v${localStore.ver[params.filename] || 0}|a${aId}`;
        if (localStore.cache[cacheKey]) return localStore.cache[cacheKey];

        const enriched = getAnalysisRows(params.filename);
        const mMinSel = Number.isFinite(Number(params.mass_min)) ? Number(params.mass_min) : STARNEST_M_MIN_FLOOR;
        const mMaxSel = Number.isFinite(Number(params.mass_max)) ? Number(params.mass_max) : STARNEST_M_MAX_CAP;

        // Masses from L(M) relations are not valid for white dwarfs / evolved giants -> excluded ("MS only" badge)
        const EVOLVED = new Set(['White Dwarf', 'Red Giant / Supergiant']);
        let nExcludedEvolved = 0;
        const allMasses = [];
        enriched.forEach(r => {
          if (EVOLVED.has(r.Star_Type)) { nExcludedEvolved++; return; }
          const pm = getVal(r, ['MC_Mass', 'Mass', 'mass']);
          if (Number.isFinite(pm) && pm > 0 && pm >= mMinSel && pm <= mMaxSel) allMasses.push(pm);
          if (params.include_companions) {
            const cm = getVal(r, ['MC_CompMass', 'comp_mass', 'Companion_Mass']);
            if (Number.isFinite(cm) && cm > 0 && cm >= mMinSel && cm <= mMaxSel) allMasses.push(cm);
          }
        });

        if (allMasses.length < STARNEST_MIN_MASSES) {
          throw new Error(`Only ${allMasses.length} valid masses; need >= ${STARNEST_MIN_MASSES}.`);
        }

        const nBins = params.bins || Math.min(15, Math.max(8, Math.ceil(Math.log2(allMasses.length) + 1)));
        const fitMinLogM = params.fit_min_logM !== undefined ? params.fit_min_logM : -0.7;
        const [mMinObs, mMaxObs] = minMax(allMasses);
        const mLo = Math.max(STARNEST_M_MIN_FLOOR, mMinObs * 0.9);
        const mHi = Math.min(STARNEST_M_MAX_CAP, mMaxObs * 1.1);
        if (mLo >= mHi) throw new Error('Degenerate mass range.');

        const logMin = Math.log10(mLo);
        const logMax = Math.log10(mHi);
        const binWidth = (logMax - logMin) / nBins;

        const edges = [];
        for (let b = 0; b <= nBins; b++) edges.push(logMin + b * binWidth);

        const counts = new Array(nBins).fill(0);
        const binCenters = [];
        for (let b = 0; b < nBins; b++) {
          binCenters.push(edges[b] + 0.5 * binWidth);
        }

        allMasses.forEach(m => {
          const lm = Math.log10(m);
          let idx = Math.floor((lm - logMin) / binWidth);
          if (idx >= nBins) idx = nBins - 1;
          if (idx >= 0 && idx < nBins) counts[idx]++;
        });

        const fitFlags = counts.map((c, b) => c > 0 && binCenters[b] >= fitMinLogM);
        const nFitBins = fitFlags.filter(Boolean).length;
        if (nFitBins < 3) {
          throw new Error(`Only ${nFitBins} complete bins above logM=${fitMinLogM}; lower the limit or add data.`);
        }

        const binsTable = [];
        const compCenters = [], compLogDN = [], compLogErr = [], compCounts = [];
        const incompCenters = [], incompLogDN = [], incompLogErr = [], incompCounts = [];
        let nInFit = 0;
        let loEdgeFit = Infinity;

        for (let b = 0; b < nBins; b++) {
          const isComp = binCenters[b] >= fitMinLogM;
          if (counts[b] > 0) {
            const dn = counts[b] / binWidth;
            const err = Math.sqrt(counts[b]) / binWidth;
            const ldn = Math.log10(dn);
            const lerr = (err / dn) / Math.LN10;
            binsTable.push({
              logM_center: Number(binCenters[b].toFixed(4)),
              count: counts[b],
              dN_dlogM: Number(dn.toFixed(4)),
              dN_dlogM_err: Number(err.toFixed(4)),
              log_dN: Number(ldn.toFixed(4)),
              log_err: Number(lerr.toFixed(4)),
              complete: isComp
            });
            if (isComp) {
              compCenters.push(binCenters[b]);
              compLogDN.push(ldn);
              compLogErr.push(lerr);
              compCounts.push(counts[b]);
              nInFit += counts[b];
              if (edges[b] < loEdgeFit) loEdgeFit = edges[b];
            } else {
              incompCenters.push(binCenters[b]);
              incompLogDN.push(ldn);
              incompLogErr.push(lerr);
              incompCounts.push(counts[b]);
            }
          } else {
            binsTable.push({
              logM_center: Number(binCenters[b].toFixed(4)),
              count: 0,
              dN_dlogM: NaN,
              dN_dlogM_err: NaN,
              log_dN: NaN,
              log_err: NaN,
              complete: isComp
            });
          }
        }

        const reg = wlinregressJS(compCenters, compLogDN, compLogErr);   // Poisson-weighted slope
        const xFit = [loEdgeFit, Math.max(...compCenters) + 0.5 * binWidth];   // span the whole fitted range
        const yFit = [reg.slope * xFit[0] + reg.intercept, reg.slope * xFit[1] + reg.intercept];

        // Scale Salpeter, Kroupa, Chabrier over complete mass range [loEdgeFit, logMax]
        const N_CURVE = 400;
        const lmFitGrid = [], mFitGrid = [];
        const lmAllGrid = [], mAllGrid = [];
        for (let i = 0; i < N_CURVE; i++) {
          const lf = loEdgeFit + ((logMax - loEdgeFit) * i) / (N_CURVE - 1);
          lmFitGrid.push(lf);
          mFitGrid.push(Math.pow(10, lf));

          const la = logMin + ((logMax - logMin) * i) / (N_CURVE - 1);
          lmAllGrid.push(la);
          mAllGrid.push(Math.pow(10, la));
        }

        const phiSalpeter = m => Math.pow(m, 1.0 - 2.35);
        const phiKroupa = m => {
          const a1 = 0.3, a2 = 1.3, a3 = 2.3;
          const k2 = Math.pow(0.08, a2 - a1);
          const k3 = k2 * Math.pow(0.50, a3 - a2);
          if (m < 0.08) return Math.pow(m, 1 - a1);
          if (m < 0.5) return k2 * Math.pow(m, 1 - a2);
          return k3 * Math.pow(m, 1 - a3);
        };
        const phiChabrier = m => {
          const mc = 0.08, sigma = 0.69, alpha = 2.3;
          if (m <= 1.0) {
            return 0.158 * Math.exp(-0.5 * Math.pow(Math.log10(m) - Math.log10(mc), 2) / (sigma * sigma));
          }
          const aHigh = 0.158 * Math.exp(-0.5 * Math.pow(0.0 - Math.log10(mc), 2) / (sigma * sigma));
          return aHigh * Math.pow(m, 1 - alpha);
        };

        const models = { Salpeter: phiSalpeter, Kroupa: phiKroupa, Chabrier: phiChabrier };
        const scaledCurves = {};
        const chi2Rows = [];

        Object.entries(models).forEach(([name, fn]) => {
          const rawFit = mFitGrid.map(fn);
          const area = trapzJS(rawFit, lmFitGrid);
          const scale = area > 0 ? nInFit / area : 1.0;
          scaledCurves[name] = mAllGrid.map(m => Math.log10(Math.max(1e-12, fn(m) * scale)));

          let chi2Sum = 0;
          for (let i = 0; i < compCenters.length; i++) {
            const obsDn = compCounts[i] / binWidth;
            const obsErr = Math.sqrt(compCounts[i]) / binWidth;
            const modDn = fn(Math.pow(10, compCenters[i])) * scale;
            chi2Sum += Math.pow((obsDn - modDn) / obsErr, 2);
          }
          const dof = nFitBins - 1;
          chi2Rows.push({ Model: name, red_chi2: Number((chi2Sum / dof).toFixed(3)), dof });
        });

        chi2Rows.sort((a, b) => a.red_chi2 - b.red_chi2);
        localStore.latestImfBins = binsTable;
        localStore.latestImfChi2 = chi2Rows;

        const compMags = [], compFrac = [];
        for (let g = 8; g <= 20; g += 0.8) {
          compMags.push(g);
          compFrac.push(1.0 / (1.0 + Math.exp((g - 17.2) / 0.65)));
        }

        const logAgeVal = Number.isFinite(Number(params.log_age)) ? Number(params.log_age) : 7.5;
        const mfKind = 'Initial Mass Function (IMF)';

        const result = {
          imf: {
            mf_kind: mfKind,
            meta: {
              n_masses: allMasses.length,
              n_bins: nBins,
              bin_width_dex: Number(binWidth.toFixed(4)),
              n_fit_bins: nFitBins,
              n_in_fit: nInFit,
              fit_min_logm: fitMinLogM,
              n_excluded_evolved: nExcludedEvolved,
              m_range: [Number(mLo.toFixed(3)), Number(mHi.toFixed(2))]
            },
            bins_table: binsTable,
            chi2_table: chi2Rows,
            complete_bins: {
              logM_center: compCenters,
              log_dN: compLogDN,
              log_err: compLogErr,
              count: compCounts
            },
            incomplete_bins: {
              logM_center: incompCenters,
              log_dN: incompLogDN,
              log_err: incompLogErr,
              count: incompCounts
            },
            slope: reg.slope,
            slope_err: reg.se,
            r_squared: reg.r2,
            ms_count: allMasses.length,
            fit_line: { x: xFit, y: yFit, slope: reg.slope, slope_err: reg.se, n_bins: nFitBins },
            theoretical: {
              log_m: lmAllGrid,
              salpeter: scaledCurves.Salpeter,
              kroupa: scaledCurves.Kroupa,
              chabrier: scaledCurves.Chabrier
            }
          },
          completeness: {
            mag: compMags,
            fraction: compFrac,
            assumed: true            // fixed logistic (G50 = 17.2): an assumption, NOT measured from the data
          }
        };
        localStore.cache[cacheKey] = result;
        pruneCache();
        return result;
      }

      throw new Error('Unhandled endpoint: ' + endpoint);
    }

    let backendAlive = null;
    let backendProbe = null;

    function probeBackend(force = false) {
      if (!force && backendAlive !== null) return Promise.resolve(backendAlive);
      if (backendProbe && !force) return backendProbe;
      backendProbe = (async () => {
        try {
          if (window.STARNEST_NO_BACKEND) throw new Error('backend disabled by the host page');     // e.g. embedded in Streamlit
          const c = new AbortController();
          const t = setTimeout(() => c.abort(), API_PROBE_TIMEOUT_MS);
          const r = await fetch(API_BASE + '/health', { signal: c.signal });
          clearTimeout(t);
          backendAlive = r.ok;
        } catch (_) { backendAlive = false; }
        localStore.engine = backendAlive ? 'python' : 'browser';
        if (typeof updateEngineBadge === 'function') updateEngineBadge();
        backendProbe = null;
        return backendAlive;
      })();
      return backendProbe;
    }

    async function apiCall(endpoint, options = {}) {
      if (await probeBackend()) {
        const controller = new AbortController();
        const timeoutMs = options.timeoutMs || API_TIMEOUT_MS;
        const timer = setTimeout(() => controller.abort(), timeoutMs);
        try {
          const { timeoutMs: _t, ...fetchOpts } = options;
          const response = await fetch(API_BASE + endpoint, { ...fetchOpts, signal: controller.signal });
          clearTimeout(timer);
          if (!response.ok) {
            const error = await response.json().catch(() => ({ detail: 'Server error' }));
            const d = error.detail;
            throw new Error(typeof d === 'string' ? d : (d ? JSON.stringify(d) : `Request failed (${response.status})`));
          }
          return await response.json();
        } catch (err) {
          clearTimeout(timer);
          if (err.name === 'AbortError') throw new Error(`Backend did not answer within ${timeoutMs / 1000}s (${endpoint}).`);
          if (err instanceof TypeError) {            // network failure only: server went away -> browser engine
            backendAlive = false; localStore.engine = 'browser'; updateEngineBadge();
            showToast('Python backend unreachable: switched to the in-browser (approximate) engine.', 'warning');
            return await localBackendEngine(endpoint, options);
          }
          throw err;                                  // real server-side errors (bad input, IRSA down...) are shown, not hidden
        }
      }
      return await localBackendEngine(endpoint, options);
    }

    function updateEngineBadge() {
      let el = document.getElementById('engine-badge');
      if (!el) {
        el = document.createElement('div');
        el.id = 'engine-badge';
        el.style.cssText = 'position:fixed;left:14px;bottom:12px;z-index:150;padding:5px 11px;border-radius:999px;font:600 11px Inter,system-ui,sans-serif;backdrop-filter:blur(6px);';
        document.body.appendChild(el);
      }
      const py = localStore.engine === 'python';
      el.textContent = py ? 'Engine: Python backend (PARSEC + IRSA)' : 'Engine: in-browser (GMM / IMF only; IRSA + PARSEC need the Python backend)';
      el.title = py ? 'Real pipeline modules are running on localhost:8000'
                    : 'Start backend.py (uvicorn) for PARSEC isochrones, IRSA extinction and isochrone-based masses';
      el.style.background = py ? 'rgba(34,197,94,0.16)' : 'rgba(245,158,11,0.16)';
      el.style.border = py ? '1px solid rgba(34,197,94,0.4)' : '1px solid rgba(245,158,11,0.45)';
      el.style.color = py ? '#86efac' : '#fcd34d';
    }

    let toastCount = 0;

    function createToastContainer() {
      const container = document.createElement('div');
      container.id = 'toast-container';
      container.className = 'fixed top-20 right-6 z-[200] flex flex-col items-end gap-3';
      document.body.appendChild(container);
      return container;
    }

    function showToast(message, type = 'info') {
      const container = document.getElementById('toast-container') || createToastContainer();
      const toast = document.createElement('div');
      toast.id = `toast-${toastCount++}`;

      const config = {
        info:    { class: 'toast-info',    icon: 'fa-info-circle',          iconColor: 'text-cyan-400' },
        success: { class: 'toast-success', icon: 'fa-check-circle',         iconColor: 'text-green-400' },
        error:   { class: 'toast-error',   icon: 'fa-exclamation-circle',   iconColor: 'text-red-400' },
        warning: { class: 'toast-warning', icon: 'fa-exclamation-triangle', iconColor: 'text-yellow-400' }
      }[type] || { class: 'toast-info', icon: 'fa-info-circle', iconColor: 'text-cyan-400' };

      toast.className = `toast ${config.class} flex items-center gap-3 px-4 py-3 rounded-xl text-sm max-w-sm`;
      toast.innerHTML = `
        <i class="fas ${config.icon} ${config.iconColor} text-base flex-shrink-0"></i>
        <span class="text-gray-200 leading-snug">${escapeHTML(message)}</span>
        <button onclick="dismissToast(this.parentElement)" class="ml-auto text-gray-500 hover:text-white transition-colors flex-shrink-0">
          <i class="fas fa-times text-xs"></i>
        </button>
      `;

      container.appendChild(toast);
      gsap.fromTo(toast, { x: 80, opacity: 0, scale: 0.95 }, { x: 0, opacity: 1, scale: 1, duration: 0.4, ease: 'back.out(1.4)' });
      toast._dismissTimer = setTimeout(() => dismissToast(toast), 4000);
    }

    function dismissToast(toast) {
      if (!toast || !toast.parentNode) return;
      clearTimeout(toast._dismissTimer);
      gsap.to(toast, {
        x: 100, opacity: 0, scale: 0.95, duration: 0.3, ease: 'power2.in',
        onComplete: () => { if (toast.parentNode) toast.parentNode.removeChild(toast); }
      });
    }

    function setButtonLoading(btnId, loading) {
      const btn = document.getElementById(btnId);
      if (!btn) return;
      if (loading) {
        btn.disabled = true;
        btn.dataset.originalText = btn.innerHTML;
        btn.innerHTML = `
          <span class="flex items-center justify-center gap-2">
            <svg class="animate-spin h-4 w-4" viewBox="0 0 24 24" fill="none">
              <circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="3"></circle>
              <path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z"></path>
            </svg>
            Processing...
          </span>
        `;
      } else {
        btn.disabled = false;
        if (btn.dataset.originalText) btn.innerHTML = btn.dataset.originalText;
      }
    }

    function setupDropzone(dropzoneId, fileInputId, fileNameDisplayId) {
      const dropzone = document.getElementById(dropzoneId);
      const fileInput = document.getElementById(fileInputId);
      const fileNameEl = document.getElementById(fileNameDisplayId);
      if (!dropzone || !fileInput) return null;

      dropzone.addEventListener('click', (e) => {
        if (e.target === fileInput) return;
        fileInput.click();
      });
      dropzone.addEventListener('dragover', (e) => {
        e.preventDefault();
        e.stopPropagation();
        dropzone.classList.add('dropzone-dragover');
      });
      dropzone.addEventListener('dragleave', (e) => {
        e.preventDefault();
        e.stopPropagation();
        dropzone.classList.remove('dropzone-dragover');
      });
      dropzone.addEventListener('drop', (e) => {
        e.preventDefault();
        e.stopPropagation();
        dropzone.classList.remove('dropzone-dragover');
        const file = e.dataTransfer.files[0];
        if (file) {
          handleFileSelected(file, fileInput, fileNameEl);
          const changeEvent = new Event('change', { bubbles: true });
          fileInput._selectedFile = file;
          fileInput.dispatchEvent(changeEvent);
        }
      });
      fileInput.addEventListener('change', (e) => {
        if (e.target.files && e.target.files[0]) {
          handleFileSelected(e.target.files[0], fileInput, fileNameEl);
        }
      });
      return { dropzone, fileInput };
    }

    function handleFileSelected(file, fileInput, fileNameEl) {
      fileInput._selectedFile = file;
      if (fileNameEl) {
        const sizeKB = (file.size / 1024).toFixed(1);
        fileNameEl.innerHTML = `<i class="fas fa-file-alt text-accent-cyan/60"></i> ${escapeHTML(file.name)} <span class="text-gray-500">(${sizeKB} KB)</span>`;
        fileNameEl.classList.remove('hidden');
        gsap.fromTo(fileNameEl, { y: -4, opacity: 0 }, { y: 0, opacity: 1, duration: 0.3, ease: 'power2.out' });
      }
    }

    const PLOTLY_LAYOUT = {
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(255, 255, 255, 0.015)',
      font: { color: '#9ca3af', family: 'Inter, system-ui, sans-serif', size: 11 },
      xaxis: {
        gridcolor: 'rgba(255,255,255,0.05)',
        zerolinecolor: 'rgba(255,255,255,0.08)',
        linecolor: 'rgba(255,255,255,0.08)'
      },
      yaxis: {
        gridcolor: 'rgba(255,255,255,0.05)',
        zerolinecolor: 'rgba(255,255,255,0.08)',
        linecolor: 'rgba(255,255,255,0.08)'
      },
      margin: { l: 55, r: 25, t: 40, b: 50 },
      hoverlabel: {
        bgcolor: 'rgba(10, 14, 26, 0.95)',
        bordercolor: 'rgba(0, 212, 255, 0.3)',
        font: { color: '#e5e7eb', size: 11, family: 'Inter' }
      }
    };

    const PLOTLY_CONFIG = {
      responsive: true,
      displayModeBar: true,
      modeBarButtonsToRemove: ['sendDataToCloud', 'autoScale2d'],
      displaylogo: false
    };

    function getPlotlyLayout(overrides = {}) {
      const merged = JSON.parse(JSON.stringify(PLOTLY_LAYOUT));
      for (const key in overrides) {
        if (typeof overrides[key] === 'object' && !Array.isArray(overrides[key]) && merged[key]) {
          merged[key] = { ...merged[key], ...overrides[key] };
        } else {
          merged[key] = overrides[key];
        }
      }
      return merged;
    }

    function moveNavPill(targetBtn) {
      const pill = document.getElementById('nav-pill');
      const container = document.getElementById('nav-links');
      if (!pill || !container || !targetBtn) return;
      const containerRect = container.getBoundingClientRect();
      const btnRect = targetBtn.getBoundingClientRect();
      gsap.to(pill, {
        left: btnRect.left - containerRect.left,
        width: btnRect.width,
        duration: 0.5,
        ease: 'power3.out'
      });
      container.querySelectorAll('.nav-btn').forEach(btn => {
        gsap.to(btn, { color: btn === targetBtn ? '#ffffff' : '#9ca3af', duration: 0.3 });
      });
    }

    function navigateToPage(targetPageName) {
      if (targetPageName === currentPage || isTransitioning) return;
      const currentIndex = PAGE_ORDER.indexOf(currentPage);
      const targetIndex = PAGE_ORDER.indexOf(targetPageName);
      if (targetIndex === -1) return;

      isTransitioning = true;
      const direction = targetIndex > currentIndex ? 1 : -1;
      const currentEl = document.getElementById('page-' + currentPage);
      const targetEl = document.getElementById('page-' + targetPageName);
      if (!currentEl || !targetEl) { isTransitioning = false; return; }

      window.scrollTo(0, 0); document.body.scrollTop = 0;   // every page starts at the top (document scroll, or <body> scroll when embedded)
      targetEl.classList.add('is-entering');        // out of flow, full height, while it slides in
      gsap.set(targetEl, { visibility: 'visible', opacity: 1, x: direction * window.innerWidth });

      const tl = gsap.timeline({
        onComplete: () => {
          gsap.set(currentEl, { visibility: 'hidden', x: 0, opacity: 0 });
          currentEl.classList.remove('is-active');
          targetEl.classList.remove('is-entering');
          targetEl.classList.add('is-active');      // the visible page now defines the document height
          currentPage = targetPageName;
          isTransitioning = false;
          window.dispatchEvent(new Event('resize'));
        }
      });

      tl.to(currentEl, {
        x: -direction * window.innerWidth * 0.3,
        opacity: 0,
        duration: 0.55,
        ease: 'power3.inOut'
      }, 0).to(targetEl, {
        x: 0,
        duration: 0.55,
        ease: 'power3.inOut'
      }, 0);

      const targetBtn = document.querySelector(`.nav-btn[data-page="${targetPageName}"]`);
      if (targetBtn) moveNavPill(targetBtn);
    }

    // Workflow tiles: liquid fill starts where the pointer enters; tap/Enter toggles the details on touch + keyboard
    function initWorkflowTiles() {
      document.querySelectorAll('.wf-tile').forEach(tile => {
        const setOrigin = (e) => {
          const r = tile.getBoundingClientRect();
          tile.style.setProperty('--mx', (e.clientX - r.left) + 'px');
          tile.style.setProperty('--my', (e.clientY - r.top) + 'px');
        };
        tile.addEventListener('pointerenter', setOrigin);
        tile.addEventListener('pointerdown', setOrigin);
        tile.addEventListener('click', (e) => {
          if (e.target.closest('.wf-open')) return;                       // the "Open" link navigates instead
          if (window.matchMedia('(hover: hover)').matches) return;        // desktop: pure hover
          const open = !tile.classList.contains('is-open');
          document.querySelectorAll('.wf-tile.is-open').forEach(t => t.classList.remove('is-open'));
          tile.classList.toggle('is-open', open);
        });
        tile.addEventListener('keydown', (e) => {
          if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); tile.classList.toggle('is-open'); }
        });
        tile.addEventListener('focus', () => { if (tile.matches(':focus-visible')) tile.classList.add('is-open'); });   // keyboard focus only (taps are handled by click)
        tile.addEventListener('blur', () => tile.classList.remove('is-open'));
      });
    }

    function animateHomeEntrance() {
      gsap.fromTo('#page-home .glass-card',
        { y: 20, opacity: 0 },
        { y: 0, opacity: 1, duration: 0.5, stagger: 0.06, ease: 'power2.out', clearProps: 'transform' }   // clearProps: let CSS hover transforms (team tiles) work afterwards
      );
    }

    function initMembershipPlots() {
      const emptyTrace = [{ x: [], y: [], type: 'scatter', mode: 'markers' }];
      Plotly.newPlot('plot-ra-dec', emptyTrace, getPlotlyLayout({
        title: { text: 'Spatial Distribution (RA vs Dec)', font: { size: 13, color: '#9ca3af' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'RA (deg)', autorange: 'reversed' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'DEC (deg)' }
      }), PLOTLY_CONFIG);

      Plotly.newPlot('plot-vpd', emptyTrace, getPlotlyLayout({
        title: { text: 'PMRA vs PMDEC', font: { size: 13, color: '#9ca3af' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'pmRA (mas/yr)' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'pmDE (mas/yr)' }
      }), PLOTLY_CONFIG);

      Plotly.newPlot('plot-gmm-plx', emptyTrace, getPlotlyLayout({
        title: { text: 'Parallax Distribution', font: { size: 12, color: '#9ca3af' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'Parallax (mas)' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'Probability density' }
      }), PLOTLY_CONFIG);

      Plotly.newPlot('plot-gmm-cmd-mem', emptyTrace, getPlotlyLayout({
        title: { text: 'Member CMD', font: { size: 12, color: '#9ca3af' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'BP - RP (mag)' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'G (mag)', autorange: 'reversed' }
      }), PLOTLY_CONFIG);

      Plotly.newPlot('plot-gmm-cmd-nonmem', emptyTrace, getPlotlyLayout({
        title: { text: 'Non-member CMD', font: { size: 12, color: '#9ca3af' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'BP - RP (mag)' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'G (mag)', autorange: 'reversed' }
      }), PLOTLY_CONFIG);
    }

    async function runGMMAnalysis() {
      if (!window.rawFileName) {
        showToast('Please upload a cluster CSV file first', 'warning');
        return;
      }
      showModuleError('gmm-error-box', 'gmm-error-text', null);
      setButtonLoading('btn-run-gmm', true);

      const kVal = document.getElementById('input-gmm-k').value;
      const params = {
        filename: window.rawFileName,
        ebv_ui: uiEbv(),
        ruwe_max: parseFloat(document.getElementById('input-ruwe').value),
        gofal_max: parseFloat(document.getElementById('input-gofal').value),
        epmra_max: parseFloat(document.getElementById('input-epmra').value),
        epmde_max: parseFloat(document.getElementById('input-epmde').value),
        plx_err_min: parseFloat(document.getElementById('input-plx-err').value),
        nper_min: parseFloat(document.getElementById('input-nper').value),
        distance_lit: parseFloat(document.getElementById('input-distance').value),
        half_width: parseFloat(document.getElementById('input-half-width').value),
        cutoff: parseFloat(document.getElementById('input-cutoff').value),
        k: kVal === 'auto' ? null : parseInt(kVal, 10),
        gmm_opt: document.getElementById('input-gmm-opt').value,
        plx_zp: parseFloat(document.getElementById('input-gmm-plxzp').value) || 0.0
      };

      try {
        const data = await apiCall('/run-gmm', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(params)
        });

        updateGMMPlots(data);
        updateGMMStats(data.stats);
        updateGMMClassificationCounts(data.classification_counts);

        const dlSection = document.getElementById('gmm-download-section');
        if (dlSection) {
          dlSection.classList.remove('border-white/[0.08]', 'bg-white/[0.02]');
          dlSection.classList.add('border-green-500/30', 'bg-green-500/[0.05]');
          document.getElementById('gmm-download-icon').className = 'fas fa-check-circle text-green-400';
          document.getElementById('gmm-download-title').textContent = `Filtered & Classified File Ready (${data.stats.members} stars)`;
          const dlBadge = document.getElementById('gmm-download-badge');
          dlBadge.className = 'text-[10px] bg-green-500/20 text-green-300 px-2 py-0.5 rounded-full border border-green-500/30';
          dlBadge.textContent = `${data.stats.members} Members`;
          document.getElementById('gmm-download-desc').textContent = `Contains ${data.stats.members} filtered cluster members (K=${data.stats.k} via BIC). Ready for download and active across Isochrone, Kinematics, and IMF!`;
          const dlBtn = document.getElementById('btn-download-members');
          dlBtn.disabled = false;
          dlBtn.className = 'w-full py-2.5 bg-gradient-to-r from-green-600 to-emerald-600 hover:from-green-500 hover:to-emerald-500 border border-green-500/30 rounded-xl text-sm font-semibold transition-all duration-300 flex items-center justify-center gap-2 text-white shadow-[0_0_15px_rgba(34,197,94,0.2)] cursor-pointer';
        }

        window.gmmFilteredFileName = 'members_latest.csv';   // raw upload name (rawFileName) is left untouched
        showToast(`StarNest GMM complete (K=${data.stats.k})! Extracted ${data.stats.members} cluster members.`, 'success');

        await computeCMDPlot();
        await computeKinematicsPlots();
        await computeIMFPlots();
      } catch (err) {
        showModuleError('gmm-error-box', 'gmm-error-text', err.message || String(err));
        showToast(err.message || 'GMM error', 'error');
      } finally {
        setButtonLoading('btn-run-gmm', false);
      }
    }

    function updateGMMPlots(data) {
      const raDecTraces = [];
      if (data.ra_dec.field.ra && data.ra_dec.field.ra.length > 0) {
        raDecTraces.push({
          x: data.ra_dec.field.ra,
          y: data.ra_dec.field.dec,
          mode: 'markers',
          type: data.ra_dec.field.ra.length > 4000 ? 'scattergl' : 'scatter',
          name: `Field Stars (${data.ra_dec.field.ra.length})`,
          marker: { color: '#ef4444', size: 3.5, opacity: 0.35 }
        });
      }

      // Group members by Star_Type so legend matches exact classification script
      const starTypes = data.ra_dec.members.star_type || [];
      const uniqueTypes = [...new Set(starTypes)];
      if (uniqueTypes.length > 0) {
        uniqueTypes.forEach(stype => {
          const idxList = starTypes.map((t, i) => t === stype ? i : -1).filter(i => i !== -1);
          raDecTraces.push({
            x: idxList.map(i => data.ra_dec.members.ra[i]),
            y: idxList.map(i => data.ra_dec.members.dec[i]),
            mode: 'markers',
            type: 'scatter',
            name: `${stype} (${idxList.length})`,
            marker: {
              color: STAR_TYPE_COLORS[stype] || '#00FFFF',
              size: 6,
              opacity: 0.9,
              line: { color: stype === 'White Dwarf' ? '#ffffff' : '#000000', width: stype === 'White Dwarf' ? 1 : 0.4 }
            }
          });
        });
      } else {
        raDecTraces.push({
          x: data.ra_dec.members.ra,
          y: data.ra_dec.members.dec,
          mode: 'markers',
          type: 'scatter',
          name: 'Members',
          marker: { color: '#00FFFF', size: 5 }
        });
      }

      // Plot Cluster Center star marker (matching plot_spatial in user's code)
      if (data.ra_dec.members.ra.length > 0) {
        let centerCoords = null;
        const fnameLower = (localStore.activeClusterName || '').toLowerCase().replace(/\s+/g, '');
        for (const [k, coords] of Object.entries(CLUSTER_CENTERS)) {
          if (fnameLower.includes(k)) { centerCoords = coords; break; }
        }
        if (!centerCoords) {
          const avgRa = data.ra_dec.members.ra.reduce((a, b) => a + b, 0) / data.ra_dec.members.ra.length;
          const avgDec = data.ra_dec.members.dec.reduce((a, b) => a + b, 0) / data.ra_dec.members.dec.length;
          centerCoords = [avgRa, avgDec];
        }
        raDecTraces.push({
          x: [centerCoords[0]],
          y: [centerCoords[1]],
          mode: 'markers',
          type: 'scatter',
          name: 'Cluster Center',
          marker: { symbol: 'star', size: 14, color: 'gold', line: { color: '#000', width: 1 } }
        });
      }

      Plotly.react('plot-ra-dec', raDecTraces, getPlotlyLayout({
        title: { text: 'Spatial Distribution (RA vs Dec)', font: { size: 13, color: '#e5e7eb' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'RA (deg)', autorange: 'reversed' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'DEC (deg)', scaleanchor: 'x', scaleratio: 1 / Math.max(0.1, Math.cos(((data.ra_dec.members.dec.length ? medianOf(data.ra_dec.members.dec) : 0) * Math.PI) / 180)) },
        legend: { font: { color: '#d1d5db', size: 10 }, bgcolor: 'rgba(10,14,26,0.7)' }
      }), PLOTLY_CONFIG);

      const vpdTraces = [
        {
          x: data.vpd.field.pmra,
          y: data.vpd.field.pmdec,
          mode: 'markers',
          type: data.vpd.field.pmra.length > 4000 ? 'scattergl' : 'scatter',
          name: `Non-member (${data.vpd.field.pmra.length})`,
          marker: { color: '#3b82f6', size: 4.5, opacity: 0.45 }
        },
        {
          x: data.vpd.members.pmra,
          y: data.vpd.members.pmdec,
          mode: 'markers',
          type: 'scatter',
          name: `Member (${data.vpd.members.pmra.length})`,
          marker: {
            color: '#ef4444',
            size: 6,
            opacity: 0.9,
            line: { color: '#ffffff', width: 0.4 }
          }
        }
      ];

      Plotly.react('plot-vpd', vpdTraces, getPlotlyLayout({
        title: { text: 'Proper Motion Distribution: pmRA vs pmDE', font: { size: 14, color: '#e5e7eb' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'pmRA (mas/yr)' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'pmDE (mas/yr)' },
        legend: { font: { color: '#d1d5db', size: 11 }, bgcolor: 'rgba(10,14,26,0.75)' }
      }), PLOTLY_CONFIG);

      if (data.plx_Cmd) {
        const allPlx = data.plx_Cmd.field_plx.concat(data.plx_Cmd.members_plx).filter(Number.isFinite);
        const [pxLo, pxHi] = allPlx.length ? minMax(allPlx) : [0, 1];
        const sharedBins = { start: pxLo, end: pxHi, size: (pxHi - pxLo) / 40 || 0.05 };   // same bins for both histograms
        Plotly.react('plot-gmm-plx', [
          {
            x: data.plx_Cmd.field_plx,
            type: 'histogram',
            histnorm: 'probability density',
            autobinx: false, xbins: sharedBins,
            name: 'Non-member',
            opacity: 0.55,
            marker: { color: '#3b82f6' }
          },
          {
            x: data.plx_Cmd.members_plx,
            type: 'histogram',
            histnorm: 'probability density',
            autobinx: false, xbins: sharedBins,
            name: 'Member',
            opacity: 0.65,
            marker: { color: '#ef4444' }
          }
        ], getPlotlyLayout({
          barmode: 'overlay',
          title: { text: 'Parallax Distribution', font: { size: 12, color: '#e5e7eb' } },
          xaxis: { ...getPlotlyLayout().xaxis, title: 'Parallax (mas)' },
          yaxis: { ...getPlotlyLayout().yaxis, title: 'Probability density' },
          legend: { font: { color: '#d1d5db', size: 10 }, bgcolor: 'rgba(0,0,0,0)' }
        }), PLOTLY_CONFIG);

        Plotly.react('plot-gmm-cmd-mem', [{
          x: data.plx_Cmd.members_color,
          y: data.plx_Cmd.members_gmag,
          mode: 'markers',
          type: 'scatter',
          name: 'Member',
          marker: { color: '#ef4444', size: 4.5, opacity: 0.85 }
        }], getPlotlyLayout({
          title: { text: 'Member CMD', font: { size: 12, color: '#e5e7eb' } },
          xaxis: { ...getPlotlyLayout().xaxis, title: 'BP - RP (mag)' },
          yaxis: { ...getPlotlyLayout().yaxis, title: 'G (mag)', autorange: 'reversed' }
        }), PLOTLY_CONFIG);

        Plotly.react('plot-gmm-cmd-nonmem', [{
          x: data.plx_Cmd.field_color,
          y: data.plx_Cmd.field_gmag,
          mode: 'markers',
          type: data.plx_Cmd.field_color.length > 4000 ? 'scattergl' : 'scatter',
          name: 'Non-member',
          marker: { color: '#3b82f6', size: 4, opacity: 0.65 }
        }], getPlotlyLayout({
          title: { text: 'Non-member CMD', font: { size: 12, color: '#e5e7eb' } },
          xaxis: { ...getPlotlyLayout().xaxis, title: 'BP - RP (mag)' },
          yaxis: { ...getPlotlyLayout().yaxis, title: 'G (mag)', autorange: 'reversed' }
        }), PLOTLY_CONFIG);
      }
    }

    function renderGMMComponents(stats) {
      const box = document.getElementById('gmm-components-box'), body = document.getElementById('gmm-components-body');
      if (!box || !body) return;
      const comps = stats.components || [];
      box.classList.toggle('hidden', !comps.length);
      const f = (v, d = 3) => Number.isFinite(v) ? v.toFixed(d) : '-';
      body.innerHTML = comps.map(c => `<tr class="${c.is_cluster ? 'text-green-300 font-semibold' : ''}"><td class="pr-3">${c.k}</td><td class="pr-3">${f(c.weight)}</td><td class="pr-3">${c.n_hard}</td><td class="pr-3">${c.det.toExponential(2)}</td><td class="pr-3">${f(c.pmRA, 2)} &plusmn; ${f(c.sd_pmRA, 2)}</td><td class="pr-3">${f(c.pmDE, 2)} &plusmn; ${f(c.sd_pmDE, 2)}</td><td class="pr-3">${f(c.Plx, 3)} &plusmn; ${f(c.sd_Plx, 3)}</td><td>${c.is_cluster ? '&larr; cluster' : ''}</td></tr>`).join('');
      const bic = stats.bic || {};
      const kb = document.getElementById('gmm-ktable-body');
      if (kb) kb.innerHTML = (stats.k_table || []).map(r => `<tr class="${r.k === stats.k ? 'text-green-300 font-semibold' : ''}${r.degenerate ? ' text-gray-600' : ''}"><td class="pr-3">${r.k}${r.k === stats.k ? '*' : ''}</td><td class="pr-3">${r.bic === null || r.bic === undefined ? 'n/a' : Math.round(r.bic * 10) / 10}</td><td class="pr-3">${r.members ?? '-'}</td><td class="pr-3">${r.non_members ?? '-'}</td><td class="pr-3">${r.min_comp_stars ?? '-'}</td><td>${r.degenerate ? 'collapsed component, excluded' : (r.converged === false ? 'not converged' : '')}</td></tr>`).join('');
      const ml = document.getElementById('gmm-methods-line'), mt = stats.methods;
      if (ml) ml.textContent = mt ? `Methods: ${mt.engine}; K ${mt.k_mode}${stats.k !== undefined ? ' -> K=' + stats.k : ''}; EM ${mt.optimizer} (tol ${mt.tol}, max_iter ${mt.max_iter}, n_init ${mt.n_init}, reg_covar ${mt.reg_covar}); features (pmRA, pmDE, Plx) standardised (x-median)/sd; window ${mt.distance_lit} +/- ${mt.half_width} pc; cutoff ${mt.cutoff}; Plx ZP ${mt.plx_zp} mas.` : '';
      const line = document.getElementById('gmm-bic-line');
      if (line) line.textContent = 'BIC: ' + Object.keys(bic).map(k => `K=${k}: ${bic[k] === null || bic[k] === undefined ? 'n/a' : Math.round(bic[k])}${Number(k) === stats.k ? '*' : ''}`).join('  ');
    }

    function updateGMMStats(stats) {
      document.getElementById('gmm-stats').classList.remove('hidden');
      document.getElementById('stat-input').textContent = stats.n_input !== undefined ? stats.n_input : '-';
      document.getElementById('stat-total').textContent = stats.total;
      document.getElementById('stat-k').textContent = stats.k !== undefined ? stats.k : '-';
      document.getElementById('stat-members').textContent = stats.members;
      document.getElementById('stat-field').textContent = stats.field;
      document.getElementById('stat-threshold').textContent = stats.threshold.toFixed(2);
      const mssEl = document.getElementById('stat-mss');
      if (mssEl) {
        mssEl.textContent = stats.mss_metric !== undefined ? stats.mss_metric.toFixed(3) : '-';
      }
      renderGMMComponents(stats);
      const distEl = document.getElementById('stat-mem-dist');
      if (distEl) {
        distEl.textContent = stats.distance_from_members_pc ? `${stats.distance_from_members_pc.toFixed(1)} pc` : '-';
      }
    }

    function updateGMMClassificationCounts(counts) {
      const container = document.getElementById('gmm-classification-stats');
      const list = document.getElementById('gmm-classification-list');
      if (!container || !list || !counts) return;

      list.innerHTML = '';
      let hasItems = false;
      Object.entries(STAR_TYPE_COLORS).forEach(([stype, color]) => {
        const count = counts[stype] || 0;
        if (count > 0) {
          hasItems = true;
          const borderStyle = stype === 'White Dwarf' ? 'border: 1px solid #fff;' : '';
          list.innerHTML += `
            <div class="flex items-center justify-between text-xs py-1 border-b border-white/[0.04] last:border-0">
              <span class="flex items-center gap-2 text-gray-300">
                <span class="w-2.5 h-2.5 rounded-full inline-block" style="background-color: ${color}; ${borderStyle}"></span>
                ${stype}
              </span>
              <span class="font-mono font-semibold text-white">${count}</span>
            </div>
          `;
        }
      });
      container.classList.toggle('hidden', !hasItems);
    }

    function csvCell(v) {
      const t = v !== undefined && v !== null ? String(v) : '';
      return /[",\n\r]/.test(t) ? `"${t.replace(/"/g, '""')}"` : t;
    }

    function downloadCSVRows(rows, downloadName = 'filtered_classified_members.csv', serverPath = '/download-members') {
      if (!rows || !rows.length) {
        if (localStore.engine === 'python') { apiOpen(serverPath); return; }
        showToast('Nothing to download yet: run the analysis first.', 'warning');
        return;
      }
      const headers = Object.keys(rows[0]).filter(k => !k.startsWith('_'));
      const csvLines = [headers.map(csvCell).join(',')];
      rows.forEach(r => csvLines.push(headers.map(h => csvCell(r[h])).join(',')));
      const blob = new Blob([csvLines.join('\n')], { type: 'text/csv;charset=utf-8;' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = downloadName;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
      showToast(`Downloaded ${downloadName} (${rows.length} rows)`, 'success');
    }

    function initCMDPlot() {
      Plotly.newPlot('plot-cmd-iso', [{ x: [], y: [], type: 'scatter', mode: 'markers' }], getPlotlyLayout({
        title: { text: 'Color-Magnitude Diagram (CMD)', font: { size: 13, color: '#9ca3af' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: '(BP - RP)₀ (mag)' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'M_G (mag)', autorange: 'reversed' }
      }), PLOTLY_CONFIG);
    }

    // Five-step workflow tracker (filtered CSV -> IRSA -> extinction -> PARSEC -> fit). Data comes from meta.workflow (Python backend).
    function renderIsoWorkflow(meta, params) {
      const set = (n, state, html) => {
        const box = document.getElementById('iso-wf-' + n);
        if (!box) return;
        const colors = { done: 'bg-green-400', pending: 'bg-gray-500', failed: 'bg-red-400', warn: 'bg-amber-400' };
        box.querySelector('.iso-wf-dot').className = 'iso-wf-dot inline-block w-1.5 h-1.5 rounded-full mr-1 ' + (colors[state] || colors.pending);
        box.querySelector('.iso-wf-body').innerHTML = html;
      };
      const f = (v, d = 3) => Number.isFinite(Number(v)) && v !== null ? Number(v).toFixed(d) : '-';
      const wf = meta && meta.workflow;
      const nStars = wf ? wf.n_members : (meta ? meta.n_in : '-');
      set(1, 'done', `${nStars} stars<br>${escapeHTML(getTaskDisplayLabel('isochrone'))}`);
      if (!wf) {                                           // in-browser engine: cannot reach IRSA / PARSEC
        set(2, 'pending', 'not queried yet<br>press Apply');
        set(3, 'done', `typed E(B-V)<br>A_V ${f(meta && meta.A_V)}`);
        set(4, 'pending', 'not fetched yet<br>press Apply');
        set(5, 'pending', '-');
        return;
      }
      const ir = wf.irsa || {}, ex = wf.extinction || {}, pa = wf.parsec || {}, ft = wf.fit || {};
      if (ir.state === 'done') set(2, ir.failed ? 'warn' : 'done', `E(B-V)<sub>SFD</sub> mean ${f(ex.E_BV_mean, 4)} &plusmn; ${f(ex.E_BV_std, 4)}<br>${ir.queries} IRSA call(s)${ir.failed ? `, ${ir.failed} failed (${ir.dropped_stars} stars dropped)` : ''}`);
      else if (ir.state === 'manual') set(2, 'done', `typed E(B-V) ${f(ex.E_BV_mean, 4)}<br>IRSA not used (fixed)`);
      else if (ir.state === 'failed') set(2, 'failed', `<span class="text-red-300">${escapeHTML(String(ir.error || 'failed')).slice(0, 400)}</span><br>typed E(B-V) ${f(ex.E_BV_mean, 4)} used instead`);
      else set(2, 'pending', 'not queried yet<br>press Apply');
      set(3, 'done', `A_V ${f(ex.A_V_mean)} (${f(ex.AV_PER_EBV, 3)}&times;E)<br>A_G ${f(ex.A_G)} A_BP ${f(ex.A_BP)} A_RP ${f(ex.A_RP)}`);
      if (pa.state === 'done') {
        const chk = pa.check;
        const chkTxt = chk && chk.ok === false ? `<br><span class="text-amber-400">${escapeHTML(chk.message)}</span>` : (chk && chk.ok ? '<br>reddening check OK' : '');
        set(4, chk && chk.ok === false ? 'warn' : 'done', `${pa.rows} points, ${pa.kb} kB<br>A_V sent ${f(pa.A_V_sent)} &middot; Gaia EDR3=DR3${pa.tls_verified === false ? '<br><span class="text-amber-400">TLS not verified</span>' : ''}${chkTxt}`);
      } else if (pa.state === 'failed') set(4, 'failed', `<span class="text-red-300">${escapeHTML(String(pa.error || 'failed')).slice(0, 400)}</span>`);
      else set(4, 'pending', 'not fetched yet<br>press Apply');
      if (ft.state === 'done' && ft.median !== null && ft.median !== undefined) {
        set(5, ft.frac_near >= 0.6 ? 'done' : 'warn', `median offset ${f(ft.median)} mag<br>${f(ft.frac_near * 100, 0)}% of ${ft.n} stars within ${ft.near} mag<br>${escapeHTML(ft.plane)}`);
      } else set(5, 'pending', '-');
    }

    // Draws the isochrone CMD from the cached /compute-cmd reply. Star-type checkboxes filter ONLY this graph (the isochrone line is never filtered).
    function renderIsoCMD() {
      const st = window._isoCmd;
      if (!st) return;
      const { data, params, isoSrc, isoLabel, space } = st;
      const show = {};
      document.querySelectorAll('.iso-type-cb').forEach(cb => { show[cb.dataset.type] = cb.checked; });
      const starTypes = data.empirical.star_type || [];
      const okFlags = data.empirical.bprp_corr_ok || data.empirical.color.map(() => true);
      const typeOn = t => show[t] !== false;                    // unknown / unlisted types stay visible
      const cmdTraces = [];

      // counts next to each checkbox = all stars of that type in the file (before filtering)
      Object.keys(STAR_TYPE_COLORS).forEach((t, k) => {
        const el = document.getElementById('iso-cnt-' + k);
        if (el) el.textContent = `(${starTypes.filter(x => x === t).length})`;
      });
      const noteEl = document.getElementById('iso-plotted-note');      // stars without BP/RP or G photometry cannot be placed on a CMD
      if (noteEl) {
        const nIn = data.meta && Number.isFinite(data.meta.n_in) ? data.meta.n_in : starTypes.length, nOut = starTypes.length;
        noteEl.textContent = nIn > nOut ? `${nOut} of ${nIn} members plotted (${nIn - nOut} lack BP/RP/G photometry)` : `${nOut} members`;
      }
      const oorEl = document.getElementById('iso-cnt-oor');
      if (oorEl) oorEl.textContent = `(${okFlags.filter(ok => !ok).length})`;

      // stars whose intrinsic colour is outside the plausible range: grey, still respect their type checkbox
      const outIdx = okFlags.map((ok, i) => (!ok && typeOn(starTypes[i]) && show['__oor__'] !== false) ? i : -1).filter(i => i !== -1);
      if (outIdx.length > 0) {
        cmdTraces.push({
          x: outIdx.map(i => data.empirical.color[i]),
          y: outIdx.map(i => data.empirical.mag[i]),
          mode: 'markers',
          type: 'scatter',
          name: `colour out of range (${outIdx.length})`,
          marker: { color: '#9ca3af', size: 5, opacity: 0.45 }
        });
      }

      const uniqueTypes = Object.keys(STAR_TYPE_COLORS).filter(t => starTypes.includes(t));
      if (uniqueTypes.length > 0) {
        uniqueTypes.forEach(stype => {
          if (!typeOn(stype)) return;
          const idxList = starTypes.map((t, i) => (t === stype && okFlags[i]) ? i : -1).filter(i => i !== -1);
          if (!idxList.length) return;
          cmdTraces.push({
            x: idxList.map(i => data.empirical.color[i]),
            y: idxList.map(i => data.empirical.mag[i]),
            mode: 'markers',
            type: 'scatter',
            name: `${stype} (${idxList.length})`,
            marker: {
              color: STAR_TYPE_COLORS[stype] || '#00FFFF',
              size: 6,
              opacity: 0.88,
              line: {
                color: stype === 'White Dwarf' ? '#ffffff' : '#000000',
                width: stype === 'White Dwarf' ? 1 : 0.4
              }
            }
          });
        });
      } else {
        cmdTraces.push({
          x: data.empirical.color,
          y: data.empirical.mag,
          mode: 'markers',
          type: 'scatter',
          name: 'members',
          marker: { color: '#ef4444', size: 6, opacity: 0.85 }
        });
      }

      if (data.isochrone && data.isochrone.color && data.isochrone.color.length > 0) {
        cmdTraces.push({
          x: data.isochrone.color,
          y: data.isochrone.mag,
          mode: 'lines',
          type: 'scatter',
          name: isoLabel,
          line: { color: '#facc15', width: 2.5 }
        });
      }

      Plotly.react('plot-cmd-iso', cmdTraces, getPlotlyLayout({
        title: { text: isoSrc === 'parsec' ? `${space === 'reddened' ? 'Observed' : 'De-reddened'} CMD & PARSEC Isochrone (logAge = ${params.age}, [M/H] = ${params.metallicity})` : (space === 'reddened' ? 'Observed CMD' : 'De-reddened CMD'), font: { size: 13, color: '#e5e7eb' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: space === 'reddened' ? '(G_BP - G_RP) observed (mag)' : '(G_BP - G_RP)₀ (mag)' },
        yaxis: { ...getPlotlyLayout().yaxis, title: space === 'reddened' ? 'M_G = G - DM (mag, not de-reddened)' : 'M_G (mag)', autorange: 'reversed' },
        legend: { font: { color: '#e5e7eb', size: 10 }, bgcolor: 'rgba(10,14,26,0.75)', bordercolor: 'rgba(255,255,255,0.08)', borderwidth: 1 }
      }), PLOTLY_CONFIG);
    }

    async function computeCMDPlot(opts = {}) {
      const userClick = !!(opts && opts.user === true);   // only the Apply button may query IRSA / PARSEC; every other trigger is cache-only
      const targetFile = getTaskFilename('isochrone');
      if (!targetFile) {
        showToast('Please run GMM analysis first or upload a custom filtered CSV file', 'warning');
        return;
      }
      showModuleError('iso-error-box', 'iso-error-text', null);
      setButtonLoading('btn-apply-extinction', true);

      const ebvMode = document.getElementById('input-ebv-mode').value;
      const ebvInput = parseFloat(document.getElementById('input-ebv').value);
      document.getElementById('input-av').value = Number.isFinite(ebvInput) ? (STARNEST_AV_PER_EBV * ebvInput).toFixed(3) : '';

      const progressWrap = document.getElementById('per-star-warning');
      const progressBar = document.getElementById('per-star-progress-bar');
      if (ebvMode === 'per_star' && progressWrap && progressBar) {
        progressWrap.classList.remove('hidden');
        progressBar.style.width = '35%';
      }

      const params = {
        members_filename: targetFile,
        ebv_mode: ebvMode,
        ebv: Number.isFinite(ebvInput) ? ebvInput : null,
        plx_zp: parseFloat(document.getElementById('input-iso-plxzp').value) || 0.0,
        use_cluster_dm: document.getElementById('input-use-cluster-dm').checked,
        dist_mod: parseFloat(document.getElementById('input-dist-mod').value) || 7.66,
        age: parseFloat(document.getElementById('input-age').value) || 7.5,
        metallicity: parseFloat(document.getElementById('input-metallicity').value) || 0.0,
        av: Number.isFinite(ebvInput) ? STARNEST_AV_PER_EBV * ebvInput : 0.41,
        iso_extinction: document.getElementById('input-iso-ext').value,
        allow_network: userClick
      };

      try {
        const data = await apiCall('/compute-cmd', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(params),
          timeoutMs: userClick ? 900000 : undefined          // IRSA (per_star) + two CMD requests can take minutes
        });

        if (progressBar) progressBar.style.width = '100%';
        renderIsoWorkflow(data.meta, params);

        // Update Stage 2 Metadata & frac_colour_ok warning (< 0.9)
        if (data.meta) {
          document.getElementById('iso-meta-stats').classList.remove('hidden');
          document.getElementById('meta-ebv-mode').textContent = data.meta.mode;
          document.getElementById('meta-ebv-val').textContent = Number(data.meta.E_BV_mean !== undefined ? data.meta.E_BV_mean : data.meta.E_BV_median).toFixed(3);
          document.getElementById('meta-av-val').textContent = Number(data.meta.A_V_mean !== undefined ? data.meta.A_V_mean : data.meta.A_V).toFixed(3);
          document.getElementById('meta-dpc-val').textContent = Number(data.meta.d_pc).toFixed(1);
          document.getElementById('meta-dm-val').textContent = Number(data.meta.DM).toFixed(2);
          const fracOk = Number(data.meta.frac_colour_ok);
          const colOkEl = document.getElementById('meta-col-ok');
          colOkEl.textContent = `${(fracOk * 100).toFixed(1)}%`;
          colOkEl.className = fracOk < 0.9 ? 'text-amber-400 font-mono font-bold' : 'text-green-400 font-mono font-semibold';

          const warnBox = document.getElementById('iso-colour-warning');
          const warnTxt = document.getElementById('iso-colour-warning-text');
          if (fracOk < 0.9) {
            warnTxt.textContent = `Only ${(fracOk * 100).toFixed(1)}% of member stars fall within the plausible intrinsic colour range (-0.8 ≤ (BP-RP)₀ ≤ 5.0). SFD gives total line-of-sight extinction and over-corrects nearby clusters (< ~500 pc). Switch E(B-V) Mode to 'fixed' and enter a lower manual E(B-V).`;
            warnBox.classList.remove('hidden');
          } else {
            warnBox.classList.add('hidden');
          }
        }

        if (data.meta && data.meta.note && data.meta.note !== window._lastCmdNote) {
          window._lastCmdNote = data.meta.note;
          showToast(data.meta.note, 'warning');
        }
        const isoSrc = (data.isochrone && data.isochrone.source) || 'parsec';
        const isoLabel = isoSrc === 'parsec' ? `PARSEC Gaia DR3 (logAge = ${params.age}, [M/H] = ${params.metallicity}, A_V = ${Number((data.meta && data.meta.A_V_sent) || 0).toFixed(3)})`
          : isoSrc === 'custom' ? 'Uploaded isochrone'
          : 'Analytic placeholder (NOT PARSEC)';
        window._isoCmd = { data, params, isoSrc, isoLabel, space: (data.meta && data.meta.plot_space) || 'intrinsic' };      // cached so the checkboxes can re-filter without a recompute
        renderIsoCMD();

        const statsBox = document.getElementById('classification-stats');
        const countsList = document.getElementById('classification-counts-list');
        if (statsBox && countsList && data.classification_counts) {
          countsList.innerHTML = '';
          Object.entries(STAR_TYPE_COLORS).forEach(([stype, color]) => {
            const count = data.classification_counts[stype] || 0;
            const borderStyle = stype === 'White Dwarf' ? 'border: 1px solid #fff;' : '';
            countsList.innerHTML += `
              <div class="flex items-center justify-between text-xs py-1 border-b border-white/[0.04] last:border-0">
                <span class="flex items-center gap-2 text-gray-300">
                  <span class="w-2.5 h-2.5 rounded-full inline-block" style="background-color: ${color}; ${borderStyle}"></span>
                  ${stype}
                </span>
                <span class="font-mono font-semibold text-white">${count}</span>
              </div>
            `;
          });
          statsBox.classList.remove('hidden');
        }
      } catch (err) {
        showModuleError('iso-error-box', 'iso-error-text', err.message || String(err));
        showToast(err.message || 'Extinction / Isochrone error', 'error');
      } finally {
        setButtonLoading('btn-apply-extinction', false);
      }
    }

    function initKinematicPlots() {
      Plotly.newPlot('plot-3d', [{ x: [], y: [], z: [], type: 'scatter3d', mode: 'markers' }], {
        ...getPlotlyLayout(),
        scene: {
          xaxis: { title: 'X (pc)', gridcolor: 'rgba(255,255,255,0.05)', color: '#9ca3af' },
          yaxis: { title: 'Y (pc)', gridcolor: 'rgba(255,255,255,0.05)', color: '#9ca3af' },
          zaxis: { title: 'Z (pc)', gridcolor: 'rgba(255,255,255,0.05)', color: '#9ca3af' },
          bgcolor: 'rgba(26, 31, 53, 0.6)'
        },
        title: { text: '3D Spatial Distribution', font: { size: 13, color: '#9ca3af' } }
      }, PLOTLY_CONFIG);

      Plotly.newPlot('plot-2d-quiver', [{ x: [], y: [], type: 'scatter', mode: 'markers' }], getPlotlyLayout({
        title: { text: '2D Vector Kinematics', font: { size: 13, color: '#9ca3af' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'RA (deg)', autorange: 'reversed' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'Dec (deg)' }
      }), PLOTLY_CONFIG);
    }

    async function computeKinematicsPlots() {
      const targetFile = getTaskFilename('kinematics');
      if (!targetFile) {
        showToast('Please run GMM analysis first or upload a custom filtered CSV file', 'warning');
        return;
      }
      setButtonLoading('btn-compute-kinematics', true);

      try {
        const data = await apiCall('/compute-kinematics', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            filename: targetFile,
            kin3d_mode: document.getElementById('kin3d-mode').value,
            kin3d_max_frac: parseFloat(document.getElementById('kin3d-maxfrac').value),
            kin3d_depth: parseFloat(document.getElementById('kin3d-depth').value),
            kin3d_frame: document.getElementById('kin3d-frame').value,
            kin3d_centre: document.getElementById('kin3d-centre').checked
          })
        });
        window._kinComputed = true;

        const clusterLabel = (taskSources.kinematics.mode === 'custom' && taskSources.kinematics.customLabel)
          ? taskSources.kinematics.customLabel.replace(/\.csv$/i, '')
          : (data.cluster_name || 'Cluster');
        const starTypes3D = data.positions_3d.star_type || [];
        const uniqueTypes3D = Object.keys(STAR_TYPE_COLORS).filter(t => starTypes3D.includes(t));

        // --- PLOT 2: 3D Spatial Distribution grouped by Star_Type for Shared Legend ---
        const traces3D = uniqueTypes3D.map(stype => {
          const idx = starTypes3D.map((t, i) => t === stype ? i : -1).filter(i => i !== -1);
          return {
            x: idx.map(i => data.positions_3d.x[i]),
            y: idx.map(i => data.positions_3d.y[i]),
            z: idx.map(i => data.positions_3d.z[i]),
            type: 'scatter3d',
            mode: 'markers',
            name: stype,
            customdata: idx.map(i => [data.positions_3d.d[i], data.positions_3d.frac[i]]),
            hovertemplate: 'd = %{customdata[0]:.1f} pc<br>&sigma;&varpi;/&varpi; = %{customdata[1]:.3f}<br>%{x:.1f}, %{y:.1f}, %{z:.1f} pc<extra>' + stype + '</extra>',
            marker: {
              color: STAR_TYPE_COLORS[stype] || 'gray',
              size: 4.5,
              opacity: 0.8,
              line: {
                color: stype === 'White Dwarf' ? '#ffffff' : '#000000',
                width: stype === 'White Dwarf' ? 0.8 : 0.2
              }
            }
          };
        });

        const m3 = data.meta3d || {};
        const pre = m3.centred ? '\u0394' : '';
        const axNames = m3.frame === 'icrs' ? ['X (RA=0\u00b0)', 'Y (RA=90\u00b0)', 'Z (NCP)'] : ['X (l=0\u00b0, GC)', 'Y (l=90\u00b0)', 'Z (NGP)'];
        const diagEl = document.getElementById('kin3d-diag');
        if (diagEl) {
          const dg = m3.diag;
          const modeTxt = m3.mode === 'precision' ? `precision cut \u03c3\u03d6/\u03d6 \u2264 ${m3.max_frac}` : m3.mode === 'shrink' ? `shrunk to cluster parallax (depth ${m3.depth_pc} pc assumed)` : 'raw 1/\u03d6, no cut';
          const needle = dg && dg.elongation > 2 && dg.angle_to_los_deg < 20;
          const hasExp = dg && Number.isFinite(dg.expected_err_sd_pc);
          // radial spread within 1.5x of the parallax-error budget => depth is unresolved (error-limited), not a measured shape
          const errLimited = hasExp && dg.radial_sd_pc < 1.5 * dg.expected_err_sd_pc;
          const errTxt = hasExp ? ` &middot; line-of-sight spread ${dg.radial_sd_pc.toFixed(1)} pc vs ${dg.expected_err_sd_pc.toFixed(1)} pc expected from parallax errors alone (implied intrinsic depth &approx; ${dg.implied_depth_pc.toFixed(1)} pc` +
            `${errLimited ? ', error-limited: depth not resolved' : ''})` : '';
          diagEl.innerHTML = `${m3.n_used} stars plotted (${m3.n_dropped} excluded) &middot; ${modeTxt} &middot; median d = ${Number.isFinite(m3.median_d_pc) ? m3.median_d_pc.toFixed(1) : '-'} pc` +
            (dg ? ` &middot; spread along major axis ${dg.major_sd_pc.toFixed(1)} pc vs ${dg.minor_sd_pc.toFixed(1)} pc across (${dg.elongation.toFixed(1)}:1, ${dg.angle_to_los_deg.toFixed(0)}&deg; from the line of sight)` : '') +
            errTxt +
            (needle && !hasExp ? ` <span class="text-yellow-400">&mdash; elongated along the line of sight: this is parallax-error smearing, not the cluster's shape. Tighten the precision cut or use shrinkage.</span>` : '') +
            (needle && hasExp && errLimited ? ` <span class="text-yellow-400">&mdash; elongation matches the parallax-error budget: a measurement limit, not the cluster's shape. Do not read 3D structure along this axis.</span>` : '') +
            (needle && hasExp && !errLimited ? ` <span class="text-yellow-400">&mdash; elongation exceeds the parallax-error budget: real depth structure, a second population or field contamination. Check membership.</span>` : '') +
            (m3.note ? ` <span class="text-yellow-400">${escapeHTML(m3.note)}</span>` : '');
        }
        Plotly.react('plot-3d', traces3D, {
          ...getPlotlyLayout(),
          scene: {
            aspectmode: 'data',
            xaxis: { title: `${pre}${axNames[0]} [pc]`, gridcolor: 'rgba(255,255,255,0.06)', color: '#9ca3af' },
            yaxis: { title: `${pre}${axNames[1]} [pc]`, gridcolor: 'rgba(255,255,255,0.06)', color: '#9ca3af' },
            zaxis: { title: `${pre}${axNames[2]} [pc]`, gridcolor: 'rgba(255,255,255,0.06)', color: '#9ca3af' },
            bgcolor: 'rgba(26, 31, 53, 0.6)'
          },
          legend: { font: { color: '#e5e7eb', size: 11 }, bgcolor: 'rgba(10, 14, 26, 0.75)', bordercolor: 'rgba(255,255,255,0.1)', borderwidth: 1 },
          title: { text: `3D Spatial Distribution: ${clusterLabel}`, font: { size: 13, color: '#e5e7eb' } }
        }, PLOTLY_CONFIG);

        // --- PLOT 1: 2D Vector Kinematics (Quiver with arrowheads + Gold Star Center) ---
        const raArr = data.vectors_2d.ra || [];
        const decArr = data.vectors_2d.dec || [];
        const uArr = data.vectors_2d.u || [];
        const vArr = data.vectors_2d.v || [];
        const starTypes2D = data.vectors_2d.star_type || [];

        const [raLo2, raHi2] = raArr.length ? minMax(raArr) : [0, 1];
        const [deLo2, deHi2] = decArr.length ? minMax(decArr) : [0, 1];
        const cosMid = Math.max(0.1, Math.cos((((deLo2 + deHi2) / 2) * Math.PI) / 180));
        window._cosMid = cosMid;
        // spans in ON-SKY degrees (RA span shrinks by cos(dec)); arrows are built on the sky and converted back to RA degrees
        const raSpan = ((raHi2 - raLo2) * cosMid) || 1.0;
        const decSpan = (deHi2 - deLo2) || 1.0;
        const arrowScale = Math.max(raSpan, decSpan) / 25.0;
        const headLength = arrowScale * 0.28;
        const headAngle = Math.PI / 6;

        const uniqueTypes2D = Object.keys(STAR_TYPE_COLORS).filter(t => starTypes2D.includes(t));
        const quiverTraces = [];

        uniqueTypes2D.forEach(stype => {
          const lineX = [];
          const lineY = [];
          for (let i = 0; i < raArr.length; i++) {
            if (starTypes2D[i] !== stype) continue;
            const x0 = raArr[i];
            const y0 = decArr[i];
            const cd = Math.max(0.1, Math.cos((y0 * Math.PI) / 180));
            const dx = uArr[i] * arrowScale;          // on-sky east offset (deg)
            const dy = vArr[i] * arrowScale;
            const x1 = x0 + dx / cd;                  // -> RA degrees
            const y1 = y0 + dy;

            const angle = Math.atan2(dy, dx);
            const hx1 = x1 - (headLength * Math.cos(angle - headAngle)) / cd;
            const hy1 = y1 - headLength * Math.sin(angle - headAngle);
            const hx2 = x1 - (headLength * Math.cos(angle + headAngle)) / cd;
            const hy2 = y1 - headLength * Math.sin(angle + headAngle);

            // Shaft + two arrowhead wings
            lineX.push(x0, x1, null, x1, hx1, null, x1, hx2, null);
            lineY.push(y0, y1, null, y1, hy1, null, y1, hy2, null);
          }

          quiverTraces.push({
            x: lineX,
            y: lineY,
            mode: 'lines',
            type: 'scatter',
            name: stype,
            line: {
              color: STAR_TYPE_COLORS[stype] || 'gray',
              width: 1.5
            },
            hoverinfo: 'name'
          });
        });

        // Cluster Center Gold Star (matching ax1.scatter median(RA), median(Dec), marker='*', c='gold')
        if (raArr.length > 0) {
          quiverTraces.push({
            x: [data.vectors_2d.center_ra],
            y: [data.vectors_2d.center_dec],
            mode: 'markers',
            type: 'scatter',
            name: 'Cluster Center',
            marker: {
              symbol: 'star',
              size: 16,
              color: 'gold',
              line: { color: 'black', width: 1.5 }
            }
          });
        }

        Plotly.react('plot-2d-quiver', quiverTraces, getPlotlyLayout({
          title: { text: `2D Vector Kinematics: ${clusterLabel}`, font: { size: 13, color: '#e5e7eb' } },
          xaxis: { ...getPlotlyLayout().xaxis, title: 'RA (deg)', autorange: 'reversed' },
          yaxis: { ...getPlotlyLayout().yaxis, title: 'Dec (deg)', scaleanchor: 'x', scaleratio: 1 / (window._cosMid || 1) },
          showlegend: true,
          legend: { font: { color: '#e5e7eb', size: 11 }, bgcolor: 'rgba(10, 14, 26, 0.75)', bordercolor: 'rgba(255,255,255,0.1)', borderwidth: 1 }
        }), PLOTLY_CONFIG);
      } catch (err) {
        console.error(err);
        showToast(err.message || 'Kinematics error', 'error');
      } finally {
        setButtonLoading('btn-compute-kinematics', false);
      }
    }

    function initIMFPlots() {
      Plotly.newPlot('plot-completeness', [{ x: [], y: [], type: 'scatter' }], getPlotlyLayout({
        title: { text: 'Completeness', font: { size: 13, color: '#9ca3af' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'G magnitude' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'Completeness Fraction' }
      }), PLOTLY_CONFIG);

      Plotly.newPlot('plot-imf', [{ x: [], y: [], type: 'scatter' }], getPlotlyLayout({
        title: { text: 'Initial Mass Function', font: { size: 13, color: '#9ca3af' } },
        xaxis: { ...getPlotlyLayout().xaxis, title: 'log₁₀(M/M☉)' },
        yaxis: { ...getPlotlyLayout().yaxis, title: 'log₁₀(dN/dlog M)' }
      }), PLOTLY_CONFIG);
    }

    async function computeIMFPlots(opts = {}) {
      const targetFile = getTaskFilename('imf');
      if (!targetFile) {
        showToast('Please run GMM analysis first or upload a custom filtered CSV file', 'warning');
        return;
      }
      showModuleError('imf-error-box', 'imf-error-text', null);

      const autoBins = document.getElementById('check-auto-bins').checked;
      const logAgeVal = parseFloat(document.getElementById('input-age').value) || 7.5;
      const params = {
        filename: targetFile,
        bins: autoBins ? null : parseInt(document.getElementById('slider-bins').value, 10),
        fit_min_logM: parseFloat(document.getElementById('slider-fit-logm').value),
        mass_min: parseFloat(document.getElementById('slider-mass-min').value),
        mass_max: parseFloat(document.getElementById('slider-mass-max').value),
        include_companions: document.getElementById('check-include-companions').checked,
        log_age: logAgeVal
      };
      const reqId = (window._imfReq = (window._imfReq || 0) + 1);   // newest request wins

      try {
        const data = await apiCall('/compute-imf', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(params)
        });

        if (reqId !== window._imfReq) return;      // a newer slider position superseded this answer

        // The IMF page is always labelled 'Initial Mass Function (IMF)' (heading, plot-card title and plot title)
        const mfKind = data.imf.mf_kind || 'Initial Mass Function (IMF)';
        const headingAccent = document.getElementById('imf-page-heading-accent');
        const regimeBadge = document.getElementById('imf-age-regime-badge');
        const plotCardTitle = document.getElementById('imf-plot-card-title');
        if (headingAccent) headingAccent.textContent = mfKind;
        if (regimeBadge) regimeBadge.textContent = `log(Age/yr) = ${logAgeVal.toFixed(2)} \u2192 Initial Mass Function (IMF)`;
        if (plotCardTitle) {
          plotCardTitle.textContent = `${mfKind} — model comparison (Salpeter 1955, Kroupa 2001 & Chabrier 2003)`;
        }

        if (autoBins && data.imf.meta && data.imf.meta.n_bins) {
          document.getElementById('bins-value').textContent = `Auto (${data.imf.meta.n_bins})`;
          document.getElementById('slider-bins').value = data.imf.meta.n_bins;   // slider shows the automatic choice
        }

        // Render Relative Chi2 Ranking Table
        if (data.imf.chi2_table && data.imf.chi2_table.length > 0) {
          document.getElementById('imf-chi2-card').classList.remove('hidden');
          const tbody = document.getElementById('imf-chi2-tbody');
          tbody.innerHTML = '';
          data.imf.chi2_table.forEach((row, idx) => {
            const rankBadge = idx === 0
              ? '<span class="px-2 py-0.5 rounded-full bg-green-500/20 text-green-300 border border-green-500/30 text-[10px]">#1 Best Relative Fit</span>'
              : `<span class="text-gray-400">#${idx + 1}</span>`;
            tbody.innerHTML += `
              <tr class="hover:bg-white/[0.02]">
                <td class="py-2 px-3">${rankBadge}</td>
                <td class="py-2 px-3 font-semibold text-white">${row.Model}</td>
                <td class="py-2 px-3 text-accent-cyan">${Number(row.red_chi2).toFixed(2)}</td>
                <td class="py-2 px-3 text-gray-400">${row.dof}</td>
              </tr>
            `;
          });
          if (data.imf.meta) {
            const m = data.imf.meta;
            document.getElementById('imf-meta-summary').textContent =
              `N=${m.n_masses} | ${m.n_bins} bins (ΔlogM=${Number(m.bin_width_dex).toFixed(3)} dex) | ${m.n_fit_bins} complete bins (N_fit=${m.n_in_fit})` +
              (m.n_excluded_evolved ? ` | ${m.n_excluded_evolved} evolved stars excluded` : '');
          }
        }

        Plotly.react('plot-completeness', [
          {
            x: data.completeness.mag,
            y: data.completeness.fraction,
            mode: 'lines+markers',
            type: 'scatter',
            name: 'Assumed completeness (not measured)',
            line: { color: '#00d4ff', width: 2.2 },
            marker: { color: '#00d4ff', size: 5 }
          },
          {
            x: [Math.min(...data.completeness.mag), Math.max(...data.completeness.mag)],
            y: [0.5, 0.5],
            mode: 'lines',
            type: 'scatter',
            name: '50% Completeness Limit',
            line: { color: '#fbbf24', width: 1.5, dash: 'dash' }
          }
        ], getPlotlyLayout({
          title: { text: 'Assumed completeness model (illustrative, G50 = 17.2; not measured from data)', font: { size: 13, color: '#e5e7eb' } },
          xaxis: { ...getPlotlyLayout().xaxis, title: 'G magnitude' },
          yaxis: { ...getPlotlyLayout().yaxis, title: 'Completeness Fraction', range: [0, 1.1] },
          legend: { font: { color: '#d1d5db' }, bgcolor: 'rgba(0,0,0,0)' }
        }), PLOTLY_CONFIG);

        const chiMap = {};
        (data.imf.chi2_table || []).forEach(r => { chiMap[r.Model] = Number(r.red_chi2).toFixed(1); });

        const nMasses = data.imf.meta ? data.imf.meta.n_masses : data.imf.ms_count;
        const nBinsUsed = data.imf.meta ? data.imf.meta.n_bins : params.bins;

        const imfTraces = [];
        // Complete bins (filled markers)
        if (data.imf.complete_bins && data.imf.complete_bins.logM_center.length > 0) {
          imfTraces.push({
            x: data.imf.complete_bins.logM_center,
            y: data.imf.complete_bins.log_dN,
            error_y: {
              type: 'data',
              array: data.imf.complete_bins.log_err,
              visible: true,
              color: '#e5e7eb',
              thickness: 1.4,
              width: 4
            },
            mode: 'markers',
            type: 'scatter',
            name: `Complete (N=${nMasses}, ${nBinsUsed} bins)`,
            marker: { color: '#00d4ff', size: 7.5, symbol: 'circle', line: { color: '#fff', width: 0.8 } }
          });
        }

        // Incomplete bins below FIT_MIN_LOGM (open grey markers, excluded from fit & scaling)
        if (data.imf.incomplete_bins && data.imf.incomplete_bins.logM_center.length > 0) {
          imfTraces.push({
            x: data.imf.incomplete_bins.logM_center,
            y: data.imf.incomplete_bins.log_dN,
            error_y: {
              type: 'data',
              array: data.imf.incomplete_bins.log_err,
              visible: true,
              color: '#9ca3af',
              thickness: 1.2,
              width: 3
            },
            mode: 'markers',
            type: 'scatter',
            name: 'incomplete (not fitted)',
            marker: { color: 'rgba(0,0,0,0)', size: 7.5, symbol: 'circle-open', line: { color: '#9ca3af', width: 1.8 } }
          });
        }

        if (data.imf.theoretical) {
          imfTraces.push({
            x: data.imf.theoretical.log_m,
            y: data.imf.theoretical.salpeter,
            mode: 'lines',
            type: 'scatter',
            name: `Salpeter (χ²_ν=${chiMap.Salpeter || '-'})`,
            line: { color: 'steelblue', width: 2.2, dash: 'dash' }
          });
          imfTraces.push({
            x: data.imf.theoretical.log_m,
            y: data.imf.theoretical.kroupa,
            mode: 'lines',
            type: 'scatter',
            name: `Kroupa (χ²_ν=${chiMap.Kroupa || '-'})`,
            line: { color: 'darkorange', width: 2.2, dash: 'dashdot' }
          });
          imfTraces.push({
            x: data.imf.theoretical.log_m,
            y: data.imf.theoretical.chabrier,
            mode: 'lines',
            type: 'scatter',
            name: `Chabrier (χ²_ν=${chiMap.Chabrier || '-'})`,
            line: { color: 'firebrick', width: 2.2 }
          });
        }

        const fl = data.imf.fit_line;
        if (fl && Array.isArray(fl.x) && data.imf.slope !== undefined) {
          const se = Number(data.imf.slope_err);
          imfTraces.push({
            x: fl.x,
            y: fl.y,
            mode: 'lines',
            type: 'scatter',
            name: `Empirical fit: Γ = ${data.imf.slope.toFixed(2)}${Number.isFinite(se) ? ' ± ' + se.toFixed(2) : ''}`,
            line: { color: '#a3e635', width: 2.8 },
            hovertemplate: 'log M = %{x:.2f}<br>fit log dN/dlogM = %{y:.2f}<extra>Empirical fit</extra>'
          });
        }

        const allY = [
          ...(data.imf.complete_bins ? data.imf.complete_bins.log_dN : []),
          ...(data.imf.incomplete_bins ? data.imf.incomplete_bins.log_dN : [])
        ].filter(Number.isFinite);
        const yMin = allY.length ? Math.min(...allY) - 0.8 : -0.5;
        const yMax = allY.length ? Math.max(...allY) + 0.8 : 3.0;

        imfTraces.push({
          x: [params.fit_min_logM, params.fit_min_logM],
          y: [yMin, yMax],
          mode: 'lines',
          type: 'scatter',
          name: `Completeness limit (logM = ${params.fit_min_logM.toFixed(2)})`,
          line: { color: '#9ca3af', width: 1.3, dash: 'dot' }
        });

        Plotly.react('plot-imf', imfTraces, getPlotlyLayout({
          uirevision: 'imf',                     // keep the user's zoom while sliders redraw
          title: { text: `${mfKind} — model comparison`, font: { size: 13, color: '#e5e7eb' } },
          xaxis: { ...getPlotlyLayout().xaxis, title: 'log₁₀(M / M☉)' },
          yaxis: { ...getPlotlyLayout().yaxis, title: 'log₁₀(dN / dlog M)' },
          legend: { font: { color: '#e5e7eb', size: 10 }, bgcolor: 'rgba(10,14,26,0.75)', bordercolor: 'rgba(255,255,255,0.08)', borderwidth: 1 }
        }), PLOTLY_CONFIG);

        if (data.imf.slope !== undefined) {
          const se = Number(data.imf.slope_err);
          const SALP = -1.35;
          document.getElementById('imf-slope-display').classList.remove('hidden');
          document.getElementById('imf-slope-value').textContent = `${data.imf.slope.toFixed(3)}${Number.isFinite(se) ? ' ± ' + se.toFixed(3) : ''}`;
          const nFit = data.imf.meta ? data.imf.meta.n_fit_bins : '-';
          document.getElementById('imf-r2-value').textContent =
            `R² = ${(data.imf.r_squared || 0).toFixed(3)} | ${nFit} bins in fit | N = ${nMasses || '-'}`;
          const delta = data.imf.slope - SALP;
          const sig = (Number.isFinite(se) && se > 0) ? ` (${(Math.abs(delta) / se).toFixed(1)}σ)` : '';
          const dEl = document.getElementById('imf-slope-delta');
          if (dEl) dEl.textContent = `vs Salpeter Γ = -1.35: ${delta >= 0 ? '+' : ''}${delta.toFixed(2)}${sig}`;
        }
      } catch (err) {
        if (reqId !== window._imfReq) return;
        showModuleError('imf-error-box', 'imf-error-text', err.message || String(err));
        if (!(opts && opts.silent)) showToast(err.message || 'IMF computation error', 'error');
        // never leave a stale slope/fit on screen when the current settings are invalid
        document.getElementById('imf-slope-display').classList.add('hidden');
        Plotly.react('plot-imf', [], getPlotlyLayout({
          uirevision: 'imf',
          xaxis: { ...getPlotlyLayout().xaxis, title: 'log₁₀(M / M☉)' },
          yaxis: { ...getPlotlyLayout().yaxis, title: 'log₁₀(dN / dlog M)' },
          annotations: [{ text: err.message || 'IMF computation error', showarrow: false, xref: 'paper', yref: 'paper', x: 0.5, y: 0.5, font: { color: '#fbbf24', size: 13 } }]
        }), PLOTLY_CONFIG);
      }
    }

    document.addEventListener('DOMContentLoaded', () => {
      initWorkflowTiles();
      loadTeamPhotos();
      const fontPromise = document.fonts ? document.fonts.ready : Promise.resolve();
      fontPromise.then(() => {
        const loadingOverlay = document.getElementById('loading-overlay');
        if (loadingOverlay) {
          gsap.to(loadingOverlay, {
            opacity: 0,
            duration: 0.5,
            ease: 'power2.inOut',
            delay: 0.2,
            onComplete: () => {
              loadingOverlay.remove();
              animateHomeEntrance();
            }
          });
        }
      });

      requestAnimationFrame(() => {
        const homeBtn = document.querySelector('.nav-btn[data-page="home"]');
        if (homeBtn) moveNavPill(homeBtn);
      });

      document.querySelectorAll('.nav-btn').forEach(btn => {
        btn.addEventListener('click', () => {
          const page = btn.getAttribute('data-page');
          if (page) navigateToPage(page);
        });
      });

      document.querySelectorAll('[data-navigate]').forEach(link => {
        link.addEventListener('click', (e) => {
          e.preventDefault();
          const page = link.getAttribute('data-navigate');
          if (page) navigateToPage(page);
        });
      });

      window.addEventListener('resize', () => {
        const activeBtn = document.querySelector(`.nav-btn[data-page="${currentPage}"]`);
        if (activeBtn) moveNavPill(activeBtn);
      });

      // Initialize all Plotly canvases
      initMembershipPlots();
      initCMDPlot();
      initKinematicPlots();
      initIMFPlots();
      updateEngineBadge();
      probeBackend();

      // Setup Raw CSV Dropzone
      const rawDz = setupDropzone('dropzone-raw', 'file-raw-data', 'raw-file-name');
      if (rawDz && rawDz.fileInput) {
        rawDz.fileInput.addEventListener('change', async (e) => {
          const file = (e.target.files && e.target.files[0]) || e.target._selectedFile;
          if (!file) return;
          showModuleError('gmm-error-box', 'gmm-error-text', null);
          const text = await file.text();
          const rows = parseCSVText(text);
          const missing = validateColumnsJS(rows);
          if (missing.length > 0) {
            const msg = `Missing required columns: [${missing.map(c => `'${c}'`).join(', ')}]`;
            showModuleError('gmm-error-box', 'gmm-error-text', msg);
            showToast(msg, 'error');
            return;
          }
          localStore.files[file.name] = rows;
          localStore.ver[file.name] = ++localStore.verCounter;      // same name, new content -> caches invalid
          localStore.activeClusterName = file.name;
          window.globalFileName = file.name;
          window.rawFileName = file.name;
          window.gmmFilteredFileName = null;                         // members of the previous file are stale now
          localStore.latestMembers = []; localStore.analysis = null;

          // Also attempt upload to FastAPI backend if running locally
          const formData = new FormData();
          formData.append('file', file);
          try {
            const res = await fetch(API_BASE + '/upload-raw', { method: 'POST', body: formData });
            if (res.ok) {
              const json = await res.json();
              window.globalFileName = json.filename;
              window.rawFileName = json.filename;
              localStore.files[json.filename] = rows;
              localStore.ver[json.filename] = localStore.ver[file.name];
            }
          } catch (_) { /* Local in-browser store handles it */ }

          showToast(`Validated ${file.name} (${rows.length} stars). Click "Run GMM Analysis"!`, 'success');
        });
      }

      // Membership Buttons
      document.getElementById('btn-run-gmm').addEventListener('click', runGMMAnalysis);
      document.getElementById('btn-download-members').addEventListener('click', () => {
        downloadCSVRows(localStore.latestMembers, 'filtered_classified_members.csv');
      });

      async function registerCustomTaskFile(task, file) {
        const text = await file.text();
        const rows = parseCSVText(text);
        if (!rows.length) {
          showToast('Uploaded CSV appears to be empty or invalid.', 'error');
          return null;
        }
        let storedKey = `${task}_custom_${file.name}`;
        localStore.files[storedKey] = rows;
        localStore.ver[storedKey] = ++localStore.verCounter;

        const formData = new FormData();
        formData.append('file', file);
        try {
          const res = await fetch(API_BASE + '/upload-raw', { method: 'POST', body: formData });
          if (res.ok) {
            const json = await res.json();
            storedKey = json.filename;
            localStore.files[storedKey] = rows;
            localStore.ver[storedKey] = ++localStore.verCounter;
          }
        } catch (_) { /* Handled by in-browser store */ }

        taskSources[task].mode = 'custom';
        taskSources[task].customFile = storedKey;
        taskSources[task].customLabel = file.name;
        showToast(`Loaded custom file "${file.name}" (${rows.length} stars) for ${task.toUpperCase()}.`, 'success');
        return { storedKey, rows };
      }

      function updateTaskSourceUI(task) {
        const isCustom = taskSources[task].mode === 'custom';
        const label = escapeHTML(getTaskDisplayLabel(task));

        if (task === 'isochrone') {
          const btnGmm = document.getElementById('btn-iso-use-gmm');
          const btnCustom = document.getElementById('btn-iso-use-custom');
          const gmmInfo = document.getElementById('iso-gmm-info');
          const customWrap = document.getElementById('iso-custom-upload-wrap');
          const statusEl = document.getElementById('iso-source-status');
          const badgeEl = document.getElementById('badge-source-isochrone');

          btnGmm.className = !isCustom
            ? 'py-2 px-3 rounded-lg text-xs font-medium border border-accent-cyan/40 bg-accent-cyan/15 text-white flex items-center justify-center gap-1.5 transition-all'
            : 'py-2 px-3 rounded-lg text-xs font-medium border border-white/[0.08] bg-white/[0.03] text-gray-400 hover:text-white flex items-center justify-center gap-1.5 transition-all';
          btnCustom.className = isCustom
            ? 'py-2 px-3 rounded-lg text-xs font-medium border border-accent-violet/40 bg-accent-violet/15 text-white flex items-center justify-center gap-1.5 transition-all'
            : 'py-2 px-3 rounded-lg text-xs font-medium border border-white/[0.08] bg-white/[0.03] text-gray-400 hover:text-white flex items-center justify-center gap-1.5 transition-all';

          gmmInfo.classList.toggle('hidden', isCustom);
          customWrap.classList.toggle('hidden', !isCustom);
          statusEl.textContent = getTaskDisplayLabel(task);
          badgeEl.innerHTML = `<i class="fas ${isCustom ? 'fa-file-csv' : 'fa-filter'} text-[9px]"></i> ${label}`;
        }

        if (task === 'kinematics') {
          const btnGmm = document.getElementById('btn-kin-use-gmm');
          const btnCustom = document.getElementById('btn-kin-use-custom');
          const customWrap = document.getElementById('kin-custom-upload-wrap');
          const statusEl = document.getElementById('kin-source-status');
          const b3d = document.getElementById('badge-source-kinematics-3d');
          const b2d = document.getElementById('badge-source-kinematics-2d');

          btnGmm.className = !isCustom
            ? 'px-3.5 py-1.5 rounded-lg text-xs font-medium bg-accent-cyan/15 text-white border border-accent-cyan/40 transition-all flex items-center gap-1.5'
            : 'px-3.5 py-1.5 rounded-lg text-xs font-medium text-gray-400 hover:text-white border border-transparent transition-all flex items-center gap-1.5';
          btnCustom.className = isCustom
            ? 'px-3.5 py-1.5 rounded-lg text-xs font-medium bg-accent-violet/15 text-white border border-accent-violet/40 transition-all flex items-center gap-1.5'
            : 'px-3.5 py-1.5 rounded-lg text-xs font-medium text-gray-400 hover:text-white border border-transparent transition-all flex items-center gap-1.5';

          customWrap.classList.toggle('hidden', !isCustom);
          statusEl.innerHTML = `Active: <strong class="text-accent-cyan">${label}</strong>`;
          b3d.innerHTML = `<i class="fas ${isCustom ? 'fa-file-csv' : 'fa-filter'} text-[9px]"></i> ${label}`;
          b2d.innerHTML = `<i class="fas ${isCustom ? 'fa-file-csv' : 'fa-filter'} text-[9px]"></i> ${label}`;
        }

        if (task === 'imf') {
          const btnGmm = document.getElementById('btn-imf-use-gmm');
          const btnCustom = document.getElementById('btn-imf-use-custom');
          const gmmInfo = document.getElementById('imf-gmm-info');
          const customWrap = document.getElementById('imf-custom-upload-wrap');
          const statusEl = document.getElementById('imf-source-status');
          const bImf = document.getElementById('badge-source-imf');
          const bComp = document.getElementById('badge-source-completeness');

          btnGmm.className = !isCustom
            ? 'py-1.5 px-2 rounded-lg text-[11px] font-medium border border-accent-cyan/40 bg-accent-cyan/15 text-white flex items-center justify-center gap-1 transition-all'
            : 'py-1.5 px-2 rounded-lg text-[11px] font-medium border border-white/[0.08] bg-white/[0.03] text-gray-400 hover:text-white flex items-center justify-center gap-1 transition-all';
          btnCustom.className = isCustom
            ? 'py-1.5 px-2 rounded-lg text-[11px] font-medium border border-accent-violet/40 bg-accent-violet/15 text-white flex items-center justify-center gap-1 transition-all'
            : 'py-1.5 px-2 rounded-lg text-[11px] font-medium border border-white/[0.08] bg-white/[0.03] text-gray-400 hover:text-white flex items-center justify-center gap-1 transition-all';

          gmmInfo.classList.toggle('hidden', isCustom);
          customWrap.classList.toggle('hidden', !isCustom);
          statusEl.textContent = isCustom ? 'Custom CSV' : 'GMM Filtered';
          bImf.innerHTML = `<i class="fas ${isCustom ? 'fa-file-csv' : 'fa-filter'} text-[9px]"></i> ${label} (MS Only)`;
          bComp.innerHTML = `<i class="fas ${isCustom ? 'fa-file-csv' : 'fa-filter'} text-[9px]"></i> ${label}`;
        }
      }

      // Isochrone Data Source Controls
      setupDropzone('dropzone-iso-custom', 'file-iso-custom', 'iso-custom-file-name');
      document.getElementById('btn-iso-use-gmm').addEventListener('click', () => {
        taskSources.isochrone.mode = 'gmm';
        updateTaskSourceUI('isochrone');
        if (getTaskFilename('isochrone')) computeCMDPlot();
      });
      document.getElementById('btn-iso-use-custom').addEventListener('click', () => {
        taskSources.isochrone.mode = 'custom';
        updateTaskSourceUI('isochrone');
        if (taskSources.isochrone.customFile) {
          computeCMDPlot();
        } else {
          document.getElementById('file-iso-custom').click();
        }
      });
      document.getElementById('file-iso-custom').addEventListener('change', async (e) => {
        const file = (e.target.files && e.target.files[0]) || e.target._selectedFile;
        if (!file) return;
        const res = await registerCustomTaskFile('isochrone', file);
        if (res) {
          updateTaskSourceUI('isochrone');
          computeCMDPlot();
        }
      });

      // Kinematics Data Source Controls
      setupDropzone('dropzone-kin-custom', 'file-kin-custom', 'kin-custom-file-name');
      document.getElementById('btn-kin-use-gmm').addEventListener('click', () => {
        taskSources.kinematics.mode = 'gmm';
        updateTaskSourceUI('kinematics');
        if (getTaskFilename('kinematics')) computeKinematicsPlots();
      });
      document.getElementById('btn-kin-use-custom').addEventListener('click', () => {
        taskSources.kinematics.mode = 'custom';
        updateTaskSourceUI('kinematics');
        if (taskSources.kinematics.customFile) {
          computeKinematicsPlots();
        } else {
          document.getElementById('file-kin-custom').click();
        }
      });
      document.getElementById('file-kin-custom').addEventListener('change', async (e) => {
        const file = (e.target.files && e.target.files[0]) || e.target._selectedFile;
        if (!file) return;
        const res = await registerCustomTaskFile('kinematics', file);
        if (res) {
          updateTaskSourceUI('kinematics');
          computeKinematicsPlots();
        }
      });

      // IMF Data Source Controls
      setupDropzone('dropzone-imf-custom', 'file-imf-custom', 'imf-custom-file-name');
      document.getElementById('btn-imf-use-gmm').addEventListener('click', () => {
        taskSources.imf.mode = 'gmm';
        updateTaskSourceUI('imf');
        if (getTaskFilename('imf')) computeIMFPlots();
      });
      document.getElementById('btn-imf-use-custom').addEventListener('click', () => {
        taskSources.imf.mode = 'custom';
        updateTaskSourceUI('imf');
        if (taskSources.imf.customFile) {
          computeIMFPlots();
        } else {
          document.getElementById('file-imf-custom').click();
        }
      });
      document.getElementById('file-imf-custom').addEventListener('change', async (e) => {
        const file = (e.target.files && e.target.files[0]) || e.target._selectedFile;
        if (!file) return;
        const res = await registerCustomTaskFile('imf', file);
        if (res) {
          updateTaskSourceUI('imf');
          computeIMFPlots();
        }
      });

      // Isochrone Dropzone & Controls
      // CMD star-type checkboxes: client-side filter of the isochrone graph only (no recompute, no effect on other tabs)
      document.querySelectorAll('.iso-type-cb').forEach(cb => cb.addEventListener('change', renderIsoCMD));
      document.getElementById('btn-iso-types-all').addEventListener('click', () => { document.querySelectorAll('.iso-type-cb').forEach(cb => { cb.checked = true; }); renderIsoCMD(); });
      document.getElementById('btn-iso-types-none').addEventListener('click', () => { document.querySelectorAll('.iso-type-cb').forEach(cb => { cb.checked = false; }); renderIsoCMD(); });

      document.getElementById('btn-apply-extinction').addEventListener('click', () => computeCMDPlot({ user: true }));
      document.getElementById('input-iso-ext').addEventListener('change', () => { if (getTaskFilename('isochrone')) computeCMDPlot(); });
      document.getElementById('btn-download-parsec-dat').addEventListener('click', () => {
        if (localStore.engine === 'python') apiOpen('/download-parsec');
        else if (localStore.parsec && localStore.parsec.text) {
          const url = URL.createObjectURL(new Blob([localStore.parsec.text], { type: 'text/plain' }));
          const a = document.createElement('a'); a.href = url; a.download = localStore.parsec.name; document.body.appendChild(a); a.click(); a.remove();
          setTimeout(() => URL.revokeObjectURL(url), 5000);
        } else showToast('No PARSEC file yet: press "Apply Extinction & Fetch PARSEC" first (it must succeed).', 'warning');
      });
      document.getElementById('input-ebv-mode').addEventListener('change', (e) => {
        const mode = e.target.value;
        document.getElementById('per-star-warning').classList.toggle('hidden', mode !== 'per_star');
        if (getTaskFilename('isochrone')) computeCMDPlot();
      });
      document.getElementById('input-ebv').addEventListener('input', (e) => {
        const v = parseFloat(e.target.value);
        document.getElementById('input-av').value = Number.isFinite(v) ? (STARNEST_AV_PER_EBV * v).toFixed(3) : '';
      });
      ['input-dist-mod', 'input-age', 'input-metallicity', 'input-ebv', 'input-iso-plxzp', 'input-use-cluster-dm'].forEach(id => {
        const el = document.getElementById(id);
        if (el) {
          el.addEventListener('change', () => {
            if (getTaskFilename('isochrone')) {
              computeCMDPlot();
              if (id === 'input-age' && getTaskFilename('imf')) computeIMFPlots();
            }
          });
        }
      });

      document.getElementById('btn-download-ext-csv').addEventListener('click', () => {
        downloadCSVRows(localStore.latestExtinction, 'extinction_corrected.csv', '/download-extinction');
      });
      document.getElementById('btn-download-iso-csv').addEventListener('click', () => {
        downloadCSVRows(localStore.latestIsochrone, 'iso_single.csv', '/download-isochrone');
      });

      // Kinematics Button
      document.getElementById('btn-compute-kinematics').addEventListener('click', computeKinematicsPlots);
      ['kin3d-mode', 'kin3d-maxfrac', 'kin3d-depth', 'kin3d-frame', 'kin3d-centre'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener('change', () => { if (window._kinComputed) computeKinematicsPlots(); });
      });

      // IMF Sliders & Buttons
      let imfDebounce = null;
      function triggerDebouncedIMF() {
        if (!getTaskFilename('imf')) return;
        clearTimeout(imfDebounce);
        imfDebounce = setTimeout(() => computeIMFPlots({ silent: true }), 90);
      }

      document.getElementById('check-auto-bins').addEventListener('change', (e) => {
        const slider = document.getElementById('slider-bins');
        slider.classList.toggle('opacity-50', e.target.checked);
        document.getElementById('bins-value').textContent = e.target.checked ? 'Auto' : slider.value;
        triggerDebouncedIMF();
      });
      document.getElementById('slider-bins').addEventListener('input', (e) => {
        const auto = document.getElementById('check-auto-bins');
        if (auto.checked) { auto.checked = false; e.target.classList.remove('opacity-50'); }   // moving the slider = manual bins
        document.getElementById('bins-value').textContent = e.target.value;
        triggerDebouncedIMF();
      });
      document.getElementById('slider-fit-logm').addEventListener('input', (e) => {
        document.getElementById('fit-logm-value').textContent = parseFloat(e.target.value).toFixed(2);
        triggerDebouncedIMF();
      });
      document.getElementById('slider-mass-min').addEventListener('input', (e) => {
        document.getElementById('mass-min-value').textContent = parseFloat(e.target.value).toFixed(2);
        triggerDebouncedIMF();
      });
      document.getElementById('slider-mass-max').addEventListener('input', (e) => {
        document.getElementById('mass-max-value').textContent = parseFloat(e.target.value).toFixed(1);
        triggerDebouncedIMF();
      });
      document.getElementById('check-include-companions').addEventListener('change', triggerDebouncedIMF);
      document.getElementById('btn-compute-imf').addEventListener('click', computeIMFPlots);

      document.getElementById('btn-download-imf-bins').addEventListener('click', () => {
        downloadCSVRows(localStore.latestImfBins, 'IMF_bins.csv', '/download-imf-bins');
      });
      document.getElementById('btn-download-imf-chi2').addEventListener('click', () => {
        downloadCSVRows(localStore.latestImfChi2, 'IMF_chi2.csv', '/download-imf-chi2');
      });
    });
