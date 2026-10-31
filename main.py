"""
main.py -- FastAPI server that runs the REAL ocmass modules for the StarNest front end.
Run:   pip install -r requirements.txt fastapi uvicorn python-multipart
       uvicorn main:app --host 127.0.0.1 --port 8000
The front end probes GET /api/health; if it answers, every analysis call goes here (real GMM-BIC, IRSA
extinction, PARSEC isochrone, isochrone-based masses). Otherwise the page falls back to its in-browser engine.
State is kept in memory PER VISITOR: the page served from '/' sets a cookie (sn_sid); every request is mapped to its own STATE.
Requests without the cookie (page opened from a file, local use) share one 'default' session. Idle sessions expire; at most SESSION_MAX are kept.

TERMINOLOGY
  members_latest.csv : GMM output (PMemb >= cutoff)   |  analysis : members + extinction/DM/colour/mass/type
  MC_Mass            : mass per star; here interpolated on the PARSEC main-sequence (M_G -> Mass), else L(M) fallback

WORKFLOW
  upload-raw -> run-gmm (cuts, BIC-GMM, classify) -> compute-cmd (E(B-V), A_V, DM, PARSEC, masses)
  -> compute-kinematics / compute-imf (both read the analysis table) -> download-* (CSV)
"""
from __future__ import annotations

import contextvars
import io
import math
import sklearn
import os
import re
import time
import uuid
from functools import lru_cache

import numpy as np
import pandas as pd
from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles

import ocmass_extinction_isochrone as EXT
import ocmass_gmm as GMM
import ocmass_imf as IMF
import ocmass_kinematics as KIN

# ===================== GLOBAL CONFIGURATION =====================
API_PREFIX = '/api'
FRONTEND_HTML = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'index.html')   # served at '/'
MEMBERS_NAME = 'members_latest.csv'
AV_PER_EBV = EXT.AV_PER_EBV
EVOLVED_TYPES = {'White Dwarf', 'Red Giant / Supergiant'}      # masses from L(M)/MS-interp invalid for these
IMF_LOG_AGE_SPLIT = 7.0                                        # < 7.0 (10 Myr) -> 'IMF', else 'Present-day MF'
COMPLETENESS_G50, COMPLETENESS_W = 17.2, 0.65                  # ASSUMED curve (flagged assumed=True)
KEEP_ROWS_MAX = 250000
SESSION_COOKIE = 'sn_sid'
SESSION_MAX = 20                                               # sessions kept in memory (oldest idle one is dropped)
SESSION_TTL_S = 6 * 3600                                       # idle time after which a session is dropped
MAX_UPLOAD_MB = 50
# =================================================================

def _new_state() -> dict:
    return dict(files={}, members=None, analysis=None, iso=None, parsec=None, imf_bins=None, imf_chi2=None, name='Cluster')


_SESSIONS: dict = {}                                           # sid -> {'s': state dict, 't': last use}
_sid_var: contextvars.ContextVar = contextvars.ContextVar('sn_sid', default='default')


class _SessionState:
    """dict-like view of the CURRENT visitor's state (the visitor is set per request by the middleware below)."""
    def _d(self) -> dict:
        sid, now = _sid_var.get(), time.time()
        for k in [k for k, v in _SESSIONS.items() if k != 'default' and now - v['t'] > SESSION_TTL_S]:
            _SESSIONS.pop(k, None)
        if sid not in _SESSIONS:
            while len(_SESSIONS) >= SESSION_MAX:
                _SESSIONS.pop(min(_SESSIONS, key=lambda k: _SESSIONS[k]['t']), None)
            _SESSIONS[sid] = dict(s=_new_state(), t=now)
        _SESSIONS[sid]['t'] = now
        return _SESSIONS[sid]['s']

    def __getitem__(self, k): return self._d()[k]
    def __setitem__(self, k, v): self._d()[k] = v
    def get(self, k, default=None): return self._d().get(k, default)
    def __contains__(self, k): return k in self._d()


STATE = _SessionState()

app = FastAPI(title='StarNest OCMass backend')
app.add_middleware(CORSMiddleware, allow_origins=['*'], allow_methods=['*'], allow_headers=['*'], expose_headers=['Content-Disposition'])


@app.middleware('http')
async def _session_middleware(request: Request, call_next):
    # visitor id: cookie (page served by this app) -> X-Session-Id header / ?sid= (page hosted elsewhere, e.g. Streamlit wrapper)
    sid = request.cookies.get(SESSION_COOKIE) or request.headers.get('x-session-id') or request.query_params.get('sid') or ''
    tok = _sid_var.set(sid if re.fullmatch(r'[0-9a-zA-Z_-]{8,64}', sid) else 'default')
    try:
        return await call_next(request)
    finally:
        _sid_var.reset(tok)


# ------------------------------- helpers -------------------------------
def _none(o):
    """NaN/inf -> None recursively, numpy -> python (strict JSON)."""
    if isinstance(o, dict):
        return {k: _none(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_none(v) for v in o]
    if isinstance(o, np.ndarray):
        return _none(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def classify_star(mg: float, bprp0: float) -> str:
    """Same heuristic thresholds as the front end (classifyStarRow)."""
    if not (np.isfinite(mg) and np.isfinite(bprp0)):
        return 'Main Sequence'
    if mg > 10.0 and bprp0 < 0.6:
        return 'White Dwarf'
    if mg < 3.0 and bprp0 > 0.8:
        return 'Red Giant / Supergiant'
    if 4.0 < mg < 3.0 * bprp0 + 2.0:
        return 'Pre-Main Sequence'
    return 'Main Sequence'


def ml_mass(lum: np.ndarray) -> np.ndarray:
    """Piecewise MS mass-luminosity: L=0.23M^2.3 (<0.43), M^4 (0.43-2), 1.4M^3.5 (>2). Fallback only."""
    lum = np.asarray(lum, float)
    out = np.full(lum.shape, np.nan)
    ok = lum > 0
    lo, mid, hi = ok & (lum < 0.033), ok & (lum >= 0.033) & (lum < 16), ok & (lum >= 16)
    out[lo] = (lum[lo] / 0.23) ** (1 / 2.3)
    out[mid] = lum[mid] ** 0.25
    out[hi] = (lum[hi] / 1.4) ** (1 / 3.5)
    return out


def approx_mass(mg, color):
    mbol = mg + (-0.012 - 0.25 * color - 0.15 * color ** 2)
    return ml_mass(10 ** (0.4 * (4.74 - mbol)))


def derive_columns(df: pd.DataFrame, mg_col='M_G_corr', col_col='bp_rp_corr', iso: pd.DataFrame | None = None) -> pd.DataFrame:
    """Star_Type (+ MC_Mass: isochrone-MS interpolation, else L(M) fallback). User-supplied values are kept."""
    d = df.copy()
    mg, c = d[mg_col].to_numpy(float), d[col_col].to_numpy(float)
    if 'Star_Type' in d and d.get('Star_Type_src', pd.Series(['user'] * len(d))).ne('auto').all() and d['Star_Type'].isin(
            ['Main Sequence', 'Pre-Main Sequence', 'White Dwarf', 'Red Giant / Supergiant']).all():
        d['Star_Type_src'] = 'user'
    else:
        d['Star_Type'] = [classify_star(a, b) for a, b in zip(mg, c)]
        d['Star_Type_src'] = 'auto'
    d['type of stars'] = d['Star_Type']
    given = next((k for k in ('MC_Mass', 'Mass', 'mass') if k in d.columns), None)
    if given is not None and not (d.get('MC_Mass_src', pd.Series(['user'] * len(d))) == 'approx_ML').any() \
            and not (d.get('MC_Mass_src', pd.Series([''] * len(d))) == 'isochrone').any():
        d['MC_Mass'] = pd.to_numeric(d[given], errors='coerce')
        d['MC_Mass_src'] = 'user'
        return d
    mass = np.full(len(d), np.nan)
    src = np.array(['approx_ML'] * len(d), dtype=object)
    if iso is not None and len(iso) and 'label' in iso:
        ms = iso[iso['label'].isin([0, 1])].sort_values('Gmag')
        if len(ms) > 3:
            inside = (mg >= ms['Gmag'].min()) & (mg <= ms['Gmag'].max())
            mass[inside] = np.interp(mg[inside], ms['Gmag'].to_numpy(), ms['Mass'].to_numpy())
            src[inside] = 'isochrone'
    miss = ~np.isfinite(mass)
    mass[miss] = approx_mass(mg[miss], c[miss])
    d['MC_Mass'] = np.round(mass, 4)
    d['MC_Mass_src'] = src
    return d


def get_rows(name: str) -> pd.DataFrame:
    if name == MEMBERS_NAME:
        if STATE['members'] is None:
            raise HTTPException(400, 'No member stars loaded. Run GMM Analysis or upload a custom CSV first.')
        return STATE['members']
    if name not in STATE['files']:
        raise HTTPException(400, f"File '{name}' not found on the server. Upload it again.")
    return STATE['files'][name]


def analysis_rows(name: str) -> pd.DataFrame:
    a = STATE['analysis']
    if a is not None and a['src'] == name:
        return a['df']
    df = get_rows(name)
    if 'M_G_corr' not in df or 'bp_rp_corr' not in df:
        d, _ = EXT.apply_extinction(df, mode='fixed', ebv_fixed=0.0)
        df = d
    return derive_columns(df)


def _csv(df, name):
    if df is None or not len(df):
        raise HTTPException(404, 'Nothing to download yet: run the analysis first.')
    return Response(df.to_csv(index=False), media_type='text/csv', headers={'Content-Disposition': f'attachment; filename={name}'})


# ------------------------------- endpoints -------------------------------
@app.get('/_stcore/health')
@app.get('/healthz')
def platform_health():
    """Health probe for hosts that check the Streamlit path (Streamlit Community Cloud runs `streamlit run main.py` and serves this ASGI app)."""
    return Response('ok', media_type='text/plain')


@app.get(API_PREFIX + '/health')
def health():
    return dict(status='ok', engine='python', modules=['ocmass_gmm', 'ocmass_extinction_isochrone', 'ocmass_kinematics', 'ocmass_imf'])


@app.post(API_PREFIX + '/upload-raw')
async def upload_raw(file: UploadFile = File(...)):
    name = re.sub(r'[^A-Za-z0-9._-]', '_', os.path.basename(file.filename or 'upload.csv'))
    raw = await file.read()
    if len(raw) > MAX_UPLOAD_MB * 1024 * 1024:
        raise HTTPException(413, f'File larger than {MAX_UPLOAD_MB} MB.')
    try:
        df = GMM.load_raw(io.BytesIO(raw))
    except Exception as e:
        raise HTTPException(400, f'Could not read CSV: {e}')
    if len(df) > KEEP_ROWS_MAX:
        raise HTTPException(413, f'File too large ({len(df)} rows > {KEEP_ROWS_MAX}).')
    STATE['files'][name] = df
    return dict(filename=name, n_rows=len(df))


@app.post(API_PREFIX + '/run-gmm')
def run_gmm(p: dict):
    raw = get_rows(p['filename'])
    miss = GMM.validate_columns(raw)
    if miss:
        raise HTTPException(400, f'Missing required columns: {miss}')
    q = dict(RUWE_MAX=p.get('ruwe_max', 1.6), E_PM_MAX=p.get('epmra_max', 0.1), E_PMRA_MAX=p.get('epmra_max', 0.1),
             E_PMDE_MAX=p.get('epmde_max', 0.1), PLX_SNR_MIN=p.get('plx_err_min', 5.0),
             GOFAL_MAX=p.get('gofal_max', 3.0), NPER_MIN=p.get('nper_min', 10.0))
    clean = GMM.apply_quality_cuts(raw, q)
    try:
        res = GMM.run_gmm(clean, distance_lit_pc=p.get('distance_lit', 341.0), half_width_pc=p.get('half_width', 195.5),
                          cutoff=p.get('cutoff', 0.5), k=p.get('k'), plx_zp=p.get('plx_zp', 0.0) or 0.0,
                          optimizer=p.get('gmm_opt', GMM.OPTIMIZER))
    except ValueError as e:
        raise HTTPException(400, str(e))
    zp = p.get('plx_zp', 0.0) or 0.0
    optimizer_used = p.get('gmm_opt', GMM.OPTIMIZER)
    dist = res.summary['distance_from_members_pc']
    dm = 5 * math.log10(dist) - 5 if math.isfinite(dist) else 7.66
    ebv = float(p.get('ebv_ui') or 0.0)

    def enrich(df):
        d = df.copy()
        av = AV_PER_EBV * ebv
        d['A_G'], d['A_BP'], d['A_RP'] = EXT.AG_AV * av, EXT.ABP_AV * av, EXT.ARP_AV * av
        d['M_G_corr'] = d[GMM.COL['gmag']] - dm - d['A_G']
        if GMM.COL['bp'] in d and GMM.COL['rp'] in d:
            d['bp_rp_corr'] = (d[GMM.COL['bp']] - d['A_BP']) - (d[GMM.COL['rp']] - d['A_RP'])
        else:
            d['bp_rp_corr'] = d[GMM.COL['color']] - (d['A_BP'] - d['A_RP'])
        d['bprp_corr_ok'] = d['bp_rp_corr'].between(*EXT.BPRP0_RANGE)
        return derive_columns(d)

    members, field = enrich(res.members), res.non_members
    STATE['members'], STATE['analysis'], STATE['iso'], STATE['parsec'] = members, None, None, None
    STATE['name'] = os.path.splitext(p['filename'])[0]
    ra, de, col, g = GMM.COL['ra'], GMM.COL['dec'], GMM.COL['color'], GMM.COL['gmag']
    colors = {'Main Sequence': '#00FFFF', 'Pre-Main Sequence': '#BA8E23', 'White Dwarf': '#000000', 'Red Giant / Supergiant': '#FF0000'}
    s = res.summary
    out = dict(
        ra_dec=dict(members=dict(ra=members.get(ra, pd.Series(dtype=float)).tolist(), dec=members.get(de, pd.Series(dtype=float)).tolist(),
                                 color=[colors[t] for t in members['Star_Type']], star_type=members['Star_Type'].tolist(),
                                 pmemb=members['PMemb'].round(6).tolist()),
                    field=dict(ra=field.get(ra, pd.Series(dtype=float)).tolist(), dec=field.get(de, pd.Series(dtype=float)).tolist())),
        vpd=dict(members=dict(pmra=members['pmRA'].tolist(), pmdec=members['pmDE'].tolist(), color=[colors[t] for t in members['Star_Type']]),
                 field=dict(pmra=field['pmRA'].tolist(), pmdec=field['pmDE'].tolist())),
        plx_Cmd=dict(members_plx=members['Plx'].tolist(), field_plx=field['Plx'].tolist(),
                     members_color=members[col].tolist(), members_gmag=members[g].tolist(),
                     field_color=field[col].tolist(), field_gmag=field[g].tolist()),
        stats=dict(n_input=len(clean), total=s['n_in_window'], n_fit=s['n_fit'], k=s['k'], bic={str(k): v for k, v in res.bic.items()},
                   cluster_component=res.cluster_component, components=res.components, k_table=res.k_table,
                   methods=dict(engine='scikit-learn ' + sklearn.__version__ + ' GaussianMixture(full)', optimizer=optimizer_used, tol=(GMM.GMM_TOL if optimizer_used == 'converged' else 1e-3), max_iter=(GMM.GMM_MAX_ITER if optimizer_used == 'converged' else 100), n_init=GMM.N_INIT, reg_covar=1e-6, random_state=GMM.RANDOM_STATE, k_mode=('fixed' if p.get('k') else 'BIC over 2..6'), cutoff=p.get('cutoff', 0.5), distance_lit=p.get('distance_lit', 341.0), half_width=p.get('half_width', 195.5), plx_zp=zp), members=s['n_members'], field=s['n_non_members'], threshold=p.get('cutoff', 0.5),
                   mss_metric=round(s['mss'], 4), median_plx=round(s['median_plx'], 4) if math.isfinite(s['median_plx']) else 0,
                   distance_from_members_pc=round(dist, 1) if math.isfinite(dist) else 0),
        classification_counts=members['Star_Type'].value_counts().to_dict())
    return _none(out)


@app.post(API_PREFIX + '/compute-cmd')
def compute_cmd(p: dict):
    """Workflow: members CSV -> IRSA E(B-V) (mean) -> reference-code extinction columns -> CMD 3.8 PARSEC (Gaia EDR3/DR3, A_V forwarded) -> fit.
    iso_extinction='reddened' (default): the MEAN A_V is sent to CMD, stars are plotted observed (G-DM, BP-RP) -> extinction applied once, in the isochrone.
    iso_extinction='intrinsic': CMD gets A_V=0 and the stars are de-reddened (reference code) -> extinction applied once, in the stars.
    allow_network=False (automatic refreshes): only cached IRSA / PARSEC values are used, nothing is requested from the shared servers."""
    members = get_rows(p['members_filename'])
    age, mh = float(p['age']), float(p.get('metallicity') or 0.0)
    mode = p.get('ebv_mode') or 'center'
    handling = p.get('iso_extinction') or 'reddened'
    allow_net = p.get('allow_network', True) is not False
    if handling not in ('reddened', 'intrinsic'):
        raise HTTPException(400, f"iso_extinction must be 'reddened' or 'intrinsic', got {handling!r}")
    if age < EXT.LOGAGE_FLOOR:
        raise HTTPException(400, f'logAge {age} below PARSEC floor {EXT.LOGAGE_FLOOR}')
    if mode == 'fixed' and p.get('ebv') is None:
        raise HTTPException(400, "mode='fixed' requires ebv_fixed")
    wf = dict(n_members=len(members), irsa=dict(), extinction=dict(), parsec=dict(), fit=dict())
    notes = []
    kw = dict(plx_zp=p.get('plx_zp') or 0.0, use_cluster_dm=p.get('use_cluster_dm') is not False)
    try:
        try:
            d, meta = EXT.apply_extinction(members, mode=mode, ebv_fixed=p.get('ebv'), cache_only=not allow_net, **kw)
            wf['irsa'] = dict(state='done' if mode != 'fixed' else 'manual', mode=mode, queries=meta['irsa_queries'], failed=meta['irsa_failed'],
                              dropped_stars=meta['irsa_dropped_stars'])
        except EXT.NetworkNotAllowed:             # automatic refresh before the user pressed Apply: use the typed E(B-V), ask nobody
            if p.get('ebv') is None:
                raise ValueError('Press "Apply Extinction & Fetch PARSEC" to query IRSA (or choose mode fixed and type an E(B-V)).')
            d, meta = EXT.apply_extinction(members, mode='fixed', ebv_fixed=p['ebv'], **kw)
            meta['mode'] = 'fixed (IRSA not queried yet)'
            wf['irsa'] = dict(state='pending', mode=mode)
            notes.append('IRSA / PARSEC have not been queried yet: showing the typed E(B-V). Press "Apply Extinction & Fetch PARSEC".')
    except RuntimeError as e:
        raise HTTPException(502, f'{e}. Switch E(B-V) mode to fixed and type a value.')
    except ValueError as e:
        raise HTTPException(400, str(e))
    av_mean = meta['A_V_mean']
    wf['extinction'] = dict(E_BV_mean=meta['E_BV_mean'], E_BV_std=meta['E_BV_std'], A_V_mean=av_mean, A_G=EXT.AG_AV * av_mean,
                            A_BP=EXT.ABP_AV * av_mean, A_RP=EXT.ARP_AV * av_mean, AV_PER_EBV=EXT.AV_PER_EBV, DM=meta['DM'], d_pc=meta['d_pc'])
    av_send = round(av_mean, 3) if handling == 'reddened' else 0.0

    iso0 = iso_plot = raw_used = None
    check = None
    try:
        r0 = EXT.get_parsec(round(age, 4), round(mh, 4), 0.0, cache_only=not allow_net)             # intrinsic: masses, de-reddened plane, shift check
        iso0 = r0['iso']
        rU = r0
        if handling == 'reddened' and av_send > 0.0005:
            rU = EXT.get_parsec(round(age, 4), round(mh, 4), av_send, cache_only=not allow_net)
            check = EXT.check_reddening_shift(iso0, rU['iso'], av_send)
            if check['ok'] is False:
                notes.append('WARNING: ' + check['message'])
        iso_plot, raw_used = rU['iso'], rU['raw']
        wf['parsec'] = dict(state='done', rows=int(len(iso_plot)), kb=round(len(raw_used) / 1024, 1), logAge=age, MH=mh, A_V_sent=av_send,
                            system='Gaia EDR3 = DR3 photometry (Vega)', service='CMD 3.8', tls_verified=bool(rU['info'].get('ssl_verified')),
                            file_url=rU['info'].get('source_url'), check=check)
    except EXT.NetworkNotAllowed:
        wf['parsec'] = dict(state='pending')
        notes.append('PARSEC not fetched yet: press "Apply Extinction & Fetch PARSEC" (one request per logAge / [M/H] / A_V).')
    except (EXT.CMDRequestError, RuntimeError, KeyError, ValueError) as e:        # CMD down / rejected: keep going, tell the user
        wf['parsec'] = dict(state='failed', error=str(e))
        notes.append(f'PARSEC download failed ({e}); isochrone not drawn and masses use the L(M) fallback.')
    STATE['iso'] = iso_plot
    STATE['parsec'] = dict(raw=raw_used, name=f'PARSEC_gaiaEDR3_logAge{age:.2f}_MH{mh:+.2f}_AV{av_send:.3f}.dat') if raw_used else None
    d = derive_columns(d, iso=iso0)
    STATE['analysis'] = dict(src=p['members_filename'], df=d)

    reddened = handling == 'reddened'
    cx, my = ('bp_rp_obs', 'M_G_obs') if reddened else ('bp_rp_corr', 'M_G_corr')
    if iso_plot is not None:
        i = EXT.iso_for_plot(iso_plot)
        icol, imag = (i['G_BPmag'] - i['G_RPmag']).to_numpy(float), i['Gmag'].to_numpy(float)
        isochrone = dict(color=icol.tolist(), mag=imag.tolist(), rows=[], source='parsec')
        okm = d['bprp_corr_ok'].to_numpy(bool) & (d['Star_Type'] != 'White Dwarf').to_numpy(bool)
        wf['fit'] = dict(state='done', plane='reddened (observed colour, G-DM)' if reddened else 'de-reddened',
                         **EXT.fit_quality(d.loc[okm, cx], d.loc[okm, my], icol, imag))
    else:
        isochrone = dict(color=[], mag=[], rows=[], source='none')
        wf['fit'] = dict(state='pending')
    meta.update(engine='python', note=' '.join(notes), workflow=wf, plot_space='reddened' if reddened else 'intrinsic', A_V_sent=av_send)
    out = dict(empirical=dict(color=d[cx].tolist(), mag=d[my].tolist(), star_type=d['Star_Type'].tolist(),
                              bprp_corr_ok=d['bprp_corr_ok'].tolist()),
               isochrone=isochrone, meta=meta, classification_counts=d['Star_Type'].value_counts().to_dict(),
               total_classified=len(d))
    return _none(out)


@app.post(API_PREFIX + '/compute-kinematics')
def compute_kinematics(p: dict):
    """3-D map + PM vectors: the maths lives in ocmass_kinematics.py (this route only reads the request and shapes the reply)."""
    base = analysis_rows(p['filename'])
    try:
        res = KIN.compute_kinematics(base, mode=p.get('kin3d_mode') or KIN.KIN3D_MODE,
                                     max_frac=float(p.get('kin3d_max_frac') or KIN.KIN3D_MAX_FRAC),
                                     depth_pc=float(p.get('kin3d_depth') or KIN.KIN3D_DEPTH_PC),
                                     frame=p.get('kin3d_frame') or KIN.KIN3D_FRAME, centre=bool(p.get('kin3d_centre', KIN.KIN3D_CENTRE)))
    except ValueError as e:
        raise HTTPException(400, str(e))
    v = res['vectors_2d']; v.pop('res_pmra'); v.pop('res_pmde')            # raw residuals are for the CSV helper, not the page
    name = re.sub(r'newcorrected.*$', '', re.sub(r'^Extinctioncorr_', '', STATE['name'], flags=re.I), flags=re.I).split('_')[0]
    return _none(dict(cluster_name=name, **res))


def _wls(x, y, sig):
    """Poisson-weighted line fit (matches the front end's wlinregressJS)."""
    x, y, sig = map(np.asarray, (x, y, sig))
    w = np.where(sig > 0, 1 / np.where(sig > 0, sig, 1) ** 2, 0.0)
    n = len(x)
    mx, my = np.sum(w * x) / w.sum(), np.sum(w * y) / w.sum()
    sxx, sxy, syy = np.sum(w * (x - mx) ** 2), np.sum(w * (x - mx) * (y - my)), np.sum(w * (y - my) ** 2)
    slope = sxy / sxx if sxx > 0 else 0.0
    icpt = my - slope * mx
    chi2r = np.sum(w * (y - icpt - slope * x) ** 2) / (n - 2) if n > 2 else 1.0
    se = math.sqrt(max(1.0, chi2r) / sxx) if sxx > 0 else 0.0
    r2 = sxy ** 2 / (sxx * syy) if sxx > 0 and syy > 0 else 0.0
    return slope, icpt, r2, se


@app.post(API_PREFIX + '/compute-imf')
def compute_imf(p: dict):
    d = analysis_rows(p['filename'])
    d = d[~d['Star_Type'].isin(EVOLVED_TYPES)] if 'Star_Type' in d else d
    n_excl = 0 if 'Star_Type' not in d else len(analysis_rows(p['filename'])) - len(d)
    lo, hi = float(p.get('mass_min') or IMF.M_MIN_FLOOR), float(p.get('mass_max') or IMF.M_MAX_CAP)
    prim = d['MC_Mass'].to_numpy(float) if 'MC_Mass' in d else np.array([])
    prim = prim[np.isfinite(prim) & (prim >= lo) & (prim <= hi)]
    comp = None
    if p.get('include_companions') and 'MC_CompMass' in d:
        c = d['MC_CompMass'].to_numpy(float)
        comp = c[np.isfinite(c) & (c >= lo) & (c <= hi)]
    try:
        r = IMF.compute_imf(prim, comp, fit_min_logm=float(p.get('fit_min_logM', IMF.FIT_MIN_LOGM)), n_bins=p.get('bins'))
    except ValueError as e:
        raise HTTPException(400, str(e))
    STATE['imf_bins'], STATE['imf_chi2'] = r.bins, r.chi2
    b = r.bins
    cb, ib = b[(b['count'] > 0) & b['complete']], b[(b['count'] > 0) & ~b['complete']]
    pack = lambda t: dict(logM_center=t['logM_center'].tolist(), log_dN=t['log_dN'].tolist(), log_err=t['log_err'].tolist(), count=t['count'].tolist())
    slope, icpt, r2, se = _wls(cb['logM_center'], cb['log_dN'], cb['log_err'])
    half = r.meta['bin_width_dex'] / 2          # line spans the whole fitted range (bin edges), like the page
    xf = [float(cb['logM_center'].min()) - half, float(cb['logM_center'].max()) + half]
    age = float(p.get('log_age') or 7.5)
    g = np.arange(8, 20.01, 0.8)
    out = dict(imf=dict(
        mf_kind='IMF' if age < IMF_LOG_AGE_SPLIT else 'Present-day MF',
        meta=dict(n_masses=r.meta['n_masses'], n_bins=r.meta['n_bins'], bin_width_dex=round(r.meta['bin_width_dex'], 4),
                  n_fit_bins=r.meta['n_fit_bins'], n_in_fit=r.meta['n_in_fit'], fit_min_logm=r.meta['fit_min_logm'],
                  n_excluded_evolved=n_excl, m_range=[round(r.meta['m_range'][0], 3), round(r.meta['m_range'][1], 2)]),
        bins_table=r.bins.to_dict('records'),
        chi2_table=r.chi2.sort_values('red_chi2').round(3).to_dict('records'),
        complete_bins=pack(cb), incomplete_bins=pack(ib), slope=slope, slope_err=se, r_squared=r2, ms_count=r.meta['n_masses'],
        fit_line=dict(x=xf, y=[slope * xf[0] + icpt, slope * xf[1] + icpt], slope=slope, slope_err=se, n_bins=r.meta['n_fit_bins']),
        theoretical=dict(log_m=r.curves['logM'], salpeter=np.log10(np.maximum(r.curves['Salpeter'], 1e-12)),
                         kroupa=np.log10(np.maximum(r.curves['Kroupa'], 1e-12)), chabrier=np.log10(np.maximum(r.curves['Chabrier'], 1e-12)))),
        completeness=dict(mag=g, fraction=1 / (1 + np.exp((g - COMPLETENESS_G50) / COMPLETENESS_W)), assumed=True))
    return _none(out)


@app.get(API_PREFIX + '/download-members')
def dl_members():
    return _csv(STATE['members'], 'filtered_classified_members.csv')


@app.get(API_PREFIX + '/download-extinction')
def dl_ext():
    return _csv(STATE['analysis']['df'] if STATE['analysis'] else None, 'extinction_corrected.csv')


@app.get(API_PREFIX + '/download-isochrone')
def dl_iso():
    return _csv(STATE['iso'], 'iso_single.csv')


@app.get(API_PREFIX + '/download-parsec')
def dl_parsec():
    """The PARSEC file exactly as CMD 3.8 returned it (the isochrone used in the plot), preceded by '#' provenance lines."""
    pz = STATE['parsec']
    if not pz or not pz.get('raw'):
        raise HTTPException(404, 'No PARSEC file yet: press Apply Extinction & Fetch PARSEC first.')
    return Response(pz['raw'], media_type='text/plain', headers={'Content-Disposition': f"attachment; filename={pz['name']}"})


@app.get(API_PREFIX + '/download-imf-bins')
def dl_imf_bins():
    return _csv(STATE['imf_bins'], 'IMF_bins.csv')


@app.get(API_PREFIX + '/download-imf-chi2')
def dl_imf_chi2():
    return _csv(STATE['imf_chi2'], 'IMF_chi2.csv')


# static files of the front end (JS / CSS / team photos) so that http://localhost:8000/ shows the complete page
for _d in ('assets', 'images'):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), _d)
    if os.path.isdir(_p):
        app.mount('/' + _d, StaticFiles(directory=_p), name=_d)


@app.get('/')
def index(request: Request):
    if os.path.exists(FRONTEND_HTML):
        resp = FileResponse(FRONTEND_HTML)
        if not request.cookies.get(SESSION_COOKIE):                # a new visitor gets a private session
            resp.set_cookie(SESSION_COOKIE, uuid.uuid4().hex, max_age=SESSION_TTL_S, httponly=True, samesite='lax')
        return resp
    raise HTTPException(404, 'Put the front end next to backend.py as index.html to serve it here.')
