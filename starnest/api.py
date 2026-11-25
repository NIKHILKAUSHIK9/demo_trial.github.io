"""
main.py -- FastAPI server that runs the REAL StarNest modules for the StarNest front end.
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
import threading
import time
import uuid

import numpy as np
import pandas as pd
from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles

from . import binaries as BIN
from . import dynamics as DYN
from . import extinction_isochrone as EXT
from . import gmm as GMM
from . import imf as IMF
from . import kinematics as KIN
from .schemas import (BinariesReq, ComputeCmdReq, DynamicsReq, GmmDiagReq, ImfReq, IsoFitReq, KinematicsReq, RunGmmReq,
                     VelocitiesReq)

# ===================== GLOBAL CONFIGURATION =====================
API_PREFIX = '/api'
ROOT = os.environ.get('STARNEST_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root: index.html, assets/, images/
FRONTEND_HTML = os.path.join(ROOT, 'index.html')                                          # served at '/'
MEMBERS_NAME = 'members_latest.csv'
AV_PER_EBV = EXT.AV_PER_EBV
EVOLVED_TYPES = {'White Dwarf', 'Red Giant / Supergiant'}      # masses from L(M)/MS-interp invalid for these
COMPLETENESS_G50, COMPLETENESS_W = 17.2, 0.65                  # ASSUMED logistic curve: fallback only, used when no measured curve exists (assumed=True)
MASS_QUALITY_OK = ('isochrone', 'user')                          # MC_Mass_src values whose masses are trusted (approx_ML = crude L(M) fallback)
KEEP_ROWS_MAX = 250000
SESSION_COOKIE = 'sn_sid'
SESSION_MAX = 20                                               # sessions kept in memory (oldest idle one is dropped)
SESSION_TTL_S = 6 * 3600                                       # idle time after which a session is dropped
MAX_UPLOAD_MB = 50
# =================================================================

def _new_state() -> dict:
    return dict(files={}, members=None, analysis=None, iso=None, parsec=None, imf_bins=None, imf_chi2=None, name='Cluster',
                completeness=None, iso0=None, ext=None, gmm_diag=None, gmm_res=None, binaries=None, velocities=None, dynamics=None, isofit=None, imf_res=None, iso_stale=None)


_SESSIONS: dict = {}                                           # sid -> {'s': state dict, 't': last use}
_SESSIONS_LOCK = threading.RLock()                             # sync routes run in a thread pool: prune / evict / create must be atomic
_sid_var: contextvars.ContextVar = contextvars.ContextVar('sn_sid', default='default')


class _SessionState:
    """dict-like view of the CURRENT visitor's state (the visitor is set per request by the middleware below)."""
    def _d(self) -> dict:
        sid, now = _sid_var.get(), time.time()
        with _SESSIONS_LOCK:                                   # one thread at a time may prune, evict or create (no KeyError / lost session)
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


def _stale_isochrone_message() -> str | None:
    """TERMINOLOGY  stale isochrone: the analysis table's masses came from an isochrone fetched for one (logAge, [M/H]), but the inputs on the
    Isochrone tab have since changed and the new isochrone has NOT been fetched (a no-network refresh cannot fetch it).
    WORKFLOW  1) read STATE['iso_stale'] (set by compute_cmd); 2) None -> fine; 3) else return the sentence every dependent route shows."""
    st = STATE['iso_stale']
    if not st:
        return None
    return (f"The isochrone masses are stale: they were assigned with logAge {st['logAge']}, [M/H] {st['mh']} but the inputs are now logAge "
            f"{st['requested_logAge']}, [M/H] {st['requested_mh']} and that isochrone has not been fetched. Press 'Apply Extinction & Fetch PARSEC' "
            f"on the Isochrone tab, then run this step again.")

app = FastAPI(title='StarNest backend')
app.add_middleware(CORSMiddleware, allow_origins=['*'], allow_methods=['*'], allow_headers=['*'], expose_headers=['Content-Disposition'])


@app.exception_handler(RequestValidationError)
async def _validation_400(request: Request, exc: RequestValidationError):
    """Bad request body -> HTTP 400 whose `detail` names each offending field (front end shows string details verbatim)."""
    parts = []
    for e in exc.errors():
        loc = [str(x) for x in e.get('loc', ()) if x not in ('body', 'query', 'path')]
        parts.append(f"{'.'.join(loc) or 'body'}: {e.get('msg', 'invalid')}")
    return JSONResponse(status_code=400, content=dict(detail='Invalid request: ' + '; '.join(parts), fields=[p.split(':')[0] for p in parts]))


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
def _pv(p: dict, key: str, default):
    """Request value with a default that also covers an explicit JSON null (the browser sends null for a blank number box: JSON.stringify(NaN))."""
    v = p.get(key)
    return default if v is None else v


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
        d['Mass_quality'] = 'ok'
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
    # quality flag: 'ok' = mass read off the PARSEC main sequence; 'approx' = L(M) fallback (outside the isochrone / no isochrone),
    # which ignores age, metallicity and PMS contraction and is NOT fit for an IMF slope
    d['Mass_quality'] = np.where(np.isin(src, MASS_QUALITY_OK), 'ok', 'approx')
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
    return dict(status='ok', engine='python', modules=['gmm', 'extinction_isochrone', 'kinematics', 'imf'])


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
def run_gmm(req: RunGmmReq):
    p = req.model_dump(exclude_unset=True)
    raw = get_rows(p['filename'])
    miss = GMM.validate_columns(raw)
    if miss:
        raise HTTPException(400, f'Missing required columns: {miss}')
    q = dict(RUWE_MAX=_pv(p, 'ruwe_max', 1.6), E_PM_MAX=_pv(p, 'epmra_max', 0.1), E_PMRA_MAX=_pv(p, 'epmra_max', 0.1),
             E_PMDE_MAX=_pv(p, 'epmde_max', 0.1), PLX_SNR_MIN=_pv(p, 'plx_err_min', 5.0),
             GOFAL_MAX=_pv(p, 'gofal_max', 3.0), NPER_MIN=_pv(p, 'nper_min', 10.0))
    clean = GMM.apply_quality_cuts(raw, q)
    try:
        res = GMM.run_gmm(clean, distance_lit_pc=_pv(p, 'distance_lit', 341.0), half_width_pc=_pv(p, 'half_width', 195.5),
                          cutoff=_pv(p, 'cutoff', 0.5), k=p.get('k'), plx_zp=p.get('plx_zp', 0.0) or 0.0,
                          optimizer=_pv(p, 'gmm_opt', GMM.OPTIMIZER), cluster_rule=p.get('cluster_rule') or GMM.CLUSTER_RULE,
                          lit_pm=p.get('lit_pm') or GMM.LIT_PM, lit_plx=p.get('lit_plx') or GMM.LIT_PLX)
    except ValueError as e:
        raise HTTPException(400, str(e))
    zp = p.get('plx_zp', 0.0) or 0.0
    optimizer_used = _pv(p, 'gmm_opt', GMM.OPTIMIZER)
    dist = res.summary['distance_from_members_pc']
    dm = 5 * math.log10(dist) - 5 if math.isfinite(dist) else 7.66
    ebv_raw = p.get('ebv_ui')
    ebv_stamp = float(ebv_raw) if ebv_raw not in (None, '') and float(ebv_raw) > 0 else None   # only an EXPLICIT value is used and recorded

    def enrich(df):
        """Stage 1 must not invent an extinction. With no explicit E(B-V): observed plane (M_G_obs, bp_rp_obs) is used only for the
        PROVISIONAL Star_Type colours, and no M_G_corr/bp_rp_corr/MC_Mass columns are written. With an explicit E(B-V) the old behaviour
        is kept AND the value is stamped in column E_BV_used so a downstream reader can tell what was assumed."""
        d = df.copy()
        if ebv_stamp is None:
            d['M_G_obs'] = d[GMM.COL['gmag']] - dm
            d['bp_rp_obs'] = d[GMM.COL['color']]
            d['Star_Type'] = [classify_star(a, b) for a, b in zip(d['M_G_obs'], d['bp_rp_obs'])]   # reddened plane: provisional only
            d['Star_Type_src'], d['Star_Type_plane'] = 'auto', 'observed (provisional)'
            d['type of stars'] = d['Star_Type']
            return d
        av = AV_PER_EBV * ebv_stamp
        d['E_BV_used'] = ebv_stamp
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
    STATE['iso0'], STATE['ext'], STATE['iso_stale'] = None, None, None
    STATE['gmm_res'], STATE['binaries'], STATE['velocities'], STATE['dynamics'], STATE['isofit'], STATE['imf_res'] = res, None, None, None, None, None
    # measured completeness: quality-cut survivors / raw stars, both inside the same distance window (see gmm.completeness_curve)
    raw_win = GMM.distance_window(raw, _pv(p, 'distance_lit', 341.0), _pv(p, 'half_width', 195.5), zp)
    clean_win = res.working
    if GMM.COMPLETE_PREFILTER and p.get('complete_prefilter', True) is not False and len(res.members) >= 5:
        # same cluster-compatible subset on both sides of the ratio: raw (before quality cuts) vs cleaned (after)
        raw_win = GMM.cluster_compatible(raw_win, res.members)
        clean_win = GMM.cluster_compatible(clean_win, res.members)
    STATE['completeness'] = GMM.completeness_curve(raw_win, clean_win, frac=p.get('complete_frac'))
    STATE['completeness']['prefiltered'] = bool(GMM.COMPLETE_PREFILTER and p.get('complete_prefilter', True) is not False)
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
                   methods=dict(engine='scikit-learn ' + sklearn.__version__ + ' GaussianMixture(full)', optimizer=optimizer_used, tol=(GMM.GMM_TOL if optimizer_used == 'converged' else 1e-3), max_iter=(GMM.GMM_MAX_ITER if optimizer_used == 'converged' else 100), n_init=GMM.N_INIT, reg_covar=1e-6, random_state=GMM.RANDOM_STATE, k_mode=('fixed' if p.get('k') else 'BIC over 2..6'), cutoff=_pv(p, 'cutoff', 0.5), distance_lit=_pv(p, 'distance_lit', 341.0), half_width=_pv(p, 'half_width', 195.5), plx_zp=zp), members=s['n_members'], field=s['n_non_members'], threshold=_pv(p, 'cutoff', 0.5),
                   mss_metric=round(s['mss'], 4), median_plx=round(s['median_plx'], 4) if math.isfinite(s['median_plx']) else 0,
                   distance_from_members_pc=round(dist, 1) if math.isfinite(dist) else 0),
        classification_counts=members['Star_Type'].value_counts().to_dict())
    return _none(out)


@app.post(API_PREFIX + '/compute-cmd')
def compute_cmd(req: ComputeCmdReq):
    """Workflow: members CSV -> IRSA E(B-V) (mean) -> reference-code extinction columns -> CMD 3.8 PARSEC (Gaia EDR3/DR3, A_V forwarded) -> fit.
    iso_extinction='reddened' (default): the MEAN A_V is sent to CMD, stars are plotted observed (G-DM, BP-RP) -> extinction applied once, in the isochrone.
    iso_extinction='intrinsic': CMD gets A_V=0 and the stars are de-reddened (reference code) -> extinction applied once, in the stars.
    allow_network=False (automatic refreshes): only cached IRSA / PARSEC values are used, nothing is requested from the shared servers."""
    p = req.model_dump(exclude_unset=True)
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
                            system='Gaia EDR3 = DR3 photometry (Vega)', service=rU['info'].get('service', 'CMD 3.8'), tls_verified=bool(rU['info'].get('ssl_verified')),
                            file_url=rU['info'].get('source_url'), check=check, source=rU['info'].get('source', 'cmd'),
                            interpolated=rU['info'].get('interpolated'), grid_source=rU['info'].get('grid_source'))
        if rU['info'].get('grid_note'):
            notes.append(f"Local PARSEC grid not used ({rU['info']['grid_note']}); CMD was queried instead.")
        if check is not None and rU['info'].get('source') == 'local_grid':
            check['independent'] = False                         # the shift was applied by us with the same coefficients: sanity check only
            check['message'] = 'reddening applied analytically from the local grid (coefficients stored with the grid); sanity check only, not independent'
    except EXT.NetworkNotAllowed:
        wf['parsec'] = dict(state='pending')
        notes.append('PARSEC not fetched yet: press "Apply Extinction & Fetch PARSEC" (one request per logAge / [M/H] / A_V).')
    except (EXT.CMDRequestError, RuntimeError, KeyError, ValueError) as e:        # CMD down / rejected: keep going, tell the user
        wf['parsec'] = dict(state='failed', error=str(e))
        notes.append(f'PARSEC download failed ({e}); isochrone not drawn and masses use the L(M) fallback.')
    d = derive_columns(d, iso=iso0)
    d['E_BV_used'] = meta['E_BV_mean']                      # the E(B-V) actually behind M_G_corr / bp_rp_corr in this table
    prior, prior_ext = STATE['analysis'], STATE['ext']
    # A refresh that could NOT fetch this (logAge, [M/H]) (no network on an automatic refresh, or CMD down) must not replace a table whose masses
    # came from a real isochrone with one that only has the L(M) fallback: keep the good table, flag it stale, and let dependent steps say so.
    keep_prior = (iso0 is None and prior is not None and prior['src'] == p['members_filename'] and prior_ext is not None
                  and 'MC_Mass_src' in prior['df'] and bool((prior['df']['MC_Mass_src'] == 'isochrone').any()))
    if keep_prior:
        STATE['iso_stale'] = dict(logAge=prior_ext.get('logAge'), mh=prior_ext.get('MH'), requested_logAge=age, requested_mh=mh)
        notes.append(f"Isochrone for logAge {age}, [M/H] {mh} is not fetched: the mass table from logAge {prior_ext.get('logAge')} is kept but marked stale. "
                     "Press \"Apply Extinction & Fetch PARSEC\" to refresh it.")
    else:
        STATE['iso_stale'] = None
        STATE['iso'] = iso_plot
        STATE['parsec'] = dict(raw=raw_used, name=f'PARSEC_gaiaEDR3_logAge{age:.2f}_MH{mh:+.2f}_AV{av_send:.3f}.dat') if raw_used else None
        STATE['analysis'] = dict(src=p['members_filename'], df=d)
        STATE['iso0'] = iso0                                # intrinsic (A_V=0) isochrone: G_lim -> mass conversion in /compute-imf
        STATE['ext'] = dict(DM=meta['DM'], A_G=EXT.AG_AV * av_mean, A_BP=EXT.ABP_AV * av_mean, A_RP=EXT.ARP_AV * av_mean, A_V=av_mean,
                            E_BV=meta['E_BV_mean'], d_pc=meta['d_pc'], logAge=age, MH=mh)

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
def compute_kinematics(req: KinematicsReq):
    """3-D map + PM vectors: the maths lives in kinematics.py (this route only reads the request and shapes the reply)."""
    p = req.model_dump(exclude_unset=True)
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


_wls = IMF.wls_slope        # kept as a module-level alias (same Poisson-weighted line fit, now defined once in imf.py)


def _fit_limit(p: dict, a_src: str | None):
    """Pick the completeness limit for the mass-function fit. Priority: 1) explicit user override (`fit_min_logM` present in the request);
    2) MEASURED limit (G_lim from the retained-fraction curve -> mass on the A_V=0 isochrone); 3) config default FIT_MIN_LOGM.
    Returns (log10 m_min, info dict). The measured limit is only trusted when the table being analysed is the stage-1 member list."""
    comp, iso0, ext = STATE['completeness'], STATE['iso0'], STATE['ext']
    info = dict(source='default', g_lim=None, m_min=None, note='')
    measured = None
    if comp is not None and iso0 is not None and ext is not None and a_src == MEMBERS_NAME and np.isfinite(comp['g_lim']):
        try:
            measured = IMF.g_limit_to_mass(comp['g_lim'], ext['DM'], ext['A_G'], iso0)
            info.update(g_lim=float(comp['g_lim']), m_min=measured['m_min'], clipped=measured['clipped'], M_G_lim=measured['M_G_lim'])
            if measured['clipped']:
                info['note'] = 'G_lim maps outside the isochrone magnitude range: the mass limit is clipped to the isochrone end (a bound, not a measurement).'
        except ValueError as e:
            info['note'] = str(e)
    if p.get('fit_min_logM') is not None:
        info['source'] = 'override'
        return float(p['fit_min_logM']), info
    if measured is not None:
        info['source'] = 'measured'
        return measured['logm'], info
    if not info['note']:
        info['note'] = 'No measured completeness available (needs stage-1 members + fetched PARSEC isochrone); using the config default.'
    return IMF.FIT_MIN_LOGM, info


@app.post(API_PREFIX + '/compute-imf')
def compute_imf(req: ImfReq):
    p = req.model_dump(exclude_unset=True)
    d_all = analysis_rows(p['filename'])
    a = STATE['analysis']
    a_src = a['src'] if a is not None and a['src'] == p['filename'] else None
    stale = _stale_isochrone_message() if a_src is not None else None
    if stale:
        raise HTTPException(400, stale)
    d = d_all[~d_all['Star_Type'].isin(EVOLVED_TYPES)] if 'Star_Type' in d_all else d_all
    n_excl = len(d_all) - len(d)
    lo, hi = float(p.get('mass_min') or IMF.M_MIN_FLOOR), float(p.get('mass_max') or IMF.M_MAX_CAP)
    # mass-quality gate: approx_ML masses (L(M) fallback) are dropped unless the request opts in
    prim, n_approx = IMF.mass_sample(d, include_approx_ml=bool(p.get('include_approx_ml', IMF.INCLUDE_APPROX_ML)))
    prim = prim[np.isfinite(prim) & (prim >= lo) & (prim <= hi)]
    comp = None
    if p.get('include_companions') and 'MC_CompMass' in d:
        c = d['MC_CompMass'].to_numpy(float)
        comp = c[np.isfinite(c) & (c >= lo) & (c <= hi)]
    logm_lim, lim = _fit_limit(p, a_src)
    method = p.get('method') or IMF.METHOD
    try:
        r = IMF.compute_imf(prim, comp, fit_min_logm=logm_lim, n_bins=p.get('bins'), method=method)
    except ValueError as e:
        extra = f' ({n_approx} approx_ML masses were excluded; fetch the PARSEC isochrone first or set include_approx_ml)' if n_approx else ''
        raise HTTPException(400, str(e) + extra)
    STATE['imf_bins'], STATE['imf_chi2'], STATE['imf_res'] = r.bins, r.chi2, r         # imf_res: full IMFResult for figure export (report.py)
    b = r.bins
    cb, ib = b[(b['count'] > 0) & b['complete']], b[(b['count'] > 0) & ~b['complete']]
    pack = lambda t: dict(logM_center=t['logM_center'].tolist(), log_dN=t['log_dN'].tolist(), log_err=t['log_err'].tolist(), count=t['count'].tolist())
    m = r.meta
    slope, se, icpt, r2 = m['fit_gamma'], m['fit_gamma_err'], m['fit_intercept'], m['binned']['r2']
    half = m['bin_width_dex'] / 2          # line spans the whole fitted range (bin edges), like the page
    xf = [float(cb['logM_center'].min()) - half, float(cb['logM_center'].max()) + half]
    age = float(p.get('log_age') or 7.5)
    comp_c = STATE['completeness'] if a_src == MEMBERS_NAME else None      # measured curve belongs to the stage-1 member list only
    if comp_c is not None:
        cv = comp_c['curve'].dropna(subset=['fraction'])
        completeness = dict(mag=cv['Gmag'].tolist(), fraction=cv['fraction'].tolist(), n_raw=cv['n_raw'].tolist(), n_clean=cv['n_clean'].tolist(),
                            assumed=False, source='measured', g_lim=lim['g_lim'], frac=comp_c['threshold'], plateau=comp_c['plateau'], rule=comp_c['rule'],
                            mode=comp_c['mode'], norm=comp_c['norm'], prefiltered=comp_c.get('prefiltered', False))
    else:
        g = np.arange(8, 20.01, 0.8)
        completeness = dict(mag=g, fraction=1 / (1 + np.exp((g - COMPLETENESS_G50) / COMPLETENESS_W)), assumed=True, source='assumed')
    out = dict(imf=dict(
        mf_kind=IMF.mf_label(age),
        meta=dict(n_masses=m['n_masses'], n_bins=m['n_bins'], bin_width_dex=round(m['bin_width_dex'], 4),
                  n_fit_bins=m['n_fit_bins'], n_in_fit=m['n_in_fit'], fit_min_logm=m['fit_min_logm'],
                  n_excluded_evolved=n_excl, n_excluded_approx_ml=n_approx, m_range=[round(m['m_range'][0], 3), round(m['m_range'][1], 2)],
                  method=m['method'], stat=m['stat'], fit_min_source=lim['source'], g_lim=lim['g_lim'], m_min=lim['m_min'], limit_note=lim['note'],
                  mle=m['mle'], binned=m['binned']),
        bins_table=r.bins.to_dict('records'),
        chi2_table=r.chi2.sort_values('red_chi2').round(3).to_dict('records'),
        complete_bins=pack(cb), incomplete_bins=pack(ib), slope=slope, slope_err=se, r_squared=r2, ms_count=m['n_masses'],
        fit_line=dict(x=xf, y=[slope * xf[0] + icpt, slope * xf[1] + icpt], slope=slope, slope_err=se, n_bins=m['n_fit_bins']),
        theoretical=dict(log_m=r.curves['logM'], salpeter=np.log10(np.maximum(r.curves['Salpeter'], 1e-12)),
                         kroupa=np.log10(np.maximum(r.curves['Kroupa'], 1e-12)), chabrier=np.log10(np.maximum(r.curves['Chabrier'], 1e-12)))),
        completeness=completeness)
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


# ============================ Phase-1 endpoints (NEW; the routes above are unchanged) ============================
def _stage2(name: str):
    """Analysis table + intrinsic isochrone + extinction meta, or a 400 that says which step is missing."""
    d = analysis_rows(name)
    stale = _stale_isochrone_message()
    if stale:
        raise HTTPException(400, stale)
    if STATE['iso0'] is None or STATE['ext'] is None:
        raise HTTPException(400, 'Run "Apply Extinction & Fetch PARSEC" first: this step needs the isochrone, DM and A_V.')
    return d, STATE['iso0'], STATE['ext']


def _clean_json(o):
    """Replace numpy arrays by lists and drop non-finite floats so the dict is JSON-safe."""
    return _none(o)


@app.post(API_PREFIX + '/gmm-diagnostics')
def gmm_diagnostics(req: GmmDiagReq):
    """HDBSCAN cross-check + PMemb resampling on the stored GMM fit of the last /run-gmm."""
    p = req.model_dump(exclude_unset=True)
    res = STATE['gmm_res']
    if res is None:
        raise HTTPException(400, 'Run the GMM first.')
    try:
        rs = GMM.resample_pmemb(res, p.get('n_resample'), p.get('seed'))
        hd = GMM.hdbscan_crosscheck(res, p.get('min_cluster_size'))
    except ValueError as e:
        raise HTTPException(400, str(e))
    mem_idx = res.members.index
    m = rs.loc[mem_idx]
    tab = res.working.loc[rs.index, ['PMemb']].join(rs).join(hd['table'])
    out = dict(n_resample=int(p.get('n_resample') or GMM.N_RESAMPLE), hdbscan=hd['stats'],
               members=dict(n=int(len(m)), mean_pmemb_std=float(m['PMemb_std'].mean()), max_pmemb_std=float(m['PMemb_std'].max()),
                            n_unstable=int((m['PMemb_frac_member'] < 0.9).sum()), n_unstable_def='member in < 90 % of the draws'),
               cluster_rule=p.get('cluster_rule') or GMM.CLUSTER_RULE, table=tab.reset_index(drop=True).to_dict('records'))
    STATE['gmm_diag'] = tab
    return _clean_json(out)


@app.post(API_PREFIX + '/compute-binaries')
def compute_binaries(req: BinariesReq):
    """Monte-Carlo masses with unresolved binaries (binaries.py) + total cluster mass."""
    p = req.model_dump(exclude_unset=True)
    d, iso0, ext = _stage2(p['filename'])
    a = STATE['analysis']
    a_src = a['src'] if a is not None and a['src'] == p['filename'] else None
    iso_ms = iso0[iso0['label'].isin([0, 1])] if 'label' in iso0 else iso0       # PMS + MS only: evolved giants are not matched to MS masses
    f_bin = float(p.get('f_bin', BIN.F_BIN))
    try:
        tab = BIN.assign_masses(d, iso_ms, ext['DM'], ext['A_G'], ext['A_BP'], ext['A_RP'], n_synth=int(p.get('n_synth', BIN.N_SYNTH)),
                                f_bin=f_bin, n_mc=int(p.get('n_mc', BIN.N_MC)), seed=int(p.get('seed', BIN.BIN_SEED)))
    except ValueError as e:
        raise HTTPException(400, str(e))
    mass = None
    note = ''
    prim = tab['MC_Mass'].to_numpy(float)
    if 'Star_Type' in d:
        prim = prim[~d['Star_Type'].isin(EVOLVED_TYPES).to_numpy()]
    prim = prim[np.isfinite(prim)]
    try:
        logm, lim = _fit_limit({}, a_src)
        mle = IMF.mle_slope(prim, 10 ** logm, seed=IMF.BOOT_SEED)
        m_to = float(iso_ms.loc[iso_ms['label'] == 1, 'Mass'].max()) if 'label' in iso_ms else float(iso_ms['Mass'].max())
        mass = BIN.cluster_mass(tab, mle, m_to, f_bin=f_bin)
        mass.update(alpha=mle['alpha'], gamma=mle['gamma'], m_min=mle['m_min'], m_turnoff=m_to, fit_min_source=lim['source'])
    except ValueError as e:
        note = f'cluster mass not computed: {e}'
    out_tab = tab.copy()
    for c in ('Source', 'Gmag', 'BPmag', 'RPmag'):
        if c in d:
            out_tab.insert(0, c, d[c].to_numpy())
    STATE['binaries'] = dict(table=out_tab, mass=mass, f_bin=f_bin)
    if a_src is not None:                                                         # additive: makes the (unchecked) companion-mass option usable
        for c in ('MC_CompMass', 'is_binary', 'q'):
            a['df'][c] = tab[c].to_numpy()
    n_un = int(tab['MC_Mass'].isna().sum())
    warns = []
    if n_un > 0.2 * len(tab):
        warns.append(f'{n_un}/{len(tab)} stars have no synthetic match within {BIN.MATCH_MAX_DIST} mag (isochrone / extinction / DM do not describe them): '
                     'masses and M_T are unreliable.')
    if mass and mass.get('alpha') is not None and min(abs(mass['alpha'] - IMF.ALPHA_BOUNDS[0]), abs(mass['alpha'] - IMF.ALPHA_BOUNDS[1])) < 1e-3:
        warns.append('The mass-function slope hit its bound, so the unseen-mass and white-dwarf extrapolations in M_T are not meaningful.')
    out = dict(n=int(len(tab)), n_binary=int(tab['is_binary'].sum()), observed_binary_fraction=float(tab['is_binary'].mean()),
               n_unmatched=n_un, mass=mass, note=note, warnings=warns,
               caveat='is_binary / binary_prob are probabilities under the assumed F_bin prior, not detections.',
               table=out_tab.reset_index(drop=True).to_dict('records'))
    return _clean_json(out)


@app.post(API_PREFIX + '/compute-velocities')
def compute_velocities(req: VelocitiesReq):
    """Tangential velocities, v_out / v_az, ML velocity dispersion (kinematics.compute_velocities)."""
    p = req.model_dump(exclude_unset=True)
    d = analysis_rows(p['filename'])
    mass = p.get('mass') or ((STATE['binaries'] or {}).get('mass') or {}).get('M_T')
    try:
        r = KIN.compute_velocities(d, eta=float(p.get('eta', KIN.ETA_VIRIAL)), plx_mode=p.get('plx_mode') or KIN.VEL_PLX_MODE,
                                   epsi_max=p.get('epsi_max', KIN.EPSI_MAX), perspective=p.get('perspective', True) is not False,
                                   mass_msun=float(mass) if mass else None, n_boot=int(p.get('n_boot', KIN.N_BOOT_KIN)))
    except ValueError as e:
        raise HTTPException(400, str(e))
    STATE['velocities'] = r
    return _clean_json(r)


@app.post(API_PREFIX + '/compute-dynamics')
def compute_dynamics(req: DynamicsReq):
    """Mass segregation (MSR), tidal radius, relaxation time, optional orbit / tail check (dynamics.py)."""
    p = req.model_dump(exclude_unset=True)
    d = analysis_rows(p['filename'])
    if 'MC_Mass' not in d:
        raise HTTPException(400, 'No MC_Mass column: run "Apply Extinction & Fetch PARSEC" first.')
    bt = STATE['binaries']
    masses = d['MC_Mass'].to_numpy(float)
    if bt is not None and len(bt['table']) == len(d):                             # prefer the binary-aware primary mass; fall back per star when unmatched
        mc = bt['table']['MC_Mass'].to_numpy(float)
        masses = np.where(np.isfinite(mc), mc, masses)
    ra, dec = d['RA_ICRS'].to_numpy(float), d['DE_ICRS'].to_numpy(float)
    if ra.max() - ra.min() > 180:
        ra = np.where(ra < 180, ra + 360, ra)
    d_pc = float(1000.0 / np.nanmedian(d['Plx']))
    tx, ty = KIN.tangent_offsets(ra, dec)
    xy = np.column_stack([tx, ty]) * d_pc
    msr = DYN.mass_segregation_ratio(xy, masses, p.get('n_massive') or DYN.N_MASSIVE, int(p.get('n_random', DYN.N_RANDOM)))
    m_t = p.get('mass') or ((bt or {}).get('mass') or {}).get('M_T') or float(np.nansum(masses))
    out = dict(msr=msr.to_dict('records'), mass_msun=float(m_t), mass_source='user' if p.get('mass') else ('binaries.M_T' if (bt or {}).get('mass') else 'sum of observed MC_Mass (lower limit)'),
               r_t_pc=DYN.tidal_radius(float(m_t)), d_pc=d_pc, notes=[])
    r_hn = KIN.half_number_radius_pc(ra, dec, d_pc)
    out['r_half_number_pc'] = r_hn
    vel = STATE['velocities']
    if vel is None and p.get('sigma_1d') is None:
        try:
            vel = KIN.compute_velocities(d, n_boot=100)
        except ValueError as e:
            out['notes'].append(f'velocity dispersion unavailable: {e}')
    sig = p.get('sigma_1d') or (vel or {}).get('sigma_1d_kms')
    if sig:
        n_tot = int(p.get('n_total') or max(len(d), round(float(m_t) / max(np.nanmean(masses), 1e-3))))
        out['relaxation'] = DYN.relaxation_time(n_tot, r_hn * 4.0 / 3.0, float(sig))
        age = (STATE['ext'] or {}).get('logAge')
        if age is not None:
            out['relaxation']['age_myr'] = float(10 ** age / 1e6)
            out['relaxation']['age_over_trelax'] = float(10 ** age / 1e6 / out['relaxation']['t_relax_myr'])
    if p.get('orbit', True) is not False:
        rv = (vel or {}).get('rv_cluster_kms')
        if rv is not None and np.isfinite(rv):
            try:
                k3 = KIN.kin3d(ra, dec, d['Plx'].to_numpy(float), d['e_Plx'].to_numpy(float))
                tail = DYN.orbit_tail_check(float(np.median(ra)), float(np.median(dec)), float(np.nanmedian(d['Plx'])),
                                            float(np.median(d['pmRA'])), float(np.median(d['pmDE'])), float(rv), k3['P'])
                tail['orbit_path_pc'] = tail['orbit_path_pc'][::10]
                tail['t_myr'] = tail['t_myr'][::10]
                out['orbit'] = tail
            except (RuntimeError, ValueError) as e:
                out['notes'].append(f'orbit / tail check skipped: {e}')
        else:
            out['notes'].append('orbit / tail check skipped: no cluster radial velocity (needs >= 3 Gaia RVs)')
    STATE['dynamics'] = out
    return _clean_json(out)


@app.post(API_PREFIX + '/fit-isochrone-likelihood')
def fit_isochrone_likelihood(req: IsoFitReq):
    """Likelihood isochrone fit (A_V, DM, binary-sequence fraction, outlier fraction) with a selectable dust source."""
    p = req.model_dump(exclude_unset=True)
    d, iso0, ext = _stage2(p['filename'])
    src = p.get('dust_source') or EXT.DUST_SOURCE
    try:
        if src == 'bayestar':
            ebv, _ = EXT.get_dust_ebv('bayestar', d, d_pc=ext['d_pc'])
            av_prior = float(EXT.AV_PER_EBV * ebv[0])
        elif src == 'free':
            av_prior = None
        elif src == 'sfd':
            av_prior = float(ext['A_V'])                                           # the A_V already used by /compute-cmd (IRSA or typed value)
        else:
            raise HTTPException(400, "dust_source must be 'sfd', 'bayestar' or 'free'")
        iso_ms = iso0[iso0['label'] <= EXT.ISO_MAX_LABEL] if 'label' in iso0 else iso0
        r = EXT.fit_isochrone_likelihood(d, iso_ms, ext['DM'], av_prior, src, fit_binaries=p.get('fit_binaries', True) is not False)
    except RuntimeError as e:
        raise HTTPException(400, str(e))
    except ValueError as e:
        raise HTTPException(400, str(e))
    pb, po = r.pop('p_binary_seq'), r.pop('p_outlier')
    r.update(n_binary_seq=int((pb > 0.5).sum()), n_outlier=int((po > 0.5).sum()), legacy_note='fit_quality (nearest-point distance) is unchanged in /compute-cmd')
    STATE['isofit'] = r
    return _clean_json(r)


@app.get(API_PREFIX + '/hunt-reffert')
def hunt_reffert(name: str):
    """Optional literature prefill (Hunt & Reffert 2023) for distance / PM / parallax. Never applied automatically."""
    try:
        return _clean_json(GMM.hunt_reffert_prefill(name))
    except RuntimeError as e:
        raise HTTPException(502, str(e))


@app.get(API_PREFIX + '/download-binaries')
def dl_binaries():
    b = STATE['binaries']
    return _csv(b['table'] if b else None, 'binaries_masses.csv')


# static files of the front end (JS / CSS / team photos) so that http://localhost:8000/ shows the complete page
for _d in ('assets', 'images'):
    _p = os.path.join(ROOT, _d)
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
