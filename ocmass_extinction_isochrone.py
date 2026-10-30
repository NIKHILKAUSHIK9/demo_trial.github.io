"""
ocmass_extinction_isochrone.py -- Stage 2: IRSA dust extinction -> PARSEC (CMD 3.8, Gaia DR3 photometry) -> isochrone fit (no UI code).

TERMINOLOGY
  E(B-V)        : colour excess from the SFD dust map (IRSA, 'ext SFD mean'). A_V = AV_PER_EBV * E(B-V)  (Schlafly & Finkbeiner 2011: 2.742).
  A_G/A_BP/A_RP : Gaia-band extinction = (A_lambda/A_V) * A_V with fixed ratios below (used for the de-reddened columns).
  DM            : distance modulus = 5*log10(d_pc) - 5.
  M_G_corr      : absolute, de-reddened G = G - DM - A_G.      bp_rp_corr : intrinsic colour (BP-RP)_0.
  M_G_obs       : absolute but NOT de-reddened G = G - DM.     bp_rp_obs  : observed BP-RP.
  CMD 3.8       : the PARSEC web service (https://stev.oapd.inaf.it/cgi-bin/cmd_3.8). It can apply an interstellar extinction A_V to the
                  isochrone itself, so the returned magnitudes are "reddened absolute magnitudes" (M_band + A_band).
  Gaia DR3      : DR3 photometry (G, BP, RP) IS the EDR3 photometry, so the CMD system 'Gaia EDR3' (Vega) is the correct one.
                  NOT 'Gaia DR2' / 'gaiaDR2maiz' / 'gaiaDR2weiler': different passbands -> wrong colours.

WORKFLOW (what the web app runs, in this order)
  1 filtered members CSV (GMM) ->
  2 IRSA Dust: E(B-V) per star (deduplicated by sky cell) or at the cluster centre -> mean E(B-V), mean A_V ->
  3 reference-code extinction columns (A_G, A_BP, A_RP, M_G_corr, bp_rp_corr) ->
  4 CMD 3.8 request: logAge, [M/H], Gaia EDR3 Vega, A_V = mean A_V ('reddened' handling) or 0 ('intrinsic' handling) -> raw PARSEC file ->
  5 overlay the isochrone on the stars in the SAME space (reddened isochrone <-> observed stars, or Av=0 isochrone <-> de-reddened stars)
    + fit-quality numbers.  Extinction is applied ONCE: never a reddened isochrone on de-reddened stars.

NETWORK POLITENESS (IRSA and CMD are shared public services)
  results are cached, requests are rate limited, deterministic errors (bad parameters) are never retried, nothing is fetched on a slider move.
"""
from __future__ import annotations

import gzip
import io
import json
import os
import re
import threading
import time
import zlib
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

# ===================== GLOBAL CONFIGURATION =====================
COL = dict(ra='RA_ICRS', dec='DE_ICRS', g='Gmag', bp='BPmag', rp='RPmag', plx='Plx')

# --- network (change here, or set the environment variable shown) ---
NET_USE_SYSTEM_PROXY = os.environ.get('OCMASS_NO_PROXY', '0') != '1'   # False: ignore HTTP(S)_PROXY / Windows proxy settings (set OCMASS_NO_PROXY=1)
NET_PROXY_URL = os.environ.get('OCMASS_PROXY') or None                 # e.g. 'http://proxy.du.ac.in:3128' : force this proxy for IRSA and CMD
NET_IPV4_ONLY = os.environ.get('OCMASS_IPV4_ONLY', '0') == '1'         # True: never use IPv6 (fixes hanging connections on broken IPv6 routes)
NET_CONNECT_TIMEOUT_S = 30                                             # TCP connect timeout (a blocked host is reported after ~this x retries)

# --- IRSA dust (E(B-V)) ---
EBV_MODE = 'center'            # 'fixed' | 'center' (ONE IRSA call at the cluster centre) | 'per_star' (one call per sky cell, then averaged)
EBV_FIXED = None               # float used when EBV_MODE == 'fixed'
EBV_KEY = 'ext SFD mean'       # IRSA ebv-section column ('ext SandF mean' for the Schlafly&Finkbeiner map)
IRSA_RETRIES, IRSA_WAIT_S = 4, 3.0           # connection failures only (flaky networks); bad answers are not retried
IRSA_ROUND_DEG = 0.05          # per_star: stars inside the same 0.05 deg cell share ONE IRSA query (SFD beam is ~6 arcmin)
IRSA_MAX_QUERIES = 150         # per_star: refuse more unique cells than this (use 'center')
IRSA_MIN_INTERVAL_S = 0.3      # pause between consecutive IRSA calls

# --- extinction law used for the de-reddened columns ---
AV_PER_EBV = 2.742             # SFD E(B-V) -> A_V (S&F 2011). Original code used 3.1 (=13% too high for SFD)
AG_AV, ABP_AV, ARP_AV = 0.789, 1.002, 0.589   # fixed Gaia A_lambda/A_V (as in the reference code)
USE_CLUSTER_DM = True          # True: one DM from median member parallax; False: per-star DM
PLX_ZP = 0.0                   # parallax zero-point (mas)
BPRP0_RANGE = (-0.8, 5.0)      # plausible intrinsic BP-RP; outside -> bprp_corr_ok = False (over-correction flag)

# --- CMD 3.8 / PARSEC ---
CMD_URL = 'https://stev.oapd.inaf.it/cgi-bin/cmd_3.8'  # the CMD 3.8 CGI (same address as its form page)
CMD_URL_FALLBACKS = {}                                   # {primary_url: [other urls]}; the generic /cgi-bin/cmd now serves CMD 3.9, so it is NOT used
CMD_SEND_VARIANTS = ('body', 'query')                      # 'body': form-encoded POST body (what the browser form sends); 'query': fields in the URL (what ezpadova sends). Tried in this order, only while no output link comes back
CMD_VERSION = '3.8'
PHOTSYS = 'gaiaEDR3'                                  # Gaia DR3 photometry = EDR3 photometry (Vega, passbands from the ESA/Gaia site)
PHOTSYS_FILE = f'YBC_tab_mag_odfnew/tab_mag_{PHOTSYS}.dat'
PHOTSYS_VERSION = 'YBCnewVega'                        # YBC bolometric corrections + Vega SED of Bohlin et al. 2020
CMD_TIMEOUT_S = 180
CMD_MIN_INTERVAL_S = 2.0       # pause between consecutive CMD requests
ISO_RETRIES, ISO_WAIT_S = 4, 5.0
LOGAGE_RANGE = (6.6, 10.13)    # PARSEC range accepted by the form
MH_RANGE = (-2.5, 0.8)         # loose sanity range for [M/H]; the server rejects anything it cannot interpolate
LOGAGE_FLOOR = LOGAGE_RANGE[0]
LOGAGE_HALF_WIDTH, DLOGAGE = 0.5, 0.02
ISO_COLS = ['Gmag', 'G_BPmag', 'G_RPmag', 'Mass', 'logL', 'logTe', 'logg']
ISO_EXTRA = ['logAge', 'MH', 'Mini', 'label']          # kept if present in the download
ISO_NAME_ALIASES = {'Gmag': ['Gmag'], 'G_BPmag': ['G_BPmag', 'G_BPmag_bright', 'BPmag'],
                    'G_RPmag': ['G_RPmag', 'RPmag'], 'Mass': ['Mass'], 'logL': ['logL'],
                    'logTe': ['logTe'], 'logg': ['logg']}
ISO_MAX_LABEL = 7              # PARSEC evolutionary labels drawn: 0 PMS, 1 MS, 2 SGB, 3 RGB, 4 core-He, 5-7 E-AGB (8, 9 = TP-AGB, post-AGB dropped)
CMD_BASE_PARAMS = {            # every field of the CMD 3.8 form (defaults of the form), so nothing is left to server-side guessing
    'cmd_version': CMD_VERSION, 'track_omegai': '0.00', 'track_parsec': 'parsec_CAF09_v1.2S',
    'track_colibri': 'parsec_CAF09_v1.2S_S_LMC_08_web', 'track_postagb': 'no', 'n_inTPC': '10', 'eta_reimers': '0.2',
    'kind_interp': '1', 'kind_postagb': '-1', 'dust_sourceM': 'dpmod60alox40', 'dust_sourceC': 'AMCSIC15',
    'kind_mag': '2', 'kind_dust': '0', 'extinction_coeff': 'constant', 'extinction_curve': 'cardelli', 'kind_LPV': '3',
    'imf_file': 'tab_imf/imf_kroupa_orig.dat', 'output_kind': '0', 'output_evstage': '1',
    'lf_maginf': '-15', 'lf_magsup': '20', 'lf_deltamag': '0.5', 'sim_mtot': '1.0e4', 'submit_form': 'Submit'}

# --- reddening self-check (reddened isochrone must be shifted by roughly the Gaia extinction coefficients) ---
CHECK_DG_RANGE = (0.60, 1.00)  # expected  dG   / A_V
CHECK_DC_RANGE = (0.25, 0.65)  # expected  d(BP-RP) / A_V
FIT_NEAR_MAG = 0.2             # a star "fits" when it lies within this many mag of the isochrone in the CMD plane
OUT_DIR = './ocmass_output'
# =================================================================


class CMDRequestError(RuntimeError):
    """Deterministic CMD failure (bad parameters / unexpected page): NOT retried, so the shared server is not hammered."""


class NetworkNotAllowed(RuntimeError):
    """Raised in cache-only mode when a value is not cached yet (nothing was fetched)."""


def _net_setup():
    """Apply the NET_* settings once: IPv4-only resolution for urllib3 (requests/astroquery)."""
    if NET_IPV4_ONLY:
        import socket
        import urllib3.util.connection as _c
        _c.allowed_gai_family = lambda: socket.AF_INET


def _proxies():
    return {'http': NET_PROXY_URL, 'https': NET_PROXY_URL} if NET_PROXY_URL else None


def _retry(fn, n, wait, what, fatal=(CMDRequestError, NetworkNotAllowed)):
    last = None
    for i in range(n):
        try:
            return fn()
        except fatal:
            raise
        except Exception as e:           # network / server errors
            last = e
            if i < n - 1:
                time.sleep(wait * (i + 1))
    raise RuntimeError(f'{what} failed after {n} attempts: {last}')


# ------------------------------ step 2: IRSA dust ------------------------------
_IRSA_CACHE: dict = {}
_IRSA_LOCK = threading.Lock()
_irsa_last = [0.0]


def _irsa_ebv_value(ra: float, dec: float) -> float:
    """THE single IRSA network call: E(B-V) from the SFD map at (ra, dec) in degrees."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.ipac.irsa.irsa_dust import IrsaDust
    wait = IRSA_MIN_INTERVAL_S - (time.time() - _irsa_last[0])
    if wait > 0:
        time.sleep(wait)
    _net_setup()
    IrsaDust.TIMEOUT = NET_CONNECT_TIMEOUT_S
    sess = getattr(IrsaDust, '_session', None)
    if sess is not None:
        sess.trust_env = NET_USE_SYSTEM_PROXY
        if NET_PROXY_URL:
            sess.proxies.update(_proxies())
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
    tab = IrsaDust.get_query_table(c, section='ebv')
    _irsa_last[0] = time.time()
    if EBV_KEY not in tab.colnames:
        raise KeyError(f"IRSA table has no column '{EBV_KEY}'. Columns: {tab.colnames}")
    return float(tab[EBV_KEY][0])


def _irsa_cell(ra: float, dec: float) -> tuple:
    g = IRSA_ROUND_DEG
    return (round(round(ra / g) * g, 4), round(round(dec / g) * g, 4))


def _irsa_cached(cell: tuple, cache_only: bool) -> float:
    with _IRSA_LOCK:
        if cell in _IRSA_CACHE:
            return _IRSA_CACHE[cell]
    if cache_only:
        raise NetworkNotAllowed('IRSA value not cached')
    v = _retry(lambda: _irsa_ebv_value(*cell), IRSA_RETRIES, IRSA_WAIT_S, 'IRSA query')
    with _IRSA_LOCK:
        _IRSA_CACHE[cell] = v
    return v


def get_ebv_info(df: pd.DataFrame, mode: str = EBV_MODE, ebv_fixed: float | None = EBV_FIXED,
                 progress=None, cache_only: bool = False):
    """E(B-V) per star + info dict(n_queries, n_failed, cells). progress: optional callable(frac) for per_star."""
    n = len(df)
    if mode == 'fixed':
        if ebv_fixed is None:
            raise ValueError("mode='fixed' requires ebv_fixed")
        return np.full(n, float(ebv_fixed)), dict(n_queries=0, n_failed=0)
    ra, dec = df[COL['ra']].to_numpy(float), df[COL['dec']].to_numpy(float)
    if mode == 'center':
        v = _irsa_cached(_irsa_cell(float(np.nanmedian(ra)), float(np.nanmedian(dec))), cache_only)
        return np.full(n, v), dict(n_queries=1, n_failed=0)
    if mode == 'per_star':
        cells = [_irsa_cell(a, b) if np.isfinite(a) and np.isfinite(b) else None for a, b in zip(ra, dec)]
        uniq = sorted({c for c in cells if c is not None})
        todo = [c for c in uniq if c not in _IRSA_CACHE]
        if cache_only and todo:
            raise NetworkNotAllowed('IRSA values not cached')
        if len(todo) > IRSA_MAX_QUERIES:
            raise ValueError(f'per_star needs {len(todo)} IRSA queries (limit {IRSA_MAX_QUERIES}, cell {IRSA_ROUND_DEG} deg). Use E(B-V) mode center or fixed.')
        vals, failed = {}, 0
        for k, c in enumerate(uniq):
            try:
                vals[c] = _irsa_cached(c, cache_only)
            except RuntimeError:
                vals[c] = np.nan
                failed += 1
            if progress and (k % 5 == 0 or k == len(uniq) - 1):
                progress((k + 1) / len(uniq))
        out = np.array([vals.get(c, np.nan) if c is not None else np.nan for c in cells], float)
        return out, dict(n_queries=len(todo), n_failed=failed, n_cells=len(uniq))
    raise ValueError(f'Unknown EBV mode: {mode}')


def get_ebv(df: pd.DataFrame, mode: str = EBV_MODE, ebv_fixed: float | None = EBV_FIXED, progress=None) -> np.ndarray:
    return get_ebv_info(df, mode, ebv_fixed, progress)[0]


# ------------------------------ step 3: extinction columns ------------------------------
def apply_extinction(members: pd.DataFrame, mode: str = EBV_MODE, ebv_fixed: float | None = EBV_FIXED,
                     plx_zp: float = PLX_ZP, use_cluster_dm: bool = USE_CLUSTER_DM, progress=None, cache_only: bool = False):
    """Returns (df_ext, meta). df_ext adds E_BV, A_G, A_BP, A_RP, distance_modulus, M_G_corr, M_BP_corr, M_RP_corr,
    bp_rp_corr, bprp_corr_ok, M_G_obs, bp_rp_obs. meta carries the MEAN E(B-V)/A_V that is forwarded to CMD."""
    need = [COL[k] for k in ('ra', 'dec', 'g', 'bp', 'rp', 'plx')]
    miss = [c for c in need if c not in members.columns]
    if miss:
        raise ValueError(f'Members file missing columns: {miss}')
    d = members.copy()
    ebv, info = get_ebv_info(d, mode, ebv_fixed, progress, cache_only)
    d['E_BV'] = ebv
    d = d.dropna(subset=['E_BV', COL['g'], COL['bp'], COL['rp']]).copy()
    if d.empty:
        raise ValueError('No stars left after E(B-V)/photometry NaN removal.')
    av = AV_PER_EBV * d['E_BV']
    d['A_G'], d['A_BP'], d['A_RP'] = AG_AV * av, ABP_AV * av, ARP_AV * av
    if use_cluster_dm:
        d_pc = 1000.0 / np.nanmedian(d[COL['plx']] - plx_zp)
        d['distance_modulus'] = 5 * np.log10(d_pc) - 5
    else:
        d['distance_modulus'] = 5 * np.log10(1000.0 / (d[COL['plx']] - plx_zp)) - 5
        d_pc = float(10 ** ((d['distance_modulus'].median() + 5) / 5))
    d['M_G_corr'] = d[COL['g']] - d['distance_modulus'] - d['A_G']
    d['M_BP_corr'] = d[COL['bp']] - d['distance_modulus'] - d['A_BP']
    d['M_RP_corr'] = d[COL['rp']] - d['distance_modulus'] - d['A_RP']
    d['bp_rp_corr'] = (d[COL['bp']] - d['A_BP']) - (d[COL['rp']] - d['A_RP'])
    d['bprp_corr_ok'] = d['bp_rp_corr'].between(*BPRP0_RANGE)
    d['M_G_obs'] = d[COL['g']] - d['distance_modulus']              # same plane as a REDDENED isochrone (no extinction removed)
    d['bp_rp_obs'] = d[COL['bp']] - d[COL['rp']]
    e_mean = float(d['E_BV'].mean())
    meta = dict(mode=mode, n_in=len(members), n_out=len(d), E_BV_median=float(d['E_BV'].median()),
                E_BV_mean=e_mean, E_BV_std=float(d['E_BV'].std(ddof=0)), A_V=float(AV_PER_EBV * d['E_BV'].median()),
                A_V_mean=float(AV_PER_EBV * e_mean), d_pc=float(d_pc),
                DM=float(d['distance_modulus'].median()), frac_colour_ok=float(d['bprp_corr_ok'].mean()),
                irsa_queries=int(info.get('n_queries', 0)), irsa_failed=int(info.get('n_failed', 0)),
                irsa_dropped_stars=int(len(members) - len(d)))
    return d, meta


# ------------------------------ step 4: CMD 3.8 -> PARSEC file ------------------------------
_cmd_last = [0.0]
_CMD_LOCK = threading.Lock()
_TMP_RE = re.compile(r'output\d+')


def build_cmd_params(logage: float, mh: float, av: float = 0.0) -> dict:
    """Form payload for ONE isochrone (single log-age, single [M/H]) in the Gaia EDR3 (= DR3) Vega system with interstellar extinction A_V=av."""
    if not (LOGAGE_RANGE[0] <= logage <= LOGAGE_RANGE[1]):
        raise ValueError(f'logAge {logage} outside the PARSEC range {LOGAGE_RANGE}')
    if not (MH_RANGE[0] <= mh <= MH_RANGE[1]):
        raise ValueError(f'[M/H] {mh} outside the accepted range {MH_RANGE}')
    if not (0.0 <= av <= 20.0):
        raise ValueError(f'A_V {av} outside 0..20 mag')
    kw = dict(CMD_BASE_PARAMS)
    kw.update(photsys_file=PHOTSYS_FILE, photsys_version=PHOTSYS_VERSION,
              isoc_isagelog='1', isoc_lagelow=f'{logage:.4f}', isoc_lageupp=f'{logage:.4f}', isoc_dlage='0',
              isoc_ismetlog='1', isoc_metlow=f'{mh:.4f}', isoc_metupp=f'{mh:.4f}', isoc_dmet='0',
              extinction_av=f'{av:.4f}')
    # defaults of the unused linear-age / Z fields (the form always sends them)
    kw.update(isoc_agelow='1.0e9', isoc_ageupp='1.0e10', isoc_dage='0.0', isoc_zlow='0.0152', isoc_zupp='0.03', isoc_dz='0.0')
    return kw


def _http(method: str, url: str, **kw):
    """requests call (NET_* settings applied): TLS verified first; if the certificate is rejected, repeat once unverified and report it (public data, no credentials)."""
    import requests
    _net_setup()
    s = requests.Session()
    s.trust_env = NET_USE_SYSTEM_PROXY
    if NET_PROXY_URL:
        s.proxies.update(_proxies())
    t = kw.pop('timeout', None)
    kw['timeout'] = (NET_CONNECT_TIMEOUT_S, t) if t is not None else NET_CONNECT_TIMEOUT_S      # (connect, read)
    try:
        return s.request(method, url, verify=True, **kw), True
    except requests.exceptions.SSLError:
        import urllib3
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        return s.request(method, url, verify=False, **kw), False


def _reply_summary(html: str, n: int = 300) -> str:
    """Short text of a CMD reply that has no output link: sentences that look like an error first, else the tail of the page (the header is boilerplate)."""
    t = _page_text(html)
    hits = re.findall(r'[^.!?]*\b(?:error|invalid|wrong|must|cannot|can not|not allowed|out of range|please|failed|warning)\b[^.!?]*[.!?]?', t, flags=re.I)
    hits = [h.strip() for h in hits if len(h.strip()) > 8 and 'sorrow' not in h]
    return ('; '.join(hits)[:n]) if hits else ('...' + t[-n:])


def _page_text(html: str) -> str:
    return re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', html)).strip()


def fetch_parsec_raw(logage: float, mh: float, av: float = 0.0, cmd_url: str | None = None):
    """POST one request to CMD 3.8 and download the resulting file. Returns (raw_bytes, params, info)."""
    url = cmd_url or CMD_URL
    params = build_cmd_params(float(logage), float(mh), float(av))

    urls = [url] + CMD_URL_FALLBACKS.get(url, [])

    def go():
        with _CMD_LOCK:                                   # one request at a time + minimum spacing: be gentle with the shared server
            seen = []
            try:
                for u, variant in [(u, v) for u in urls for v in CMD_SEND_VARIANTS]:
                    wait = CMD_MIN_INTERVAL_S - (time.time() - _cmd_last[0])
                    if wait > 0:
                        time.sleep(wait)
                    kw = dict(data=params) if variant == 'body' else dict(params=params)
                    r, ssl_ok = _http('POST', u, timeout=CMD_TIMEOUT_S, allow_redirects=True, **kw)
                    _cmd_last[0] = time.time()
                    if r.status_code != 200:
                        raise RuntimeError(f'CMD server answered HTTP {r.status_code}')
                    m = _TMP_RE.findall(r.text)
                    if m:
                        break
                    seen.append(f'[{variant}] {u}: ' + _reply_summary(r.text))     # a form page came back instead of a file link
                    try:                                                          # keep the raw reply for debugging
                        os.makedirs(OUT_DIR, exist_ok=True)
                        with open(os.path.join(OUT_DIR, f'cmd_last_response_{variant}.html'), 'w', encoding='utf-8') as fh:
                            fh.write(r.text)
                    except OSError:
                        pass
                else:
                    raise CMDRequestError('CMD did not return an output file. Server said: ' + ' | '.join(seen)
                                          + f'  (raw replies saved in {OUT_DIR}/cmd_last_response_*.html)')
                dom = '/'.join(u.split('/')[:3])
                src = f'{dom}/tmp/{m[0]}.dat'
                r2, ssl2 = _http('GET', src, timeout=CMD_TIMEOUT_S)
                if r2.status_code != 200:
                    raise RuntimeError(f'CMD output file {src} answered HTTP {r2.status_code}')
                data = r2.content
            finally:
                _cmd_last[0] = time.time()
        if data[:2] == b'\x1f\x8b':
            data = gzip.decompress(data)
        elif data[:1] == b'\x78':
            try:
                data = zlib.decompress(data)
            except zlib.error:
                pass
        if not data.strip():
            raise RuntimeError('CMD output file is empty')
        return data, dict(source_url=src, ssl_verified=bool(ssl_ok and ssl2), send_variant=variant)

    data, info = _retry(go, ISO_RETRIES, ISO_WAIT_S, 'PARSEC download')
    return data, params, info


def parse_parsec(raw: bytes | str) -> pd.DataFrame:
    """Parse a CMD isochrone table: header = the last '#' line before the first data row; data rows are whitespace separated."""
    text = raw.decode('utf-8', errors='replace') if isinstance(raw, (bytes, bytearray)) else raw
    header, rows = None, []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith('#'):
            if not rows and 'logAge' in s and ('Gmag' in s or 'mbolmag' in s):
                header = s.lstrip('#').split()
            continue
        rows.append(s.split())
    if header is None or not rows:
        raise CMDRequestError('Could not find the column header / data rows in the PARSEC file (is it an isochrone table in a Gaia system?).')
    width = len(header)
    rows = [r[:width] for r in rows if len(r) >= width]
    df = pd.DataFrame(rows, columns=header).apply(pd.to_numeric, errors='coerce')
    return df.dropna(how='all').reset_index(drop=True)


def select_iso_columns(raw: pd.DataFrame) -> pd.DataFrame:
    """Pick ISO_COLS by NAME (never by position); raise with the available names if one is missing."""
    cols, rename = [], {}
    for std in ISO_COLS:
        found = next((a for a in ISO_NAME_ALIASES[std] if a in raw.columns), None)
        if found is None:
            raise KeyError(f"PARSEC output lacks '{std}'. Available: {list(raw.columns)}")
        cols.append(found); rename[found] = std
    extra = [c for c in ISO_EXTRA if c in raw.columns]
    return raw[cols + extra].rename(columns=rename).round(6).copy()


def validate_parsec(iso: pd.DataFrame, logage: float, mh: float):
    """Fail loudly instead of plotting the wrong isochrone."""
    if len(iso) < 20:
        raise CMDRequestError(f'PARSEC file has only {len(iso)} rows')
    if not np.isfinite(iso[['Gmag', 'G_BPmag', 'G_RPmag']].to_numpy(float)).all(axis=1).mean() > 0.9:
        raise CMDRequestError('PARSEC file has too many non-finite magnitudes')
    if 'logAge' in iso and abs(float(iso['logAge'].median()) - logage) > 0.01:
        raise CMDRequestError(f"PARSEC returned logAge {float(iso['logAge'].median()):.3f}, requested {logage:.3f}")
    if 'MH' in iso and abs(float(iso['MH'].median()) - mh) > 0.05:
        raise CMDRequestError(f"PARSEC returned [M/H] {float(iso['MH'].median()):.3f}, requested {mh:.3f}")


def with_provenance(raw: bytes, params: dict, info: dict, note: str = '') -> bytes:
    """Prefix the untouched CMD file with '#' comment lines recording exactly what was requested (still parseable by parse_parsec)."""
    head = ['# ---- OCMass web app: PARSEC file provenance ----',
            f'# retrieved_utc   : {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}',
            f'# service         : CMD {CMD_VERSION}  (https://stev.oapd.inaf.it/cgi-bin/cmd_{CMD_VERSION})',
            f'# photometric sys : Gaia EDR3 = Gaia DR3 photometry, Vega  ({params.get("photsys_file")}, {params.get("photsys_version")})',
            f'# logAge / [M/H]  : {params.get("isoc_lagelow")} / {params.get("isoc_metlow")}',
            f'# interstellar A_V: {params.get("extinction_av")} mag ({params.get("extinction_curve")}, coefficients: {params.get("extinction_coeff")})',
            f'# tracks          : {params.get("track_parsec")} + {params.get("track_colibri")}',
            f'# tls verified    : {info.get("ssl_verified")}',
            f'# file url        : {info.get("source_url")}',
            f'# request sent as : {info.get("send_variant")}']
    if note:
        head.append(f'# note            : {note}')
    head.append('# request params  : ' + json.dumps(params, sort_keys=True))
    head.append('# ---- original CMD output follows unchanged ----')
    return ('\n'.join(head) + '\n').encode('utf-8') + raw


_PARSEC_CACHE: dict = {}


def get_parsec(logage: float, mh: float, av: float = 0.0, cache_only: bool = False, cmd_url: str | None = None) -> dict:
    """Cached: dict(iso (named columns), raw (bytes incl. provenance), params, info, key). One CMD request per distinct (logAge, [M/H], A_V)."""
    key = (round(float(logage), 4), round(float(mh), 4), round(float(av), 3), cmd_url or CMD_URL)
    if key in _PARSEC_CACHE:
        return _PARSEC_CACHE[key]
    if cache_only:
        raise NetworkNotAllowed('PARSEC file not cached')
    raw, params, info = fetch_parsec_raw(key[0], key[1], key[2], cmd_url)
    iso = select_iso_columns(parse_parsec(raw))
    validate_parsec(iso, key[0], key[1])
    out = dict(iso=iso, raw=with_provenance(raw, params, info), params=params, info=info, key=key)
    _PARSEC_CACHE[key] = out
    if len(_PARSEC_CACHE) > 24:
        _PARSEC_CACHE.pop(next(iter(_PARSEC_CACHE)))
    return out


def fetch_isochrone(logage: float, mh: float) -> pd.DataFrame:
    """Intrinsic (A_V=0) single-age isochrone, named columns."""
    if logage < LOGAGE_FLOOR:
        raise ValueError(f'logAge {logage} below PARSEC floor {LOGAGE_FLOOR}')
    return get_parsec(logage, mh, 0.0)['iso']


def fetch_isochrone_grid(logage: float, mh: float, half_width: float = LOGAGE_HALF_WIDTH,
                         dlogage: float = DLOGAGE) -> pd.DataFrame:
    """Grid logAge +/- half_width (lower edge clamped to LOGAGE_FLOOR): ONE request with a step."""
    lo = max(logage - half_width, LOGAGE_FLOOR)
    hi = logage + half_width
    if lo >= hi:
        raise ValueError(f'Empty age grid: {lo}..{hi}')
    p = build_cmd_params(lo, mh, 0.0)
    p.update(isoc_lageupp=f'{hi:.4f}', isoc_dlage=f'{dlogage:.4f}')
    n = int(round((hi - lo) / dlogage)) + 1
    if n > 400:
        raise ValueError(f'Age grid of {n} isochrones is too large for the shared CMD server')
    import requests  # noqa: F401  (clear ImportError early)
    raw_params = p

    def go():
        with _CMD_LOCK:
            r, _ = _http('POST', CMD_URL, params=raw_params, timeout=CMD_TIMEOUT_S, allow_redirects=True)
            m = _TMP_RE.findall(r.text)
            if r.status_code != 200 or not m:
                raise CMDRequestError('CMD grid request failed: ' + _page_text(r.text)[:300])
            r2, _ = _http('GET', '/'.join(CMD_URL.split('/')[:3]) + f'/tmp/{m[0]}.dat', timeout=CMD_TIMEOUT_S)
            data = r2.content
            return gzip.decompress(data) if data[:2] == b'\x1f\x8b' else data
    return select_iso_columns(parse_parsec(_retry(go, ISO_RETRIES, ISO_WAIT_S, 'PARSEC grid download')))


# ------------------------------ step 5: fit ------------------------------
def iso_for_plot(iso: pd.DataFrame, max_label: int = ISO_MAX_LABEL) -> pd.DataFrame:
    return iso[iso['label'] <= max_label] if 'label' in iso else iso


def check_reddening_shift(iso0: pd.DataFrame, iso_av: pd.DataFrame, av: float) -> dict:
    """Compare the A_V>0 isochrone with the A_V=0 one at equal initial mass: the shift must look like Gaia extinction.
    Guards against assuming CMD applied the extinction when it did not (or applied it twice)."""
    out = dict(av=float(av), dG=None, dBPRP=None, ok=None, message='')
    if av <= 0.0005:
        out.update(ok=True, message='A_V = 0: nothing to check')
        return out
    key = 'Mini' if ('Mini' in iso0 and 'Mini' in iso_av) else 'Mass'
    a, b = iso0.sort_values(key), iso_av.sort_values(key)
    a, b = a[a['label'] == 1] if 'label' in a else a, b[b['label'] == 1] if 'label' in b else b     # main sequence only
    if len(a) < 5 or len(b) < 5:
        out.update(ok=None, message='too few main-sequence points to verify the reddening shift')
        return out
    x = b[key].to_numpy(float)
    g0 = np.interp(x, a[key].to_numpy(float), a['Gmag'].to_numpy(float))
    c0 = np.interp(x, a[key].to_numpy(float), (a['G_BPmag'] - a['G_RPmag']).to_numpy(float))
    dG = float(np.median(b['Gmag'].to_numpy(float) - g0))
    dC = float(np.median((b['G_BPmag'] - b['G_RPmag']).to_numpy(float) - c0))
    ok = (CHECK_DG_RANGE[0] * av <= dG <= CHECK_DG_RANGE[1] * av) and (CHECK_DC_RANGE[0] * av <= dC <= CHECK_DC_RANGE[1] * av)
    out.update(dG=dG, dBPRP=dC, A_G_over_A_V=dG / av, E_BPRP_over_A_V=dC / av, ok=bool(ok),
               message=('reddening applied by CMD as expected' if ok else
                        f'UNEXPECTED reddening shift (dG={dG:.3f}, d(BP-RP)={dC:.3f} for A_V={av:.3f}): do not trust this fit'))
    return out


def _densify(x: np.ndarray, y: np.ndarray, n: int = 3000):
    seg = np.hypot(np.diff(x), np.diff(y))
    s = np.r_[0, np.cumsum(seg)]
    if s[-1] <= 0:
        return x, y
    t = np.linspace(0, s[-1], n)
    return np.interp(t, s, x), np.interp(t, s, y)


def fit_quality(color, mag, iso_color, iso_mag, near: float = FIT_NEAR_MAG) -> dict:
    """Distance (mag, in the CMD plane) from each star to the nearest point of the polyline isochrone."""
    c, m = np.asarray(color, float), np.asarray(mag, float)
    ok = np.isfinite(c) & np.isfinite(m)
    ic, im = np.asarray(iso_color, float), np.asarray(iso_mag, float)
    g = np.isfinite(ic) & np.isfinite(im)
    if ok.sum() == 0 or g.sum() < 3:
        return dict(n=int(ok.sum()), median=None, rms=None, frac_near=None, near=near)
    dc, dm = _densify(ic[g], im[g])
    d = np.sqrt((c[ok, None] - dc[None, :]) ** 2 + (m[ok, None] - dm[None, :]) ** 2).min(axis=1)
    return dict(n=int(ok.sum()), median=float(np.median(d)), rms=float(np.sqrt(np.mean(d ** 2))), frac_near=float((d <= near).mean()), near=near)


def plot_cmd(df_ext: pd.DataFrame, iso: pd.DataFrame | None = None, title: str = '', reddened: bool = False) -> Figure:
    """CMD with optional isochrone overlay. reddened=False: de-reddened absolute CMD; True: observed colour, G-DM (use the A_V>0 isochrone)."""
    f = Figure(figsize=(7, 7)); ax = f.subplots()
    cx, my = ('bp_rp_obs', 'M_G_obs') if reddened else ('bp_rp_corr', 'M_G_corr')
    ok = df_ext['bprp_corr_ok'] if 'bprp_corr_ok' in df_ext else np.ones(len(df_ext), bool)
    ax.scatter(df_ext.loc[~ok, cx], df_ext.loc[~ok, my], s=8, c='lightgrey', label='colour out of range')
    ax.scatter(df_ext.loc[ok, cx], df_ext.loc[ok, my], s=8, c='tab:red', label='members')
    if iso is not None and len(iso):
        i = iso_for_plot(iso)
        ax.plot(i['G_BPmag'] - i['G_RPmag'], i['Gmag'], c='k', lw=1.3, label='PARSEC')
    ax.invert_yaxis()
    ax.set(xlabel=r'$(G_{BP}-G_{RP})$' + ('' if reddened else r'$_0$'), ylabel=r'$M_G$', title=title or ('CMD (reddened isochrone)' if reddened else 'De-reddened CMD'))
    ax.legend()
    return f


def save_outputs(tag: str, df_ext: pd.DataFrame, iso_single: pd.DataFrame | None = None,
                 iso_grid: pd.DataFrame | None = None, out_dir: str = OUT_DIR, parsec_raw: bytes | None = None) -> dict:
    """CSV exports, index=False (+ the raw PARSEC .dat when given)."""
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    p = os.path.join(out_dir, f'{tag}_extinction_corrected.csv'); df_ext.to_csv(p, index=False); paths['extinction'] = p
    if iso_single is not None:
        p = os.path.join(out_dir, f'{tag}_iso_single.csv'); iso_single.to_csv(p, index=False); paths['iso_single'] = p
    if iso_grid is not None:
        p = os.path.join(out_dir, f'{tag}_grid_iso_aligned.csv'); iso_grid.to_csv(p, index=False); paths['iso_grid'] = p
    if parsec_raw is not None:
        p = os.path.join(out_dir, f'{tag}_PARSEC_gaiaEDR3.dat'); open(p, 'wb').write(parsec_raw); paths['parsec_raw'] = p
    return paths
