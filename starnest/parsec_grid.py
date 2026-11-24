"""
parsec_grid.py -- local PARSEC isochrone grid: read, interpolate in (logAge, [M/H]), shift by A_V (Phase 2.2).

TERMINOLOGY
  Grid node          : one isochrone stored at an exact (logAge, [M/H]) pair; the grid is A_V = 0 (intrinsic) only.
  Evolutionary label : PARSEC phase flag (0 PMS, 1 MS, 2 SGB, 3 RGB, 4 core-He, 5-7 E-AGB, 8-9 TP-AGB / post-AGB).
  Arc coordinate u   : 0..1 position along ONE label segment, measured as normalised path length in the (logTe, logL) plane. Isochrones of
                       different ages have different numbers of points, so they can only be blended at equal u (an EEP-like coordinate).
  Bilinear blend     : value = sum over the 4 surrounding nodes of weight * value(u); weights are linear in logAge and [M/H].
  A_V shift          : reddening applied analytically: mag_band += (A_band / A_V) * A_V, with the constant coefficients CMD 3.8 itself uses
                       (extinction_coeff='constant', Cardelli curve). The coefficients are MEASURED from CMD at build time and stored in the grid.

WORKFLOW (interpolate_isochrone)
  1) locate the 4 bracketing nodes and the linear weights; an exact node is returned untouched;
  2) per label present in the corners: resample every corner segment onto a common u grid (length = longest segment);
  3) blend the columns (Mini and Mass in log10) with the bilinear weights; a label missing from some corners is blended from the corners that
     have it and kept only if those corners carry >= LABEL_MIN_WEIGHT of the total weight;
  4) concatenate the segments in label order and stamp logAge / MH / label.
"""
from __future__ import annotations

import json
import os
import threading
from datetime import datetime, timezone

import numpy as np
import pandas as pd

# ===================== GLOBAL CONFIGURATION =====================
_DIR = os.environ.get('STARNEST_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root (data/ lives next to the package)
GRID_PATH = os.environ.get('STARNEST_PARSEC_GRID') or os.path.join(_DIR, 'data', 'parsec_grid.parquet')
ALLOW_MOCK = os.environ.get('STARNEST_ALLOW_MOCK_GRID', '0') == '1'     # a grid whose metadata says source='MOCK' is refused unless this is set
GRID_NODE_COLS = ['logAge_node', 'MH_node']                              # exact node coordinates stored with every row
DATA_COLS = ['Mini', 'Mass', 'logL', 'logTe', 'logg', 'Gmag', 'G_BPmag', 'G_RPmag']   # blended columns
LOG_COLS = ('Mini', 'Mass')                                              # blended in log10 (positive, span >1 dex)
OUT_COLS = ['logAge', 'MH'] + DATA_COLS[:1] + DATA_COLS[1:] + ['label']  # columns of the returned isochrone (EXT.select_iso_columns + extras)
INTERP_MAX_POINTS = 400                                                  # points per label segment after resampling
LABEL_MIN_WEIGHT = 0.5                                                   # a label missing from some corners is kept if >= this weight has it
EDGE_TOL = 1e-9                                                          # weight below this = exact node
META_KEY = b'starnest_grid'
# =================================================================


class GridUnavailable(Exception):
    """No usable local grid for this request (file missing, request outside the grid, mock refused). The caller falls back to CMD."""


class Grid:
    """In-memory grid: axes + {(i_age, j_mh): {column: ndarray}} + metadata."""
    def __init__(self, ages, mhs, tables, meta, path):
        self.ages, self.mhs, self.tables, self.meta, self.path = ages, mhs, tables, meta, path


_CACHE: dict = {}
_LOCK = threading.Lock()


def grid_available(path: str | None = None) -> bool:
    return os.path.isfile(path or GRID_PATH)


def load_grid(path: str | None = None, allow_mock: bool | None = None) -> Grid:
    """Read the parquet grid once per file version (path + mtime) and index it by node."""
    path = path or GRID_PATH
    allow_mock = ALLOW_MOCK if allow_mock is None else allow_mock
    if not os.path.isfile(path):
        raise GridUnavailable(f'grid file not found: {path}')
    key = (path, os.path.getmtime(path), allow_mock)
    with _LOCK:
        if key in _CACHE:
            return _CACHE[key]
        import pyarrow.parquet as pq
        t = pq.read_table(path)
        raw_meta = (t.schema.metadata or {}).get(META_KEY)
        meta = json.loads(raw_meta.decode()) if raw_meta else {}
        df = t.to_pandas()
        miss = [c for c in GRID_NODE_COLS + DATA_COLS + ['label'] if c not in df.columns]
        if miss:
            raise GridUnavailable(f'grid file lacks columns {miss}')
        if str(meta.get('source', '')).upper() == 'MOCK' and not allow_mock:
            raise GridUnavailable('grid is a MOCK test grid (set STARNEST_ALLOW_MOCK_GRID=1 only for tests)')
        ages = np.array(sorted(df['logAge_node'].round(4).unique()), float)
        mhs = np.array(sorted(df['MH_node'].round(4).unique()), float)
        ai = {round(a, 4): i for i, a in enumerate(ages)}
        mi = {round(m, 4): j for j, m in enumerate(mhs)}
        tables = {}
        for (a, m), g in df.groupby([df['logAge_node'].round(4), df['MH_node'].round(4)], sort=False):
            tables[(ai[a], mi[m])] = {c: g[c].to_numpy(float) for c in DATA_COLS} | {'label': g['label'].to_numpy(int)}
        _CACHE.clear()                                          # keep one grid in memory
        _CACHE[key] = Grid(ages, mhs, tables, meta, path)
        return _CACHE[key]


def _bracket(axis: np.ndarray, x: float):
    """(lower index, upper index, weight of upper) of x on a sorted axis; raises GridUnavailable outside it."""
    if x < axis[0] - EDGE_TOL or x > axis[-1] + EDGE_TOL:
        raise GridUnavailable(f'{x:.4f} outside the grid range {axis[0]:.3f}..{axis[-1]:.3f}')
    x = float(np.clip(x, axis[0], axis[-1]))
    if len(axis) == 1:
        return 0, 0, 0.0
    i1 = int(np.clip(np.searchsorted(axis, x, side='right'), 1, len(axis) - 1))
    i0 = i1 - 1
    w = (x - axis[i0]) / (axis[i1] - axis[i0])
    return i0, i1, float(w)


def _arc(logte: np.ndarray, logl: np.ndarray) -> np.ndarray:
    """Normalised path length 0..1 along a segment in the (logTe, logL) plane (equal spacing if the segment is a single point)."""
    n = len(logte)
    if n < 2:
        return np.zeros(n)
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(logte), np.diff(logl)))])
    return s / s[-1] if s[-1] > 0 else np.linspace(0, 1, n)


def _segment(tab: dict, label: int):
    m = tab['label'] == label
    return None if not m.any() else {c: tab[c][m] for c in DATA_COLS}


def interpolate_isochrone(grid: Grid, logage: float, mh: float):
    """Bilinear isochrone at (logage, mh). Returns (DataFrame with OUT_COLS, info dict)."""
    i0, i1, wa = _bracket(grid.ages, logage)
    j0, j1, wm = _bracket(grid.mhs, mh)
    corners = [((i, j), (1 - wa if i == i0 else wa) * (1 - wm if j == j0 else wm))
               for i in sorted({i0, i1}) for j in sorted({j0, j1})]
    if i0 == i1:                                                                    # degenerate axis: weights collapse onto the single node
        corners = [((i0, j), (1 - wm if j == j0 else wm)) for j in sorted({j0, j1})]
    corners = [(k, w) for k, w in corners if w > EDGE_TOL]
    info = dict(interpolated=len(corners) > 1, corners=[dict(logAge=float(grid.ages[k[0]]), MH=float(grid.mhs[k[1]]), weight=round(w, 6))
                                                       for k, w in corners], labels_partial=[])
    for k, _ in corners:
        if k not in grid.tables:
            raise GridUnavailable(f'grid node {k} missing')
    if len(corners) == 1:
        tab = grid.tables[corners[0][0]]
        out = pd.DataFrame({c: tab[c] for c in DATA_COLS})
        out['label'] = tab['label']
    else:
        labels = sorted({int(l) for k, _ in corners for l in np.unique(grid.tables[k]['label'])})
        parts = []
        for lab in labels:
            segs = [(_segment(grid.tables[k], lab), w) for k, w in corners]
            have = [(s, w) for s, w in segs if s is not None]
            wsum = sum(w for _, w in have)
            if wsum < LABEL_MIN_WEIGHT:                                              # phase mostly absent at this age / metallicity
                info['labels_partial'].append(lab)
                continue
            if len(have) < len(segs):
                info['labels_partial'].append(lab)
            n = int(np.clip(max(len(s['logTe']) for s, _ in have), 2, INTERP_MAX_POINTS))
            u = np.linspace(0.0, 1.0, n)
            acc = {c: np.zeros(n) for c in DATA_COLS}
            for s, w in have:
                su = _arc(s['logTe'], s['logL'])
                for c in DATA_COLS:
                    v = np.log10(np.clip(s[c], 1e-12, None)) if c in LOG_COLS else s[c]
                    acc[c] += (w / wsum) * np.interp(u, su, v)                       # weights renormalised over the corners that have this label
            p = pd.DataFrame({c: (10 ** acc[c] if c in LOG_COLS else acc[c]) for c in DATA_COLS})
            p['label'] = lab
            parts.append(p)
        if not parts:
            raise GridUnavailable('no evolutionary phase survived the interpolation')
        out = pd.concat(parts, ignore_index=True)
    out.insert(0, 'MH', float(mh))
    out.insert(0, 'logAge', float(logage))
    return out[OUT_COLS].round(6), info


def apply_av_shift(iso: pd.DataFrame, av: float, coeffs: dict) -> pd.DataFrame:
    """Add interstellar extinction to an A_V = 0 isochrone: each Gaia band moves by (A_band / A_V) * A_V (constant coefficients, as CMD 3.8)."""
    out = iso.copy()
    if av:
        out['Gmag'] = out['Gmag'] + coeffs['AG_AV'] * av
        out['G_BPmag'] = out['G_BPmag'] + coeffs['ABP_AV'] * av
        out['G_RPmag'] = out['G_RPmag'] + coeffs['ARP_AV'] * av
    return out.round(6)


def _provenance_text(iso: pd.DataFrame, logage: float, mh: float, av: float, info: dict, meta: dict, coeffs: dict) -> bytes:
    """Re-parseable '#'-commented file (extinction_isochrone.parse_parsec reads it) recording exactly how this isochrone was produced."""
    head = ['# ---- StarNest web app: PARSEC file provenance (LOCAL GRID) ----',
            f'# produced_utc    : {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}',
            f'# service         : local PARSEC grid built from {meta.get("service", "CMD 3.8")} on {meta.get("built_utc", "?")}',
            f'# photometric sys : {meta.get("photsys", "Gaia EDR3 = Gaia DR3 photometry, Vega")}',
            f'# logAge / [M/H]  : {logage:.4f} / {mh:+.4f}  ({"bilinear blend of 4 nodes" if info["interpolated"] else "exact grid node"})',
            f'# interstellar A_V: {av:.4f} mag, applied analytically, A_G/A_V={coeffs["AG_AV"]:.4f} A_BP/A_V={coeffs["ABP_AV"]:.4f} A_RP/A_V={coeffs["ARP_AV"]:.4f}'
            f' ({coeffs.get("source", "?")})',
            f'# nodes used      : {json.dumps(info["corners"])}',
            f'# grid source     : {meta.get("source", "?")}']
    cols = list(iso.columns)
    body = '\n'.join(' '.join(f'{v:.6g}' if isinstance(v, float) else str(v) for v in row) for row in iso.itertuples(index=False))
    return ('\n'.join(head) + '\n# ' + ' '.join(cols) + '\n' + body + '\n').encode()


def get_parsec_from_grid(logage: float, mh: float, av: float, default_coeffs: dict, path: str | None = None) -> dict:
    """Same return shape as extinction_isochrone.get_parsec (iso, raw, params, info) but served from the local grid. Raises GridUnavailable."""
    g = load_grid(path)
    iso0, info = interpolate_isochrone(g, float(logage), float(mh))
    cm = g.meta.get('ext_coeff')
    coeffs = dict(AG_AV=cm['AG_AV'], ABP_AV=cm['ABP_AV'], ARP_AV=cm['ARP_AV'], source='measured from CMD 3.8 at build time') if cm else \
        dict(default_coeffs, source='StarNest defaults (grid has no measured coefficients)')
    iso = apply_av_shift(iso0, float(av), coeffs)
    info.update(source='local_grid', service='local PARSEC grid (CMD 3.8)', ssl_verified=None, source_url=None, grid_path=g.path,
                grid_source=g.meta.get('source'), grid_built_utc=g.meta.get('built_utc'), av_shift='analytic', ext_coeff=coeffs)
    params = dict(isoc_lagelow=f'{logage:.4f}', isoc_metlow=f'{mh:.4f}', extinction_av=f'{av:.4f}', extinction_curve='cardelli',
                  extinction_coeff='constant (analytic shift)', photsys_file=g.meta.get('photsys_file'), photsys_version=g.meta.get('photsys_version'))
    return dict(iso=iso, raw=_provenance_text(iso, logage, mh, av, info, g.meta, coeffs), params=params, info=info)
