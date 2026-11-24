"""
build_grid.py -- download the local PARSEC grid once, so the web app never needs CMD at run time (Phase 2.2).

USAGE   python build_grid.py   (or: python -m starnest.build_grid)                      # real build from CMD 3.8 (needs internet; ~11 requests, resumable)
        python build_grid.py --out data/parsec_grid.parquet --resume
        python build_grid.py --mock --out /tmp/mock_grid.parquet     # synthetic grid for tests / CI ONLY (metadata source='MOCK', refused at run time)

TERMINOLOGY
  Node           : one isochrone at an exact (logAge, [M/H]).  Grid = all nodes at A_V = 0 in the Gaia EDR3 (= DR3) Vega system.
  Reference run  : two single isochrones requested WITH extinction (A_V = REF_AV) at (REF_LOGAGE, REF_MH). Matching them to the A_V = 0 one at
                   equal initial mass gives the per-band shift A_band / A_V that CMD applies (extinction_coeff='constant'); stored in the file.
  Interp check   : an OFF-grid isochrone is requested from CMD and compared with the bilinear interpolation of the stored nodes (CMD-plane
                   distance and mass error on the main sequence); the numbers are written into the grid metadata.

WORKFLOW
  1) for each [M/H] node: ONE CMD request returns every logAge node (cached on disk -> --resume continues after an interruption);
  2) parse + validate (rows, finite magnitudes, logAge / [M/H] echo) and stack into one table with exact node coordinates;
  3) measure the extinction coefficients from the reference run (and their linearity at 3 x A_V);
  4) interpolation check against direct CMD isochrones at INTERP_CHECK_POINTS;
  5) write parquet (zstd) with the metadata JSON embedded under schema key 'starnest_grid'.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

# ===================== GLOBAL CONFIGURATION =====================
_ROOT = os.environ.get('STARNEST_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_PATH = os.path.join(_ROOT, 'data', 'parsec_grid.parquet')
CACHE_DIR = os.path.join(_ROOT, 'data', '_grid_cache')   # raw per-[M/H] downloads (resume)
AGE_LO, AGE_HI, DLOGAGE = 6.6, 10.1, 0.05          # log10(age/yr) nodes
MH_LO, MH_HI, DMH = -0.5, 0.5, 0.1                 # [M/H] nodes
REF_LOGAGE, REF_MH, REF_AV = 8.0, 0.0, 1.0         # reference isochrone for the extinction coefficients (also checked at 3 x REF_AV)
INTERP_CHECK_POINTS = [(8.025, 0.05), (7.125, -0.25), (9.575, 0.35)]    # off-grid (logAge, [M/H]) requested from CMD for the interpolation check
CHECK_MASS_RANGE = (0.3, 1.5)                      # Msun: main-sequence window used for the mass-error statistic
CHECK_TOL_MAG = 0.03                               # warn if the median CMD-plane error of the interpolation exceeds this (mag)
FLOAT_DTYPE = 'float32'
# =================================================================

from . import extinction_isochrone as EXT   # noqa: E402  (after config so a user can edit the block first)
from . import parsec_grid as PG             # noqa: E402


# ------------------------------ mock CMD (tests / CI only) ------------------------------
MOCK_COEFF = dict(AG_AV=0.836, ABP_AV=1.084, ARP_AV=0.634)   # deliberately different from the StarNest defaults: proves the measured values are used


def mock_isochrone(logage: float, mh: float, av: float = 0.0) -> pd.DataFrame:
    """Smooth analytic isochrone (NOT PARSEC): main sequence up to a turn-off mass that falls with age, then SGB + RGB; PMS only for young ages.
    Point counts depend on age on purpose, so the interpolation is tested with unequal segment lengths."""
    mto = 10 ** ((10.0 - logage) / 2.5)                                  # t_MS ~ M^-2.5 (1 Msun at 10 Gyr)
    n_ms = int(150 + 30 * (logage - 6.6))
    parts = []
    m_zams = 0.15 + 0.25 * np.clip((7.8 - logage) / 1.2, 0, 1)             # PMS -> MS joint mass moves continuously with age (PMS gone by logAge 7.8)
    if m_zams > 0.1501:
        parts.append((np.geomspace(0.15, m_zams, 25, endpoint=False), 0))
    m = np.geomspace(m_zams, mto, n_ms)
    parts.append((m, 1))
    rows = []
    for m, lab in parts:
        lm = np.log10(m)
        logl = 3.8 * lm + 0.1 * mh
        logte = 3.55 + 0.30 * lm - 0.05 * lm ** 2
        rows.append(pd.DataFrame(dict(Mini=m, label=lab, logL=logl, logTe=logte)))
    ms = rows[-1]
    t = np.linspace(0, 1, 40)[1:]
    sgb = pd.DataFrame(dict(Mini=mto * (1 + 0.005 * t), label=2, logL=ms['logL'].iloc[-1] + 0.3 * t, logTe=ms['logTe'].iloc[-1] - 0.06 * t))
    t = np.linspace(0, 1, 80)[1:]
    rgb = pd.DataFrame(dict(Mini=mto * (1.005 + 0.002 * t), label=3, logL=sgb['logL'].iloc[-1] + 1.8 * t, logTe=sgb['logTe'].iloc[-1] - 0.08 * t))
    d = pd.concat(rows + [sgb, rgb], ignore_index=True)
    d['Mass'] = d['Mini'] * 0.99
    d['logg'] = 4.4 - 0.6 * (d['logL'] - 3.8 * np.log10(d['Mini'])) - 2 * (d['logTe'] - 3.7)
    mbol = 4.74 - 2.5 * d['logL']
    bc_g = -0.1 + 4.0 * (d['logTe'] - 3.75) ** 2
    d['Gmag'] = mbol - bc_g + 0.05 * mh + MOCK_COEFF['AG_AV'] * av
    col = 0.9 - 4.5 * (d['logTe'] - 3.68) + 0.06 * mh
    d['G_BPmag'] = d['Gmag'] + 0.4 * col + (MOCK_COEFF['ABP_AV'] - MOCK_COEFF['AG_AV']) * av
    d['G_RPmag'] = d['Gmag'] - 0.6 * col + (MOCK_COEFF['ARP_AV'] - MOCK_COEFF['AG_AV']) * av
    d.insert(0, 'MH', mh)
    d.insert(0, 'logAge', logage)
    return d[['logAge', 'MH', 'Mini', 'Mass', 'logL', 'logTe', 'logg', 'label', 'Gmag', 'G_BPmag', 'G_RPmag']].round(6)


def _as_cmd_text(df: pd.DataFrame) -> bytes:
    head = '# ' + ' '.join(df.columns)
    return (head + '\n' + '\n'.join(' '.join(f'{v:.6g}' for v in r) for r in df.itertuples(index=False)) + '\n').encode()


class MockFetcher:
    """Same two calls as the real fetcher, answered analytically (no network)."""
    service = 'MOCK (analytic isochrones, not PARSEC)'

    def age_grid(self, lo, hi, dl, mh):
        ages = np.round(np.arange(lo, hi + dl / 2, dl), 4)
        return _as_cmd_text(pd.concat([mock_isochrone(a, mh) for a in ages], ignore_index=True))

    def single(self, logage, mh, av):
        return _as_cmd_text(mock_isochrone(logage, mh, av))


class CMDFetcher:
    """Real CMD 3.8 requests through extinction_isochrone (TLS / proxy / retry settings of that module apply)."""
    service = f'CMD {EXT.CMD_VERSION}'

    def age_grid(self, lo, hi, dl, mh):
        return EXT.fetch_age_grid_raw(lo, hi, dl, mh)

    def single(self, logage, mh, av):
        return EXT.fetch_parsec_raw(logage, mh, av)[0]


# ------------------------------ build steps ------------------------------
def _parse(raw: bytes, mh_req: float) -> pd.DataFrame:
    iso = EXT.select_iso_columns(EXT.parse_parsec(raw))
    if len(iso) < 20 or not np.isfinite(iso[['Gmag', 'G_BPmag', 'G_RPmag']].to_numpy(float)).all():
        raise RuntimeError(f'[M/H]={mh_req}: table has {len(iso)} rows or non-finite magnitudes')
    if abs(float(iso['MH'].median()) - mh_req) > 0.15:               # echo check: the server answered the metallicity we asked for
        raise RuntimeError(f"[M/H] requested {mh_req}, table says {float(iso['MH'].median()):.3f}")
    return iso


def download_nodes(fetch, ages, mhs, cache_dir: str, resume: bool) -> pd.DataFrame:
    """One request per [M/H] node (all ages); raw files cached so an interrupted build resumes."""
    os.makedirs(cache_dir, exist_ok=True)
    frames = []
    for mh in mhs:
        f = os.path.join(cache_dir, f'mh_{mh:+.2f}.txt')
        if resume and os.path.isfile(f):
            raw = open(f, 'rb').read()
            print(f'[M/H]={mh:+.2f}: cached', file=sys.stderr)
        else:
            print(f'[M/H]={mh:+.2f}: requesting {len(ages)} isochrones ...', file=sys.stderr)
            raw = fetch.age_grid(float(ages[0]), float(ages[-1]), DLOGAGE, float(mh))
            open(f, 'wb').write(raw)
        iso = _parse(raw, float(mh))
        iso['logAge_node'] = (np.round(iso['logAge'] / DLOGAGE) * DLOGAGE).round(4)       # snap printed ages (6.6, 6.65, ...) to the node
        iso['MH_node'] = round(float(mh), 4)
        got = set(iso['logAge_node'].unique())
        missing = [a for a in ages if round(float(a), 4) not in got]
        if missing:
            raise RuntimeError(f'[M/H]={mh:+.2f}: ages missing from the CMD answer: {missing[:5]}...')
        frames.append(iso)
    return pd.concat(frames, ignore_index=True)


def _ms_match(iso0: pd.DataFrame, iso_av: pd.DataFrame):
    """Main-sequence points of both isochrones paired by initial mass (interpolating the A_V=0 one at the A_V>0 masses)."""
    a, b = iso0[iso0['label'] == 1].sort_values('Mini'), iso_av[iso_av['label'] == 1].sort_values('Mini')
    x = b['Mini'].to_numpy(float)
    ok = (x >= a['Mini'].min()) & (x <= a['Mini'].max())
    return a, b, x, ok


def measure_ext_coeff(fetch) -> dict:
    """A_band / A_V from CMD itself: median magnitude shift of the matched main sequence for A_V = REF_AV and 3 x REF_AV."""
    iso0 = _parse(fetch.single(REF_LOGAGE, REF_MH, 0.0), REF_MH)
    res = {}
    for av in (REF_AV, 3 * REF_AV):
        isa = _parse(fetch.single(REF_LOGAGE, REF_MH, av), REF_MH)
        a, b, x, ok = _ms_match(iso0, isa)
        res[av] = {k: float(np.median(b[col].to_numpy(float)[ok] - np.interp(x[ok], a['Mini'], a[col]))) / av
                   for k, col in (('AG_AV', 'Gmag'), ('ABP_AV', 'G_BPmag'), ('ARP_AV', 'G_RPmag'))}
    lo, hi = res[REF_AV], res[3 * REF_AV]
    out = dict(lo)
    out.update(ref_logAge=REF_LOGAGE, ref_MH=REF_MH, ref_AV=REF_AV, linearity_max_dev=max(abs(lo[k] - hi[k]) for k in lo),
               method='median matched-main-sequence shift (equal initial mass) / A_V')
    return out


def interp_check(grid: PG.Grid, fetch) -> list:
    """Interpolated vs direct-CMD isochrone at off-grid points. Error = CMD-plane distance (mag) from each interpolated main-sequence point to the
    direct polyline, and relative mass error at equal G within CHECK_MASS_RANGE."""
    out = []
    for age, mh in INTERP_CHECK_POINTS:
        try:
            ip, _ = PG.interpolate_isochrone(grid, age, mh)
        except PG.GridUnavailable as e:
            out.append(dict(logAge=age, MH=mh, error=str(e)))
            continue
        tr = _parse(fetch.single(age, mh, 0.0), mh)
        t1, i1 = tr[tr['label'] == 1].sort_values('Mini'), ip[ip['label'] == 1].sort_values('Mini')
        tc, tg = (t1['G_BPmag'] - t1['G_RPmag']).to_numpy(float), t1['Gmag'].to_numpy(float)
        ex, ey = EXT._densify(tc, tg)
        ic, ig = (i1['G_BPmag'] - i1['G_RPmag']).to_numpy(float), i1['Gmag'].to_numpy(float)
        dist = np.array([np.min(np.hypot(ex - c, ey - g)) for c, g in zip(ic, ig)])
        sel = (i1['Mass'] >= CHECK_MASS_RANGE[0]) & (i1['Mass'] <= CHECK_MASS_RANGE[1])
        order = np.argsort(tg)
        mt = np.interp(i1['Gmag'][sel], tg[order], t1['Mass'].to_numpy(float)[order])
        rel = np.abs(i1['Mass'][sel].to_numpy(float) / mt - 1) if sel.any() else np.array([np.nan])
        out.append(dict(logAge=age, MH=mh, median_cmd_err_mag=float(np.median(dist)), max_cmd_err_mag=float(dist.max()),
                        median_mass_err_rel=float(np.nanmedian(rel)), max_mass_err_rel=float(np.nanmax(rel)), n=int(len(dist))))
    return out


def write_grid(df: pd.DataFrame, meta: dict, out: str):
    import pyarrow as pa
    import pyarrow.parquet as pq
    cols = PG.GRID_NODE_COLS + PG.DATA_COLS + ['label']
    d = df[cols].copy()
    d[PG.DATA_COLS] = d[PG.DATA_COLS].astype(FLOAT_DTYPE)
    d['label'] = d['label'].astype('int8')
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    t = pa.Table.from_pandas(d, preserve_index=False).replace_schema_metadata({PG.META_KEY: json.dumps(meta).encode()})
    pq.write_table(t, out, compression='zstd')


def build(fetch, out: str = OUT_PATH, cache_dir: str = CACHE_DIR, resume: bool = True, check: bool = True,
          age_lo=AGE_LO, age_hi=AGE_HI, mh_lo=MH_LO, mh_hi=MH_HI) -> dict:
    ages = np.round(np.arange(age_lo, age_hi + DLOGAGE / 2, DLOGAGE), 4)
    mhs = np.round(np.arange(mh_lo, mh_hi + DMH / 2, DMH), 4)
    df = download_nodes(fetch, ages, mhs, cache_dir, resume)
    meta = dict(source='MOCK' if isinstance(fetch, MockFetcher) else 'CMD', service=fetch.service, cmd_version=EXT.CMD_VERSION,
                photsys='Gaia EDR3 = Gaia DR3 photometry, Vega (YBCnewVega)', photsys_file=EXT.PHOTSYS_FILE, photsys_version=EXT.PHOTSYS_VERSION,
                tracks=f"{EXT.CMD_BASE_PARAMS['track_parsec']} + {EXT.CMD_BASE_PARAMS['track_colibri']}",
                logAge=[float(ages[0]), float(ages[-1]), DLOGAGE], MH=[float(mhs[0]), float(mhs[-1]), DMH], av=0.0, n_rows=int(len(df)),
                built_utc=datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S'))
    meta['ext_coeff'] = measure_ext_coeff(fetch)
    write_grid(df, meta, out)
    if check:                                                       # needs the file on disk (reads it back exactly as the app will)
        g = PG.load_grid(out, allow_mock=True)       # the builder may read its own mock output
        meta['interp_check'] = interp_check(g, fetch)
        bad = [c for c in meta['interp_check'] if c.get('median_cmd_err_mag', 0) > CHECK_TOL_MAG]
        meta['interp_check_ok'] = not bad
        write_grid(df, meta, out)
        if bad:
            print(f'WARNING: interpolation error above {CHECK_TOL_MAG} mag at {[(c["logAge"], c["MH"]) for c in bad]}; consider a finer grid', file=sys.stderr)
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser(description='Build the local PARSEC grid (A_V = 0, Gaia EDR3 Vega).')
    ap.add_argument('--out', default=OUT_PATH)
    ap.add_argument('--cache-dir', default=CACHE_DIR)
    ap.add_argument('--resume', action='store_true', help='reuse raw per-[M/H] downloads already in the cache dir')
    ap.add_argument('--mock', action='store_true', help='synthetic grid for tests / CI (never use for science)')
    ap.add_argument('--no-check', action='store_true', help='skip the off-grid interpolation check (saves 3 CMD requests)')
    a = ap.parse_args(argv)
    meta = build(MockFetcher() if a.mock else CMDFetcher(), a.out, a.cache_dir, a.resume, not a.no_check)
    print(json.dumps({k: v for k, v in meta.items() if k != 'interp_check'}, indent=1))
    for c in meta.get('interp_check', []):
        print('interp check:', c)


if __name__ == '__main__':
    main()
