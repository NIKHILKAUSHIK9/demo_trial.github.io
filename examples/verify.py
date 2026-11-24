"""
verify.py -- compare a StarNest batch run with an INDEPENDENT literature table and flag out-of-tolerance values.

TERMINOLOGY
  summary.csv      : written by `python -m starnest run` (one row per cluster).
  literature table : a CSV YOU supply, one row per cluster, with a non-empty `source` (paper / catalogue name). It must not be built from
                     StarNest's own earlier output: agreement with yourself verifies nothing.
  tolerance        : the largest |pipeline - literature| still counted as agreement. Fixed in GLOBAL CONFIGURATION below BEFORE the run;
                     if you change a tolerance after seeing the result, say so where you report it.
  PASS / FAIL / MISSING / SKIP : inside tolerance / outside / literature has a value but the pipeline has none (counts as a failure) /
                     the literature table has no value (not checked).

WORKFLOW
  1) read summary.csv and the literature table (exit 2 on a missing file, column or `source`);
  2) match rows by cluster name (case, spaces and punctuation ignored: "UPK 545" == "UPK_545");
  3) summary-level checks: N, d, E(B-V); member-level checks (median of <run dir>/<cluster slug>/members.csv): pmRA, pmDE, Plx;
  4) print one line per check with its tolerance, optionally write them to --csv; exit 0 if nothing failed, 1 otherwise.

NOT CHECKED, ON PURPOSE: logAge is the age you typed into the config (an input, not a result); Gamma, masses, MSR and velocities have no
literature value that is comparable without matching the method, and masses are meaningless when the run used a MOCK isochrone grid.
E(B-V) from the IRSA/SFD map is the total line-of-sight dust, an upper limit for nearby clusters, so a FAIL there for a close cluster is a
real finding, not a bug.

USAGE   python examples/verify.py starnest_output/summary.csv my_literature_table.csv [--csv checks.csv]
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys

import numpy as np
import pandas as pd

# ===================== GLOBAL CONFIGURATION =====================
D_TOL_FRAC = 0.05               # |d - d_lit| / d_lit; both are 1000/parallax-type distances, so a few per cent is the realistic floor
N_TOL_FRAC = 0.30               # |N - N_lit| / N_lit; member counts depend on the membership method, hence wide
EBV_TOL_MAG = 0.05              # |E(B-V) - E(B-V)_lit| in mag
PM_NSIGMA = 3.0                 # proper-motion tolerance = max(PM_NSIGMA * pm_err_masyr, PM_TOL_FLOOR_MASYR)
PM_TOL_FLOOR_MASYR = 0.15       # used when the literature gives no error (or a smaller one)
PLX_NSIGMA = 3.0                # parallax tolerance = max(PLX_NSIGMA * plx_err_mas, PLX_TOL_FLOOR_MAS)
PLX_TOL_FLOOR_MAS = 0.03        # about the Gaia DR3 zero-point scale; the pipeline default PLX_ZP = 0 is not corrected for
MEMBERS_FILENAME = 'members.csv'                  # per-cluster member list written by the CLI (starnest/cli.py DOWNLOADS)
SLUG_PATTERN = r'[^A-Za-z0-9_.-]+'                # same rule the CLI uses for the per-cluster folder name
LIT_REQUIRED = ('cluster', 'source')     # optional value columns: d_pc, plx_mas, plx_err_mas, pmra_masyr, pmde_masyr, pm_err_masyr, ebv, n_members
# =================================================================

CHECK_COLUMNS = ['cluster', 'quantity', 'pipeline', 'literature', 'tolerance', 'status', 'source']


def _key(name: str) -> str:
    """'UPK 545' / 'UPK_545' / 'upk-545' -> 'upk545'."""
    return re.sub(r'[^a-z0-9]', '', str(name).lower())


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float('nan')
    return v if math.isfinite(v) else float('nan')


def _judge(pipe: float, lit: float, tol: float) -> str:
    """SKIP if the literature has no value, MISSING if only the pipeline lacks one, else PASS / FAIL on |pipe - lit| <= tol."""
    if not math.isfinite(lit):
        return 'SKIP'
    if not math.isfinite(pipe):
        return 'MISSING'
    return 'PASS' if abs(pipe - lit) <= tol else 'FAIL'


def load_tables(summary_path: str, lit_path: str):
    """Read both CSVs; SystemExit(2) with a readable message if anything the comparison needs is absent."""
    for p in (summary_path, lit_path):
        if not os.path.isfile(p):
            raise SystemExit(f'verify: file not found: {p}')
    summary, lit = pd.read_csv(summary_path), pd.read_csv(lit_path)
    miss = [c for c in LIT_REQUIRED if c not in lit.columns]
    if miss:
        raise SystemExit(f'verify: literature table needs columns {list(LIT_REQUIRED)}; missing {miss}')
    if lit['source'].isna().any() or (lit['source'].astype(str).str.strip() == '').any():
        raise SystemExit('verify: every literature row needs a non-empty `source` (paper or catalogue). An unsourced reference verifies nothing.')
    if 'cluster' not in summary.columns:
        raise SystemExit(f'verify: {summary_path} has no `cluster` column (is it a starnest summary.csv?)')
    if lit['cluster'].map(_key).duplicated().any():
        raise SystemExit('verify: duplicate cluster names in the literature table')
    return summary, lit


def member_medians(run_dir: str, cluster: str) -> dict:
    """Median pmRA, pmDE, Plx of the cluster's members.csv (NaN when the file or a column is absent)."""
    path = os.path.join(run_dir, re.sub(SLUG_PATTERN, '_', str(cluster)), MEMBERS_FILENAME)
    out = {'pmRA': float('nan'), 'pmDE': float('nan'), 'Plx': float('nan')}
    if not os.path.isfile(path):
        return out
    m = pd.read_csv(path)
    for c in out:
        if c in m.columns:
            out[c] = float(np.nanmedian(pd.to_numeric(m[c], errors='coerce'))) if len(m) else float('nan')
    return out


def run_checks(summary: pd.DataFrame, lit: pd.DataFrame, run_dir: str) -> pd.DataFrame:
    """One row per (cluster, quantity). A literature cluster absent from the summary is reported as MISSING for every value it gives."""
    by_key = {_key(r['cluster']): r for _, r in summary.iterrows()}
    rows = []
    for _, L in lit.iterrows():
        S = by_key.get(_key(L['cluster']))
        src = str(L['source'])
        ok_run = S is not None and str(S.get('status', 'ok')) == 'ok'
        med = member_medians(run_dir, S['cluster']) if S is not None else {}
        d_lit, n_lit, ebv_lit = _num(L.get('d_pc')), _num(L.get('n_members')), _num(L.get('ebv'))
        pm_tol = max(PM_NSIGMA * (_num(L.get('pm_err_masyr')) if math.isfinite(_num(L.get('pm_err_masyr'))) else 0.0), PM_TOL_FLOOR_MASYR)
        plx_err = _num(L.get('plx_err_mas'))
        plx_tol = max(PLX_NSIGMA * (plx_err if math.isfinite(plx_err) else 0.0), PLX_TOL_FLOOR_MAS)
        spec = [                                                    # (quantity, pipeline value, literature value, tolerance)
            ('N',     _num(S.get('N')) if S is not None else float('nan'),         n_lit,                    N_TOL_FRAC * n_lit if math.isfinite(n_lit) else float('nan')),
            ('d_pc',  _num(S.get('d')) if S is not None else float('nan'),         d_lit,                    D_TOL_FRAC * d_lit if math.isfinite(d_lit) else float('nan')),
            ('E(B-V)', _num(S.get('E(B-V)')) if S is not None else float('nan'),   ebv_lit,                  EBV_TOL_MAG),
            ('pmRA',  med.get('pmRA', float('nan')),                               _num(L.get('pmra_masyr')), pm_tol),
            ('pmDE',  med.get('pmDE', float('nan')),                               _num(L.get('pmde_masyr')), pm_tol),
            ('Plx',   med.get('Plx', float('nan')),                                _num(L.get('plx_mas')),   plx_tol),
        ]
        for q, pipe, lv, tol in spec:
            status = _judge(pipe, lv, tol)
            if status == 'PASS' and not ok_run:
                status = 'FAIL'                                      # a partial / failed run must not pass on numbers it may not have produced
            rows.append(dict(cluster=L['cluster'], quantity=q, pipeline=pipe, literature=lv, tolerance=tol, status=status, source=src))
    return pd.DataFrame(rows, columns=CHECK_COLUMNS)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description='Compare a StarNest summary.csv with an independent literature table.')
    ap.add_argument('summary'); ap.add_argument('literature'); ap.add_argument('--csv', help='also write the checks to this CSV (index=False)')
    a = ap.parse_args(argv)
    summary, lit = load_tables(a.summary, a.literature)
    checks = run_checks(summary, lit, os.path.dirname(os.path.abspath(a.summary)))
    print(f'tolerances: N +-{N_TOL_FRAC:.0%}, d +-{D_TOL_FRAC:.0%}, E(B-V) +-{EBV_TOL_MAG} mag, PM max({PM_NSIGMA} sigma, {PM_TOL_FLOOR_MASYR} mas/yr), '
          f'Plx max({PLX_NSIGMA} sigma, {PLX_TOL_FLOOR_MAS} mas)', file=sys.stderr)
    if 'isochrone_source' in summary.columns and summary['isochrone_source'].astype(str).str.upper().str.contains('MOCK').any():
        print('WARNING: at least one run used a MOCK isochrone grid: masses, Gamma and A_V from that run are invalid (not checked here).', file=sys.stderr)
    with pd.option_context('display.width', 200, 'display.max_columns', None, 'display.float_format', '{:.4g}'.format):
        print(checks.drop(columns=['source']).to_string(index=False))
    n = checks['status'].value_counts().to_dict()
    print(f"\n{n.get('PASS', 0)} pass, {n.get('FAIL', 0)} fail, {n.get('MISSING', 0)} missing, {n.get('SKIP', 0)} skipped", file=sys.stderr)
    if a.csv:
        checks.to_csv(a.csv, index=False)
    return 1 if (n.get('FAIL', 0) or n.get('MISSING', 0)) else 0


if __name__ == '__main__':
    sys.exit(main())
