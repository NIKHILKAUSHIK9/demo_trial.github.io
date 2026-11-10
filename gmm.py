"""
gmm.py -- Stage 1: Gaia quality cuts + GMM membership (web-app-ready, no UI code).

TERMINOLOGY
  RUWE / gofAL / Nper : Gaia astrometric quality flags (lower RUWE/gofAL = cleaner single-star fit).
  GMM                 : Gaussian Mixture Model; fits K Gaussians in (pmRA, pmDE, Plx) space.
  PMemb               : posterior probability that a star belongs to the CLUSTER component.
  Cluster component   : the Gaussian with the smallest det(covariance) = tightest blob in PM+parallax.
  BIC                 : Bayesian Information Criterion; lower = better fit penalised for extra components.
  MSS                 : Mean Separation Score; mean over features of (std_nonmember - std_member)/max(std).

WORKFLOW
  load_raw -> validate_columns -> apply_quality_cuts -> run_gmm (distance window -> standardise ->
  fit K=2..6 by BIC -> pick tightest component -> PMemb -> members/non-members -> MSS) -> make_gmm_figures
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from sklearn.mixture import GaussianMixture

# ===================== GLOBAL CONFIGURATION =====================
COL = dict(                                   # VizieR Gaia DR3 column names (do not rename)
    ra='RA_ICRS', dec='DE_ICRS', plx='Plx', e_plx='e_Plx',
    pmra='pmRA', pmde='pmDE', e_pmra='e_pmRA', e_pmde='e_pmDE',
    gmag='Gmag', bp='BPmag', rp='RPmag', color='BP-RP',
    ruwe='RUWE', gofal='gofAL', nper='Nper',
)
QUALITY = dict(RUWE_MAX=1.6, E_PM_MAX=0.1, PLX_SNR_MIN=5, GOFAL_MAX=3, NPER_MIN=10)   # strict </> as in original
FEATURES = ['pmRA', 'pmDE', 'Plx']
DIST_LIT_PC = 341.0          # literature distance (pc)
HALF_WIDTH_PC = 195.5        # distance window half-width (pc)
CUTOFF = 0.5                 # member: PMemb >= CUTOFF ; non-member: PMemb <= 1 - CUTOFF
K_DEFAULT = None             # None -> choose K by BIC over K_RANGE (report the K-robustness table) ; 2 = Colab script (cluster + field, no BIC) ; other int -> force that K
K_RANGE = (2, 6)
N_INIT = 5
OPTIMIZER = 'converged'          # 'colab' = sklearn defaults tol=1e-3, max_iter=100 (reproduces the Colab script) ; 'converged' = GMM_TOL / GMM_MAX_ITER below
GMM_TOL = 1e-6               # 'converged' mode: EM convergence tol on mean log-lik. sklearn default 1e-3 stops EM early on overlapping components
GMM_MAX_ITER = 1000          # sklearn default 100 often ends unconverged (ConvergenceWarning)
RANDOM_STATE = 42
PLX_ZP = 0.0                 # parallax zero-point (mas); distance = 1000/(Plx - PLX_ZP). DR3 global ~ -0.017
MIN_STARS_PER_COMPONENT = 10
MIN_COMPONENT_STARS = 5      # a K whose fitted component holds fewer stars is a collapsed (singular) likelihood spike -> excluded from BIC choice
REQUIRED_COLUMNS = [COL[k] for k in ('plx', 'e_plx', 'pmra', 'pmde', 'e_pmra', 'e_pmde',
                                     'gmag', 'color', 'ruwe', 'gofal', 'nper')]
OUT_DIR = './starnest_output'  # only used by save_* helpers
# =================================================================


@dataclass
class GMMResult:
    working: pd.DataFrame                  # distance-window stars incl. PMemb (NaN where features missing)
    members: pd.DataFrame
    non_members: pd.DataFrame
    k: int
    bic: dict = field(default_factory=dict)
    cluster_component: int = -1
    mss: float = float('nan')
    summary: dict = field(default_factory=dict)
    components: list = field(default_factory=list)   # per-component diagnostics for the chosen K
    k_table: list = field(default_factory=list)      # per-K robustness rows: k, bic, members, non_members, min_comp_stars


def validate_columns(df: pd.DataFrame, required=REQUIRED_COLUMNS) -> list:
    """Return list of missing required columns ([] = OK). Call before anything else on uploads."""
    return [c for c in required if c not in df.columns]


def load_raw(file_or_path) -> pd.DataFrame:
    """Read a VizieR CSV (path or uploaded file-like). latin-1 avoids UnicodeDecodeError; all cols -> numeric."""
    if hasattr(file_or_path, 'seek'):
        file_or_path.seek(0)
    df = pd.read_csv(file_or_path, encoding='latin-1')
    return df.apply(pd.to_numeric, errors='coerce')


def apply_quality_cuts(df: pd.DataFrame, quality: dict = QUALITY) -> pd.DataFrame:
    """RUWE<, e_pm<, Plx/e_Plx>, gofAL<, Nper> cuts (NaN rows are dropped, as in the original)."""
    miss = validate_columns(df)
    if miss:
        raise ValueError(f'Missing required columns: {miss}')
    c = COL
    keep = ((df[c['ruwe']] < quality['RUWE_MAX'])
            & (df[c['e_pmra']] < quality.get('E_PMRA_MAX', quality['E_PM_MAX']))
            & (df[c['plx']] / df[c['e_plx']] > quality['PLX_SNR_MIN'])
            & (df[c['e_pmde']] < quality.get('E_PMDE_MAX', quality['E_PM_MAX']))
            & (df[c['gofal']] < quality['GOFAL_MAX'])
            & (df[c['nper']] > quality['NPER_MIN']))
    return df.loc[keep].copy()


def _standardise(X: pd.DataFrame) -> pd.DataFrame:
    """(x - median)/std per column, as in the original."""
    sd = X.std(ddof=0).replace(0, np.nan)
    if sd.isna().any():
        raise ValueError('Zero variance in a GMM feature column.')
    return (X - X.median()) / sd


def _mss(members, non_members, features, eps=1e-7) -> float:
    if len(members) < 2 or len(non_members) < 2:
        return 0.0
    m = []
    for f in features:
        sm, sn = np.std(members[f]), np.std(non_members[f])
        m.append((sn - sm) / max(sm, sn, eps))
    return float(np.mean(m))


def run_gmm(df_clean: pd.DataFrame, distance_lit_pc: float = DIST_LIT_PC, half_width_pc: float = HALF_WIDTH_PC,
            cutoff: float = CUTOFF, k: int | None = K_DEFAULT, k_range=K_RANGE, features=FEATURES,
            plx_zp: float = PLX_ZP, random_state: int = RANDOM_STATE, optimizer: str = OPTIMIZER) -> GMMResult:
    """Distance window -> standardise -> GMM (K by BIC unless k given) -> tightest component = cluster."""
    plx = COL['plx']
    d = df_clean.copy()
    plx_corr = d[plx] - plx_zp
    d = d.loc[plx_corr > 0].copy()
    d['distance_pc'] = 1000.0 / (d[plx] - plx_zp)
    wd = d.loc[(d['distance_pc'] - distance_lit_pc).abs() < half_width_pc].copy()
    X = wd[list(features)].dropna()
    k_lo, k_hi = (k, k) if k is not None else k_range
    if len(X) < MIN_STARS_PER_COMPONENT * k_hi:
        raise ValueError(f'Only {len(X)} stars in distance window; need >= {MIN_STARS_PER_COMPONENT * k_hi}. '
                         'Check distance / half-width / quality cuts.')
    Z = _standardise(X)

    tol, max_iter = (GMM_TOL, GMM_MAX_ITER) if optimizer == 'converged' else (1e-3, 100)
    models, bic, k_table = {}, {}, []
    for kk in range(k_lo, k_hi + 1):
        g = GaussianMixture(n_components=kk, covariance_type='full', n_init=N_INIT,
                            tol=tol, max_iter=max_iter, random_state=random_state).fit(Z)
        if g.weights_.min() * len(Z) < MIN_COMPONENT_STARS:      # degenerate: skip this K
            k_table.append(dict(k=kk, bic=None, members=None, non_members=None, min_comp_stars=int(round(g.weights_.min() * len(Z))), converged=bool(g.converged_), degenerate=True))
            continue
        models[kk], bic[kk] = g, float(g.bic(Z))
        pm_k = g.predict_proba(Z)[:, int(np.argmin([np.linalg.det(s) for s in g.covariances_]))]
        k_table.append(dict(k=kk, bic=bic[kk], members=int((pm_k >= cutoff).sum()), non_members=int((pm_k <= 1 - cutoff).sum()),
                            min_comp_stars=int(round(g.weights_.min() * len(Z))), converged=bool(g.converged_), degenerate=False))
    if not bic:
        raise ValueError(f'Every K in {k_lo}..{k_hi} gave a component with < {MIN_COMPONENT_STARS} stars '
                         f'({len(X)} stars in window). Widen half-width or relax quality cuts.')
    k_best = min(bic, key=bic.get)
    gmm = models[k_best]

    c_idx = int(np.argmin([np.linalg.det(s) for s in gmm.covariances_]))   # tightest blob = cluster
    wd['PMemb'] = np.nan
    wd.loc[X.index, 'PMemb'] = gmm.predict_proba(Z)[:, c_idx]
    wd['gmm_label'] = np.nan
    wd.loc[X.index, 'gmm_label'] = (wd.loc[X.index, 'PMemb'] >= cutoff).astype(int)

    members = wd.loc[wd['PMemb'] >= cutoff].copy()
    non_members = wd.loc[wd['PMemb'] <= (1 - cutoff)].copy()
    mss = _mss(members, non_members, list(features))

    summary = dict(
        n_input=len(df_clean), n_in_window=len(wd), n_fit=len(X), k=k_best, n_members=len(members),
        n_non_members=len(non_members), mss=mss,
        median_plx=float(members[plx].median()) if len(members) else float('nan'),
        median_pmRA=float(members[COL['pmra']].median()) if len(members) else float('nan'),
        median_pmDE=float(members[COL['pmde']].median()) if len(members) else float('nan'),
        distance_from_members_pc=float(1000.0 / (members[plx].median() - plx_zp)) if len(members) else float('nan'),
    )
    mu = gmm.means_ * X.std(ddof=0).values + X.median().values          # back to pmRA/pmDE/Plx units
    comps = [dict(k=i, weight=float(gmm.weights_[i]), n_hard=int((gmm.predict(Z) == i).sum()), det=float(np.linalg.det(gmm.covariances_[i])),
                  pmRA=float(mu[i][0]), pmDE=float(mu[i][1]), Plx=float(mu[i][2]),
                  sd_pmRA=float(np.sqrt(gmm.covariances_[i][0, 0]) * X.std(ddof=0).values[0]),
                  sd_pmDE=float(np.sqrt(gmm.covariances_[i][1, 1]) * X.std(ddof=0).values[1]),
                  sd_Plx=float(np.sqrt(gmm.covariances_[i][2, 2]) * X.std(ddof=0).values[2]), is_cluster=(i == c_idx))
             for i in range(k_best)]
    return GMMResult(wd, members, non_members, k_best, bic, c_idx, mss, summary, comps, k_table)


def make_gmm_figures(res: GMMResult, tag: str = '') -> dict:
    """Return {'pm','plx','cmd_member','cmd_nonmember'} -> matplotlib Figure (no pyplot state, no show())."""
    c = COL
    t = f'{tag} ' if tag else ''
    figs = {}

    f = Figure(figsize=(7, 6)); ax = f.subplots()
    ax.scatter(res.non_members[c['pmra']], res.non_members[c['pmde']], s=4, c='tab:blue', label='Non-member')
    ax.scatter(res.members[c['pmra']], res.members[c['pmde']], s=6, c='tab:red', label='Member')
    ax.set(xlabel='pmRA (mas/yr)', ylabel='pmDE (mas/yr)', title=f'{t}Proper motion'); ax.legend()
    figs['pm'] = f

    f = Figure(figsize=(7, 6)); ax = f.subplots()
    ax.hist(res.non_members[c['plx']].dropna(), bins=40, density=True, alpha=0.5, color='tab:blue', label='Non-member')
    ax.hist(res.members[c['plx']].dropna(), bins=40, density=True, alpha=0.5, color='tab:red', label='Member')
    ax.set(xlabel='Parallax (mas)', ylabel='Probability density', title=f'{t}Parallax'); ax.legend()
    figs['plx'] = f

    for key, df, col, name in (('cmd_member', res.members, 'tab:red', 'Member CMD'),
                               ('cmd_nonmember', res.non_members, 'tab:blue', 'Non-member CMD')):
        f = Figure(figsize=(7, 6)); ax = f.subplots()
        ax.scatter(df[c['color']], df[c['gmag']], s=6, c=col)
        ax.invert_yaxis(); ax.set(xlabel='BP - RP (mag)', ylabel='G (mag)', title=f'{t}{name}')
        figs[key] = f
    return figs


def save_members(res: GMMResult, tag: str, out_dir: str = OUT_DIR) -> str:
    """Write members CSV (index=False); returns path."""
    os.makedirs(out_dir, exist_ok=True)
    p = os.path.join(out_dir, f'{tag}_members.csv')
    res.members.to_csv(p, index=False)
    return p
