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
COMPLETE_FRAC = 0.9          # completeness: retained fraction (quality-cut survivors / raw window stars) that still counts as "complete"
G_BINS = (6.0, 20.0, 0.5)    # (G_min, G_max, step) in mag for the completeness curve
G_BIN_MIN_RAW = 5            # a G bin with fewer raw stars than this is too noisy to judge -> treated as not complete
G_LIM_RULE = 'contiguous'    # 'contiguous': G_lim = faintest bin reached before the fraction first drops below COMPLETE_FRAC (robust to a noisy faint tail)
                             # 'faintest'  : G_lim = faintest bin anywhere with fraction >= COMPLETE_FRAC
COMPLETE_PREFILTER = True    # measure completeness on raw stars that are PM/parallax-compatible with the cluster (else on every window star)
G_FRACTION_MODE = 'cumulative'   # 'cumulative': f(<G) = N_clean(<G)/N_raw(<G) (smooth, works for ~30-star clusters) | 'differential': per-bin ratio
COMPLETE_NORM = 'plateau'    # 'plateau': threshold = COMPLETE_FRAC x bright-end plateau (loss that is independent of magnitude -- RUWE / binaries -- does
                             #            not move the limit; only the magnitude-dependent drop does) | 'absolute': threshold = COMPLETE_FRAC as is
PLATEAU_MIN_RAW = 10         # a bin enters the plateau estimate only once this many raw stars lie brighter than it
PLATEAU_N_BINS = 6           # plateau = median fraction over the brightest this-many eligible bins
COMPLETE_PM_CHI2 = 11.83     # chi2 (2 dof) limit on the proper-motion offset, ~3 sigma
COMPLETE_PLX_NSIG = 3.0      # parallax offset limit in sigma
CLUSTER_RULE = 'smallest_det'    # which Gaussian is "the cluster": 'smallest_det' (tightest blob) | 'closest_to_literature' (needs LIT_PM + LIT_PLX)
LIT_PM = None                    # literature cluster proper motion (pmRA, pmDE) mas/yr, for CLUSTER_RULE='closest_to_literature'
LIT_PLX = None                   # literature cluster parallax (mas)
N_RESAMPLE = 100                 # PMemb Monte-Carlo draws within the Gaia astrometric errors
RESAMPLE_SEED = 42
HDBSCAN_MIN_CLUSTER_SIZE = 10    # smallest group HDBSCAN may call a cluster
HDBSCAN_MIN_SAMPLES = 5          # core-distance neighbours; the sklearn default (= min_cluster_size) returned 'all noise' on a clean 150-star cluster
HR_VIZIER_CATALOG = 'J/A+A/673/A114'   # Hunt & Reffert (2023) cluster catalogue on VizieR (column names are looked up defensively, see hunt_reffert_prefill)
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
    model: object = None                             # fitted GaussianMixture of the chosen K (used by resample_pmemb)
    x_median: object = None                          # feature medians / stds used to standardise (pd.Series)
    x_std: object = None
    cutoff: float = CUTOFF
    features: list = field(default_factory=lambda: list(FEATURES))


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


def distance_window(df: pd.DataFrame, distance_lit_pc: float = DIST_LIT_PC, half_width_pc: float = HALF_WIDTH_PC,
                    plx_zp: float = PLX_ZP) -> pd.DataFrame:
    """Stars whose parallax distance lies within the literature window.
    TERMINOLOGY: distance_pc = 1000/(Plx - plx_zp) [pc, Plx in mas]; window = |distance_pc - distance_lit_pc| < half_width_pc.
    WORKFLOW: 1) drop stars with Plx - zp <= 0 (no finite distance); 2) compute distance_pc; 3) keep stars inside the window.
    Shared by run_gmm and completeness_curve so both count exactly the same stars."""
    plx = COL['plx']
    d = df.loc[(df[plx] - plx_zp) > 0].copy()
    d['distance_pc'] = 1000.0 / (d[plx] - plx_zp)
    return d.loc[(d['distance_pc'] - distance_lit_pc).abs() < half_width_pc].copy()


def cluster_compatible(raw_window: pd.DataFrame, members: pd.DataFrame, chi2_pm: float | None = None,
                       nsig_plx: float | None = None) -> pd.DataFrame:
    """Raw window stars whose astrometry is consistent with the cluster, judged WITH their own errors.
    TERMINOLOGY: mu_c, sigma_c = median and std of the member pmRA / pmDE / Plx (the intrinsic cluster spread); e = a star's Gaia error;
    chi2_pm = sum over (pmRA, pmDE) of (pm - mu_c)^2 / (sigma_c^2 + e_pm^2) (2 dof); dPlx = |Plx - mu_c| / sqrt(sigma_c^2 + e_Plx^2).
    Why: the plain distance window is ~80 % field stars, whose own cut failures (RUWE, gofAL) say nothing about how completely the CLUSTER
    survives. Inflating sigma by the star's own error keeps noisy cluster members in the denominator instead of silently dropping them.
    WORKFLOW: 1) member medians and stds; 2) per raw star compute chi2_pm and dPlx; 3) keep chi2_pm < limit and dPlx < nsig.
    Stars with a missing proper motion (2-parameter solutions) cannot be tested and are left out -> faint-end fraction is mildly optimistic."""
    chi2_pm = COMPLETE_PM_CHI2 if chi2_pm is None else chi2_pm
    nsig_plx = COMPLETE_PLX_NSIG if nsig_plx is None else nsig_plx
    c = COL
    mu = {k: float(members[c[k]].median()) for k in ('pmra', 'pmde', 'plx')}
    sd = {k: float(members[c[k]].std(ddof=0)) for k in ('pmra', 'pmde', 'plx')}
    def e(col):                                       # missing error -> 0 (star then has to match on the intrinsic spread alone)
        return raw_window[col].fillna(0.0) if col in raw_window else 0.0
    chi2 = ((raw_window[c['pmra']] - mu['pmra']) ** 2 / (sd['pmra'] ** 2 + e(c['e_pmra']) ** 2)
            + (raw_window[c['pmde']] - mu['pmde']) ** 2 / (sd['pmde'] ** 2 + e(c['e_pmde']) ** 2))
    dplx = (raw_window[c['plx']] - mu['plx']).abs() / np.sqrt(sd['plx'] ** 2 + e(c['e_plx']) ** 2)
    return raw_window.loc[(chi2 < chi2_pm) & (dplx < nsig_plx)].copy()


def completeness_curve(raw_window: pd.DataFrame, clean_window: pd.DataFrame, bins=G_BINS, frac: float | None = None,
                       rule: str | None = None, min_raw: int | None = None, mode: str | None = None, norm: str | None = None) -> dict:
    """Measured sample completeness vs apparent G (replaces the assumed logistic curve).
    TERMINOLOGY: raw_window = window stars BEFORE the Gaia quality cuts; clean_window = the same stars AFTER; f(G) = N_clean/N_raw
    ('cumulative' mode counts everything brighter than G, 'differential' mode counts one bin); plateau = bright-end value of f;
    threshold = frac (norm 'absolute') or frac x plateau (norm 'plateau'); G_lim = faintest G before f first falls below the threshold.
    Faint stars fail the RUWE / e_pm / Plx-SNR / Nper cuts first, so f falls with G: the cuts, not the detector, set the usable depth.
    NOTE: this measures completeness of the QUALITY-SELECTED sample relative to what Gaia lists; it cannot see stars Gaia itself missed.
    WORKFLOW: 1) histogram both samples in G bins; 2) cumulate if requested; 3) f = clean/raw (NaN where raw < min_raw); 4) plateau from
    the brightest eligible bins; 5) G_lim by `rule` ('contiguous' walks bright -> faint and stops at the first failing bin);
    6) return the curve table, G_lim, threshold and plateau."""
    frac, rule, min_raw = (COMPLETE_FRAC if frac is None else frac, G_LIM_RULE if rule is None else rule,
                           G_BIN_MIN_RAW if min_raw is None else min_raw)       # config read at call time so edits take effect
    mode, norm = G_FRACTION_MODE if mode is None else mode, COMPLETE_NORM if norm is None else norm
    g = COL['gmag']
    lo, hi, step = bins
    edges = np.arange(lo, hi + step / 2, step)
    n_raw, _ = np.histogram(raw_window[g].dropna(), bins=edges)
    n_cln, _ = np.histogram(clean_window[g].dropna(), bins=edges)
    if mode == 'cumulative':
        n_raw, n_cln = np.cumsum(n_raw), np.cumsum(n_cln)
        ctr = edges[1:]                                     # a cumulative bin is labelled by its faint edge ("brighter than G")
    else:
        ctr = edges[:-1] + step / 2
    has = n_raw >= min_raw
    with np.errstate(invalid='ignore', divide='ignore'):
        f = np.where(has, n_cln / np.where(n_raw > 0, n_raw, 1), np.nan)
    curve = pd.DataFrame({'Gmag': ctr, 'n_raw': n_raw, 'n_clean': n_cln, 'fraction': f})
    plateau = np.nan
    if norm == 'plateau':
        elig = np.where(n_raw >= PLATEAU_MIN_RAW)[0][:PLATEAU_N_BINS]
        plateau = float(np.nanmedian(f[elig])) if len(elig) else np.nan
    thr = frac * plateau if norm == 'plateau' and np.isfinite(plateau) else frac
    ok = np.isfinite(f) & (f >= thr)
    if rule == 'faintest':
        g_lim = float(ctr[ok].max()) if ok.any() else np.nan
    else:                                                   # contiguous from the bright end, ignoring bins with too few stars
        g_lim = np.nan
        for i in range(len(ctr)):
            if not has[i]:
                continue
            if not ok[i]:
                break
            g_lim = float(ctr[i])
    return dict(curve=curve, g_lim=g_lim, frac=frac, threshold=float(thr), plateau=plateau, rule=rule, mode=mode, norm=norm,
                n_raw=int(n_raw[-1]), n_clean=int(n_cln[-1]))


def _cluster_component(g: GaussianMixture, x_med: pd.Series, x_sd: pd.Series, rule: str, lit_pm, lit_plx, features) -> int:
    """Index of the Gaussian that is "the cluster".
    TERMINOLOGY: smallest_det = component whose covariance has the smallest determinant (the tightest blob in PM + parallax);
    closest_to_literature = component whose MEAN is nearest the literature (pmRA, pmDE, Plx), in standardised units (each axis divided by the
    sample std so a mas/yr and a mas are comparable).
    Why a second rule: with K > 2 a tight FIELD sub-structure can beat the real cluster on determinant; anchoring on a published
    PM + parallax removes that ambiguity (at the price of trusting the literature value).
    WORKFLOW: 1) smallest_det -> argmin det; 2) closest_to_literature -> standardise the literature point with the sample median / std,
    take the nearest component mean."""
    if rule == 'smallest_det':
        return int(np.argmin([np.linalg.det(s) for s in g.covariances_]))
    if rule != 'closest_to_literature':
        raise ValueError(f"cluster_rule must be 'smallest_det' or 'closest_to_literature', got {rule!r}")
    if lit_pm is None or lit_plx is None:
        raise ValueError("cluster_rule='closest_to_literature' needs lit_pm=(pmRA, pmDE) and lit_plx.")
    point = {COL['pmra']: lit_pm[0], COL['pmde']: lit_pm[1], COL['plx']: lit_plx}
    try:
        lit = np.array([(point[f] - x_med[f]) / x_sd[f] for f in features], float)
    except KeyError as e:
        raise ValueError(f'closest_to_literature needs features among pmRA, pmDE, Plx; missing {e}')
    return int(np.argmin(np.linalg.norm(g.means_ - lit, axis=1)))


def resample_pmemb(res: GMMResult, n_resample: int | None = None, seed: int | None = None) -> pd.DataFrame:
    """Per-star membership uncertainty from the Gaia astrometric errors.
    TERMINOLOGY: each draw perturbs (pmRA, pmDE, Plx) by N(0, e_pmRA), N(0, e_pmDE), N(0, e_Plx) (errors treated as independent: the
    VizieR table has no correlation columns, so this slightly UNDER-states the scatter); PMemb_mean / PMemb_std = mean and std of the
    cluster posterior over the draws; PMemb_frac_member = fraction of draws with PMemb >= cutoff.
    WORKFLOW: 1) standardise with the ORIGINAL median/std; 2) draw n_resample noisy copies; 3) evaluate the already-fitted mixture's
    posterior for the cluster component (no refit: this isolates the astrometric error from the fit variance); 4) summarise per star.
    A star with PMemb ~0.55 but frac_member ~0.5 is a coin flip, not a member."""
    n_resample = N_RESAMPLE if n_resample is None else n_resample
    rng = np.random.default_rng(RESAMPLE_SEED if seed is None else seed)
    if res.model is None:
        raise ValueError('GMMResult carries no fitted model (was it produced by an older run_gmm?).')
    c, feats = COL, res.features
    err_col = {c['pmra']: c['e_pmra'], c['pmde']: c['e_pmde'], c['plx']: c['e_plx']}
    w = res.working.dropna(subset=feats)
    X = w[feats].to_numpy(float)
    E = np.column_stack([w[err_col[f]].fillna(0.0).to_numpy(float) for f in feats])
    draws = X[None, :, :] + rng.normal(size=(n_resample,) + X.shape) * E[None, :, :]
    Z = (draws - res.x_median[feats].to_numpy(float)) / res.x_std[feats].to_numpy(float)
    post = res.model.predict_proba(pd.DataFrame(Z.reshape(-1, X.shape[1]), columns=feats))[:, res.cluster_component].reshape(n_resample, -1)
    return pd.DataFrame({'PMemb_mean': post.mean(0), 'PMemb_std': post.std(0, ddof=1),
                         'PMemb_frac_member': (post >= res.cutoff).mean(0)}, index=w.index)


def hdbscan_crosscheck(res: GMMResult, min_cluster_size: int | None = None, min_samples: int | None = None) -> dict:
    """Independent density-based membership to cross-check the GMM.
    TERMINOLOGY: HDBSCAN = hierarchical density clustering (no assumed Gaussian shape, no K); label -1 = noise; HDBSCAN_in_cluster = star lies
    in the HDBSCAN group that overlaps the GMM members most; jaccard = |A and B| / |A or B| between that group and the GMM members.
    WORKFLOW: 1) standardise the window features exactly as the GMM did; 2) run sklearn HDBSCAN; 3) choose the group sharing the most
    stars with the GMM members; 4) report the overlap. A low jaccard means the GMM and a shape-free method disagree about who the cluster is."""
    from sklearn.cluster import HDBSCAN
    mcs = HDBSCAN_MIN_CLUSTER_SIZE if min_cluster_size is None else min_cluster_size
    w = res.working.dropna(subset=res.features)
    Z = ((w[res.features] - res.x_median[res.features]) / res.x_std[res.features]).to_numpy(float)
    lab = HDBSCAN(min_cluster_size=mcs, min_samples=HDBSCAN_MIN_SAMPLES if min_samples is None else min_samples, copy=True).fit_predict(Z)
    is_gmm = (w['PMemb'] >= res.cutoff).to_numpy()
    groups = [g for g in np.unique(lab) if g >= 0]
    best, ov = -1, 0
    for g in groups:
        o = int(((lab == g) & is_gmm).sum())
        if o > ov:
            best, ov = g, o
    in_h = (lab == best) if best >= 0 else np.zeros(len(lab), bool)
    union = int((in_h | is_gmm).sum())
    out = pd.DataFrame({'HDBSCAN_label': lab, 'HDBSCAN_in_cluster': in_h}, index=w.index)
    stats = dict(n_groups=len(groups), n_hdbscan=int(in_h.sum()), n_gmm=int(is_gmm.sum()), n_both=int((in_h & is_gmm).sum()),
                 jaccard=(float((in_h & is_gmm).sum() / union) if union else float('nan')), min_cluster_size=mcs)
    return dict(table=out, stats=stats)


def hunt_reffert_prefill(name: str, catalog: str = HR_VIZIER_CATALOG, timeout_s: int = 30) -> dict:
    """Look a cluster up in the Hunt & Reffert (2023) catalogue to prefill literature PM / parallax / distance.
    TERMINOLOGY: prefill = values offered to the user as starting inputs (distance_lit, lit_pm, lit_plx), never applied silently.
    WORKFLOW: 1) astroquery Vizier query by cluster name; 2) read the first row, trying the likely column names; 3) return a dict with
    whichever of pmRA, pmDE, plx, distance_pc were found. Needs internet + astroquery; raises RuntimeError otherwise.
    NOTE: the catalogue's column names were NOT verified against a live query here; unknown columns are reported in 'missing'."""
    try:
        from astroquery.vizier import Vizier
    except ImportError as e:
        raise RuntimeError(f'astroquery is not installed: {e}')
    v = Vizier(columns=['**'], row_limit=5)
    v.TIMEOUT = timeout_s
    try:
        tabs = v.query_constraints(catalog=catalog, Name=name)
    except Exception as e:                                  # network / service errors differ widely between astroquery versions
        raise RuntimeError(f'Hunt & Reffert lookup failed: {e}')
    if not tabs or len(tabs[0]) == 0:
        raise RuntimeError(f"No row for '{name}' in {catalog}.")
    row = tabs[0][0]
    names = set(tabs[0].colnames)
    pick = lambda *cands: next((float(row[c]) for c in cands if c in names and np.isfinite(float(row[c]))), None)
    out = dict(pmRA=pick('pmRA', 'pmRA_', 'pmra'), pmDE=pick('pmDE', 'pmDE_', 'pmdec'), plx=pick('plx', 'Plx', 'parallax'),
               distance_pc=pick('Dist', 'dist50', 'Dist50', 'distance'))
    out['missing'] = [k for k, v_ in out.items() if v_ is None]
    return out


def run_gmm(df_clean: pd.DataFrame, distance_lit_pc: float = DIST_LIT_PC, half_width_pc: float = HALF_WIDTH_PC,
            cutoff: float = CUTOFF, k: int | None = K_DEFAULT, k_range=K_RANGE, features=FEATURES,
            plx_zp: float = PLX_ZP, random_state: int = RANDOM_STATE, optimizer: str = OPTIMIZER,
            cluster_rule: str = CLUSTER_RULE, lit_pm=LIT_PM, lit_plx=LIT_PLX) -> GMMResult:
    """Distance window -> standardise -> GMM (K by BIC unless k given) -> cluster component by `cluster_rule` (default: tightest, smallest det)."""
    plx = COL['plx']
    wd = distance_window(df_clean, distance_lit_pc, half_width_pc, plx_zp)
    X = wd[list(features)].dropna()
    k_lo, k_hi = (k, k) if k is not None else k_range
    if len(X) < MIN_STARS_PER_COMPONENT * k_hi:
        raise ValueError(f'Only {len(X)} stars in distance window; need >= {MIN_STARS_PER_COMPONENT * k_hi}. '
                         'Check distance / half-width / quality cuts.')
    Z = _standardise(X)
    x_med, x_sd = X.median(), X.std(ddof=0)
    pick = lambda g: _cluster_component(g, x_med, x_sd, cluster_rule, lit_pm, lit_plx, list(features))

    tol, max_iter = (GMM_TOL, GMM_MAX_ITER) if optimizer == 'converged' else (1e-3, 100)
    models, bic, k_table = {}, {}, []
    for kk in range(k_lo, k_hi + 1):
        g = GaussianMixture(n_components=kk, covariance_type='full', n_init=N_INIT,
                            tol=tol, max_iter=max_iter, random_state=random_state).fit(Z)
        if g.weights_.min() * len(Z) < MIN_COMPONENT_STARS:      # degenerate: skip this K
            k_table.append(dict(k=kk, bic=None, members=None, non_members=None, min_comp_stars=int(round(g.weights_.min() * len(Z))), converged=bool(g.converged_), degenerate=True))
            continue
        models[kk], bic[kk] = g, float(g.bic(Z))
        pm_k = g.predict_proba(Z)[:, pick(g)]
        k_table.append(dict(k=kk, bic=bic[kk], members=int((pm_k >= cutoff).sum()), non_members=int((pm_k <= 1 - cutoff).sum()),
                            min_comp_stars=int(round(g.weights_.min() * len(Z))), converged=bool(g.converged_), degenerate=False))
    if not bic:
        raise ValueError(f'Every K in {k_lo}..{k_hi} gave a component with < {MIN_COMPONENT_STARS} stars '
                         f'({len(X)} stars in window). Widen half-width or relax quality cuts.')
    k_best = min(bic, key=bic.get)
    gmm = models[k_best]

    c_idx = pick(gmm)                                                      # tightest blob (or the one nearest the literature) = cluster
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
    return GMMResult(wd, members, non_members, k_best, bic, c_idx, mss, summary, comps, k_table,
                     model=gmm, x_median=x_med, x_std=x_sd, cutoff=cutoff, features=list(features))


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
