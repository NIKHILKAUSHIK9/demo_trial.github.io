"""
imf.py -- Stage 3: mass-function binning + Salpeter / Kroupa / Chabrier comparison (no UI code).

TERMINOLOGY
  dN/dlogM       : number of stars per dex of mass. Salpeter slope in log-log = 1 - alpha = -1.35.
  Completeness   : below FIT_MIN_LOGM bins are incomplete -> shown (open markers) but excluded from fit & scaling.
  Scaling        : each model is scaled so its integral over the COMPLETE mass range equals the observed count there.
  Reduced chi2   : sum(((obs - model)/sigma)^2)/(n_bins - 1); sigma = sqrt(count)/bin_width. Compare models RELATIVELY.
  IMF label      : the result is presented as the Initial Mass Function. For clusters older than ~10 Myr dynamical and stellar evolution can
                   modify the observed mass function, so interpret the slope with that in mind.
CAUTION: do not tune N_BINS or FIT_MIN_LOGM to make a model fit -- that is circular.

WORKFLOW
  masses (primaries + companions) -> clean -> log-mass bins -> mark complete bins -> Poisson errors ->
  scale 3 models -> reduced chi2 table -> plot_imf
"""
from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from scipy.integrate import trapezoid

# ===================== GLOBAL CONFIGURATION =====================
FIT_MIN_LOGM = -0.7            # completeness limit, log10(M/Msun) (~0.2 Msun)
N_BINS = None                  # None -> min(15, max(8, ceil(log2(N)+1)))
M_MIN_FLOOR, M_MAX_CAP = 0.08, 50.0
MIN_MASSES = 20
MF_LABEL_INITIAL = 'Initial Mass Function (IMF)'
MF_LABEL_PRESENT = 'Present-Day Mass Function'
IMF_LOGAGE_SPLIT = 7.5         # log(age/yr) above which stellar evolution + dynamics have modified the MF -> "Present-Day"
MF_KIND = MF_LABEL_INITIAL     # default label used in plot titles (kept for the plot_imf API)
N_CURVE_POINTS = 800
METHOD = 'mle'                 # 'mle' (truncated power law, default) | 'binned' (weighted line fit in log-log)
ALPHA_BOUNDS = (0.5, 4.5)      # search interval for the MLE slope alpha (dN/dM ~ M^-alpha); Gamma = 1 - alpha
N_BOOT = 1000                  # bootstrap resamples for sigma(Gamma)
BOOT_SEED = 42
MIN_MASSES_MLE = 10            # smallest sample above the completeness limit for which the MLE is attempted
CASH_MIN_COUNT = 5             # any complete bin with n < this -> Cash statistic instead of Pearson chi2
N_SUBGRID = 40                 # trapezoid points per bin when integrating a model over a bin
INCLUDE_APPROX_ML = False      # True -> keep masses with MC_Mass_src == 'approx_ML' (L(M) fallback, not isochrone-based)
OUT_DIR = './starnest_output'
# =================================================================


def phi_salpeter(m, alpha=2.35):
    """Salpeter (1955): dN/dlogM ~ M^(1-alpha)."""
    return np.asarray(m, float) ** (1.0 - alpha)


def phi_kroupa(m):
    """Kroupa (2001) broken power law in dN/dlogM form (exponents +1 vs dN/dM); continuous at 0.08 and 0.5."""
    m = np.atleast_1d(np.asarray(m, float))
    a1, a2, a3 = 0.3, 1.3, 2.3
    k2 = 0.08 ** (a2 - a1)
    k3 = k2 * 0.50 ** (a3 - a2)
    return np.where(m < 0.08, m ** (1 - a1), np.where(m < 0.5, k2 * m ** (1 - a2), k3 * m ** (1 - a3)))


def phi_chabrier(m, alpha=2.3):
    """Chabrier (2003): log-normal (mc=0.08, sigma=0.69 dex) below 1 Msun, power law above, continuous at 1."""
    m = np.atleast_1d(np.asarray(m, float))
    mc, sigma = 0.08, 0.69
    low = 0.158 * np.exp(-0.5 * (np.log10(np.minimum(m, 1.0)) - np.log10(mc)) ** 2 / sigma ** 2)
    a_high = 0.158 * np.exp(-0.5 * (0.0 - np.log10(mc)) ** 2 / sigma ** 2)
    return np.where(m <= 1.0, low, a_high * m ** (1 - alpha))


MODELS = {'Salpeter': phi_salpeter, 'Kroupa': phi_kroupa, 'Chabrier': phi_chabrier}


@dataclass
class IMFResult:
    bins: pd.DataFrame            # per-bin table (all bins; 'complete' flags fit bins)
    chi2: pd.DataFrame            # Model, red_chi2, dof
    curves: dict                  # {'logM': array, 'Salpeter': array, ...} scaled dN/dlogM
    meta: dict


def build_mass_sample(primary_masses, companion_masses=None) -> np.ndarray:
    """Concatenate primaries + companions (>0), drop non-finite / non-positive."""
    p = np.asarray(primary_masses, float)
    if companion_masses is not None:
        c = np.asarray(companion_masses, float)
        p = np.concatenate([p, c[np.isfinite(c) & (c > 0)]])
    return p[np.isfinite(p) & (p > 0)]


def mf_label(log_age: float | None) -> str:
    """Page / plot label. TERMINOLOGY: log_age = log10(age/yr). WORKFLOW: 1) no age -> IMF; 2) age above IMF_LOGAGE_SPLIT
    -> 'Present-Day Mass Function' (high-mass stars have evolved off, dynamics have evaporated low-mass stars); 3) else IMF."""
    try:
        age = float(log_age)
    except (TypeError, ValueError):
        return MF_LABEL_INITIAL
    return MF_LABEL_PRESENT if age > IMF_LOGAGE_SPLIT else MF_LABEL_INITIAL


def g_limit_to_mass(g_lim: float, dm: float, a_g: float, ms_iso: pd.DataFrame) -> dict:
    """Convert an observed completeness limit G_lim into a mass limit on the A_V=0 isochrone.
    TERMINOLOGY: g_lim = faintest complete apparent G (mag); dm = distance modulus; a_g = extinction in G; ms_iso = intrinsic PARSEC table
    with Gmag (absolute), Mass and optional label (0 PMS, 1 MS). M_G = g_lim - dm - a_g is the absolute magnitude of the faintest complete star.
    WORKFLOW: 1) keep label 0/1 rows; 2) sort by Gmag; 3) m_min = interp(M_G, Gmag, Mass); 4) flag if M_G lies outside the tabulated range
    (np.interp then clamps to the nearest isochrone end, so the limit is only a bound)."""
    ms = ms_iso[ms_iso['label'].isin([0, 1])] if 'label' in ms_iso else ms_iso
    ms = ms.sort_values('Gmag')
    if len(ms) < 4:
        raise ValueError('Isochrone has too few main-sequence rows to convert G_lim to a mass.')
    xs, ys = ms['Gmag'].to_numpy(float), ms['Mass'].to_numpy(float)
    m_g = float(g_lim - dm - a_g)
    m = float(np.interp(m_g, xs, ys))
    return dict(m_min=m, logm=float(np.log10(m)), M_G_lim=m_g, clipped=bool(m_g < xs.min() or m_g > xs.max()))


def wls_slope(x, y, sig):
    """Poisson-weighted line fit y = icpt + slope*x. Returns (slope, icpt, r2, se); se is inflated by sqrt(chi2_nu) when chi2_nu > 1."""
    x, y, sig = map(np.asarray, (x, y, sig))
    w = np.where(sig > 0, 1 / np.where(sig > 0, sig, 1) ** 2, 0.0)
    n = len(x)
    mx, my = np.sum(w * x) / w.sum(), np.sum(w * y) / w.sum()
    sxx, sxy, syy = np.sum(w * (x - mx) ** 2), np.sum(w * (x - mx) * (y - my)), np.sum(w * (y - my) ** 2)
    slope = sxy / sxx if sxx > 0 else 0.0
    icpt = my - slope * mx
    chi2r = np.sum(w * (y - icpt - slope * x) ** 2) / (n - 2) if n > 2 else 1.0
    se = float(np.sqrt(max(1.0, chi2r) / sxx)) if sxx > 0 else 0.0
    r2 = sxy ** 2 / (sxx * syy) if sxx > 0 and syy > 0 else 0.0
    return float(slope), float(icpt), float(r2), se


def _log_norm(alpha, lo, hi):
    """ln of I(alpha) = integral_lo^hi m^-alpha dm = (hi^b - lo^b)/b, b = 1-alpha (b -> 0 limit: ln(hi/lo))."""
    b = 1.0 - np.asarray(alpha, float)
    big_l = np.log(hi / lo)
    small = np.abs(b) < 1e-8
    bs = np.where(small, 1.0, b)
    return np.where(small, np.log(big_l), b * np.log(lo) + np.log(np.expm1(bs * big_l) / bs))


def _mle_alpha(n: int, sum_ln_m: float, lo: float, hi: float, bounds=ALPHA_BOUNDS) -> float:
    """Truncated power law dN/dM = A m^-alpha on [lo, hi]. NLL(alpha) = n*ln I(alpha) + alpha*sum(ln m); depends on the data only through
    (n, sum ln m), which makes the bootstrap cheap."""
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda a: n * float(_log_norm(a, lo, hi)) + a * sum_ln_m, bounds=bounds, method='bounded',
                        options=dict(xatol=1e-8))
    return float(r.x)


def mle_slope(masses, m_min: float, m_max: float | None = None, bounds=ALPHA_BOUNDS, n_boot: int = N_BOOT,
              seed: int = BOOT_SEED) -> dict:
    """Maximum-likelihood slope of the mass function above the completeness limit.
    TERMINOLOGY: alpha = index of dN/dM ~ M^-alpha (Salpeter 2.35); Gamma = 1 - alpha = slope of log dN/dlogM (Salpeter -1.35);
    [m_min, m_max] = truncation interval (m_min = completeness limit, m_max = most massive star unless given); sigma = bootstrap std.
    WORKFLOW: 1) keep masses in [m_min, m_max]; 2) minimise the truncated-power-law NLL over alpha in `bounds`; 3) Gamma = 1 - alpha;
    4) resample the masses with replacement n_boot times (seeded), refit, sigma = std(Gamma); 5) intercept so that the fitted
    dN/dlogM = ln10 * A * M^Gamma integrates to the observed count (used to draw the line)."""
    m = np.asarray(masses, float)
    m = m[np.isfinite(m) & (m > 0)]
    hi = float(m.max() if m_max is None else m_max)
    m = m[(m >= m_min) & (m <= hi)]
    n = len(m)
    if n < MIN_MASSES_MLE or hi <= m_min:
        raise ValueError(f'MLE needs >= {MIN_MASSES_MLE} masses in [{m_min:.3f}, {hi:.3f}] Msun; found {n}.')
    ln_m = np.log(m)
    alpha = _mle_alpha(n, float(ln_m.sum()), m_min, hi, bounds)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot = np.array([_mle_alpha(n, float(ln_m[i].sum()), m_min, hi, bounds) for i in idx])
    gam_b = 1.0 - boot
    a_norm = n / np.exp(_log_norm(alpha, m_min, hi))                      # A: integral of A m^-alpha over [lo, hi] = n
    return dict(method='mle', alpha=alpha, gamma=1.0 - alpha, gamma_err=float(np.std(gam_b, ddof=1)),
                intercept=float(np.log10(np.log(10.0) * a_norm)), m_min=float(m_min), m_max=hi, n=int(n), n_boot=int(n_boot),
                at_bound=bool(min(abs(alpha - bounds[0]), abs(alpha - bounds[1])) < 1e-4))


def _bin_integral(f, edges, n_sub: int = N_SUBGRID):
    """Integral of a model dN/dlogM shape over each log-mass bin (trapezoid on n_sub points per bin)."""
    x = np.linspace(edges[:-1], edges[1:], n_sub, axis=1)
    return trapezoid(f(10 ** x.ravel()).reshape(x.shape), x, axis=1)


def compute_imf(primary_masses, companion_masses=None, fit_min_logm: float = FIT_MIN_LOGM,
                n_bins: int | None = N_BINS, method: str = METHOD, m_max: float | None = None,
                n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> IMFResult:
    """Mass function above the completeness limit and its comparison with Salpeter / Kroupa / Chabrier.
    TERMINOLOGY: method 'mle' = truncated power-law maximum likelihood (Gamma +- bootstrap sigma); 'binned' = Poisson-weighted line fit to
    log dN/dlogM; expected bin counts mu_b = scale * integral of the model over bin b; Cash C = 2 sum(mu - n + n ln(n/mu)).
    WORKFLOW: 1) bin log M; 2) mark bins complete above fit_min_logm; 3) model counts per bin by integration, scaled to the complete-range
    count; 4) Pearson chi2 (sigma = sqrt(n)) if every complete bin has n >= CASH_MIN_COUNT, else Cash; 5) slope by `method`
    (the binned slope is always reported in meta['binned'])."""
    if method not in ('mle', 'binned'):
        raise ValueError(f"method must be 'mle' or 'binned', got {method!r}")
    masses = build_mass_sample(primary_masses, companion_masses)
    if len(masses) < MIN_MASSES:
        raise ValueError(f'Only {len(masses)} valid masses; need >= {MIN_MASSES}.')
    nb = n_bins or min(15, max(8, int(np.ceil(np.log2(len(masses)) + 1))))
    m_lo, m_hi = max(M_MIN_FLOOR, masses.min() * 0.9), min(M_MAX_CAP, masses.max() * 1.1)
    if m_lo >= m_hi:
        raise ValueError('Degenerate mass range.')
    edges = np.linspace(np.log10(m_lo), np.log10(m_hi), nb + 1)
    w = float(np.diff(edges)[0])
    ctr = edges[:-1] + w / 2
    cnt, _ = np.histogram(np.log10(masses), bins=edges)
    valid = cnt > 0
    complete = ctr >= fit_min_logm
    fit = valid & complete
    if fit.sum() < 3:
        raise ValueError(f'Only {fit.sum()} complete bins above logM={fit_min_logm}; lower the limit or add data.')

    dn = np.where(valid, cnt / w, np.nan)
    err = np.where(valid, np.sqrt(cnt) / w, np.nan)
    bins = pd.DataFrame({'logM_center': ctr, 'count': cnt, 'dN_dlogM': dn, 'dN_dlogM_err': err,
                         'log_dN': np.log10(dn), 'log_err': err / dn / np.log(10), 'complete': complete})

    n_fit = int(cnt[fit].sum())
    lm_all = np.linspace(edges[0], edges[-1], N_CURVE_POINTS)
    use_cash = bool((cnt[complete] < CASH_MIN_COUNT).any())
    n_c = cnt[complete].astype(float)
    curves, rows = {'logM': lm_all}, []
    for name, f in MODELS.items():
        integ = _bin_integral(f, edges)                          # model counts are INTEGRATED over each bin, not sampled at its centre
        scale = n_fit / integ[complete].sum()                    # fix the total on the complete range; shape is what is tested
        curves[name] = f(10 ** lm_all) * scale
        mu = integ[complete] * scale
        if use_cash:                                             # Poisson likelihood statistic, valid for small / empty bins
            nz = n_c > 0
            stat = 2.0 * float(np.sum(mu - n_c) + np.sum(n_c[nz] * np.log(n_c[nz] / mu[nz])))
        else:
            stat = float(np.sum((n_c - mu) ** 2 / n_c))
        dof = int(complete.sum() - 1)
        rows.append(dict(Model=name, red_chi2=stat / dof, dof=dof, stat='Cash' if use_cash else 'chi2'))

    cb = bins[valid & complete]
    b_slope, b_icpt, b_r2, b_se = wls_slope(cb['logM_center'], cb['log_dN'], cb['log_err'])
    binned = dict(gamma=b_slope, gamma_err=b_se, intercept=b_icpt, r2=b_r2)
    if method == 'mle':
        mle = mle_slope(masses, 10 ** fit_min_logm, m_max, n_boot=n_boot, seed=seed)
        main_fit = dict(gamma=mle['gamma'], gamma_err=mle['gamma_err'], intercept=mle['intercept'])
    else:
        mle, main_fit = None, dict(gamma=b_slope, gamma_err=b_se, intercept=b_icpt)
    meta = dict(n_masses=len(masses), n_bins=nb, bin_width_dex=w, n_fit_bins=int(fit.sum()),
                n_in_fit=n_fit, fit_min_logm=fit_min_logm, m_range=(float(m_lo), float(m_hi)),
                method=method, stat='Cash' if use_cash else 'chi2', binned=binned, mle=mle, **{f'fit_{k}': v for k, v in main_fit.items()})
    return IMFResult(bins, pd.DataFrame(rows), curves, meta)


def mass_sample(df: pd.DataFrame, include_approx_ml: bool = INCLUDE_APPROX_ML, mass_col: str = 'MC_Mass'):
    """Masses usable for the mass function. TERMINOLOGY: MC_Mass_src 'isochrone' = interpolated on the PARSEC main sequence,
    'user' = supplied by the user, 'approx_ML' = crude L(M) fallback (no isochrone) -> excluded unless include_approx_ml.
    WORKFLOW: 1) read MC_Mass; 2) drop approx_ML rows unless allowed; 3) return (masses, n_dropped_approx)."""
    if mass_col not in df:
        return np.array([]), 0
    m = df[mass_col].to_numpy(float)
    approx = (df['MC_Mass_src'] == 'approx_ML').to_numpy() if 'MC_Mass_src' in df else np.zeros(len(df), bool)
    if include_approx_ml:
        return m, 0
    return m[~approx], int((approx & np.isfinite(m)).sum())


def plot_imf(res: IMFResult, cluster_name: str = '', mf_kind: str = MF_KIND) -> Figure:
    f = Figure(figsize=(7.5, 6)); ax = f.subplots()
    b = res.bins[res.bins['count'] > 0]
    c, ic = b[b['complete']], b[~b['complete']]
    ax.errorbar(c['logM_center'], c['log_dN'], yerr=c['log_err'], fmt='o', ms=6, capsize=4, color='black',
                label=f"{cluster_name} (N={res.meta['n_masses']}, {res.meta['n_bins']} bins)", zorder=5)
    if len(ic):
        ax.errorbar(ic['logM_center'], ic['log_dN'], yerr=ic['log_err'], fmt='o', ms=6, capsize=4, mfc='white',
                    color='grey', label='incomplete (not fitted)', zorder=4)
    ax.axvline(res.meta['fit_min_logm'], color='grey', ls=':', lw=1)
    style = {'Salpeter': ('steelblue', '--'), 'Kroupa': ('darkorange', '-.'), 'Chabrier': ('firebrick', '-')}
    chi = res.chi2.set_index('Model')['red_chi2']
    for name, (col, ls) in style.items():
        ax.plot(res.curves['logM'], np.log10(res.curves[name]), color=col, ls=ls, lw=2,
                label=f'{name} ($\\chi^2_\\nu$={chi[name]:.1f})')
    ax.set(xlabel=r'$\log_{10}(M/M_\odot)$', ylabel=r'$\log_{10}(dN/d\log M)$',
           title=f'{mf_kind} - model comparison\n{cluster_name}')
    ax.tick_params(direction='in', top=True, right=True, which='both')
    ax.legend(frameon=False, fontsize=9, loc='upper right')
    return f


def save_imf(res: IMFResult, tag: str, out_dir: str = OUT_DIR) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    p1 = os.path.join(out_dir, f'{tag}_IMF_bins.csv'); res.bins.to_csv(p1, index=False)
    p2 = os.path.join(out_dir, f'{tag}_IMF_chi2.csv'); res.chi2.to_csv(p2, index=False)
    return {'bins': p1, 'chi2': p2}
