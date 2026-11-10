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
MF_KIND = 'Initial Mass Function (IMF)'   # label used in plot titles
N_CURVE_POINTS = 800
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


def compute_imf(primary_masses, companion_masses=None, fit_min_logm: float = FIT_MIN_LOGM,
                n_bins: int | None = N_BINS) -> IMFResult:
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
    lo_edge = float(edges[:-1][fit].min())
    lm = np.linspace(lo_edge, edges[-1], N_CURVE_POINTS)
    m_th = 10 ** lm
    # curves over the full plotted range, scaled on the complete range only
    lm_all = np.linspace(edges[0], edges[-1], N_CURVE_POINTS)
    curves, rows = {'logM': lm_all}, []
    for name, f in MODELS.items():
        scale = n_fit / trapezoid(f(m_th), lm)
        curves[name] = f(10 ** lm_all) * scale
        mod = f(10 ** ctr[fit]) * scale
        chi2 = float(np.sum(((dn[fit] - mod) / err[fit]) ** 2) / (fit.sum() - 1))
        rows.append(dict(Model=name, red_chi2=chi2, dof=int(fit.sum() - 1)))
    meta = dict(n_masses=len(masses), n_bins=nb, bin_width_dex=w, n_fit_bins=int(fit.sum()),
                n_in_fit=n_fit, fit_min_logm=fit_min_logm, m_range=(float(m_lo), float(m_hi)))
    return IMFResult(bins, pd.DataFrame(rows), curves, meta)


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
