"""
binaries.py -- Stage 5: Monte-Carlo masses with unresolved binaries + total cluster mass (no UI code).

TERMINOLOGY
  Synthetic cluster : N_SYNTH artificial systems drawn from a Chabrier IMF; a fraction F_BIN get a companion with mass U(0, M_primary).
  q                 : mass ratio M_companion / M_primary (0 for singles).
  System photometry : single -> isochrone (G, BP, RP) of M1; binary -> fluxes of M1 and M2 added, then converted back to magnitudes.
                      A binary is brighter and (slightly) different in colour than its primary, so a naive CMD mass overestimates M1.
  Nearest match     : for each observed star, the synthetic system closest in (G, BP, RP) (Euclidean, mag). Its M1 / M2 are assigned.
  N_MC              : the whole procedure is repeated with fresh synthetic clusters; MC_Mass = median M1, error = std.
  is_binary         : True when more than BINARY_VOTE_FRAC of the matches had a companion.
  M_T               : M_obs + M_unseen (1 + F_bin q_mean) + M_WD   (see cluster_mass)
  IFMR              : initial-final mass relation (Cummings et al. 2018): white-dwarf mass from the progenitor mass.

WORKFLOW
  A_V=0 isochrone + (DM, A_G, A_BP, A_RP) -> synthetic systems in the OBSERVED plane -> nearest-neighbour match per star (KD-tree)
  -> repeat N_MC times -> median / std per star -> assign_masses() table (MC_Mass, MC_CompMass, is_binary, q, errors)
  -> cluster_mass(): observed + unseen (fitted MF extrapolated to M_UNSEEN_FLOOR, with companions) + white dwarfs.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from . import imf as IMF

# ===================== GLOBAL CONFIGURATION =====================
N_SYNTH = 10000                 # synthetic systems per Monte-Carlo realisation
F_BIN = 0.5                     # assumed binary fraction of the synthetic population
N_MC = 50                       # Monte-Carlo realisations (fresh synthetic cluster each time)
BIN_SEED = 42
M_SYNTH_RANGE = None            # (lo, hi) Msun for primaries; None -> the isochrone's own mass range
BINARY_VOTE_FRAC = 0.5          # is_binary if more than this fraction of matches carried a companion
MATCH_MAX_DIST = 0.5            # mag; a star whose best match is farther than this (median over N_MC) gets NaN masses
BAND_COLUMNS = ('Gmag', 'BPmag', 'RPmag')                 # observed photometry
ISO_COLUMNS = ('Gmag', 'G_BPmag', 'G_RPmag')              # PARSEC absolute magnitudes (CMD 3.8, Gaia EDR3 Vega)
M_UNSEEN_FLOOR = 0.09           # Msun: hydrogen-burning limit; the fitted MF is integrated down to here
WD_PROGENITOR_MAX = 8.2         # Msun: above this the remnant is not a WD (Cummings et al. 2018 range)
# Cummings et al. (2018, ApJ 866, 21) piecewise-linear IFMR: (M_i_lo, M_i_hi, slope, intercept) -> M_f = slope*M_i + intercept
IFMR_PIECES = ((0.87, 2.8, 0.0873, 0.476), (2.8, 3.65, 0.187, 0.184), (3.65, 8.2, 0.107, 0.471))
N_INTEGRATION = 2000            # points for the numerical mass integrals
# =================================================================


def _chabrier_inverse_sampler(lo: float, hi: float, n_grid: int = 4000):
    """Inverse-CDF table of the Chabrier (2003) IMF in log mass: returns f(u) -> mass for u~U(0,1).
    Why log-mass: imf.phi_chabrier is dN/dlogM, so the CDF must be taken in log10 M, not in M."""
    lm = np.linspace(np.log10(lo), np.log10(hi), n_grid)
    cdf = np.cumsum(IMF.phi_chabrier(10 ** lm))
    cdf = (cdf - cdf[0]) / (cdf[-1] - cdf[0])
    return lambda u: 10 ** np.interp(u, cdf, lm)


def _iso_tables(iso: pd.DataFrame):
    """Sorted (mass, G, BP, RP) arrays of the intrinsic isochrone for mass -> absolute magnitude interpolation."""
    need = ['Mass', *ISO_COLUMNS]
    miss = [c for c in need if c not in iso]
    if miss:
        raise ValueError(f'Isochrone is missing columns: {miss}')
    t = iso[need].dropna().sort_values('Mass')
    if len(t) < 4:
        raise ValueError('Isochrone has too few rows to build a synthetic cluster.')
    return t['Mass'].to_numpy(float), [t[c].to_numpy(float) for c in ISO_COLUMNS]


def build_synthetic(iso: pd.DataFrame, dm: float, a_g: float, a_bp: float, a_rp: float, n: int = N_SYNTH, f_bin: float = F_BIN,
                    rng: np.random.Generator | None = None, m_range=M_SYNTH_RANGE) -> pd.DataFrame:
    """One synthetic cluster in the OBSERVED (G, BP, RP) plane.
    TERMINOLOGY: dm = distance modulus; a_g / a_bp / a_rp = extinction per band (mag); m1, m2 = primary / companion mass; q = m2/m1.
    WORKFLOW: 1) draw m1 from the Chabrier IMF over the isochrone mass range; 2) a fraction f_bin gets m2 ~ U(0, m1);
    3) interpolate absolute magnitudes of m1 and m2 (companions below the isochrone's lowest mass give NO light -- they are too faint to
    matter and extrapolating would invent flux); 4) add fluxes, back to magnitudes; 5) add DM and per-band extinction."""
    rng = rng or np.random.default_rng(BIN_SEED)
    mass, mags = _iso_tables(iso)
    lo, hi = m_range if m_range else (float(mass.min()), float(mass.max()))
    m1 = _chabrier_inverse_sampler(lo, hi)(rng.random(n))
    is_bin = rng.random(n) < f_bin
    m2 = np.where(is_bin, rng.random(n) * m1, 0.0)
    shifts = (dm + a_g, dm + a_bp, dm + a_rp)
    out = {'m1': m1, 'm2': m2, 'q': np.where(m1 > 0, m2 / m1, 0.0)}
    for name, col, sh in zip(BAND_COLUMNS, mags, shifts):
        f1 = 10 ** (-0.4 * np.interp(m1, mass, col))
        f2 = np.where(m2 >= mass.min(), 10 ** (-0.4 * np.interp(np.maximum(m2, mass.min()), mass, col)), 0.0)
        out[name] = -2.5 * np.log10(f1 + f2) + sh                # flux addition: two stars unresolved = one brighter point
    return pd.DataFrame(out)


def assign_masses(obs: pd.DataFrame, iso: pd.DataFrame, dm: float, a_g: float, a_bp: float, a_rp: float, n_synth: int = N_SYNTH,
                  f_bin: float = F_BIN, n_mc: int = N_MC, seed: int = BIN_SEED) -> pd.DataFrame:
    """Per-star primary mass, companion mass, binary flag and mass ratio.
    TERMINOLOGY: see module header. Output columns: MC_Mass, MC_Mass_err, MC_CompMass, MC_CompMass_err, is_binary, q, binary_prob
    (fraction of the n_mc matches that were binaries), match_dist (median matching distance, mag).
    WORKFLOW: 1) for each of n_mc realisations build a synthetic cluster and a KD-tree on (G, BP, RP); 2) query every observed star;
    3) store M1, M2, q, distance; 4) median / std over realisations; 5) is_binary by majority vote; 6) NaN where match_dist > MATCH_MAX_DIST."""
    miss = [c for c in BAND_COLUMNS if c not in obs]
    if miss:
        raise ValueError(f'Members file missing photometry columns: {miss}')
    X = obs[list(BAND_COLUMNS)].to_numpy(float)
    ok = np.isfinite(X).all(axis=1)
    rng = np.random.default_rng(seed)
    m1s, m2s, qs, ds = (np.full((n_mc, len(obs)), np.nan) for _ in range(4))
    for k in range(n_mc):
        syn = build_synthetic(iso, dm, a_g, a_bp, a_rp, n_synth, f_bin, rng)
        d, j = cKDTree(syn[list(BAND_COLUMNS)].to_numpy(float)).query(X[ok])
        m1s[k, ok], m2s[k, ok], qs[k, ok], ds[k, ok] = syn['m1'].to_numpy()[j], syn['m2'].to_numpy()[j], syn['q'].to_numpy()[j], d
    with np.errstate(all='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)          # all-NaN columns (unmatched stars) are expected and handled below
        prob = np.nanmean(m2s > 0, axis=0)
        q_bin = np.where(m2s > 0, qs, np.nan)
        out = pd.DataFrame({
            'MC_Mass': np.nanmedian(m1s, 0), 'MC_Mass_err': np.nanstd(m1s, 0, ddof=1),
            'MC_CompMass': np.nanmedian(m2s, 0), 'MC_CompMass_err': np.nanstd(m2s, 0, ddof=1),
            'is_binary': prob > BINARY_VOTE_FRAC, 'q': np.nanmedian(q_bin, 0), 'binary_prob': prob,
            'match_dist': np.nanmedian(ds, 0)}, index=obs.index)
    out['q'] = out['q'].fillna(0.0)
    far = ~(out['match_dist'] <= MATCH_MAX_DIST)
    out.loc[far, ['MC_Mass', 'MC_Mass_err', 'MC_CompMass', 'MC_CompMass_err']] = np.nan   # no synthetic system resembles this star
    return out


def ifmr_cummings2018(m_init):
    """White-dwarf final mass from the progenitor mass (Cummings et al. 2018). Zero outside the tabulated range (no WD remnant)."""
    m = np.atleast_1d(np.asarray(m_init, float))
    out = np.zeros_like(m)
    for lo, hi, a, b in IFMR_PIECES:
        sel = (m >= lo) & (m < hi)
        out[sel] = a * m[sel] + b
    return out


def _powerlaw_norm(alpha: float, n: int, lo: float, hi: float) -> float:
    """A in dN/dM = A m^-alpha so that the integral over [lo, hi] equals n (same normalisation as imf.mle_slope)."""
    return float(n / np.exp(IMF._log_norm(alpha, lo, hi)))


def unseen_mass(alpha: float, n_fit: int, m_min: float, m_max: float, floor: float = M_UNSEEN_FLOOR) -> float:
    """Mass of single stars below the completeness limit, from the fitted power law extrapolated down to `floor`.
    TERMINOLOGY: alpha = fitted index of dN/dM; n_fit = stars in [m_min, m_max] used by the fit; floor = lower integration bound.
    WORKFLOW: 1) A from n_fit; 2) M_unseen = integral_{floor}^{m_min} m A m^-alpha dm.
    CAVEAT: real mass functions flatten below ~0.5 Msun (Kroupa / Chabrier), so extrapolating a steeper fitted slope OVERestimates this."""
    if m_min <= floor:
        return 0.0
    a = _powerlaw_norm(alpha, n_fit, m_min, m_max)
    m = np.geomspace(floor, m_min, N_INTEGRATION)
    return float(np.trapezoid(m * a * m ** (-alpha), m))


def wd_mass(alpha: float, n_fit: int, m_min: float, m_max: float, m_turnoff: float, m_prog_max: float = WD_PROGENITOR_MAX) -> float:
    """Mass locked in white dwarfs: stars born between the turn-off and m_prog_max have died; each leaves IFMR(m_i).
    WORKFLOW: 1) A from the fit; 2) integrate A m^-alpha over [m_turnoff, m_prog_max]; 3) weight by IFMR. Needs m_turnoff < m_prog_max."""
    if m_turnoff >= m_prog_max:
        return 0.0
    a = _powerlaw_norm(alpha, n_fit, m_min, m_max)
    m = np.linspace(max(m_turnoff, 1e-3), m_prog_max, N_INTEGRATION)
    return float(np.trapezoid(a * m ** (-alpha) * ifmr_cummings2018(m), m))


def cluster_mass(table: pd.DataFrame, mle: dict, m_turnoff: float, f_bin: float = F_BIN, floor: float = M_UNSEEN_FLOOR) -> dict:
    """Total cluster mass M_T = M_obs + M_unseen (1 + f_bin q_mean) + M_WD.
    TERMINOLOGY: table = assign_masses output (observed members); mle = imf.mle_slope dict (alpha, m_min, m_max, n); m_turnoff = most massive
    star still on the main sequence (Msun); q_mean = mean mass ratio of the observed binaries; M_obs = sum of primary AND companion masses
    of the observed systems (companions are real mass even though unresolved).
    WORKFLOW: 1) M_obs; 2) M_unseen from the fitted MF; 3) companions of unseen stars: x(1 + f_bin q_mean); 4) M_WD; 5) sum.
    The unseen and WD terms are model extrapolations: treat M_T as an estimate with a large systematic error."""
    m1, m2 = table['MC_Mass'].fillna(0.0), table['MC_CompMass'].fillna(0.0)
    m_obs = float((m1 + m2).sum())
    q_mean = float(table.loc[table['is_binary'], 'q'].mean()) if table['is_binary'].any() else 0.0
    m_un = unseen_mass(mle['alpha'], mle['n'], mle['m_min'], mle['m_max'], floor)
    m_wd = wd_mass(mle['alpha'], mle['n'], mle['m_min'], mle['m_max'], m_turnoff)
    m_t = m_obs + m_un * (1.0 + f_bin * q_mean) + m_wd
    return dict(M_T=float(m_t), M_obs=m_obs, M_unseen=m_un, M_WD=m_wd, q_mean=q_mean, f_bin=f_bin,
                n_binary=int(table['is_binary'].sum()), n_stars=int(len(table)),
                observed_binary_fraction=float(table['is_binary'].mean()) if len(table) else float('nan'))


def save_binaries(table: pd.DataFrame, tag: str, out_dir: str = './starnest_output') -> str:
    """Write the per-star table (index=False)."""
    import os
    os.makedirs(out_dir, exist_ok=True)
    f = os.path.join(out_dir, f'{tag}_binaries.csv')
    table.to_csv(f, index=False)
    return f
