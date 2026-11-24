"""
kinematics.py -- Stage 4: 3-D cluster map + proper-motion (residual) vectors (no UI code).

TERMINOLOGY
  Plx / e_Plx       : Gaia parallax and its error (mas). Distance d = 1000 / Plx (pc).
  sigma_d           : per-star distance error d*e_Plx/Plx. It is purely radial (sky positions are exact to ~mas), so a round
                      cluster is smeared into a "needle" along the line of sight.
  precision mode    : keep only stars with e_Plx/Plx <= KIN3D_MAX_FRAC (drops the stars that smear the map).
  shrink mode       : empirical-Bayes parallax: each parallax is pulled toward the weighted cluster parallax plx_c; the prior
                      width is tau = plx_c^2 * KIN3D_DEPTH_PC / 1000 (mas), weights are 1/e_Plx^2 (data) and 1/tau^2 (prior).
  raw mode          : d = 1000/Plx for every star (no filter, no shrinkage).
  Frame             : 'galactic' rotates ICRS (x, y, z) into Galactic axes (ICRS_TO_GAL); 'icrs' keeps equatorial axes.
  Centre            : subtract the median position so the cluster sits at the origin.
  Residual PM vector: (pmRA - median pmRA, pmDE - median pmDE), normalised to unit length for the 2-D arrow map
                      (direction only; magnitudes are the columns res_pmRA / res_pmDE in the CSV output).
  Elongation        : sqrt(largest / mean-of-other-two eigenvalues) of the 3-D position covariance; ~1 = round.
  implied_depth_pc  : sqrt(observed radial variance - variance expected from parallax errors alone) = true line-of-sight depth.

WORKFLOW
  member table -> drop rows without RA/Dec/Plx/pmRA/pmDE or Plx <= 0 -> unwrap RA across 0/360 -> residual PM vectors ->
  kin3d: (precision cut | shrinkage | raw) -> d = 1000/Plx -> (x, y, z) -> optional Galactic rotation -> shape diagnostics ->
  optional centring -> compute_kinematics() dict (the web app serialises it) / save_kinematics() CSVs
"""
from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd

# ===================== GLOBAL CONFIGURATION =====================
KIN3D_MODE = 'precision'        # 'precision' | 'shrink' | 'raw'
KIN3D_MAX_FRAC = 0.05           # precision mode: keep e_Plx/Plx <= this
KIN3D_DEPTH_PC = 5.0            # shrink mode: assumed cluster depth (pc) that sets the prior width tau
KIN3D_FRAME = 'galactic'        # 'galactic' | 'icrs'
KIN3D_CENTRE = True             # subtract the median position
RA_COLUMNS = ('RA_ICRS', 'RAJ2000', '_RAJ2000', 'ra', 'RA')               # first one found is used
DEC_COLUMNS = ('DE_ICRS', 'DEJ2000', '_DEJ2000', 'dec', 'DEC')
EPLX_COLUMNS = ('e_Plx', 'e_plx', 'parallax_error', 'plx_error')
STAR_TYPE_COLORS = {'Main Sequence': '#00FFFF', 'Pre-Main Sequence': '#BA8E23',
                    'White Dwarf': '#000000', 'Red Giant / Supergiant': '#FF0000'}
MIN_STARS_3D = 3                # minimum stars that must pass the distance filter
# --- velocities and dispersion (compute_velocities) ---
PM_TO_KMS = 4.74047             # km/s per (mas/yr * kpc) i.e. v = 4.74047 * mu[mas/yr] / plx[mas]
G_PC_KMS2_MSUN = 4.30091e-3     # gravitational constant in pc (km/s)^2 / Msun
VEL_PLX_MODE = 'star'           # 'star': v_i = k mu_i / plx_i (adds each star's parallax noise, as specified) | 'cluster': every star at the cluster parallax
EPSI_MAX = 1.0                  # keep astrometric_excess_noise (column 'epsi') < this (mas); skipped when the column is absent
EPSI_COLUMNS = ('epsi', 'astrometric_excess_noise')
RV_COLUMNS = ('RV', 'radial_velocity', 'RVDR3')       # Gaia RV (km/s) used for the perspective-expansion correction
MIN_RV_STARS = 3                # need at least this many RVs to trust the cluster radial velocity
PERSPECTIVE_CORRECTION = True
N_BOOT_KIN = 500                # bootstrap resamples for median v_out / v_az and for sigma_1D
KIN_BOOT_SEED = 42
MIN_STARS_DISP = 8              # fewer stars than this -> no dispersion fit
DELTA_BIC_MIN = 10.0            # outlier mixture accepted only if BIC(single) - BIC(mixture) exceeds this
ETA_VIRIAL = 10.0               # sigma_vir^2 = G M / (eta r_hm); eta ~ 10 for a Plummer-like cluster (Spitzer 1987; Fleck et al. 2006)
OUT_DIR = './starnest_output'     # only used by save_* helpers
# =================================================================

# Rotation matrix ICRS -> Galactic (Hipparcos/Gaia definition; ESA 1997, Hipparcos Catalogue Vol. 1, Eq. 1.5.11)
ICRS_TO_GAL = np.array([[-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
                        [0.4941094278755837, -0.4448296299600112, 0.7469822444972189],
                        [-0.8676661490190047, -0.1980763734312015, 0.4559837761750669]])


def kin3d(ra, dec, plx, eplx, mode=KIN3D_MODE, max_frac=KIN3D_MAX_FRAC, depth_pc=KIN3D_DEPTH_PC, frame=KIN3D_FRAME, centre=KIN3D_CENTRE):
    """3-D cluster map. sigma_d = d*e_Plx/Plx is purely radial (sky position is exact to ~mas) -> a round cluster becomes a needle
    along the line of sight. precision: keep e_Plx/Plx <= max_frac | shrink: empirical-Bayes parallax toward the cluster value
    (tau = plx_c^2*depth/1000) | raw: 1000/Plx for all. Mirrors kin3dCompute() in index.html."""
    n = len(ra)
    has_err = bool(np.isfinite(eplx).any() and (eplx > 0).any())
    ok_err = np.isfinite(eplx) & (eplx > 0)
    keep, plx_use, note = np.ones(n, bool), plx.copy(), ''
    if mode == 'precision':
        if not has_err:
            note = 'No e_Plx column in this file: precision cut skipped (raw 1/Plx).'
        else:
            keep = ok_err & (eplx / plx <= max_frac)
    elif mode == 'shrink':
        if not has_err:
            note = 'No e_Plx column in this file: shrinkage skipped (raw 1/Plx).'
        else:
            w = np.where(ok_err, 1.0 / np.where(ok_err, eplx, 1.0) ** 2, 0.0)
            pc = np.sum(w * plx) / w.sum()
            wt = 1.0 / (pc ** 2 * depth_pc / 1000.0) ** 2
            plx_use = np.where(ok_err, (plx * w + pc * wt) / (w + wt), plx)
    if keep.sum() < MIN_STARS_3D:
        raise ValueError(f'3D map: only {int(keep.sum())} stars pass the distance filter (mode {mode}). Loosen it or use another mode.')
    d = 1000.0 / plx_use[keep]
    a, b = np.radians(ra[keep]), np.radians(dec[keep])
    P = np.c_[d * np.cos(b) * np.cos(a), d * np.cos(b) * np.sin(a), d * np.sin(b)]
    if frame == 'galactic':
        P = P @ ICRS_TO_GAL.T
    med = np.median(P, axis=0)
    diag = None
    if len(P) >= 4:
        w_, v_ = np.linalg.eigh(np.cov((P - med).T, bias=True))
        los = med / np.linalg.norm(med)
        rest = max((w_[0] + w_[1]) / 2, 1e-12)
        diag = dict(major_sd_pc=math.sqrt(w_[2]), minor_sd_pc=math.sqrt(rest), elongation=math.sqrt(w_[2] / rest),
                    angle_to_los_deg=math.degrees(math.acos(min(1.0, abs(float(v_[:, 2] @ los))))))
        # Observed radial spread vs the spread parallax errors alone would produce (sigma_d = d*e_Plx/Plx, per star, RMS)
        C = np.cov(P.T, bias=True)
        rad = math.sqrt(max(float(los @ C @ los), 0.0))
        diag['radial_sd_pc'] = rad
        diag['transverse_sd_pc'] = math.sqrt(max((float(np.trace(C)) - rad ** 2) / 2, 0.0))
        if mode != 'shrink' and has_err:
            e_ok = ok_err[keep]
            if e_ok.any():
                sd_d = d[e_ok] * eplx[keep][e_ok] / plx_use[keep][e_ok]
                diag['expected_err_sd_pc'] = float(np.sqrt(np.mean(sd_d ** 2)))
                diag['implied_depth_pc'] = math.sqrt(max(rad ** 2 - diag['expected_err_sd_pc'] ** 2, 0.0))
    if centre:
        P = P - med
    return dict(keep=keep, P=P, d=d, frac=(eplx[keep] / plx[keep]),
                meta=dict(mode=mode, frame=frame, centred=bool(centre), max_frac=max_frac, depth_pc=depth_pc, n_used=int(keep.sum()),
                          n_dropped=int(n - keep.sum()), median_d_pc=float(np.median(d)), note=note, diag=diag))


def compute_kinematics(base: pd.DataFrame, mode=KIN3D_MODE, max_frac=KIN3D_MAX_FRAC, depth_pc=KIN3D_DEPTH_PC,
                       frame=KIN3D_FRAME, centre=KIN3D_CENTRE) -> dict:
    """Member table (needs RA, Dec, Plx, pmRA, pmDE; optional e_Plx, Star_Type) -> dict with 'meta', 'positions_3d', 'meta3d',
    'vectors_2d' (the structure the web app plots). Raises ValueError when the table cannot be used."""
    ra_k = next((k for k in RA_COLUMNS if k in base), None)
    de_k = next((k for k in DEC_COLUMNS if k in base), None)
    if ra_k is None or de_k is None:
        raise ValueError('No RA/Dec columns in the file.')
    d = base.dropna(subset=[ra_k, de_k, 'Plx', 'pmRA', 'pmDE'])
    d = d[d['Plx'] > 0]
    if d.empty:
        raise ValueError('No stars with valid RA/Dec, parallax and proper motion for the kinematics plots.')
    ra = d[ra_k].to_numpy(float)
    if ra.max() - ra.min() > 180:                       # cluster straddles RA = 0/360
        ra = np.where(ra < 180, ra + 360, ra)
    dec = d[de_k].to_numpy(float)
    eplx_k = next((k for k in EPLX_COLUMNS if k in d), None)
    eplx = d[eplx_k].to_numpy(float) if eplx_k else np.full(len(d), np.nan)
    res_ra, res_de = d['pmRA'].to_numpy(float) - np.median(d['pmRA']), d['pmDE'].to_numpy(float) - np.median(d['pmDE'])
    vn = np.hypot(res_ra, res_de)
    vn = np.where(vn < 1e-8, 1.0, vn)
    st = d['Star_Type'].tolist() if 'Star_Type' in d else ['Main Sequence'] * len(d)
    k3 = kin3d(ra, dec, d['Plx'].to_numpy(float), eplx, mode=mode, max_frac=max_frac, depth_pc=depth_pc, frame=frame, centre=centre)
    st3 = [st[i] for i in np.flatnonzero(k3['keep'])]
    col = lambda types: [STAR_TYPE_COLORS.get(t, 'gray') for t in types]
    return dict(
        meta=dict(n_in=len(base), n_used=len(d), n_dropped=len(base) - len(d)),
        positions_3d=dict(x=k3['P'][:, 0], y=k3['P'][:, 1], z=k3['P'][:, 2], d=k3['d'], frac=k3['frac'], color=col(st3), star_type=st3),
        meta3d=k3['meta'],
        vectors_2d=dict(ra=ra, dec=dec, u=res_ra / vn, v=res_de / vn, color=col(st), star_type=st,
                        center_ra=float(np.median(ra)), center_dec=float(np.median(dec)),
                        res_pmra=res_ra, res_pmde=res_de))


def tangent_offsets(ra, dec):
    """Offsets (radians) of each star from the cluster centre on the tangent plane: theta_x = d(RA) cos(dec0) (east), theta_y = d(Dec) (north)."""
    ra0, dec0 = np.median(ra), np.median(dec)
    return np.radians((ra - ra0) * np.cos(np.radians(dec0))), np.radians(dec - dec0)


def perspective_correction(ra, dec, pmra, pmde, plx, rv_cluster):
    """Remove the apparent expansion / contraction that a cluster's radial motion imprints on its proper motions.
    TERMINOLOGY: rv_cluster = mean radial velocity (km/s, + = receding); theta = angular offset from the centre (rad); plx = parallax (mas).
    A star at offset theta from the centre direction moving with the common space velocity V has mu = (V_t - v_r*theta)/d, so the
    cluster looks like it is shrinking (v_r > 0) or expanding (v_r < 0) even if nothing moves internally. Delta_mu = -(v_r*plx/k)*theta.
    WORKFLOW: 1) tangent-plane offsets; 2) Delta_mu = -(rv * plx / PM_TO_KMS) * theta per axis; 3) corrected mu = mu - Delta_mu.
    Returns (pmra_corr, pmde_corr, delta_pmra, delta_pmde)."""
    tx, ty = tangent_offsets(ra, dec)
    d_x, d_y = -(rv_cluster * plx / PM_TO_KMS) * tx, -(rv_cluster * plx / PM_TO_KMS) * ty
    return pmra - d_x, pmde - d_y, d_x, d_y


def decompose_out_az(tx, ty, vx, vy):
    """Split residual velocities into the component pointing away from the centre (v_out) and the counter-clockwise tangential one (v_az).
    TERMINOLOGY: r_hat = offset direction, phi_hat = r_hat rotated +90 deg; v_out = v . r_hat (expansion > 0), v_az = v . phi_hat (rotation).
    WORKFLOW: 1) unit vectors from the offsets (stars at the exact centre get 0); 2) dot products."""
    r = np.hypot(tx, ty)
    r = np.where(r < 1e-12, np.inf, r)
    rx, ry = tx / r, ty / r
    return vx * rx + vy * ry, -vx * ry + vy * rx


def _boot_median(x, n_boot, rng):
    x = np.asarray(x, float)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    return float(np.median(x)), float(np.std(np.median(x[idx], axis=1), ddof=1))


def _nll_disp(theta, V, C, mixture):
    """Negative log-likelihood of 2-D velocities with per-star covariance (Walker et al. 2006, extended to 2-D, isotropic intrinsic sigma).
    params: mu_x, mu_y, ln sigma [, logit f_out, ln sigma_out]. Each star: N2(v - mu; sigma^2 I + C_i); the outlier component is a broad
    N2 with sigma_out^2 I + C_i around the same mean."""
    mx, my, ls = theta[0], theta[1], theta[2]
    dv = V - np.array([mx, my])

    def lnN(sig2):
        a, b, c = C[:, 0, 0] + sig2, C[:, 0, 1], C[:, 1, 1] + sig2
        det = a * c - b * b
        q = (c * dv[:, 0] ** 2 - 2 * b * dv[:, 0] * dv[:, 1] + a * dv[:, 1] ** 2) / det
        return -0.5 * q - 0.5 * np.log(det) - np.log(2 * np.pi)
    ll = lnN(np.exp(2 * ls))
    if mixture:
        f = 1.0 / (1.0 + np.exp(-theta[3]))
        ll = np.logaddexp(np.log1p(-f) + ll, np.log(f) + lnN(np.exp(2 * theta[4])))
    return -float(np.sum(ll))


def fit_dispersion(V, C, mixture: bool = False, x0_sigma: float | None = None) -> dict:
    """Maximum-likelihood mean velocity and intrinsic 1-D dispersion.
    TERMINOLOGY: V = (n,2) velocities (km/s); C = (n,2,2) per-star measurement covariance; sigma = intrinsic 1-D dispersion (the same
    along both axes); mixture adds a broad outlier population (fraction f, width sigma_out).
    WORKFLOW: 1) start from the measured scatter; 2) Nelder-Mead on the NLL; 3) report mean, sigma, (f, sigma_out), NLL and BIC."""
    from scipy.optimize import minimize
    n = len(V)
    s0 = x0_sigma if x0_sigma is not None else max(0.05, float(np.sqrt(max(np.var(V, axis=0).mean() - np.mean(C[:, [0, 1], [0, 1]]), 0.01))))
    x0 = [float(np.median(V[:, 0])), float(np.median(V[:, 1])), np.log(s0)] + ([-2.0, np.log(5 * s0 + 1.0)] if mixture else [])
    r = minimize(_nll_disp, x0, args=(V, C, mixture), method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-8, maxiter=4000))
    k = len(x0)
    out = dict(mu_x=float(r.x[0]), mu_y=float(r.x[1]), sigma=float(np.exp(r.x[2])), nll=float(r.fun), n_par=k,
               bic=float(2 * r.fun + k * np.log(n)), converged=bool(r.success))
    if mixture:
        out.update(f_out=float(1.0 / (1.0 + np.exp(-r.x[3]))), sigma_out=float(np.exp(r.x[4])))
    return out


def virial_sigma(mass_msun: float, r_hm_pc: float, eta: float = ETA_VIRIAL) -> float:
    """Virial 1-D velocity dispersion sigma = sqrt(G M / (eta r_hm)).
    TERMINOLOGY: M = cluster mass (Msun, e.g. M_T from binaries.cluster_mass); r_hm = half-mass radius (pc); eta ~ 10 (structure factor, configurable).
    Comparing it with the measured sigma_1D tells whether the cluster is bound (measured ~ virial) or expanding / supervirial (measured >> virial)."""
    if mass_msun <= 0 or r_hm_pc <= 0 or eta <= 0:
        raise ValueError('mass, r_hm and eta must be positive.')
    return float(np.sqrt(G_PC_KMS2_MSUN * mass_msun / (eta * r_hm_pc)))


def half_number_radius_pc(ra, dec, d_pc) -> float:
    """Projected half-number radius (pc): median angular distance from the centre x distance. A proxy for r_hm when masses are not weighted
    (3-D r_hm ~ 4/3 x projected for a spherical cluster; the caller decides)."""
    tx, ty = tangent_offsets(ra, dec)
    return float(np.median(np.hypot(tx, ty)) * d_pc)


def compute_velocities(base: pd.DataFrame, eta: float = ETA_VIRIAL, plx_mode: str = VEL_PLX_MODE, epsi_max: float | None = EPSI_MAX,
                       perspective: bool = PERSPECTIVE_CORRECTION, mass_msun: float | None = None, n_boot: int = N_BOOT_KIN,
                       seed: int = KIN_BOOT_SEED) -> dict:
    """Tangential velocities, outward / azimuthal motion and the intrinsic velocity dispersion of a member table.
    TERMINOLOGY: v = k mu / plx (km/s); v_out / v_az = residual velocity along / around the direction from the cluster centre; sigma_1D =
    intrinsic 1-D dispersion by maximum likelihood with each star's own error covariance; mixture = broad outlier component accepted when it
    lowers BIC by > DELTA_BIC_MIN (binaries and non-members inflate sigma if forced into one Gaussian).
    WORKFLOW: 1) drop rows lacking astrometry; apply epsi < epsi_max when the column exists; 2) cluster RV from Gaia RVs -> perspective
    correction of the PMs (skipped and reported when fewer than MIN_RV_STARS RVs); 3) v = k mu / plx (star or cluster parallax);
    4) per-star 2x2 velocity covariance from e_pm and e_plx; 5) median v_out and v_az of the residual velocities with bootstrap errors;
    6) ML fit with one Gaussian and with the outlier mixture, choose by delta-BIC; 7) sigma error by bootstrap of the chosen model;
    8) optional virial sigma if mass_msun is given."""
    d0 = base.dropna(subset=['RA_ICRS', 'DE_ICRS', 'Plx', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE', 'e_Plx'])
    d0 = d0[d0['Plx'] > 0]
    notes = []
    epsi_k = next((k for k in EPSI_COLUMNS if k in d0), None)
    if epsi_max is not None and epsi_k is not None:
        d0 = d0[d0[epsi_k] < epsi_max]
    elif epsi_max is not None:
        notes.append('no astrometric_excess_noise (epsi) column: EPSI_MAX cut skipped')
    if len(d0) < MIN_STARS_DISP:
        raise ValueError(f'Only {len(d0)} usable stars; need >= {MIN_STARS_DISP} for a velocity dispersion.')
    ra, dec = d0['RA_ICRS'].to_numpy(float), d0['DE_ICRS'].to_numpy(float)
    if ra.max() - ra.min() > 180:
        ra = np.where(ra < 180, ra + 360, ra)
    plx, e_plx = d0['Plx'].to_numpy(float), d0['e_Plx'].to_numpy(float)
    pmra, pmde = d0['pmRA'].to_numpy(float), d0['pmDE'].to_numpy(float)
    e_ra, e_de = d0['e_pmRA'].to_numpy(float), d0['e_pmDE'].to_numpy(float)

    rv_k = next((k for k in RV_COLUMNS if k in d0), None)
    rv_c, n_rv = float('nan'), 0
    if perspective and rv_k is not None:
        rv = d0[rv_k].to_numpy(float)
        ok = np.isfinite(rv)
        n_rv = int(ok.sum())
        if n_rv >= MIN_RV_STARS:
            rv_c = float(np.median(rv[ok]))                       # median: a few binaries with wild RVs must not drive the correction
            pmra, pmde, _, _ = perspective_correction(ra, dec, pmra, pmde, plx, rv_c)
        else:
            notes.append(f'only {n_rv} Gaia RVs (< {MIN_RV_STARS}): perspective correction skipped')
    elif perspective:
        notes.append('no RV column: perspective correction skipped')

    plx_c = float(np.median(plx))
    plx_use = plx if plx_mode == 'star' else np.full_like(plx, plx_c)
    vx, vy = PM_TO_KMS * pmra / plx_use, PM_TO_KMS * pmde / plx_use
    C = np.zeros((len(d0), 2, 2))
    kk = PM_TO_KMS / plx_use
    C[:, 0, 0], C[:, 1, 1] = (kk * e_ra) ** 2, (kk * e_de) ** 2
    if plx_mode == 'star':                                        # parallax error moves both components together -> off-diagonal covariance
        C[:, 0, 0] += (vx / plx_use * e_plx) ** 2
        C[:, 1, 1] += (vy / plx_use * e_plx) ** 2
        C[:, 0, 1] = C[:, 1, 0] = vx * vy / plx_use ** 2 * e_plx ** 2
    V = np.column_stack([vx, vy])

    rng = np.random.default_rng(seed)
    tx, ty = tangent_offsets(ra, dec)
    res_x, res_y = vx - np.median(vx), vy - np.median(vy)
    v_out, v_az = decompose_out_az(tx, ty, res_x, res_y)
    m_out, e_out = _boot_median(v_out, n_boot, rng)
    m_az, e_az = _boot_median(v_az, n_boot, rng)

    f1 = fit_dispersion(V, C, mixture=False)
    f2 = fit_dispersion(V, C, mixture=True)
    use_mix = (f1['bic'] - f2['bic']) > DELTA_BIC_MIN and f2['converged']
    best, key = (f2, True) if use_mix else (f1, False)
    idx = rng.integers(0, len(V), size=(max(50, n_boot // 5), len(V)))
    boots = []
    for i in idx:                                                 # bootstrap the chosen model (fewer draws: each is an optimisation)
        try:
            boots.append(fit_dispersion(V[i], C[i], mixture=key, x0_sigma=best['sigma'])['sigma'])
        except Exception:
            continue
    sig_err = float(np.std(boots, ddof=1)) if len(boots) > 5 else float('nan')
    out = dict(n=len(d0), rv_cluster_kms=rv_c, n_rv=n_rv, plx_mode=plx_mode, plx_cluster_mas=plx_c, notes=notes,
               v_out_kms=m_out, v_out_err=e_out, v_az_kms=m_az, v_az_err=e_az,
               sigma_1d_kms=best['sigma'], sigma_1d_err=sig_err, mixture_used=bool(use_mix), delta_bic=float(f1['bic'] - f2['bic']),
               single=f1, mixture=f2, mean_vx_kms=best['mu_x'], mean_vy_kms=best['mu_y'],
               median_meas_err_kms=float(np.median(np.sqrt(0.5 * (C[:, 0, 0] + C[:, 1, 1])))))
    if mass_msun:
        r_hn = half_number_radius_pc(ra, dec, 1000.0 / plx_c)
        out.update(r_hn_pc=r_hn, sigma_vir_kms=virial_sigma(mass_msun, r_hn * 4.0 / 3.0, eta), mass_msun=float(mass_msun), eta=eta)
    return out


def save_kinematics(res: dict, tag: str, out_dir: str = OUT_DIR) -> dict:
    """Write the 3-D positions and 2-D residual vectors as CSV (index=False)."""
    os.makedirs(out_dir, exist_ok=True)
    p3, v2 = res['positions_3d'], res['vectors_2d']
    f3 = os.path.join(out_dir, f'{tag}_kinematics_3d.csv')
    pd.DataFrame({'x_pc': p3['x'], 'y_pc': p3['y'], 'z_pc': p3['z'], 'd_pc': p3['d'], 'eplx_over_plx': p3['frac'],
                  'Star_Type': p3['star_type']}).to_csv(f3, index=False)
    f2 = os.path.join(out_dir, f'{tag}_kinematics_pm_vectors.csv')
    pd.DataFrame({'RA': v2['ra'], 'Dec': v2['dec'], 'res_pmRA': v2['res_pmra'], 'res_pmDE': v2['res_pmde'],
                  'u_unit': v2['u'], 'v_unit': v2['v'], 'Star_Type': v2['star_type']}).to_csv(f2, index=False)
    return {'positions_3d': f3, 'pm_vectors': f2}
