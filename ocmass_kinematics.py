"""
ocmass_kinematics.py -- Stage 4: 3-D cluster map + proper-motion (residual) vectors (no UI code).

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
OUT_DIR = './ocmass_output'     # only used by save_* helpers
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
