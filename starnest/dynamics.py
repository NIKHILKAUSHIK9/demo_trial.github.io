"""
dynamics.py -- Stage 6: mass segregation, tidal radius, relaxation time, optional orbit / tail check (no UI code).

TERMINOLOGY
  MST / MSR : minimum spanning tree = shortest set of edges connecting all points; MSR (mass-segregation ratio, Allison et al. 2009)
              Lambda = <l_random> / l_massive, where l is the MST length of the N most massive stars and <l_random> the mean over
              random N-star subsets. Lambda ~ 1: no segregation; Lambda > 1: massive stars are more concentrated than average.
  r_t       : tidal radius, r_t^3 = G M / (2 (A - B)^2) (Oort constants A, B; the Galactic tide in the disc), pc.
  t_relax   : (two-body) relaxation time ~ N / (8 ln N) * t_cross, t_cross = D / sigma_v  (Binney & Tremaine 2008). Compare with the
              cluster age: age >> t_relax means the cluster has had time to segregate dynamically (rather than being born so).
  tail check: do the member positions stretch along the orbit? (heuristic: angle between their principal axis and the orbit tangent)

WORKFLOW
  members + masses -> tangent-plane positions (pc) -> MSR sweep over N_MASSIVE -> r_t from M_T -> t_relax from (N, D, sigma_1D)
  -> optionally galpy MWPotential2014 orbit -> principal-axis / orbit angle (tail-like or not).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform

# ===================== GLOBAL CONFIGURATION =====================
N_MASSIVE = (5, 10, 20, 40)     # sizes of the "most massive" subsets in the MSR sweep
N_RANDOM = 100                  # random subsets per N (Allison et al. 2009)
MSR_SEED = 42
OORT_A, OORT_B = 15.3, -11.9    # km/s/kpc (Bovy 2017)
G_PC_KMS2_MSUN = 4.30091e-3     # pc (km/s)^2 / Msun
PC_PER_KMS_IN_MYR = 0.9778      # 1 pc / (km/s) = 0.9778 Myr
D_FACTOR = 2.0                  # cluster diameter D = D_FACTOR x r_hm for the crossing time
ORBIT_T_MYR = 50.0              # orbit integrated this far back and forward (tail check)
ORBIT_N_STEPS = 501
R0_KPC, V0_KMS = 8.0, 220.0     # galpy MWPotential2014 natural units
TAIL_ELONGATION_MIN = 1.5       # principal-axis ratio needed to call a member distribution elongated
TAIL_ANGLE_MAX_DEG = 30.0       # ... and aligned with the orbit within this angle
# =================================================================


def mst_length(xy: np.ndarray) -> float:
    """Total edge length of the Euclidean minimum spanning tree of an (n, 2) array (n >= 2)."""
    return float(minimum_spanning_tree(squareform(pdist(xy))).sum())


def mass_segregation_ratio(xy: np.ndarray, mass: np.ndarray, n_massive=N_MASSIVE, n_random: int = N_RANDOM, seed: int = MSR_SEED) -> pd.DataFrame:
    """Allison et al. (2009) MSR for several subset sizes.
    TERMINOLOGY: xy = positions (any linear unit; the ratio is scale-free); l_massive = MST length of the N heaviest stars;
    l_rand = MST length of N randomly chosen stars; Lambda = mean(l_rand) / l_massive; Lambda_err_lo / hi = (mean -+ std)(l_rand)/l_massive.
    WORKFLOW: 1) rank by mass; 2) for each N: MST of the top N; 3) n_random random N-subsets -> mean and std of their MST; 4) ratio.
    N >= number of stars gives NaN (the 'subset' would be the whole cluster). Ties in mass are broken by array order."""
    xy, mass = np.asarray(xy, float), np.asarray(mass, float)
    ok = np.isfinite(mass) & np.isfinite(xy).all(axis=1)
    xy, mass = xy[ok], mass[ok]
    n = len(mass)
    rng = np.random.default_rng(seed)
    order = np.argsort(-mass, kind='stable')
    rows = []
    for nm in n_massive:
        if nm < 2 or nm >= n:
            rows.append(dict(N=nm, l_massive=np.nan, l_random=np.nan, l_random_std=np.nan, Lambda=np.nan, Lambda_lo=np.nan, Lambda_hi=np.nan))
            continue
        l_m = mst_length(xy[order[:nm]])
        l_r = np.array([mst_length(xy[rng.choice(n, nm, replace=False)]) for _ in range(n_random)])
        mean, sd = float(l_r.mean()), float(l_r.std(ddof=1))
        rows.append(dict(N=nm, l_massive=l_m, l_random=mean, l_random_std=sd, Lambda=mean / l_m, Lambda_lo=(mean - sd) / l_m, Lambda_hi=(mean + sd) / l_m))
    return pd.DataFrame(rows)


def tidal_radius(mass_msun: float, a: float = OORT_A, b: float = OORT_B) -> float:
    """Tidal radius (pc): r_t = (G M / (2 (A - B)^2))^(1/3).
    TERMINOLOGY: A, B = Oort constants (km/s/kpc); (A - B) = angular speed of the local circular orbit; the factor 2 is the disc-tide
    approximation 4A(A-B) - kappa^2 ... reduced to the spec's form. Stars beyond r_t are lost to the Galaxy.
    WORKFLOW: 1) (A - B) in km/s/pc; 2) cube root of G M / (2 (A-B)^2)."""
    if mass_msun <= 0:
        raise ValueError('mass must be positive.')
    ab = (a - b) * 1e-3                                      # km/s/kpc -> km/s/pc
    return float((G_PC_KMS2_MSUN * mass_msun / (2.0 * ab ** 2)) ** (1.0 / 3.0))


def relaxation_time(n_stars: int, r_hm_pc: float, sigma_kms: float, d_factor: float = D_FACTOR) -> dict:
    """Relaxation time in Myr: t_E = N / (8 ln N) * D / sigma_v.
    TERMINOLOGY: N = number of stars (use the TOTAL, e.g. M_T / <m>, not just the observed members, when known); D = d_factor x r_hm (pc);
    sigma_v = velocity dispersion (km/s); t_cross = D / sigma_v.
    WORKFLOW: 1) t_cross = D / sigma in Myr; 2) multiply by N / (8 ln N)."""
    if n_stars < 2 or r_hm_pc <= 0 or sigma_kms <= 0:
        raise ValueError('need n_stars >= 2, r_hm > 0, sigma > 0.')
    t_cross = d_factor * r_hm_pc / sigma_kms * PC_PER_KMS_IN_MYR
    return dict(t_cross_myr=float(t_cross), t_relax_myr=float(n_stars / (8.0 * np.log(n_stars)) * t_cross), n_stars=int(n_stars))


def orbit_tail_check(ra, dec, plx_mas, pmra, pmde, rv_kms, xyz_members_pc: np.ndarray, t_myr: float = ORBIT_T_MYR, n_steps: int = ORBIT_N_STEPS) -> dict:
    """Integrate the cluster orbit in MWPotential2014 and test whether the members are stretched along it.
    TERMINOLOGY: ra, dec (deg), plx_mas, pmra, pmde (mas/yr), rv_kms = the cluster's mean 6-D phase-space point; xyz_members_pc = member
    positions in heliocentric Galactic Cartesian axes (e.g. kinematics.kin3d(...)['P']); elongation = sqrt(largest / mean of the other two
    covariance eigenvalues); angle = angle between the principal axis and the orbit tangent (0 or 180 = aligned).
    WORKFLOW: 1) galpy Orbit from the mean point; 2) integrate -t..+t Myr; 3) tangent at t=0 in heliocentric Galactic axes (galpy helioX/Y/Z);
    4) PCA of the member positions; 5) tail-like if elongation >= TAIL_ELONGATION_MIN and angle <= TAIL_ANGLE_MAX_DEG.
    Heuristic only: a real tail test needs the extended field, which the quality-cut member list does not contain.
    Requires galpy; raises RuntimeError when it is missing."""
    try:
        import astropy.units as u
        from galpy.orbit import Orbit
        from galpy.potential import MWPotential2014
    except ImportError as e:
        raise RuntimeError(f'galpy is not installed ({e}); orbit / tail check skipped.')
    o = Orbit([float(ra), float(dec), 1.0 / float(plx_mas), float(pmra), float(pmde), float(rv_kms)], radec=True, ro=R0_KPC, vo=V0_KMS)
    ts = np.linspace(-t_myr, t_myr, n_steps) * u.Myr
    o.integrate(ts, MWPotential2014)
    hx, hy, hz = (np.asarray(f(ts)) * 1000.0 for f in (o.helioX, o.helioY, o.helioZ))       # kpc -> pc
    path = np.column_stack([hx, hy, hz])
    mid = n_steps // 2
    tangent = path[mid + 1] - path[mid - 1]
    tangent = tangent / np.linalg.norm(tangent)
    P = np.asarray(xyz_members_pc, float)
    P = P[np.isfinite(P).all(axis=1)]
    if len(P) < 5:
        raise ValueError('need >= 5 member positions for the tail check.')
    w, v = np.linalg.eigh(np.cov((P - P.mean(0)).T))
    elong = float(np.sqrt(w[2] / max(np.mean(w[:2]), 1e-12)))
    cosang = abs(float(np.dot(v[:, 2], tangent)))
    ang = float(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
    return dict(elongation=elong, angle_deg=ang, tail_like=bool(elong >= TAIL_ELONGATION_MIN and ang <= TAIL_ANGLE_MAX_DEG),
                orbit_path_pc=path, t_myr=np.asarray(ts.value), orbit_tangent=tangent, principal_axis=v[:, 2])
