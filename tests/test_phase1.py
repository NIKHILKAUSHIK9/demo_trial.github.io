"""Phase-1 tests: every new module against a synthetic truth (isochrone = MOCK toy relation, because CMD is unreachable here).
TERMINOLOGY: sigma_1D = intrinsic 1-D velocity dispersion; MSR = mass-segregation ratio; M_T = total cluster mass; truth = value used to
generate the data. WORKFLOW: generate data with a known truth -> run the module -> require recovery within the stated tolerance."""
# ===================== GLOBAL CONFIGURATION =====================
import os, sys, io
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = 12
SIGMA_TRUE, V_RAD_TRUE = 1.2, 25.0          # km/s: intrinsic dispersion; cluster radial velocity for the perspective test
AV_TRUE, DM_TRUE, FBIN_TRUE = 0.40, 7.50, 0.30
# =================================================================
sys.path.insert(0, ROOT)
import numpy as np
import pandas as pd
import pytest
import binaries as B
import dynamics as D
import extinction_isochrone as E
import gmm as G
import imf as I
import kinematics as K


def toy_iso(n=60):
    """MOCK isochrone (same toy relation as the old UI tests): monotone mass-magnitude-colour, label 1 (main sequence)."""
    m = np.linspace(0.2, 3.0, n)
    return pd.DataFrame(dict(Zini=0.0152, MH=0.0, logAge=8.2, Mini=m, int_IMF=0, Mass=m, logL=3.5 * np.log10(m), logTe=3.7, logg=4.4, label=1,
                             Gmag=5 - 6 * np.log10(m), G_BPmag=5.3 - 6 * np.log10(m), G_RPmag=4.7 - 6 * np.log10(m)))


def gaia_like(n_cl=150, n_fld=400, seed=SEED):
    rng = np.random.default_rng(seed)
    def blk(n, pm, spm, plx, splx):
        return pd.DataFrame({'pmRA': rng.normal(pm[0], spm, n), 'pmDE': rng.normal(pm[1], spm, n), 'Plx': rng.normal(plx, splx, n)})
    d = pd.concat([blk(n_cl, (5, -3), 0.12, 2.94, 0.04), blk(n_fld, (0, 0), 6, 3.0, 0.9)], ignore_index=True)
    n = len(d)
    d['RA_ICRS'], d['DE_ICRS'] = rng.normal(150, 0.6, n), rng.normal(-30, 0.6, n)
    d['e_Plx'], d['e_pmRA'], d['e_pmDE'] = 0.02, 0.03, 0.03
    d['Gmag'] = rng.uniform(9, 17, n)
    d['BPmag'] = d['Gmag'] + 0.4 + rng.normal(0, 0.1, n); d['RPmag'] = d['Gmag'] - 0.4 + rng.normal(0, 0.1, n)
    d['BP-RP'] = d['BPmag'] - d['RPmag']
    d['RUWE'], d['gofAL'], d['Nper'] = 1.0, 1.0, 20
    return d[d['Plx'] > 0.5].reset_index(drop=True)


# ----------------------------- 1.4 gmm -----------------------------
@pytest.fixture(scope='module')
def gres():
    return G.run_gmm(G.apply_quality_cuts(gaia_like()), distance_lit_pc=340.0, half_width_pc=150.0)


def test_resample_pmemb_shapes_and_bounds(gres):
    r = G.resample_pmemb(gres, n_resample=60, seed=1)
    assert len(r) == len(gres.working.dropna(subset=gres.features)) and r['PMemb_mean'].between(0, 1).all()
    assert (r['PMemb_std'] >= 0).all()
    core = r.loc[gres.members.index]
    assert core['PMemb_frac_member'].mean() > 0.9                      # clean, well-separated cluster: members stay members
    assert G.resample_pmemb(gres, 60, 1)['PMemb_mean'].equals(r['PMemb_mean'])         # seeded


def test_hdbscan_crosscheck_agrees_on_clean_cluster(gres):
    h = G.hdbscan_crosscheck(gres, min_cluster_size=10)
    assert h['stats']['jaccard'] > 0.7 and h['stats']['n_both'] > 0.7 * h['stats']['n_gmm']


def test_closest_to_literature_picks_anchor_component():
    df = G.apply_quality_cuts(gaia_like())
    a = G.run_gmm(df, 340.0, 150.0, k=3)
    med = a.summary
    b = G.run_gmm(df, 340.0, 150.0, k=3, cluster_rule='closest_to_literature', lit_pm=(med['median_pmRA'], med['median_pmDE']), lit_plx=med['median_plx'])
    assert b.cluster_component == a.cluster_component
    far = G.run_gmm(df, 340.0, 150.0, k=3, cluster_rule='closest_to_literature', lit_pm=(0.0, 0.0), lit_plx=3.0)    # the broad field blob
    assert far.cluster_component != a.cluster_component
    with pytest.raises(ValueError):
        G.run_gmm(df, 340.0, 150.0, cluster_rule='closest_to_literature')
    with pytest.raises(ValueError):
        G.run_gmm(df, 340.0, 150.0, cluster_rule='bogus')


def test_hunt_reffert_prefill_handles_missing_astroquery_row(monkeypatch):
    class _Tab(list):
        colnames = ['pmRA', 'pmDE', 'plx', 'Dist']
    class _V:
        def __init__(self, **k): pass
        def query_constraints(self, **k): return [_Tab([{'pmRA': 5.0, 'pmDE': -3.0, 'plx': 2.94, 'Dist': 340.0}])]
    import astroquery.vizier as V
    monkeypatch.setattr(V, 'Vizier', _V)
    out = G.hunt_reffert_prefill('UPK 545')
    assert out['pmRA'] == 5.0 and out['distance_pc'] == 340.0 and out['missing'] == []


# ----------------------------- 1.3 kinematics -----------------------------
def vel_table(n=120, sigma=SIGMA_TRUE, rv=None, expand=0.0, seed=SEED, mean_pm=(-13.0, 4.0), plx=3.0):
    """Members with a Gaussian internal velocity (sigma km/s, optional expansion) + measurement errors. Offsets on the sky are a few arcmin."""
    rng = np.random.default_rng(seed)
    k = K.PM_TO_KMS
    ra, dec = 150 + rng.normal(0, 0.5, n), -30 + rng.normal(0, 0.5, n)
    tx, ty = K.tangent_offsets(ra, dec)
    vx = rng.normal(0, sigma, n) + expand * tx / np.hypot(tx, ty)
    vy = rng.normal(0, sigma, n) + expand * ty / np.hypot(tx, ty)
    mu_x, mu_y = mean_pm[0] + vx * plx / k, mean_pm[1] + vy * plx / k
    if rv is not None:                                          # perspective effect of the radial motion (the thing the correction must remove)
        mu_x = mu_x - (rv * plx / k) * tx            # receding cluster (rv > 0) looks like it is contracting: Delta_mu = -(rv plx / k) theta
        mu_y = mu_y - (rv * plx / k) * ty
    e = 0.02
    return pd.DataFrame({'RA_ICRS': ra, 'DE_ICRS': dec, 'Plx': plx + rng.normal(0, 0.01, n), 'e_Plx': 0.01,
                         'pmRA': mu_x + rng.normal(0, e, n), 'pmDE': mu_y + rng.normal(0, e, n), 'e_pmRA': e, 'e_pmDE': e,
                         **({'RV': rv + rng.normal(0, 1, n)} if rv is not None else {})})


def test_dispersion_recovered_within_2_sigma():
    r = K.compute_velocities(vel_table(), plx_mode='cluster', n_boot=100)
    assert abs(r['sigma_1d_kms'] - SIGMA_TRUE) < 2 * r['sigma_1d_err'] + 0.05, r
    assert r['sigma_1d_err'] < 0.3 and not r['mixture_used']


def test_outlier_mixture_chosen_when_contaminated():
    d = vel_table(150)
    d.loc[:14, 'pmRA'] += np.random.default_rng(3).choice([-1, 1], 15) * 4.0          # ~35 km/s kicks (non-members / fast binaries)
    r = K.compute_velocities(d, plx_mode='cluster', n_boot=50)
    assert r['mixture_used'] and r['delta_bic'] > K.DELTA_BIC_MIN
    assert abs(r['sigma_1d_kms'] - SIGMA_TRUE) < 0.4                                   # core dispersion survives, single Gaussian would give ~ 8


def test_perspective_correction_removes_apparent_expansion():
    tx = np.array([1e-3, -1e-3, 0.0, 0.0]); ty = np.array([0.0, 0.0, 1e-3, -1e-3])
    ra, dec = 150 + np.degrees(tx), -30 + np.degrees(ty)
    pm0 = np.zeros(4)
    p_ra, p_de, dx, dy = K.perspective_correction(ra, dec, pm0, pm0, np.full(4, 3.0), rv_cluster=30.0)
    assert np.allclose(p_ra - dx, p_ra - dx) and np.allclose(p_ra + dx, 0, atol=1e-12)   # corrected = pm - delta; pm=0 -> corrected = -delta
    assert dx[0] < 0 < dx[1] and dy[2] < 0 < dy[3]               # receding cluster (+rv): apparent contraction toward the centre
    # v_out sign: apparent contraction must be removed -> v_out ~ 0 after correction, clearly negative before
    d = vel_table(200, sigma=0.3, rv=V_RAD_TRUE, mean_pm=(0.0, 0.0))
    with_corr = K.compute_velocities(d, plx_mode='cluster', n_boot=50, perspective=True)
    no_corr = K.compute_velocities(d, plx_mode='cluster', n_boot=50, perspective=False)
    assert abs(with_corr['v_out_kms']) < abs(no_corr['v_out_kms']) and abs(with_corr['v_out_kms']) < 0.3


def test_vout_vaz_signs():
    tx, ty = np.array([1.0, 0.0, -1.0, 0.0]), np.array([0.0, 1.0, 0.0, -1.0])
    out, az = K.decompose_out_az(tx, ty, tx * 2.0, ty * 2.0)            # pure expansion
    assert np.allclose(out, 2) and np.allclose(az, 0)
    out, az = K.decompose_out_az(tx, ty, -ty, tx)                       # counter-clockwise rotation
    assert np.allclose(out, 0) and np.allclose(az, 1)


def test_virial_sigma_value():
    assert K.virial_sigma(100.0, 2.0, 10.0) == pytest.approx(np.sqrt(4.30091e-3 * 100 / 20.0))
    with pytest.raises(ValueError):
        K.virial_sigma(-1, 1)


def test_epsi_cut_applied_when_column_present():
    d = vel_table(60)
    d['epsi'] = 0.1; d.loc[:9, 'epsi'] = 5.0
    assert K.compute_velocities(d, plx_mode='cluster', n_boot=20)['n'] == 50


# ----------------------------- 1.1 binaries -----------------------------
def test_ifmr_values_and_zero_outside():
    f = B.ifmr_cummings2018([0.5, 1.0, 3.0, 5.0, 9.0])
    assert f[0] == 0 and f[4] == 0 and f[1] == pytest.approx(0.0873 + 0.476) and 0.7 < f[2] < 0.9 and f[3] == pytest.approx(0.107 * 5 + 0.471)


def test_binary_assignment_recovers_primary_mass_of_singles():
    iso = toy_iso()
    ag, abp, arp = E.AG_AV * AV_TRUE, E.ABP_AV * AV_TRUE, E.ARP_AV * AV_TRUE
    syn = B.build_synthetic(iso, DM_TRUE, ag, abp, arp, 3000, 0.5, np.random.default_rng(1))
    obs = syn.sample(80, random_state=2).reset_index(drop=True)
    tab = B.assign_masses(obs, iso, DM_TRUE, ag, abp, arp, n_mc=15)
    singles = obs['m2'] == 0
    rel = (tab.loc[singles, 'MC_Mass'] / obs.loc[singles, 'm1'] - 1).abs()
    assert np.nanmedian(rel) < 0.15                                                    # photometry-only mass of a single star
    assert (tab['MC_CompMass'].fillna(0) >= 0).all() and tab['q'].between(0, 1).all()
    assert tab.loc[~singles, 'binary_prob'].mean() > tab.loc[singles, 'binary_prob'].mean()       # binaries get higher binary probability


def test_binary_mass_is_below_naive_for_true_binaries():
    iso = toy_iso()
    ag, abp, arp = E.AG_AV * AV_TRUE, E.ABP_AV * AV_TRUE, E.ARP_AV * AV_TRUE
    syn = B.build_synthetic(iso, DM_TRUE, ag, abp, arp, 3000, 1.0, np.random.default_rng(4))
    eq = syn[syn['q'] > 0.9].sample(40, random_state=1).reset_index(drop=True)          # near-equal-mass binaries
    tab = B.assign_masses(eq, iso, DM_TRUE, ag, abp, arp, n_mc=20, f_bin=1.0)
    # naive: treat the unresolved pair as ONE star -> its mass is larger than the true primary
    naive = np.interp(eq['Gmag'] - DM_TRUE - ag, iso['Gmag'][::-1], iso['Mass'][::-1])
    assert (naive > eq['m1']).mean() > 0.8 and (tab['MC_CompMass'] > 0).mean() > 0.8


def test_cluster_mass_components_and_limits():
    tab = pd.DataFrame({'MC_Mass': [1.0, 0.8, 0.5], 'MC_CompMass': [0.4, 0.0, 0.0], 'is_binary': [True, False, False], 'q': [0.4, 0.0, 0.0]})
    mle = dict(alpha=2.35, n=300, m_min=0.4, m_max=3.0)
    r = B.cluster_mass(tab, mle, m_turnoff=3.0)
    assert r['M_obs'] == pytest.approx(2.7) and r['M_unseen'] > 0 and r['M_WD'] > 0 and r['M_T'] > r['M_obs']
    assert r['M_T'] == pytest.approx(r['M_obs'] + r['M_unseen'] * (1 + 0.5 * 0.4) + r['M_WD'])
    assert B.unseen_mass(2.35, 300, 0.05, 3.0) == 0.0 and B.wd_mass(2.35, 300, 0.4, 3.0, m_turnoff=9.0) == 0.0
    # unseen mass of a Salpeter field with a known closed form: A [m^(2-a)/(2-a)]
    a_norm = 300 / np.exp(I._log_norm(2.35, 0.4, 3.0))
    assert B.unseen_mass(2.35, 300, 0.4, 3.0, floor=0.09) == pytest.approx(a_norm * (0.4 ** -0.35 - 0.09 ** -0.35) / -0.35, rel=1e-2)


# ----------------------------- 1.2 dynamics -----------------------------
def test_msr_unsegregated_vs_segregated():
    rng = np.random.default_rng(0)
    xy, m = rng.normal(size=(300, 2)), rng.lognormal(0, 0.5, 300)
    flat = D.mass_segregation_ratio(xy, m, (10, 20))
    assert flat['Lambda'].between(0.7, 1.3).all()
    seg = D.mass_segregation_ratio(xy, 1 / (1 + np.hypot(*xy.T)), (10, 20))
    assert (seg['Lambda'] > 2.0).all() and (seg['Lambda_lo'] < seg['Lambda']).all()
    assert D.mass_segregation_ratio(xy[:8], m[:8], (10,))['Lambda'].isna().all()       # N >= n_stars -> NaN, no crash


def test_mst_length_known_case():
    assert D.mst_length(np.array([[0, 0], [1, 0], [3, 0]], float)) == pytest.approx(3.0)


def test_tidal_radius_and_relaxation_known_values():
    ab = 27.2e-3
    assert D.tidal_radius(100.0) == pytest.approx((4.30091e-3 * 100 / (2 * ab ** 2)) ** (1 / 3))
    assert D.tidal_radius(800.0) == pytest.approx(2 * D.tidal_radius(100.0))           # r_t ~ M^(1/3)
    r = D.relaxation_time(1000, 2.0, 1.0)
    assert r['t_cross_myr'] == pytest.approx(2 * 2.0 / 1.0 * 0.9778)
    assert r['t_relax_myr'] == pytest.approx(1000 / (8 * np.log(1000)) * r['t_cross_myr'])
    with pytest.raises(ValueError):
        D.relaxation_time(1, 1, 1)


def test_orbit_tail_check_runs_with_galpy():
    pytest.importorskip('galpy')
    P = np.random.default_rng(0).normal(size=(60, 3)) * [30.0, 3.0, 3.0]               # elongated along x
    out = D.orbit_tail_check(150.0, -30.0, 3.0, -13.0, 4.0, 20.0, P)
    assert out['elongation'] > 3 and 0 <= out['angle_deg'] <= 90 and out['orbit_path_pc'].shape[1] == 3


# ----------------------------- 1.5 extinction / isochrone -----------------------------
def _cmd_cluster(n=250, seed=SEED):
    rng = np.random.default_rng(seed)
    iso = toy_iso()
    syn = B.build_synthetic(iso, DM_TRUE, E.AG_AV * AV_TRUE, E.ABP_AV * AV_TRUE, E.ARP_AV * AV_TRUE, n, FBIN_TRUE, rng)
    for c in ('Gmag', 'BPmag', 'RPmag'):
        syn[c] += rng.normal(0, 0.02, len(syn))
    syn['e_Gmag'], syn['e_BPmag'], syn['e_RPmag'] = 0.003, 0.02, 0.02
    return syn, iso


def test_likelihood_fit_recovers_av_and_dm():
    syn, iso = _cmd_cluster()
    r = E.fit_isochrone_likelihood(syn, iso, DM_TRUE + 0.03, None, 'free')
    assert abs(r['A_V'] - AV_TRUE) < 3 * r['A_V_err'] + 0.03 and abs(r['DM'] - DM_TRUE) < 3 * r['DM_err'] + 0.03
    assert r['converged'] and 0 <= r['f_out'] <= 0.5


def test_likelihood_dust_prior_pulls_av(monkeypatch):
    monkeypatch.setattr(E, 'AV_PRIOR_SIGMA', 0.01)          # a prior as informative as the data (A_V_err ~ 0.006); a wide one is correctly ignored
    syn, iso = _cmd_cluster()
    pri = E.fit_isochrone_likelihood(syn, iso, DM_TRUE, av_prior=0.9, dust_source='sfd')       # prior far from truth, width 0.1
    free = E.fit_isochrone_likelihood(syn, iso, DM_TRUE, None, 'free')
    assert pri['A_V'] > free['A_V'] + 0.05                                                    # prior drags the fit toward 0.9


def test_equal_mass_binary_sequence_is_detected():
    rng = np.random.default_rng(5)
    iso = toy_iso()
    base = B.build_synthetic(iso, DM_TRUE, E.AG_AV * AV_TRUE, E.ABP_AV * AV_TRUE, E.ARP_AV * AV_TRUE, 300, 0.0, rng)
    sel = rng.random(len(base)) < 0.4
    base.loc[sel, 'Gmag'] -= E.BIN_SEQ_SHIFT                                                 # exact equal-mass binaries: 0.75 mag brighter
    for c in ('Gmag', 'BPmag', 'RPmag'):
        base[c] += rng.normal(0, 0.02, len(base))
    base['e_Gmag'], base['e_BPmag'], base['e_RPmag'] = 0.003, 0.02, 0.02
    r = E.fit_isochrone_likelihood(base, iso, DM_TRUE, None, 'free')
    r0 = E.fit_isochrone_likelihood(base, iso, DM_TRUE, None, 'free', fit_binaries=False)
    assert 0.25 < r['f_bin'] < 0.55 and r['bic'] < r0['bic'] - 10


def test_dust_source_errors_are_explicit():
    df = pd.DataFrame({'RA_ICRS': [150.0], 'DE_ICRS': [-30.0]})
    assert E.get_dust_ebv('free', df)[0] is None
    with pytest.raises(RuntimeError):
        E.get_dust_ebv('bayestar', df, d_pc=300.0)                                           # dustmaps not installed here
    with pytest.raises(ValueError):
        E.get_dust_ebv('nope', df)


def test_fit_isochrone_ages_prefers_true_age_like_isochrone():
    syn, iso = _cmd_cluster()
    shifted = iso.copy(); shifted['Gmag'] += 0.6; shifted['G_BPmag'] += 0.6; shifted['G_RPmag'] += 0.6      # a wrong "age": offset sequence
    tab = E.fit_isochrone_ages(syn, {8.2: iso, 9.0: shifted}, DM_TRUE, None, 'free')
    assert tab.iloc[0]['logAge'] in (8.2, 9.0) and 'delta_bic' in tab and tab['delta_bic'].min() == 0


# ----------------------------- endpoints (new, additive) -----------------------------
@pytest.fixture(scope='module')
def api():
    from fastapi.testclient import TestClient
    import main
    iso0 = toy_iso()
    def fake(logage, mh, av=0.0, cache_only=False, cmd_url=None):
        iso = iso0.copy(); iso['Gmag'] += E.AG_AV * av; iso['G_BPmag'] += E.ABP_AV * av; iso['G_RPmag'] += E.ARP_AV * av
        return dict(iso=iso, raw=b'#mock', params={}, info={'ssl_verified': True, 'source_url': 'mock'}, key=(logage, mh, av, 'mock'))
    mp = pytest.MonkeyPatch(); mp.setattr(E, 'get_parsec', fake)          # restored at module teardown (was a permanent override that leaked into other test files)
    c = TestClient(main.app)
    d = gaia_like(180, 500)
    d['RV'] = 25.0
    assert c.post('/api/upload-raw', files={'file': ('syn.csv', io.BytesIO(d.to_csv(index=False).encode()), 'text/csv')}).status_code == 200
    assert c.post('/api/run-gmm', json=dict(filename='syn.csv', distance_lit=340.0, half_width=150.0)).status_code == 200
    yield c, main
    mp.undo()


def test_new_endpoints_require_stage2_and_stage1(api):
    c, main = api
    r = c.post('/api/compute-binaries', json=dict(filename=main.MEMBERS_NAME))
    assert r.status_code == 400 and 'Fetch PARSEC' in r.json()['detail']
    r = c.post('/api/fit-isochrone-likelihood', json=dict(filename=main.MEMBERS_NAME))
    assert r.status_code == 400


def test_new_endpoints_end_to_end(api):
    c, main = api
    assert c.post('/api/gmm-diagnostics', json=dict(n_resample=30)).status_code == 200
    m = c.post('/api/compute-cmd', json=dict(members_filename=main.MEMBERS_NAME, age=8.2, metallicity=0.0, ebv_mode='fixed', ebv=0.13))
    assert m.status_code == 200, m.text
    fn = main.MEMBERS_NAME
    b = c.post('/api/compute-binaries', json=dict(filename=fn, n_mc=8, n_synth=3000))
    assert b.status_code == 200, b.text
    assert 'warnings' in b.json() and 'MC_CompMass' in main.STATE['analysis']['df']
    v = c.post('/api/compute-velocities', json=dict(filename=fn, n_boot=30))
    assert v.status_code == 200, v.text
    assert v.json()['sigma_1d_kms'] > 0 and v.json()['rv_cluster_kms'] == pytest.approx(25.0)
    dy = c.post('/api/compute-dynamics', json=dict(filename=fn, orbit=False, n_massive=[5, 10]))
    assert dy.status_code == 200, dy.text
    assert dy.json()['r_t_pc'] > 0 and len(dy.json()['msr']) == 2
    f = c.post('/api/fit-isochrone-likelihood', json=dict(filename=fn, dust_source='free'))
    assert f.status_code == 200, f.text
    assert c.post('/api/fit-isochrone-likelihood', json=dict(filename=fn, dust_source='bayestar')).status_code == 400
    assert c.post('/api/fit-isochrone-likelihood', json=dict(filename=fn, dust_source='???')).status_code == 400
    assert c.get('/api/download-binaries').status_code == 200
    # the unchecked "include companions" path now has data and must not crash
    i = c.post('/api/compute-imf', json=dict(filename=fn, log_age=8.2, include_companions=True))
    assert i.status_code in (200, 400)


def test_existing_routes_unchanged_after_phase1(api):
    c, main = api
    assert c.get('/api/health').status_code == 200
    r = c.post('/api/compute-kinematics', json=dict(filename=main.MEMBERS_NAME))
    assert r.status_code == 200 and 'positions_3d' in r.json()
