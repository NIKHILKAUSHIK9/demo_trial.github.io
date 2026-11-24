"""Phase-0 acceptance tests: synthetic cluster with a known mass-function slope + guards for 0.1-0.6.
TERMINOLOGY: Gamma = slope of log dN/dlogM (Salpeter -1.35); alpha = 1 - Gamma; sigma = bootstrap std of Gamma.
WORKFLOW: draw masses from a truncated Salpeter power law -> run IMF.compute_imf (MLE) -> require |Gamma - true| < 2 sigma;
then check label rule, G_lim -> mass, approx_ML exclusion, completeness curve, and that stage 1 applies no UI E(B-V)."""
# ===================== GLOBAL CONFIGURATION =====================
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRUE_GAMMA, N_STARS, M_LO, M_HI, SEED = -1.35, 500, 0.2, 5.0, 7
N_REPEAT, MIN_COVER = 40, 0.80          # coverage test: fraction of repeats whose |dGamma| < 2 sigma (nominal 95 %)
# =================================================================
sys.path.insert(0, ROOT)
import numpy as np
import pandas as pd
import pytest
import gmm as G
import imf as I


def draw_masses(n, gamma, lo, hi, rng):
    """Inverse-CDF sampling of dN/dM ~ M^-(1-gamma) on [lo, hi]."""
    b = gamma                                   # 1 - alpha = gamma
    u = rng.random(n)
    return (lo ** b + u * (hi ** b - lo ** b)) ** (1 / b)


def test_mle_recovers_gamma_within_2_sigma():
    m = draw_masses(N_STARS, TRUE_GAMMA, M_LO, M_HI, np.random.default_rng(SEED))
    r = I.compute_imf(m, fit_min_logm=np.log10(M_LO))
    g, s = r.meta['fit_gamma'], r.meta['fit_gamma_err']
    assert r.meta['method'] == 'mle' and 0 < s < 0.2
    assert abs(g - TRUE_GAMMA) < 2 * s, (g, s)


def test_mle_bootstrap_coverage():
    hits = 0
    for i in range(N_REPEAT):
        m = draw_masses(N_STARS, TRUE_GAMMA, M_LO, M_HI, np.random.default_rng(100 + i))
        r = I.compute_imf(m, fit_min_logm=np.log10(M_LO), n_boot=200)
        hits += abs(r.meta['fit_gamma'] - TRUE_GAMMA) < 2 * r.meta['fit_gamma_err']
    assert hits / N_REPEAT >= MIN_COVER


def test_binned_method_kept_and_mle_is_seeded():
    m = draw_masses(N_STARS, TRUE_GAMMA, M_LO, M_HI, np.random.default_rng(SEED))
    a = I.compute_imf(m, fit_min_logm=np.log10(M_LO), method='binned')
    b = I.compute_imf(m, fit_min_logm=np.log10(M_LO))
    c = I.compute_imf(m, fit_min_logm=np.log10(M_LO))
    assert a.meta['method'] == 'binned' and a.meta['mle'] is None
    assert b.meta['fit_gamma'] == c.meta['fit_gamma'] and b.meta['fit_gamma_err'] == c.meta['fit_gamma_err']
    with pytest.raises(ValueError):
        I.compute_imf(m, method='nope')


def test_cash_used_for_sparse_bins_and_pearson_for_dense():
    sparse = I.compute_imf(draw_masses(60, TRUE_GAMMA, M_LO, M_HI, np.random.default_rng(1)), fit_min_logm=np.log10(M_LO))
    assert (sparse.chi2['stat'] == 'Cash').all()
    dense = I.compute_imf(draw_masses(20000, TRUE_GAMMA, M_LO, M_HI, np.random.default_rng(1)), fit_min_logm=np.log10(M_LO), n_bins=8)
    assert (dense.chi2['stat'] == 'chi2').all()
    assert dense.chi2.set_index('Model')['red_chi2'].idxmin() == 'Salpeter'       # integrated-bin model must prefer the true shape


def test_mf_label_rule():
    assert I.mf_label(7.0) == I.MF_LABEL_INITIAL and I.mf_label(7.5) == I.MF_LABEL_INITIAL
    assert I.mf_label(7.51) == I.MF_LABEL_PRESENT and I.mf_label(8.2) == I.MF_LABEL_PRESENT
    assert I.mf_label(None) == I.MF_LABEL_INITIAL


def test_g_limit_to_mass_monotone_and_clipped():
    iso = pd.DataFrame({'Gmag': np.linspace(0, 12, 50), 'Mass': np.linspace(3.0, 0.2, 50), 'label': 1})
    a = I.g_limit_to_mass(14.0, 7.5, 0.3, iso)        # M_G = 6.2
    b = I.g_limit_to_mass(16.0, 7.5, 0.3, iso)        # fainter limit -> lower mass
    assert b['m_min'] < a['m_min'] and not a['clipped']
    assert I.g_limit_to_mass(30.0, 7.5, 0.3, iso)['clipped']


def test_mass_sample_excludes_approx_ml_by_default():
    df = pd.DataFrame({'MC_Mass': [1.0, 0.5, 0.3], 'MC_Mass_src': ['isochrone', 'approx_ML', 'user']})
    assert list(I.mass_sample(df)[0]) == [1.0, 0.3] and I.mass_sample(df)[1] == 1
    assert len(I.mass_sample(df, include_approx_ml=True)[0]) == 3


def test_completeness_curve_recovers_known_limit():
    rng = np.random.default_rng(3)
    g = rng.uniform(8, 18, 20000)
    keep = rng.random(g.size) < np.where(g < 14, 0.99, 0.3)        # retention collapses beyond G = 14
    raw, cln = pd.DataFrame({'Gmag': g}), pd.DataFrame({'Gmag': g[keep]})
    kw = dict(frac=0.9, mode='differential', norm='absolute')
    c = G.completeness_curve(raw, cln, **kw)
    assert 13.5 <= c['g_lim'] <= 14.5, c['g_lim']
    assert np.isnan(G.completeness_curve(raw, cln.iloc[:0], **kw)['g_lim'])
    cc = G.completeness_curve(raw, cln, frac=0.9)             # default: cumulative + plateau-normalised
    assert 13.5 <= cc['g_lim'] <= 15.0, cc['g_lim']


def test_plateau_norm_ignores_flat_astrometric_loss():
    """A magnitude-independent 20 % loss (RUWE / binaries) must not move G_lim; only the magnitude-dependent drop does."""
    rng = np.random.default_rng(4)
    g = rng.uniform(8, 18, 30000)
    keep = rng.random(g.size) < 0.8 * np.where(g < 15, 1.0, 0.2)
    raw, cln = pd.DataFrame({'Gmag': g}), pd.DataFrame({'Gmag': g[keep]})
    assert np.isnan(G.completeness_curve(raw, cln, norm='absolute')['g_lim'])        # 0.8 < 0.9 everywhere -> never "complete"
    assert 14.5 <= G.completeness_curve(raw, cln, norm='plateau')['g_lim'] <= 16.0


def test_distance_window_excludes_nonpositive_and_outside():
    df = pd.DataFrame({'Plx': [3.0, 2.9, 10.0, -1.0, 0.0, 3.1]})        # d = 333, 345, 100, n/a, n/a, 323 pc
    w = G.distance_window(df, 340.0, 20.0)
    assert sorted(w['Plx']) == [2.9, 3.0, 3.1]


@pytest.fixture(scope='module')
def client():
    from fastapi.testclient import TestClient
    import main
    return TestClient(main.app), main


def _synthetic_gaia(n_cl=120, n_fld=400, seed=5):
    """Cluster (tight PM/parallax blob at 340 pc) + broad field; quality columns clean."""
    rng = np.random.default_rng(seed)
    def block(n, pm, spm, plx, splx):
        return pd.DataFrame({'pmRA': rng.normal(pm[0], spm, n), 'pmDE': rng.normal(pm[1], spm, n), 'Plx': rng.normal(plx, splx, n)})
    d = pd.concat([block(n_cl, (5, -3), 0.15, 2.94, 0.05), block(n_fld, (0, 0), 6, 3.0, 0.9)], ignore_index=True)
    n = len(d)
    d['RA_ICRS'], d['DE_ICRS'] = rng.normal(150, 1, n), rng.normal(-30, 1, n)
    d['e_Plx'], d['e_pmRA'], d['e_pmDE'] = 0.02, 0.03, 0.03
    d['Gmag'] = rng.uniform(9, 17, n)
    d['BPmag'] = d['Gmag'] + 0.4 + rng.normal(0, 0.1, n); d['RPmag'] = d['Gmag'] - 0.4 + rng.normal(0, 0.1, n)
    d['BP-RP'] = d['BPmag'] - d['RPmag']
    d['RUWE'], d['gofAL'], d['Nper'] = 1.0, 1.0, 20
    return d[d['Plx'] > 0.5].reset_index(drop=True)


def test_stage1_applies_no_ui_ebv_and_stamps_when_explicit(client):
    c, main = client
    import io
    buf = io.BytesIO(_synthetic_gaia().to_csv(index=False).encode())
    assert c.post('/api/upload-raw', files={'file': ('syn.csv', buf, 'text/csv')}).status_code == 200
    r = c.post('/api/run-gmm', json=dict(filename='syn.csv', distance_lit=340.0, half_width=150.0, ebv_ui=0.3))
    assert r.status_code == 200, r.text
    assert 'M_G_corr' not in main.STATE['members'].columns or 'E_BV_used' in main.STATE['members'].columns
    r = c.post('/api/run-gmm', json=dict(filename='syn.csv', distance_lit=340.0, half_width=150.0))     # no explicit value
    m = main.STATE['members']
    assert 'M_G_corr' not in m.columns and 'bp_rp_corr' not in m.columns and 'E_BV_used' not in m.columns
    assert 'MC_Mass' not in m.columns and 'Star_Type' in m.columns
    r = c.post('/api/run-gmm', json=dict(filename='syn.csv', distance_lit=340.0, half_width=150.0, ebv_ui=0.3))
    m = main.STATE['members']
    assert (m['E_BV_used'] == 0.3).all() and 'M_G_corr' in m.columns


def test_derive_columns_mass_quality(client):
    _, main = client
    df = pd.DataFrame({'M_G_corr': [3.0, 4.0], 'bp_rp_corr': [0.5, 0.8]})
    iso = pd.DataFrame({'Gmag': [2.0, 3.0, 4.0, 4.5, 5.0], 'Mass': [2.0, 1.5, 1.2, 1.0, 0.8], 'label': 1})
    d = main.derive_columns(df, iso=iso)
    assert list(d['Mass_quality']) == ['ok', 'ok'] and set(d['MC_Mass_src']) == {'isochrone'}
    d2 = main.derive_columns(df)
    assert set(d2['Mass_quality']) == {'approx'}
