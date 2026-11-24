"""Phase 2.2 tests: local PARSEC grid (MOCK analytic isochrones -- CMD is unreachable in the sandbox, so real-PARSEC accuracy is NOT tested here).
TERMINOLOGY: node = exact grid isochrone; off-grid = (logAge, [M/H]) between nodes; truth = the analytic mock evaluated directly.
WORKFLOW: build a small mock grid once -> exact nodes must round-trip -> off-grid points must match truth -> A_V shift must use the MEASURED
coefficients -> get_parsec must use the grid, fall back to CMD when the grid is absent / out of range, and refuse a MOCK grid by default."""
# ===================== GLOBAL CONFIGURATION =====================
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AGE_BOX, MH_BOX = (7.0, 8.6), (-0.2, 0.2)       # small grid keeps the test fast
TOL_CMD_MAG, TOL_MASS_REL = 0.01, 0.01
# =================================================================
sys.path.insert(0, ROOT)
import numpy as np
import pytest
import build_grid as BG
import extinction_isochrone as EXT
import parsec_grid as PG


@pytest.fixture(scope='module')
def grid_path(tmp_path_factory):
    d = tmp_path_factory.mktemp('grid')
    out = str(d / 'g.parquet')
    BG.build(BG.MockFetcher(), out, str(d / 'cache'), resume=False, check=True, age_lo=AGE_BOX[0], age_hi=AGE_BOX[1], mh_lo=MH_BOX[0], mh_hi=MH_BOX[1])
    return out


@pytest.fixture()
def allow_mock(monkeypatch, grid_path):
    monkeypatch.setattr(PG, 'ALLOW_MOCK', True)
    monkeypatch.setattr(PG, 'GRID_PATH', grid_path)
    monkeypatch.setattr(EXT, 'PARSEC_SOURCE', 'auto')
    EXT._PARSEC_CACHE.clear()
    yield
    EXT._PARSEC_CACHE.clear()


def test_mock_grid_is_refused_by_default(grid_path, monkeypatch):
    monkeypatch.setattr(PG, 'ALLOW_MOCK', False)
    with pytest.raises(PG.GridUnavailable, match='MOCK'):
        PG.load_grid(grid_path)


def test_metadata_embedded_and_coefficients_measured(grid_path):
    g = PG.load_grid(grid_path, allow_mock=True)
    c = g.meta['ext_coeff']
    for k, v in BG.MOCK_COEFF.items():
        assert c[k] == pytest.approx(v, abs=1e-6)
    assert g.meta['source'] == 'MOCK' and g.meta['interp_check_ok'] is True
    assert c['linearity_max_dev'] < 1e-6


def test_exact_node_roundtrips(grid_path):
    g = PG.load_grid(grid_path, allow_mock=True)
    iso, info = PG.interpolate_isochrone(g, 7.6, 0.0)
    truth = BG.mock_isochrone(7.6, 0.0)
    assert not info['interpolated'] and len(iso) == len(truth)
    np.testing.assert_allclose(iso['Gmag'], truth['Gmag'], atol=1e-4)           # float32 storage
    np.testing.assert_allclose(iso['Mass'], truth['Mass'], rtol=1e-5)


@pytest.mark.parametrize('age,mh', [(7.625, 0.05), (8.175, -0.15), (7.78, 0.12)])
def test_offgrid_matches_truth_in_cmd_plane_and_mass(grid_path, age, mh):
    g = PG.load_grid(grid_path, allow_mock=True)
    iso, info = PG.interpolate_isochrone(g, age, mh)
    assert info['interpolated'] and len(info['corners']) == 4 and abs(sum(c['weight'] for c in info['corners']) - 1) < 1e-4
    tr = BG.mock_isochrone(age, mh)
    t1, i1 = tr[tr['label'] == 1], iso[iso['label'] == 1]
    ex, ey = EXT._densify((t1['G_BPmag'] - t1['G_RPmag']).to_numpy(), t1['Gmag'].to_numpy())
    dist = [np.min(np.hypot(ex - (b - r), ey - gg)) for b, r, gg in zip(i1['G_BPmag'], i1['G_RPmag'], i1['Gmag'])]
    assert np.max(dist) < TOL_CMD_MAG
    o = np.argsort(t1['Gmag'].to_numpy())
    sel = (i1['Mass'] > 0.3) & (i1['Mass'] < 1.5)
    mt = np.interp(i1['Gmag'][sel], t1['Gmag'].to_numpy()[o], t1['Mass'].to_numpy()[o])
    assert np.max(np.abs(i1['Mass'][sel] / mt - 1)) < TOL_MASS_REL


def test_interpolated_turnoff_mass_between_neighbours(grid_path):
    g = PG.load_grid(grid_path, allow_mock=True)
    lo, _ = PG.interpolate_isochrone(g, 7.6, 0.0)
    hi, _ = PG.interpolate_isochrone(g, 7.7, 0.0)
    mid, _ = PG.interpolate_isochrone(g, 7.65, 0.0)
    to = lambda d: d.loc[d['label'] == 1, 'Mini'].max()
    assert to(hi) < to(mid) < to(lo)                                           # older cluster -> lower turn-off mass


def test_outside_grid_raises(grid_path):
    g = PG.load_grid(grid_path, allow_mock=True)
    for a, m in [(6.6, 0.0), (7.5, 0.9), (9.9, 0.0)]:
        with pytest.raises(PG.GridUnavailable):
            PG.interpolate_isochrone(g, a, m)


def test_av_shift_uses_measured_coefficients_and_passes_selfcheck(allow_mock):
    r0 = EXT.get_parsec(7.65, 0.05, 0.0)
    r1 = EXT.get_parsec(7.65, 0.05, 0.8)
    assert r0['info']['source'] == r1['info']['source'] == 'local_grid'
    d = r1['iso']['Gmag'].to_numpy() - r0['iso']['Gmag'].to_numpy()
    np.testing.assert_allclose(d, BG.MOCK_COEFF['AG_AV'] * 0.8, atol=2e-6)
    dc = (r1['iso']['G_BPmag'] - r1['iso']['G_RPmag']).to_numpy() - (r0['iso']['G_BPmag'] - r0['iso']['G_RPmag']).to_numpy()
    np.testing.assert_allclose(dc, (BG.MOCK_COEFF['ABP_AV'] - BG.MOCK_COEFF['ARP_AV']) * 0.8, atol=2e-6)
    chk = EXT.check_reddening_shift(r0['iso'], r1['iso'], 0.8)
    assert chk['ok'] is True
    # A deliberately wrong coefficient must trip the self-check (guards against a corrupted grid header)
    bad = r1['iso'].copy(); bad['Gmag'] += 1.0
    assert EXT.check_reddening_shift(r0['iso'], bad, 0.8)['ok'] is False


def test_get_parsec_serves_grid_offline_and_raw_is_reparseable(allow_mock):
    r = EXT.get_parsec(8.12, 0.0, 0.3, cache_only=True)                        # cache_only: no network allowed -- the grid needs none
    iso = EXT.select_iso_columns(EXT.parse_parsec(r['raw']))
    assert len(iso) == len(r['iso']) and b'LOCAL GRID' in r['raw']
    EXT.validate_parsec(r['iso'], 8.12, 0.0)


def test_fallback_to_cmd_when_grid_missing_or_out_of_range(monkeypatch, tmp_path, allow_mock):
    called = {}

    def fake_fetch(a, m, v, url=None):
        called['args'] = (a, m, v)
        raise EXT.CMDRequestError('offline in test')
    monkeypatch.setattr(EXT, 'fetch_parsec_raw', fake_fetch)
    with pytest.raises(EXT.CMDRequestError, match='offline'):
        EXT.get_parsec(9.5, 0.0, 0.0)                                           # outside the small test grid -> CMD tried
    assert called['args'][:2] == (9.5, 0.0)
    monkeypatch.setattr(PG, 'GRID_PATH', str(tmp_path / 'nope.parquet'))
    with pytest.raises(EXT.CMDRequestError, match='offline'):
        EXT.get_parsec(7.6, 0.0, 0.0)                                           # grid file absent -> CMD tried
    monkeypatch.setattr(EXT, 'PARSEC_SOURCE', 'grid')
    with pytest.raises(EXT.CMDRequestError, match='cannot serve'):
        EXT.get_parsec(7.61, 0.0, 0.0)                                          # 'grid' mode never falls back


def test_cmd_only_mode_ignores_grid(monkeypatch, allow_mock):
    monkeypatch.setattr(EXT, 'PARSEC_SOURCE', 'cmd')
    monkeypatch.setattr(EXT, 'fetch_parsec_raw', lambda *a, **k: (_ for _ in ()).throw(EXT.CMDRequestError('cmd path used')))
    with pytest.raises(EXT.CMDRequestError, match='cmd path used'):
        EXT.get_parsec(7.6, 0.0, 0.0)


def test_compute_cmd_route_uses_grid(allow_mock, monkeypatch):
    """End to end: /compute-cmd must run with NO network when the grid covers the request."""
    import io
    from fastapi.testclient import TestClient
    import main
    sys.path.insert(0, os.path.join(ROOT, 'tests'))
    from test_phase1 import gaia_like
    c = TestClient(main.app)
    d = gaia_like()
    assert c.post('/api/upload-raw', files={'file': ('syn.csv', io.BytesIO(d.to_csv(index=False).encode()), 'text/csv')}).status_code == 200
    r = c.post('/api/run-gmm', json=dict(filename='syn.csv', distance_lit=340, half_width=150))
    assert r.status_code == 200, r.text
    r = c.post('/api/compute-cmd', json=dict(members_filename='members_latest.csv', ebv_mode='fixed', ebv=0.15, age=7.65, metallicity=0.05,
                                              iso_extinction='reddened', allow_network=True))
    assert r.status_code == 200, r.text
    wf = r.json()['meta']['workflow']['parsec']
    assert wf['state'] == 'done' and wf['source'] == 'local_grid' and wf['interpolated'] is True
    assert wf['check']['ok'] is True and wf['check']['independent'] is False
    assert len(r.json()['isochrone']['mag']) > 50
