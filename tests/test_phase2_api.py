"""Phase 2.1 tests: session lock + request validation.
TERMINOLOGY: session = per-visitor state in main._SESSIONS; 400 = client error whose `detail` must name the offending fields.
WORKFLOW: hammer the session store from many threads -> require no exception and the SESSION_MAX bound; post bad bodies -> require 400 + field
names; post the exact payloads the front end sends -> require they still validate."""
# ===================== GLOBAL CONFIGURATION =====================
import os, sys, io, threading
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N_THREADS, N_ITER = 16, 400
# =================================================================
sys.path.insert(0, ROOT)
import pytest
from fastapi.testclient import TestClient
import main
sys.path.insert(0, os.path.join(ROOT, 'tests'))
from test_phase1 import gaia_like

c = TestClient(main.app)


def test_sessions_thread_safe_and_bounded(monkeypatch):
    monkeypatch.setattr(main, 'SESSION_MAX', 6)
    main._SESSIONS.clear()
    errs, over = [], []
    old_iv = sys.getswitchinterval(); sys.setswitchinterval(1e-6)        # force frequent thread switches so an unlocked store would break

    def worker(i):
        try:
            for j in range(N_ITER):
                tok = main._sid_var.set(f'sid{(i * 7 + j) % 40:04d}xxxx')
                try:
                    main.STATE['name'] = f'c{i}-{j}'
                    len(main._SESSIONS) > 6 and over.append(1)               # bound must hold at every instant, not just at the end
                    assert isinstance(main.STATE['name'], str)           # value may be another thread's (shared sid) or a fresh state: only integrity matters
                finally:
                    main._sid_var.reset(tok)
        except Exception as e:                        # any KeyError / RuntimeError('dictionary changed size') lands here
            errs.append(repr(e))

    ts = [threading.Thread(target=worker, args=(i,)) for i in range(N_THREADS)]
    [t.start() for t in ts]; [t.join() for t in ts]
    sys.setswitchinterval(old_iv)
    assert not errs, errs[:3]
    assert not over, f'SESSION_MAX exceeded {len(over)} times'
    assert len(main._SESSIONS) <= 6


@pytest.mark.parametrize('route,body,fields', [
    ('/api/run-gmm', {'filename': 'a.csv', 'cutoff': 1.5}, ['cutoff']),
    ('/api/run-gmm', {'filename': 'a.csv', 'gmm_opt': 'fast', 'k': 1}, ['gmm_opt', 'k']),
    ('/api/run-gmm', {'ruwe_max': 'abc'}, ['filename', 'ruwe_max']),
    ('/api/compute-cmd', {'members_filename': 'm.csv'}, ['age']),
    ('/api/compute-cmd', {'members_filename': 'm.csv', 'age': 8.2, 'ebv_mode': 'magic'}, ['ebv_mode']),
    ('/api/compute-imf', {'filename': 'm.csv', 'method': 'wls', 'bins': 0}, ['method', 'bins']),
    ('/api/fit-isochrone-likelihood', {'filename': 'm.csv', 'dust_source': 'planck'}, ['dust_source']),
    ('/api/compute-velocities', {'filename': 'm.csv', 'plx_mode': 'x', 'n_boot': -3}, ['plx_mode', 'n_boot']),
    ('/api/compute-binaries', {'filename': 'm.csv', 'f_bin': 2}, ['f_bin']),
    ('/api/gmm-diagnostics', {'n_resample': 1}, ['n_resample']),
])
def test_bad_body_is_400_naming_fields(route, body, fields):
    r = c.post(route, json=body)
    assert r.status_code == 400, r.text
    j = r.json()
    assert isinstance(j['detail'], str)                      # front end prints string details verbatim
    for f in fields:
        assert f in j['detail'] and f in j['fields'], (f, j)


def test_legacy_frontend_payloads_still_validate():
    """Exact shapes produced by assets/js/app.js (NaN -> null, extra keys such as dist_mod / av)."""
    d = gaia_like()
    assert c.post('/api/upload-raw', files={'file': ('syn.csv', io.BytesIO(d.to_csv(index=False).encode()), 'text/csv')}).status_code == 200
    body = dict(filename='syn.csv', ruwe_max=1.6, gofal_max=3, epmra_max=0.1, epmde_max=0.1, plx_err_min=5, nper_min=10,
                distance_lit=340, half_width=150, cutoff=0.5, k=None, gmm_opt='converged', plx_zp=0)
    r = c.post('/api/run-gmm', json=body)
    assert r.status_code == 200, r.text
    # nulls for every optional numeric field must be accepted (JSON.stringify(NaN) === 'null')
    r = c.post('/api/run-gmm', json={**body, 'ruwe_max': None, 'plx_zp': None, 'lit_pm': None, 'extra_legacy_key': 1})
    assert r.status_code == 200, r.text
    r = c.post('/api/compute-cmd', json=dict(members_filename='members_latest.csv', ebv_mode='fixed', ebv=None, plx_zp=0, use_cluster_dm=True,
                                              dist_mod=7.6, age=8.2, metallicity=0, av=0.41, iso_extinction='reddened', allow_network=False))
    assert r.status_code == 400 and 'ebv' in r.json()['detail']     # route-level rule (fixed needs ebv) is untouched by validation


def test_unset_keys_stay_absent_so_route_defaults_apply():
    from schemas import RunGmmReq
    assert RunGmmReq(filename='a').model_dump(exclude_unset=True) == {'filename': 'a'}
    assert RunGmmReq(filename='a', k=None).model_dump(exclude_unset=True) == {'filename': 'a', 'k': None}
