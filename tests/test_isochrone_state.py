"""Isochrone-state tests: a refresh that cannot fetch an isochrone must not silently turn good masses into the L(M) fallback.
TERMINOLOGY
  user Apply      : /compute-cmd with allow_network=True (the "Apply Extinction & Fetch PARSEC" button; may fetch).
  auto refresh    : /compute-cmd with allow_network=False (any other input change; may only use what is cached).
  stale isochrone : the mass table was built from isochrone (logAge, [M/H]) but the inputs now differ and the new isochrone is not fetched.
  approx_ML       : MC_Mass_src of the crude L(M) fallback (no isochrone); excluded from the mass-function fit by default.
WORKFLOW
  1) a fake PARSEC source that honours cache_only exactly like extinction_isochrone.get_parsec (cached -> served, else NetworkNotAllowed);
  2) replay the sequence that produced "Only 0 valid masses; need >= 20" in the browser: Apply -> change age / [M/H] -> run the IMF;
  3) require the IMF (and the Phase-1 steps) to say the table is STALE and to work again after the next successful Apply;
  4) static checks that the three front-end triggers in app.js are wired the way the fix needs."""
# ===================== GLOBAL CONFIGURATION =====================
import io
import os
import re
import sys
SID = 'sn-test-iso-state-0001'                    # own visitor id: this file never shares state with another test module
AGE, MH, EBV = 7.5, 0.0, 0.13                      # the fetched isochrone
AGE_NEW, MH_NEW, EBV_NEW = 7.6, 0.1, 0.20          # changed inputs whose isochrone is NOT fetched
# =================================================================
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tests'))
import pytest
from fastapi.testclient import TestClient
import main
import starnest.extinction_isochrone as E
from test_phase1 import gaia_like, toy_iso

APP_JS = open(os.path.join(ROOT, 'assets', 'js', 'app.js'), encoding='utf-8').read()


class FakeParsec:
    """Stand-in for extinction_isochrone.get_parsec with the same cache / cache_only contract; `fail_fetch` simulates CMD being down."""
    def __init__(self):
        self.cache, self.fail_fetch, self.iso0 = {}, False, toy_iso()

    def __call__(self, logage, mh, av=0.0, cache_only=False, cmd_url=None):
        key = (round(logage, 4), round(mh, 4), round(av, 4))
        if key not in self.cache:
            if cache_only:
                raise E.NetworkNotAllowed('PARSEC file not cached')
            if self.fail_fetch:
                raise E.CMDRequestError('CMD unreachable (test)')
            iso = self.iso0.copy()
            iso['Gmag'] += E.AG_AV * av; iso['G_BPmag'] += E.ABP_AV * av; iso['G_RPmag'] += E.ARP_AV * av
            self.cache[key] = dict(iso=iso, raw=b'#mock', params={}, info={'ssl_verified': True, 'source_url': 'mock'}, key=key)
        return self.cache[key]


@pytest.fixture()
def api(monkeypatch):
    fake = FakeParsec()
    monkeypatch.setattr(E, 'get_parsec', fake)
    c = TestClient(main.app, headers={'X-Session-Id': SID})
    d = gaia_like(180, 500); d['RV'] = 25.0
    assert c.post('/api/upload-raw', files={'file': ('syn.csv', io.BytesIO(d.to_csv(index=False).encode()), 'text/csv')}).status_code == 200
    assert c.post('/api/run-gmm', json=dict(filename='syn.csv', distance_lit=340.0, half_width=150.0)).status_code == 200
    return c, fake


def cmd(c, user, **kw):
    body = dict(members_filename=main.MEMBERS_NAME, age=AGE, metallicity=MH, ebv_mode='fixed', ebv=EBV, iso_extinction='reddened', allow_network=user)
    body.update(kw)
    r = c.post('/api/compute-cmd', json=body)
    assert r.status_code == 200, r.text
    return r.json()['meta']


def imf(c):
    return c.post('/api/compute-imf', json=dict(filename=main.MEMBERS_NAME, log_age=AGE))


def test_imf_before_any_isochrone_names_the_fix(api):
    c, _ = api
    r = imf(c)
    assert r.status_code == 400 and 'approx_ML' in r.json()['detail'] and 'fetch the PARSEC isochrone' in r.json()['detail']
    assert 'stale' not in r.json()['detail']                                   # no isochrone was ever fetched: not "stale", just missing


def test_apply_then_imf_works_and_same_inputs_refresh_keeps_it(api):
    c, _ = api
    assert cmd(c, True)['workflow']['parsec']['state'] == 'done'
    assert imf(c).status_code == 200
    assert cmd(c, False)['workflow']['parsec']['state'] == 'done'              # cached: an automatic refresh with the same inputs is harmless
    assert imf(c).status_code == 200


@pytest.mark.parametrize('change', [dict(age=AGE_NEW), dict(metallicity=MH_NEW)])
def test_auto_refresh_with_unfetched_isochrone_keeps_masses_and_marks_them_stale(api, change):
    """The reported failure: Apply works, then changing age or [M/H] used to replace the mass table by the approx_ML one."""
    c, _ = api
    cmd(c, True)
    assert imf(c).status_code == 200
    meta = cmd(c, False, **change)
    assert meta['workflow']['parsec']['state'] == 'pending' and 'stale' in meta['note']
    r = imf(c)
    assert r.status_code == 400 and 'stale' in r.json()['detail'] and "Apply Extinction & Fetch PARSEC" in r.json()['detail']
    assert f"logAge {AGE}" in r.json()['detail']                               # says which isochrone the masses came from
    b = c.post('/api/compute-binaries', json=dict(filename=main.MEMBERS_NAME, n_mc=2, n_synth=500))
    assert b.status_code == 400 and 'stale' in b.json()['detail']              # Phase-1 steps that need the isochrone refuse too
    assert cmd(c, True, **change)['workflow']['parsec']['state'] == 'done'     # the next real Apply repairs it
    assert imf(c).status_code == 200


def test_returning_to_the_fetched_inputs_clears_stale(api):
    c, _ = api
    cmd(c, True)
    cmd(c, False, age=AGE_NEW)
    assert imf(c).status_code == 400
    cmd(c, False)                                                              # inputs back to the fetched ones: cache hit
    assert imf(c).status_code == 200


def test_changing_ebv_alone_is_not_stale(api):
    c, _ = api
    cmd(c, True)
    cmd(c, False, ebv=EBV_NEW)                                                 # the A_V = 0 isochrone is still cached; masses do not depend on E(B-V)
    assert imf(c).status_code == 200


def test_failed_user_fetch_keeps_the_old_table_but_marks_it_stale(api):
    c, fake = api
    cmd(c, True)
    fake.fail_fetch = True
    meta = cmd(c, True, age=AGE_NEW)
    assert meta['workflow']['parsec']['state'] == 'failed' and 'stale' in meta['note']
    r = imf(c)
    assert r.status_code == 400 and 'stale' in r.json()['detail']
    fake.fail_fetch = False
    assert cmd(c, True, age=AGE_NEW)['workflow']['parsec']['state'] == 'done'
    assert imf(c).status_code == 200


def test_new_gmm_run_resets_the_stale_flag(api):
    c, _ = api
    cmd(c, True)
    cmd(c, False, age=AGE_NEW)
    assert c.post('/api/run-gmm', json=dict(filename='syn.csv', distance_lit=340.0, half_width=150.0)).status_code == 200
    r = imf(c)
    assert r.status_code == 400 and 'stale' not in r.json()['detail']          # a fresh member list has no isochrone at all (plain "missing")


# ---------------------------- front end wiring (static) ----------------------------
def test_gmm_chain_does_not_run_imf_before_an_isochrone_exists():
    chain = APP_JS[APP_JS.index("window._lastCmdMeta = null;"):]
    chain = chain[:chain.index("await computeIMFPlots();") + 40]
    assert "backendAlive === true" in chain and "state === 'done'" in chain and "showToast('Mass function needs the PARSEC isochrone" in chain
    assert not re.search(r"await computeKinematicsPlots\(\);\s*await computeIMFPlots\(\);", APP_JS), 'IMF must not be unconditional after the GMM chain'


def test_successful_apply_reruns_the_imf():
    body = APP_JS[APP_JS.index('async function computeCMDPlot'):]
    body = body[:body.index('function initKinematicPlots')]
    assert re.search(r"if \(userClick && data\.meta && data\.meta\.workflow && data\.meta\.workflow\.parsec && data\.meta\.workflow\.parsec\.state === 'done' "
                     r"&& getTaskFilename\('imf'\)\) \{\s*computeIMFPlots\(\{ silent: true \}\);", body)


def test_age_change_waits_for_the_cmd_refresh_before_the_imf():
    m = re.search(r"el\.addEventListener\('change', async \(\) => \{\s*if \(getTaskFilename\('isochrone'\)\) \{\s*await computeCMDPlot\(\);.*?"
                  r"if \(id === 'input-age' && getTaskFilename\('imf'\)\) computeIMFPlots\(\);", APP_JS, re.S)
    assert m, 'input-age handler must await computeCMDPlot() before computeIMFPlots()'
