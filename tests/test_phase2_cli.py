"""Phase 2.4 / 2.5 tests: batch CLI, manifest, figures, one-page report (synthetic cluster + MOCK grid; no network).
TERMINOLOGY: summary.csv = one row per cluster; manifest = run_manifest.json (parameters, versions, input SHA-256).
WORKFLOW: build a small mock grid -> write config + synthetic raw CSV -> run the CLI -> check every promised output exists and is consistent with
the web-app numbers; then check config errors (exit 2), a failing cluster (exit 1, others still processed) and the skip-on-failure logic."""
# ===================== GLOBAL CONFIGURATION =====================
import hashlib, json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AGE_BOX, MH_BOX = (7.0, 8.6), (-0.2, 0.2)
CLUSTER_AGE = 7.65
# =================================================================
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tests'))
import pandas as pd
import pytest
import build_grid as BG
import extinction_isochrone as EXT
import parsec_grid as PG
from starnest import cli
from test_phase1 import gaia_like


@pytest.fixture(scope='module')
def env(tmp_path_factory):
    d = tmp_path_factory.mktemp('cli')
    grid = str(d / 'g.parquet')
    BG.build(BG.MockFetcher(), grid, str(d / 'cache'), resume=False, check=False, age_lo=AGE_BOX[0], age_hi=AGE_BOX[1], mh_lo=MH_BOX[0], mh_hi=MH_BOX[1])
    csv = d / 'syn.csv'
    gaia_like(160, 400).to_csv(csv, index=False)
    mp = pytest.MonkeyPatch()
    mp.setattr(PG, 'ALLOW_MOCK', True); mp.setattr(PG, 'GRID_PATH', grid); mp.setattr(EXT, 'PARSEC_SOURCE', 'auto')
    EXT._PARSEC_CACHE.clear()
    yield d, csv
    mp.undo(); EXT._PARSEC_CACHE.clear()


def _cfg(d, csv, extra_clusters=(), steps=None):
    import yaml
    c = dict(output_dir=str(d / 'out'), steps=steps or ['gmm', 'cmd', 'imf', 'binaries', 'velocities', 'dynamics'],
             defaults=dict(gmm=dict(distance_lit=340, half_width=150), cmd=dict(ebv_mode='fixed', ebv=0.15, allow_network=False, metallicity=0.05),
                           binaries=dict(n_mc=5, n_synth=2000), velocities=dict(n_boot=100), dynamics=dict(n_random=30, orbit=False)),
             clusters=[dict(name='SYN 1', csv=str(csv), cmd=dict(age=CLUSTER_AGE)), *extra_clusters])
    p = d / 'cfg.yaml'; p.write_text(yaml.safe_dump(c)); return str(p)


def test_full_run_outputs_summary_manifest_figures(env):
    d, csv = env
    rc = cli.run(_cfg(d, csv), argv=['python', '-m', 'starnest', 'run', 'cfg.yaml'])
    assert rc == 0
    out = d / 'out'
    s = pd.read_csv(out / 'summary.csv')
    assert list(s.columns[:11]) == ['cluster', 'status', 'N', 'd', 'E(B-V)', 'logAge', 'Gamma', 'sigma_Gamma', 'MSR', 'v_out', 'sigma_1D']
    assert not s.columns.str.startswith('Unnamed').any()                                   # index=False
    r = s.iloc[0]
    assert r['status'] == 'ok' and r['N'] > 100 and 330 < r['d'] < 350 and r['E(B-V)'] == pytest.approx(0.15) and r['logAge'] == CLUSTER_AGE
    assert pd.notna(r['Gamma']) and pd.notna(r['MSR']) and pd.notna(r['sigma_1D']) and pd.notna(r['v_out'])
    assert r['isochrone_source'] == 'local_grid:MOCK' and r['n_warnings'] >= 1              # mock grid must be loud in the summary
    cd = out / 'SYN_1'
    for f in ['members.csv', 'extinction_corrected.csv', 'isochrone.csv', 'imf_bins.csv', 'binaries_masses.csv', 'result.json', 'SYN_1_report.pdf',
              'SYN_1_pm.pdf', 'SYN_1_pm.png', 'SYN_1_cmd.pdf', 'SYN_1_imf.png', 'SYN_1_msr.pdf']:
        assert (cd / f).is_file() and (cd / f).stat().st_size > 500, f
    assert (cd / 'SYN_1_report.pdf').read_bytes()[:4] == b'%PDF'
    m = pd.read_csv(cd / 'members.csv')
    for col in ['RA_ICRS', 'DE_ICRS', 'Plx', 'pmRA', 'pmDE', 'PMemb']:                      # pipeline column names preserved
        assert col in m.columns
    man = json.load(open(out / 'run_manifest.json'))
    assert man['clusters'][0]['input_sha256'] == hashlib.sha256(open(csv, 'rb').read()).hexdigest()
    assert man['packages']['numpy'] and man['packages']['starnest'] == '1.0.0' and man['parsec_grid']['source'] == 'MOCK' and len(man['parsec_grid']['sha256']) == 64
    assert man['clusters'][0]['parameters']['cmd']['age'] == CLUSTER_AGE and 'started_utc' in man


def test_summary_matches_web_api_numbers(env):
    """The CLI must reproduce what the page shows: same app, same defaults."""
    import io
    from fastapi.testclient import TestClient
    import main
    d, csv = env
    c = TestClient(main.app)
    c.post('/api/upload-raw', files={'file': ('syn.csv', io.BytesIO(open(csv, 'rb').read()), 'text/csv')})
    g = c.post('/api/run-gmm', json=dict(filename='syn.csv', distance_lit=340, half_width=150)).json()
    c.post('/api/compute-cmd', json=dict(members_filename='members_latest.csv', ebv_mode='fixed', ebv=0.15, age=CLUSTER_AGE, metallicity=0.05, allow_network=False))
    im = c.post('/api/compute-imf', json=dict(filename='members_latest.csv', log_age=CLUSTER_AGE)).json()['imf']
    r = pd.read_csv(d / 'out' / 'summary.csv').iloc[0]
    assert r['N'] == g['stats']['members']
    assert r['Gamma'] == pytest.approx(im['slope'], rel=1e-9) and r['sigma_Gamma'] == pytest.approx(im['slope_err'], rel=1e-9)   # seeded bootstrap -> identical


def test_config_errors_exit_2(env, tmp_path):
    d, csv = env
    bad = tmp_path / 'b.yaml'
    for text, frag in [('clusters: []', 'non-empty'), (f'clusters:\n - {{name: a, csv: {csv}}}\n', 'cmd.age'), ('clusters:\n - {name: a, csv: /nope.csv}\n', 'not found'),
                       (f'steps: [gmm, magic]\nclusters:\n - {{name: a, csv: {csv}}}\n', 'unknown steps'), (f'colour: red\nclusters:\n - {{name: a, csv: {csv}}}\n', 'unknown top-level')]:
        bad.write_text(text)
        with pytest.raises(cli.ConfigError, match=frag):
            cli.load_config(str(bad))
        assert cli.run(str(bad), out_dir=str(tmp_path / 'o')) == 2


def test_failed_step_is_recorded_and_downstream_skipped(env):
    d, csv = env
    # cluster 2: age outside the mock grid and no network -> cmd fails; imf/binaries/... must be skipped, gmm outputs kept, exit code 1
    cfg = _cfg(d, csv, extra_clusters=[dict(name='NO ISO', csv=str(csv), cmd=dict(age=9.9))])
    rc = cli.run(cfg, out_dir=str(d / 'out2'))
    assert rc == 1
    s = pd.read_csv(d / 'out2' / 'summary.csv').set_index('cluster')
    assert s.loc['SYN 1', 'status'] == 'ok' and s.loc['NO ISO', 'status'] == 'partial'
    assert pd.notna(s.loc['NO ISO', 'N']) and pd.isna(s.loc['NO ISO', 'Gamma']) and pd.isna(s.loc['NO ISO', 'sigma_1D'])
    assert s.loc['NO ISO', 'failed_steps'].split(',')[0] == 'cmd' and 'imf' in s.loc['NO ISO', 'failed_steps']
    res = json.load(open(d / 'out2' / 'NO_ISO' / 'result.json'))
    assert 'skipped' in res['errors']['imf'] and (d / 'out2' / 'NO_ISO' / 'NO_ISO_report.pdf').is_file()     # report still produced, with placeholders


def test_only_filter(env):
    d, csv = env
    assert cli.run(_cfg(d, csv), out_dir=str(d / 'out3'), only=['nope']) == 2
