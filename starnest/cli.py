"""
cli.py -- batch runner: `python -m starnest run config.yaml` (Phase 2.4).

TERMINOLOGY
  Config      : YAML file: output_dir, steps, defaults (one dict per step) and clusters (name, csv, optional per-step overrides).
  Step        : one stage of the pipeline = one API route of the web app (gmm, cmd, imf, binaries, velocities, dynamics, isofit).
                The CLI calls the SAME FastAPI app in-process (TestClient), so batch numbers are identical to what the web page shows.
  Manifest    : run_manifest.json = resolved parameters + package versions + SHA-256 of every input file (+ local PARSEC grid), for reproducibility.
  Summary     : summary.csv, one row per cluster (N, d, E(B-V), logAge, Gamma, sigma_Gamma, MSR, v_out, sigma_1D + status columns).

WORKFLOW (per cluster)
  1) upload the raw Gaia CSV; 2) gmm -> cmd (extinction + isochrone) -> imf -> binaries -> velocities -> dynamics -> isofit, each with
  defaults overridden by the cluster's own keys; a failed step is recorded and later steps that need it are skipped;
  3) write members / extinction / isochrone / IMF / binaries CSVs (index=False), figures (PDF + PNG), the one-page report;
  4) after all clusters: summary.csv and run_manifest.json. Exit code 0 = all ok, 1 = at least one cluster/step failed, 2 = bad config.

EXAMPLE config.yaml
  output_dir: starnest_output
  defaults:
    gmm:  {ruwe_max: 1.6, cutoff: 0.5}
    cmd:  {ebv_mode: center, allow_network: true, iso_extinction: reddened}
    imf:  {method: mle}
  clusters:
    - name: UPK_545
      csv: examples/UPK_545_raw.csv
      gmm: {distance_lit: 325, half_width: 150}
      cmd: {age: 8.2, metallicity: 0.0}
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import sys
import uuid
from datetime import datetime, timezone

import pandas as pd

# ===================== GLOBAL CONFIGURATION =====================
DEFAULT_OUT_DIR = 'starnest_output'
STEPS = ('gmm', 'cmd', 'imf', 'binaries', 'velocities', 'dynamics', 'isofit')          # pipeline order
STEP_ENDPOINT = dict(gmm='/run-gmm', cmd='/compute-cmd', imf='/compute-imf', binaries='/compute-binaries', velocities='/compute-velocities',
                     dynamics='/compute-dynamics', isofit='/fit-isochrone-likelihood')
STEP_NEEDS = dict(gmm=(), cmd=('gmm',), imf=('cmd',), binaries=('cmd',), velocities=('cmd',), dynamics=('cmd',), isofit=('cmd',))   # a step is skipped if these failed
MSR_REPORT_N = 10                                    # summary 'MSR' = Lambda_MSR at this many heaviest stars (nearest available N)
MEMBERS_FILE = 'members_latest.csv'                  # name the web app gives the stage-1 member list
WRITE_FIGURES = True
WRITE_REPORT = True
SUMMARY_COLUMNS = ['cluster', 'status', 'N', 'd', 'E(B-V)', 'logAge', 'Gamma', 'sigma_Gamma', 'MSR', 'v_out', 'sigma_1D',
                   'v_out_err', 'sigma_1D_err', 'MSR_N', 'MF_kind', 'IMF_method', 'isochrone_source', 'n_warnings', 'failed_steps']
MANIFEST_PACKAGES = ('numpy', 'pandas', 'scipy', 'scikit-learn', 'matplotlib', 'astropy', 'fastapi', 'pydantic', 'pyarrow', 'galpy', 'dustmaps', 'PyYAML')
DOWNLOADS = {'members.csv': '/download-members', 'extinction_corrected.csv': '/download-extinction', 'isochrone.csv': '/download-isochrone',
             'parsec_isochrone.txt': '/download-parsec', 'imf_bins.csv': '/download-imf-bins', 'imf_chi2.csv': '/download-imf-chi2',
             'binaries_masses.csv': '/download-binaries'}
API = '/api'
# =================================================================


class ConfigError(Exception):
    pass


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load_config(path: str) -> dict:
    """Parse + validate the YAML. Paths in `csv` are resolved relative to the config file. Raises ConfigError with a readable message."""
    import yaml
    try:
        cfg = yaml.safe_load(open(path, encoding='utf-8'))
    except (OSError, yaml.YAMLError) as e:
        raise ConfigError(f'cannot read {path}: {e}')
    if not isinstance(cfg, dict) or not isinstance(cfg.get('clusters'), list) or not cfg['clusters']:
        raise ConfigError("config needs a non-empty 'clusters:' list")
    unknown = set(cfg) - {'output_dir', 'steps', 'defaults', 'clusters'}
    if unknown:
        raise ConfigError(f'unknown top-level keys: {sorted(unknown)}')
    steps = cfg.get('steps') or list(STEPS)
    bad = [s for s in steps if s not in STEPS]
    if bad:
        raise ConfigError(f'unknown steps {bad}; allowed: {list(STEPS)}')
    cfg['steps'] = [s for s in STEPS if s in steps]
    defaults = cfg.get('defaults') or {}
    if set(defaults) - set(STEPS):
        raise ConfigError(f'defaults may only contain step dicts {list(STEPS)}; got {sorted(set(defaults) - set(STEPS))}')
    base, seen = os.path.dirname(os.path.abspath(path)), set()
    for c in cfg['clusters']:
        if not isinstance(c, dict) or not c.get('name') or not c.get('csv'):
            raise ConfigError("every cluster needs 'name' and 'csv'")
        if set(c) - {'name', 'csv'} - set(STEPS):
            raise ConfigError(f"cluster {c['name']}: unknown keys {sorted(set(c) - {'name', 'csv'} - set(STEPS))}")
        c['slug'] = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(c['name']))
        if c['slug'] in seen:
            raise ConfigError(f"duplicate cluster name {c['name']!r}")
        seen.add(c['slug'])
        c['csv_path'] = c['csv'] if os.path.isabs(c['csv']) else os.path.join(base, c['csv'])
        if not os.path.isfile(c['csv_path']):
            raise ConfigError(f"cluster {c['name']}: csv not found: {c['csv_path']}")
        merged = {s: {**(defaults.get(s) or {}), **(c.get(s) or {})} for s in STEPS}
        if 'cmd' in cfg['steps'] and merged['cmd'].get('age') is None:
            raise ConfigError(f"cluster {c['name']}: cmd.age (log10 age/yr) is required for the isochrone step")
        c['params'] = merged
    cfg['output_dir_resolved'] = os.path.abspath(cfg.get('output_dir') or DEFAULT_OUT_DIR)
    return cfg


def _client():
    try:
        from fastapi.testclient import TestClient
    except Exception as e:                                           # httpx missing
        raise SystemExit(f"the CLI needs the 'api' extra (fastapi + httpx): pip install 'starnest[api]'  ({e})")
    from . import api
    return TestClient(api.app, base_url='http://starnest.local', raise_server_exceptions=False), api


def _msr_at(msr: list, n: int):
    ok = [x for x in msr or [] if x.get('Lambda') is not None and math.isfinite(x['Lambda'])]
    if not ok:
        return None, None
    x = min(ok, key=lambda r: abs(r['N'] - n))
    return float(x['Lambda']), int(x['N'])


def run_cluster(client, api, c: dict, steps: list, out_root: str) -> dict:
    """Run one cluster; returns dict(row=summary row, files=[...], errors={step: msg}, results={step: compact reply})."""
    sid = uuid.uuid4().hex
    hdr = {'X-Session-Id': sid}
    out_dir = os.path.join(out_root, c['slug']); os.makedirs(out_dir, exist_ok=True)
    res, errors, files = {}, {}, []
    fname = os.path.basename(c['csv_path'])

    def post(endpoint, body):
        r = client.post(API + endpoint, json=body, headers=hdr)
        if r.status_code != 200:
            try:
                detail = r.json().get('detail', r.text)
            except Exception:
                detail = r.text
            raise RuntimeError(f'HTTP {r.status_code}: {detail}')
        return r.json()

    up = client.post(API + '/upload-raw', files={'file': (fname, open(c['csv_path'], 'rb'), 'text/csv')}, headers=hdr)
    if up.status_code != 200:
        errors['upload'] = f'HTTP {up.status_code}: {up.text[:300]}'
        steps = []
    for s in steps:
        if any(n in errors or n not in res for n in STEP_NEEDS[s] if n in steps or n in errors):
            errors[s] = 'skipped: ' + ', '.join(n for n in STEP_NEEDS[s] if n in errors) + ' failed'
            continue
        p = dict(c['params'][s])
        try:
            if s == 'gmm':
                p['filename'] = fname
            elif s == 'cmd':
                p['members_filename'] = MEMBERS_FILE
            elif s in ('velocities', 'dynamics') and (res.get('binaries') or {}).get('mass') and 'mass' not in p:
                p['filename'] = MEMBERS_FILE; p['mass'] = res['binaries']['mass']['M_T']
            if s not in ('gmm', 'cmd'):
                p.setdefault('filename', MEMBERS_FILE)
            if s == 'imf' and 'log_age' not in p:
                p['log_age'] = c['params']['cmd'].get('age')
            res[s] = post(STEP_ENDPOINT[s], p)
            if s == 'cmd':
                pz = ((res[s].get('meta') or {}).get('workflow') or {}).get('parsec') or {}
                if pz.get('state') != 'done':                     # 200 reply but no isochrone (network off / CMD down / outside the grid): every mass-based step would fail
                    errors[s] = f"isochrone not obtained (state={pz.get('state')}): {(res[s].get('meta') or {}).get('note') or 'enable cmd.allow_network or build the local grid'}"[:500]
        except Exception as e:
            errors[s] = str(e)[:500]
    st = api._SESSIONS.get(sid, {}).get('s', {})
    for fn, ep in DOWNLOADS.items():                                   # CSV exports (index=False inside the API)
        r = client.get(API + ep, headers=hdr)
        if r.status_code == 200 and len(r.content) > 0:
            path = os.path.join(out_dir, fn); open(path, 'wb').write(r.content); files.append(path)

    g = res.get('gmm', {}).get('stats', {}); meta = (res.get('cmd') or {}).get('meta', {}); imfr = (res.get('imf') or {}).get('imf', {}); imf = imfr.get('meta', {})
    vel, dyn = res.get('velocities') or {}, res.get('dynamics') or {}
    msr, msr_n = _msr_at(dyn.get('msr'), MSR_REPORT_N)
    parsec = (meta.get('workflow') or {}).get('parsec') or {}
    warnings = list(meta.get('note') and [meta['note']] or []) + [w for k in ('binaries',) for w in ((res.get(k) or {}).get('warnings') or [])]
    warnings += list(vel.get('notes') or []) + list(dyn.get('notes') or [])
    if str(parsec.get('grid_source', '')).upper() == 'MOCK':
        warnings.insert(0, 'ISOCHRONE GRID IS A MOCK TEST GRID, NOT PARSEC: all masses / A_V / DM are invalid')
    if (imf.get('mle') or {}).get('at_bound'):
        warnings.append('IMF slope sits at the edge of its allowed range: Gamma is not constrained by this sample')
    if (res.get('isofit') or {}).get('converged') is False:
        warnings.append('likelihood isochrone fit did not converge')
    if imf.get('limit_note'):
        warnings.append(imf['limit_note'])
    if (res.get('dynamics') or {}).get('mass_source', '').endswith('(lower limit)'):
        warnings.append('Dynamics mass is the sum of observed stars (lower limit): tidal radius and relaxation time are underestimated')
    warnings += [f'step {k} failed: {v}' for k, v in errors.items()]
    age = c['params']['cmd'].get('age')
    row = {'cluster': c['name'], 'status': 'ok' if not errors else ('partial' if res.get('gmm') else 'failed'),
           'N': g.get('members'), 'd': meta.get('d_pc'), 'E(B-V)': meta.get('E_BV_median'), 'logAge': age,
           'Gamma': imfr.get('slope'), 'sigma_Gamma': imfr.get('slope_err'), 'MSR': msr, 'v_out': vel.get('v_out_kms'), 'sigma_1D': vel.get('sigma_1d_kms'),
           'v_out_err': vel.get('v_out_err'), 'sigma_1D_err': vel.get('sigma_1d_err'), 'MSR_N': msr_n, 'MF_kind': (res.get('imf') or {}).get('imf', {}).get('mf_kind'),
           'IMF_method': imf.get('method'), 'isochrone_source': parsec.get('source') if not parsec.get('grid_source') else f"{parsec.get('source')}:{parsec.get('grid_source')}",
           'n_warnings': len(warnings), 'failed_steps': ','.join(errors)}
    if WRITE_FIGURES or WRITE_REPORT:
        from . import report
        bundle = dict(name=c['name'], gmm_res=st.get('gmm_res'), analysis_df=(st.get('analysis') or {}).get('df'), iso=st.get('iso'), reddened=meta.get('plot_space') == 'reddened',
                      imf_res=st.get('imf_res'), mf_kind=(res.get('imf') or {}).get('imf', {}).get('mf_kind', 'Mass Function'), msr=dyn.get('msr'))
        try:
            figs, skipped = report.collect_figures(bundle)
            if WRITE_FIGURES:
                files += report.export_all(figs, out_dir, c['slug'])
            if WRITE_REPORT:
                prov = f"isochrone: {parsec.get('service', 'n/a')} ({row['isochrone_source']}); logAge {age}; input {fname}"
                files.append(report.one_page_report(os.path.join(out_dir, f"{c['slug']}_report.pdf"), c['name'], row, figs, skipped, warnings, prov))
        except Exception as e:
            errors['report'] = f'figure/report export failed: {e}'
            row['status'] = 'partial' if res.get('gmm') else 'failed'
    compact = {k: {kk: vv for kk, vv in v.items() if not isinstance(vv, (list, dict))} for k, v in res.items()}
    json.dump(dict(cluster=c['name'], errors=errors, warnings=warnings, scalars=compact), open(os.path.join(out_dir, 'result.json'), 'w'), indent=1, default=str)
    api._SESSIONS.pop(sid, None)
    return dict(row=row, files=files, errors=errors)


def _versions() -> dict:
    from importlib import metadata
    out = {}
    for p in MANIFEST_PACKAGES:
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            out[p] = None
    from . import __version__
    out['starnest'] = __version__
    return out


def build_manifest(cfg: dict, results: list, started: str, argv: list) -> dict:
    from . import extinction_isochrone as EXT
    from . import parsec_grid as PG
    grid = None
    if PG.grid_available():
        grid = dict(path=PG.GRID_PATH, sha256=sha256_file(PG.GRID_PATH), parsec_source_setting=EXT.PARSEC_SOURCE)
        try:
            gm = PG.load_grid(allow_mock=True).meta
            grid.update(source=gm.get('source'), built_utc=gm.get('built_utc'), logAge=gm.get('logAge'), MH=gm.get('MH'), ext_coeff=gm.get('ext_coeff'))
        except Exception as e:
            grid['error'] = str(e)
    return dict(started_utc=started, finished_utc=datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), command=argv, python=sys.version.split()[0],
                platform=platform.platform(), packages=_versions(), steps=cfg['steps'], parsec_grid=grid, parsec_source_setting=EXT.PARSEC_SOURCE,
                clusters=[dict(name=c['name'], input_csv=c['csv_path'], input_sha256=sha256_file(c['csv_path']), input_bytes=os.path.getsize(c['csv_path']),
                               parameters=c['params'], status=r['row']['status'], errors=r['errors'], outputs=[os.path.relpath(f, cfg['output_dir_resolved']) for f in r['files']])
                          for c, r in zip(cfg['clusters'], results)])


def run(config_path: str, out_dir: str | None = None, only: list | None = None, argv: list | None = None) -> int:
    started = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    try:
        cfg = load_config(config_path)
    except ConfigError as e:
        print(f'config error: {e}', file=sys.stderr)
        return 2
    if out_dir:
        cfg['output_dir_resolved'] = os.path.abspath(out_dir)
    if only:
        miss = set(only) - {c['name'] for c in cfg['clusters']}
        if miss:
            print(f'config error: --only names not in config: {sorted(miss)}', file=sys.stderr)
            return 2
        cfg['clusters'] = [c for c in cfg['clusters'] if c['name'] in only]
    out = cfg['output_dir_resolved']; os.makedirs(out, exist_ok=True)
    client, api = _client()
    results = []
    for c in cfg['clusters']:
        print(f"[{c['name']}] running {', '.join(cfg['steps'])}", file=sys.stderr)
        results.append(run_cluster(client, api, c, cfg['steps'], out))
        print(f"[{c['name']}] {results[-1]['row']['status']}" + (f" ({', '.join(results[-1]['errors'])})" if results[-1]['errors'] else ''), file=sys.stderr)
    summary = pd.DataFrame([r['row'] for r in results]).reindex(columns=SUMMARY_COLUMNS)
    summary.to_csv(os.path.join(out, 'summary.csv'), index=False)
    json.dump(build_manifest(cfg, results, started, argv or sys.argv), open(os.path.join(out, 'run_manifest.json'), 'w'), indent=1, default=str)
    print(f'wrote {os.path.join(out, "summary.csv")} and run_manifest.json', file=sys.stderr)
    return 0 if all(r['row']['status'] == 'ok' for r in results) else 1


def main(argv: list | None = None) -> int:
    ap = argparse.ArgumentParser(prog='python -m starnest', description='StarNest batch pipeline')
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run', help='run the pipeline for every cluster in a YAML config')
    r.add_argument('config')
    r.add_argument('--out', default=None, help=f'output directory (default: config output_dir or ./{DEFAULT_OUT_DIR})')
    r.add_argument('--only', nargs='*', default=None, help='run only these cluster names')
    a = ap.parse_args(argv)
    return run(a.config, a.out, a.only, ['python', '-m', 'starnest'] + (argv if argv is not None else sys.argv[1:]))
