"""Phase 3 (credibility) tests: citation metadata, placeholders, engine banner wiring, examples/verify.py.
TERMINOLOGY: CFF = CITATION.cff (GitHub "Cite this repository"); Zenodo JSON = .zenodo.json; footer BibTeX = the constants in assets/js/footer.js.
WORKFLOW: the software title / authors / version and the paper title / authors must be identical in every place they are written; no TODO or dead
link may ship; verify.py must pass, fail, skip and refuse exactly as documented. The verify.py numbers below are SYNTHETIC test inputs (cluster name
TESTCLUSTER), not astrophysical results."""
# ===================== GLOBAL CONFIGURATION =====================
import importlib.util, json, os, re, shutil
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLACEHOLDER_SCAN = ('README.md', 'index.html', 'CITATION.cff', '.zenodo.json', 'examples/README.md', 'examples/verify.py', 'assets/js/app.js',
                    'assets/js/panels.js', 'assets/js/footer.js', 'assets/css/styles.css')
# =================================================================
import numpy as np
import pandas as pd
import pytest
import yaml


def _read(rel):
    return open(os.path.join(ROOT, rel), encoding='utf-8').read()


CFF = yaml.safe_load(_read('CITATION.cff'))
ZEN = json.loads(_read('.zenodo.json'))
FOOTER = _read('assets/js/footer.js')


def _bib_field(block, name):
    m = re.search(rf'{name}\s*=\s*\{{(.*)\}},?\s*$', block, re.M)
    assert m, name
    return m.group(1)


def _footer_block(const):
    m = re.search(rf'const {const} = `(.*?)`;', FOOTER, re.S)
    assert m, const
    return m.group(1)


def _bib_authors(block):
    out = []
    for a in _bib_field(block, 'author').split(' and '):
        a = a.strip()
        out.append(a[1:-1] if a.startswith('{') else a)                       # '{Monu}' -> 'Monu'
    return out


def _cff_names(authors):
    return [f"{a['family-names']}, {a['given-names']}" if 'given-names' in a else a['family-names'] for a in authors]   # Monu: single name


def test_cff_structure_and_unpublished_paper():
    for k in ('cff-version', 'message', 'title', 'authors', 'version', 'license', 'repository-code'):
        assert k in CFF, k
    assert CFF['type'] == 'software'
    paper = next(r for r in CFF['references'] if r['type'] == 'unpublished')
    assert paper['status'] == 'in-preparation'
    for forbidden in ('journal', 'volume', 'issue', 'start', 'end', 'doi', 'pages', 'identifiers'):
        assert forbidden not in paper, f'an unpublished paper must not carry {forbidden}'
    assert 'preferred-citation' not in CFF, 'the unpublished paper must not become the preferred citation'


def test_cff_references_cite_the_verified_papers():
    refs = {r['title']: r for r in CFF['references']}
    a = refs['Revisiting the mass of open clusters with Gaia data']
    assert (a['volume'], a['start'], a['year']) == (525, 2315, 2023)
    k = refs['Astrophysical parameters of Galactic open clusters']
    assert (k['volume'], k['start'], k['end'], k['year']) == (438, 1163, 1173, 2005)


def test_version_is_identical_everywhere():
    init = _read('starnest/__init__.py')
    v = re.search(r"__version__\s*=\s*'([^']+)'", init).group(1)
    assert CFF['version'] == ZEN['version'] == v
    assert re.search(rf'version\s*=\s*"{re.escape(v)}"', _read('pyproject.toml'))
    assert re.search(rf'version\s*=\s*\{{{re.escape(v)}\}}', _footer_block('BIBTEX_REPO'))


def test_software_title_and_authors_identical_in_cff_zenodo_footer():
    repo = _footer_block('BIBTEX_REPO')
    title = _bib_field(repo, 'title')
    assert title.startswith('{') and title.endswith('}')
    title = title[1:-1]
    assert CFF['title'] == ZEN['title'] == title
    cff_names = _cff_names(CFF['authors'])
    assert cff_names == [c['name'] for c in ZEN['creators']] == _bib_authors(repo)


def test_paper_title_and_author_order_identical_in_cff_and_footer():
    paper = next(r for r in CFF['references'] if r['type'] == 'unpublished')
    bib = _footer_block('BIBTEX_PAPER')
    assert bib.startswith('@unpublished{'), 'an unpublished paper is @unpublished, not @article'
    assert 'journal' not in bib and 'doi' not in bib.lower()
    assert _bib_field(bib, 'title') == paper['title']
    assert _bib_authors(bib) == _cff_names(paper['authors'])
    assert str(paper['year']) in _bib_field(bib, 'year')
    bibitem = _footer_block('BIBITEM_PAPER')
    assert paper['title'] in bibitem and 'in preparation' in bibitem


def test_zenodo_json_fields():
    assert ZEN['upload_type'] == 'software' and ZEN['access_right'] == 'open' and ZEN['license'] == 'MIT'
    assert ZEN['creators'] and all(c['name'] for c in ZEN['creators'])
    assert 'publication_date' not in ZEN, 'Zenodo sets the date at release time'


def test_no_placeholders_ship():
    bad = re.compile(r'TODO|FIXME')
    for rel in PLACEHOLDER_SCAN:
        assert not bad.search(_read(rel)), f'{rel} still contains a TODO / FIXME'
    assert "PAPER_DOI_URL     = '#'" not in FOOTER and 'journal      = {In preparation}' not in FOOTER
    assert re.search(r"const PAPER_DOI_URL\s*= '';", FOOTER), 'empty = the Paper / DOI button is hidden'
    assert 'coming soon' not in FOOTER


def test_paper_button_is_hidden_not_dead():
    assert "if (!PAPER_DOI_URL)" in FOOTER and "a.style.display = 'none'" in FOOTER
    assert 'data-paper-link' in _read('index.html')


def test_engine_banner_wiring():
    html, app = _read('index.html'), _read('assets/js/app.js')
    assert html.count('id="engine-approx-banner"') == 1
    assert 'display:none' in re.search(r'<div id="engine-approx-banner"[^>]*>', html).group(0), 'hidden until the probe fails'
    assert html.index('id="engine-approx-banner"') < html.index('<main id="pages-container"')
    assert 'approximate, not the published method' in html
    assert "backendAlive === false" in app and "getElementById('engine-approx-banner')" in app


def test_readme_attribution_and_references():
    r = _read('README.md')
    assert 'not a strict replication' in r and 'inspired by' in r
    assert 'MNRAS 525, 2315' in r
    assert re.search(r'^20\. Kharchenko', r, re.M) and re.search(r'^21\. Cummings', r, re.M)
    assert 'COCD' in r and 'in preparation' in r and 'unpublished working paper' in r
    assert 'examples/verify.py' in r and 'Relation to Almeida et al. (2023)' in r


def test_readme_difference_table_matches_code():
    """The 'what differs' claims are about the code: if the code changes, the README must too."""
    b = _read('starnest/binaries.py')
    assert 'F_BIN = 0.5' in b and 'N_SYNTH = 10000' in b
    assert not re.search(r'rng\.normal|photometric err', b, re.I), 'README says no photometric errors are injected into the synthetic cluster'
    assert 'IFMR_PIECES = ((0.87, 2.8' in b and '8.2, 0.107, 0.471))' in b          # README: 0.87 to 8.2 Msun
    api = _read('starnest/api.py')
    assert "f_bin = float(p.get('f_bin', BIN.F_BIN))" in api and 'BIN.cluster_mass(tab, mle, m_to, f_bin=f_bin)' in api   # README: assumed F_BIN, not measured


@pytest.mark.skipif(shutil.which('cffconvert') is None, reason='cffconvert not installed')
def test_cff_validates_with_cffconvert():
    import subprocess
    r = subprocess.run(['cffconvert', '--validate', '-i', os.path.join(ROOT, 'CITATION.cff')], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


# ---------------------------------------------------------------- examples/verify.py (SYNTHETIC inputs)
def _load_verify():
    spec = importlib.util.spec_from_file_location('verify_example', os.path.join(ROOT, 'examples', 'verify.py'))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


V = _load_verify()


def _run_dir(tmp_path, name='TESTCLUSTER', n=100, d=300.0, ebv=0.05, pmra=10.0, pmde=-5.0, plx=3.3, status='ok', members=True, iso='local_grid:PARSEC'):
    pd.DataFrame([dict(cluster=name, status=status, N=n, d=d, **{'E(B-V)': ebv}, isochrone_source=iso)]).to_csv(tmp_path / 'summary.csv', index=False)
    if members:
        (tmp_path / re.sub(V.SLUG_PATTERN, '_', name)).mkdir(exist_ok=True)
        rng = np.random.default_rng(1)
        pd.DataFrame(dict(pmRA=pmra + rng.normal(0, 0.05, 101), pmDE=pmde + rng.normal(0, 0.05, 101), Plx=plx + rng.normal(0, 0.01, 101))
                     ).to_csv(tmp_path / re.sub(V.SLUG_PATTERN, '_', name) / V.MEMBERS_FILENAME, index=False)
    return str(tmp_path / 'summary.csv')


def _lit(tmp_path, **kw):
    row = dict(cluster='TESTCLUSTER', source='synthetic test table', d_pc=300.0, plx_mas=3.3, plx_err_mas=0.01, pmra_masyr=10.0, pmde_masyr=-5.0, pm_err_masyr=0.05,
               ebv=0.05, n_members=100)
    row.update(kw)
    p = tmp_path / 'lit.csv'; pd.DataFrame([row]).to_csv(p, index=False); return str(p)


def _status(tmp_path, **kw):
    s, l = V.load_tables(str(tmp_path / 'summary.csv'), str(tmp_path / 'lit.csv'))
    return V.run_checks(s, l, str(tmp_path)).set_index('quantity')['status'].to_dict()


def test_verify_all_pass_and_exit_zero(tmp_path):
    s = _run_dir(tmp_path); l = _lit(tmp_path)
    assert set(_status(tmp_path).values()) == {'PASS'}
    assert V.main([s, l]) == 0


def test_verify_distance_outside_tolerance_fails(tmp_path):
    s = _run_dir(tmp_path, d=300.0 * (1 + V.D_TOL_FRAC + 0.02)); l = _lit(tmp_path)
    st = _status(tmp_path)
    assert st['d_pc'] == 'FAIL' and st['N'] == 'PASS'
    assert V.main([s, l]) == 1


def test_verify_boundary_is_inclusive():
    assert V._judge(1.0, 0.0, 1.0) == 'PASS' and V._judge(1.0 + 1e-9, 0.0, 1.0) == 'FAIL'
    assert V._judge(float('nan'), 0.0, 1.0) == 'MISSING' and V._judge(1.0, float('nan'), 1.0) == 'SKIP'


def test_verify_ebv_outside_tolerance_fails(tmp_path):
    _run_dir(tmp_path, ebv=0.05 + V.EBV_TOL_MAG + 0.02); _lit(tmp_path)
    assert _status(tmp_path)['E(B-V)'] == 'FAIL'


def test_verify_member_count_uses_the_wide_tolerance(tmp_path):
    _run_dir(tmp_path, n=int(100 * (1 + V.N_TOL_FRAC) - 1)); _lit(tmp_path)
    assert _status(tmp_path)['N'] == 'PASS'
    _run_dir(tmp_path, n=int(100 * (1 + V.N_TOL_FRAC) + 2))
    assert _status(tmp_path)['N'] == 'FAIL'


def test_verify_pm_tolerance_scales_with_literature_error_with_a_floor(tmp_path):
    _run_dir(tmp_path, pmra=10.0 + 0.30); _lit(tmp_path, pm_err_masyr=0.05)                 # 3 sigma = 0.15 -> floor 0.15 < 0.30 -> FAIL
    assert _status(tmp_path)['pmRA'] == 'FAIL'
    _lit(tmp_path, pm_err_masyr=0.20)                                                        # 3 sigma = 0.60 > 0.30 -> PASS
    assert _status(tmp_path)['pmRA'] == 'PASS'


def test_verify_missing_pipeline_value_is_a_failure_and_blank_literature_is_skipped(tmp_path):
    s = _run_dir(tmp_path, d=float('nan'), members=False); l = _lit(tmp_path, ebv=np.nan)
    st = _status(tmp_path)
    assert st['d_pc'] == 'MISSING' and st['pmRA'] == 'MISSING' and st['E(B-V)'] == 'SKIP'
    assert V.main([s, l]) == 1


def test_verify_all_skipped_is_not_a_pass_of_anything(tmp_path):
    _run_dir(tmp_path); _lit(tmp_path, d_pc=np.nan, plx_mas=np.nan, pmra_masyr=np.nan, pmde_masyr=np.nan, ebv=np.nan, n_members=np.nan)
    assert set(_status(tmp_path).values()) == {'SKIP'}


def test_verify_partial_run_cannot_pass(tmp_path):
    _run_dir(tmp_path, status='partial'); _lit(tmp_path)
    assert 'PASS' not in set(_status(tmp_path).values())


def test_verify_cluster_name_matching_ignores_case_space_and_underscore(tmp_path):
    _run_dir(tmp_path, name='UPK_545'); _lit(tmp_path, cluster='upk 545')
    assert set(_status(tmp_path).values()) == {'PASS'}


def test_verify_cluster_absent_from_summary_is_missing(tmp_path):
    _run_dir(tmp_path, name='OTHER'); _lit(tmp_path, cluster='TESTCLUSTER')
    assert set(_status(tmp_path).values()) == {'MISSING'}


@pytest.mark.parametrize('bad', ['no_source_column', 'blank_source', 'dup'])
def test_verify_refuses_unsourced_or_duplicate_literature(tmp_path, bad):
    _run_dir(tmp_path)
    df = pd.DataFrame([dict(cluster='TESTCLUSTER', source='x', d_pc=300.0)])
    if bad == 'no_source_column':
        df = df.drop(columns='source')
    elif bad == 'blank_source':
        df['source'] = '  '
    else:
        df = pd.concat([df, df.assign(cluster='testcluster')])
    df.to_csv(tmp_path / 'lit.csv', index=False)
    with pytest.raises(SystemExit) as e:
        V.load_tables(str(tmp_path / 'summary.csv'), str(tmp_path / 'lit.csv'))
    assert 'verify:' in str(e.value)


def test_verify_missing_file_exits_with_message(tmp_path):
    with pytest.raises(SystemExit) as e:
        V.load_tables(str(tmp_path / 'nope.csv'), str(tmp_path / 'nope2.csv'))
    assert 'not found' in str(e.value)


def test_verify_writes_csv_without_index_and_warns_on_mock(tmp_path, capsys):
    s = _run_dir(tmp_path, iso='local_grid:MOCK'); l = _lit(tmp_path)
    out = tmp_path / 'checks.csv'
    assert V.main([s, l, '--csv', str(out)]) == 0
    assert list(pd.read_csv(out).columns) == V.CHECK_COLUMNS
    assert 'MOCK' in capsys.readouterr().err


def test_literature_template_header_matches_verify():
    head = _read('examples/literature_table.template.csv').strip().split(',')
    assert head[:2] == list(V.LIT_REQUIRED)
    assert set(head) >= {'d_pc', 'plx_mas', 'plx_err_mas', 'pmra_masyr', 'pmde_masyr', 'pm_err_masyr', 'ebv', 'n_members'}
    assert len(_read('examples/literature_table.template.csv').strip().splitlines()) == 1, 'the template must hold no invented values'
