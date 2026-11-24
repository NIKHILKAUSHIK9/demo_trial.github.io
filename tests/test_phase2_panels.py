"""Phase 2 panel wiring tests (static: no browser needed). The behavioural check ran in Playwright against a live backend (see the Phase 2 report).
TERMINOLOGY: mount = <div id="panel-..."> placeholder in index.html that panels.js fills. WORKFLOW: every mount id used by panels.js must exist
once in index.html, panels.js must load after app.js, the Streamlit inliner must include it, and the file must parse."""
# ===================== GLOBAL CONFIGURATION =====================
import os, re, shutil, subprocess
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# =================================================================
import pytest

html = open(os.path.join(ROOT, 'index.html'), encoding='utf-8').read()
js = open(os.path.join(ROOT, 'assets/js/panels.js'), encoding='utf-8').read()


def test_every_mount_exists_exactly_once():
    ids = re.findall(r"'(panel-[a-z-]+)': mount", js)
    assert len(ids) == 6
    for i in ids:
        assert html.count(f'id="{i}"') == 1, i


def test_script_order_and_streamlit_inliner():
    assert html.index('assets/js/app.js') < html.index('assets/js/panels.js') < html.index('assets/js/footer.js')
    assert '"panels"' in open(os.path.join(ROOT, 'streamlit_app/streamlit_app.py'), encoding='utf-8').read()


def test_app_js_hooks_present():
    app = open(os.path.join(ROOT, 'assets/js/app.js'), encoding='utf-8').read()
    assert 'starnestClusterRuleParams' in app and 'window._lastCmdMeta = data.meta' in app


@pytest.mark.skipif(shutil.which('node') is None, reason='node not installed')
def test_panels_js_parses():
    r = subprocess.run(['node', '--check', os.path.join(ROOT, 'assets/js/panels.js')], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_backend_only_panels_never_call_the_browser_engine():
    # every panel request goes through call() which probes first; a direct apiCall elsewhere would fall into localBackendEngine('Unhandled endpoint')
    assert len(re.findall(r'\bapiCall\(', js)) == 1
