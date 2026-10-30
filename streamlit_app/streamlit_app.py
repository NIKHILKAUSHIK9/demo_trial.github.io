# ===== GLOBAL CONFIGURATION =====
PAGE_TITLE = "StarNest - Open Star Cluster Analysis"
FRAME_HEIGHT = 1000                 # height of the embedded page in px
PHOTO_EXTS = (".jpg", ".jpeg", ".png", ".webp")
BACKEND_URL = ""                    # optional external backend (Render / HF), e.g. "https://x.onrender.com"; or secret/env STARNEST_BACKEND_URL
PROBE_TIMEOUT_MS = 20000            # first call wakes the Streamlit<->page bridge; give it time
# ================================
# Streamlit Community Cloud front end for StarNest. ONE file to deploy; it picks the engine by itself:
#   1. BACKEND_URL / STARNEST_BACKEND_URL set  -> the page calls that external FastAPI server (CORS, X-Session-Id).
#   2. else, if main.py (FastAPI backend) imports -> BRIDGE mode: Streamlit runs the backend IN-PROCESS. The page is a
#      Streamlit component; each API call goes page -> setComponentValue -> this script -> FastAPI app (TestClient) ->
#      component args -> page. No second port or server is needed (Community Cloud exposes only one port).
#   3. else -> browser-only (in-browser GMM / kinematics / IMF; no IRSA / PARSEC).
import base64, json, os, pathlib, sys, tempfile
import streamlit as st
import streamlit.components.v1 as components

ROOT = pathlib.Path(__file__).resolve().parents[1]
st.set_page_config(page_title=PAGE_TITLE, layout="wide", initial_sidebar_state="collapsed")
st.markdown("<style>header,footer{visibility:hidden;height:0} .block-container{padding:0!important;max-width:100%!important}</style>", unsafe_allow_html=True)


def _secret() -> str:
    try:
        return str(st.secrets.get("STARNEST_BACKEND_URL", ""))
    except Exception:
        return ""


@st.cache_resource(show_spinner=False)
def load_backend():
    """FastAPI app of main.py wrapped in a TestClient, or None if it cannot be imported."""
    try:
        sys.path.insert(0, str(ROOT))
        os.chdir(ROOT)
        import main as sn_main
        from fastapi.testclient import TestClient
        return TestClient(sn_main.app, base_url="http://starnest.local", raise_server_exceptions=False)
    except Exception as e:  # missing dependency etc.
        print("StarNest backend unavailable:", repr(e), file=sys.stderr)
        return None


def inline_page(flags: str, extra_js: str = "") -> str:
    html = (ROOT / "index.html").read_text(encoding="utf-8")
    css = (ROOT / "assets/css/styles.css").read_text(encoding="utf-8")
    html = html.replace('<link rel="stylesheet" href="assets/css/styles.css">', f"<style>{css}</style>")
    for name in ("tailwind-config", "app"):
        js = (ROOT / f"assets/js/{name}.js").read_text(encoding="utf-8").replace("</script>", "<\\/script>")
        tag = f'<script src="assets/js/{name}.js"></script>'
        pre = f"<script>{extra_js}</script>" if (name == "app" and extra_js) else ""
        html = html.replace(tag, pre + "<script>__JS__</script>".replace("__JS__", js))
    photos = {}
    for f in sorted((ROOT / "images" / "team").glob("*")):
        if f.suffix.lower() in PHOTO_EXTS:
            mime = "image/jpeg" if f.suffix.lower() in (".jpg", ".jpeg") else f"image/{f.suffix.lower()[1:]}"
            photos[f.stem] = f"data:{mime};base64," + base64.b64encode(f.read_bytes()).decode()
    head = f"<script>{flags}window.STARNEST_PHOTOS={json.dumps(photos)};</script>"
    return html.replace('<script src="https://cdn.tailwindcss.com"></script>', head + '<script src="https://cdn.tailwindcss.com"></script>', 1)


@st.cache_resource(show_spinner=False)
def bridge_component():
    """Component whose front end is the whole StarNest page + bridge.js (written once to a temp folder)."""
    d = pathlib.Path(tempfile.mkdtemp(prefix="starnest_bridge_"))
    flags = ("window.STARNEST_API_BASE='http://starnest.bridge/api';window.STARNEST_BRIDGE=true;"
             f"window.STARNEST_PROBE_TIMEOUT_MS={PROBE_TIMEOUT_MS};window.STARNEST_FRAME_HEIGHT={FRAME_HEIGHT};")
    bridge = (pathlib.Path(__file__).parent / "bridge.js").read_text(encoding="utf-8").replace("</script>", "<\\/script>")
    (d / "index.html").write_text(inline_page(flags, bridge), encoding="utf-8")
    return components.declare_component("starnest_bridge", path=str(d))


def serve(req: dict) -> dict:
    """Run one page request on the in-process FastAPI app."""
    client = load_backend()
    hdr = {k: v for k, v in (req.get("headers") or {}).items() if k.lower() not in ("host", "content-length", "connection")}
    body = base64.b64decode(req.get("body") or "")
    try:
        r = client.request(req.get("method", "GET"), "/api" + req["path"], headers=hdr, content=body)
        keep = {k: v for k, v in r.headers.items() if k.lower() in ("content-type", "content-disposition")}
        return {"id": req["id"], "status": r.status_code, "headers": keep, "body": base64.b64encode(r.content).decode()}
    except Exception as e:
        msg = json.dumps({"detail": f"backend error: {e}"}).encode()
        return {"id": req["id"], "status": 500, "headers": {"content-type": "application/json"}, "body": base64.b64encode(msg).decode()}


url = (BACKEND_URL or os.environ.get("STARNEST_BACKEND_URL") or _secret()).rstrip("/")
if url:                                                    # mode 1: external backend
    components.html(inline_page(f"window.STARNEST_API_BASE={json.dumps(url + '/api')};window.STARNEST_PROBE_TIMEOUT_MS=60000;"), height=FRAME_HEIGHT, scrolling=True)
elif load_backend() is not None:                           # mode 2: in-process backend through the bridge
    comp = bridge_component()
    resp = st.session_state.get("sn_resp")
    val = comp(resp=resp, key="starnest", default=None)
    key = f'{val["id"]}:{val.get("part", 0)}' if val else None
    if val and key != st.session_state.get("sn_last_id"):
        st.session_state["sn_last_id"] = key
        buf = st.session_state.setdefault("sn_buf", {})
        buf.setdefault(val["id"], []).append(val.get("body") or "")
        if val.get("part", 0) + 1 < val.get("parts", 1):         # more parts to come: acknowledge this one
            st.session_state["sn_resp"] = {"id": val["id"], "ack": val.get("part", 0)}
        else:
            val["body"] = "".join(buf.pop(val["id"]))
            st.session_state["sn_resp"] = serve(val)
        st.rerun()
else:                                                      # mode 3: browser only
    components.html(inline_page("window.STARNEST_NO_BACKEND=true;"), height=FRAME_HEIGHT, scrolling=True)
