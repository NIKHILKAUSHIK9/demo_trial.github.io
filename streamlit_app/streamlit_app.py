# ===== GLOBAL CONFIGURATION =====
PAGE_TITLE = "StarNest - Open Star Cluster Analysis"
FRAME_HEIGHT = 1000                 # height of the embedded page in px
PHOTO_EXTS = (".jpg", ".jpeg", ".png", ".webp")
# ================================
# What this file does (Streamlit Community Cloud wrapper, BROWSER-ONLY mode):
#   1. read index.html and inline assets/css + assets/js (an embedded page cannot load relative files)
#   2. inline the team photos from images/team as data URIs
#   3. tell the page there is no backend (it then uses its in-browser GMM / kinematics / IMF engine)
#   4. show it with components.html
# The IRSA + PARSEC step needs the FastAPI backend (see docs/DEPLOYMENT.md); it is not available here.
import base64, json, pathlib
import streamlit as st
import streamlit.components.v1 as components

ROOT = pathlib.Path(__file__).resolve().parents[1]
st.set_page_config(page_title=PAGE_TITLE, layout="wide", initial_sidebar_state="collapsed")
st.markdown("<style>header,footer{visibility:hidden;height:0} .block-container{padding:0!important;max-width:100%!important}</style>", unsafe_allow_html=True)


@st.cache_data(show_spinner=False)
def build_page() -> str:
    html = (ROOT / "index.html").read_text(encoding="utf-8")
    css = (ROOT / "assets/css/styles.css").read_text(encoding="utf-8")
    html = html.replace('<link rel="stylesheet" href="assets/css/styles.css">', f"<style>{css}</style>")
    for name in ("tailwind-config", "app"):
        js = (ROOT / f"assets/js/{name}.js").read_text(encoding="utf-8").replace("</script>", "<\\/script>")
        tag = f'<script src="assets/js/{name}.js"></script>'
        html = html.replace(tag, "<script>__JS__</script>".replace("__JS__", js))
    photos = {}
    for f in sorted((ROOT / "images" / "team").glob("*")):
        if f.suffix.lower() in PHOTO_EXTS:
            mime = "image/jpeg" if f.suffix.lower() in (".jpg", ".jpeg") else f"image/{f.suffix.lower()[1:]}"
            photos[f.stem] = f"data:{mime};base64," + base64.b64encode(f.read_bytes()).decode()
    flags = f"<script>window.STARNEST_NO_BACKEND=true;window.STARNEST_PHOTOS={json.dumps(photos)};</script>"
    return html.replace("<script src=\"https://cdn.tailwindcss.com\"></script>", flags + "<script src=\"https://cdn.tailwindcss.com\"></script>", 1)


components.html(build_page(), height=FRAME_HEIGHT, scrolling=True)
