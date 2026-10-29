# ===== GLOBAL CONFIGURATION =====
PAGE_TITLE = "StarNest - Open Star Cluster Analysis"
PAGE_URL = "/starnest/"            # where main.py mounts the StarNest page + API under Streamlit
FRAME_HEIGHT = 1000                # px
# ================================
# Streamlit page used ONLY when the app runs under `streamlit run main.py` (Streamlit Community Cloud).
# Workflow: Streamlit loads this script (so the platform's health/handshake checks succeed) -> it embeds the
# full StarNest page (served by the FastAPI app mounted at PAGE_URL) in a full-width frame.
import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(page_title=PAGE_TITLE, layout="wide", initial_sidebar_state="collapsed")
st.markdown("<style>header,footer{visibility:hidden;height:0} .block-container{padding:0!important;max-width:100%!important}</style>", unsafe_allow_html=True)
components.iframe(PAGE_URL, height=FRAME_HEIGHT, scrolling=True)
