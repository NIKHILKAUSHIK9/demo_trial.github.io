# ===== GLOBAL CONFIGURATION =====
PAGE_TITLE = "StarNest - Open Star Cluster Analysis"
PAGE_URL = "/starnest/"            # FastAPI page + API mounted by streamlitapp.py
FRAME_HEIGHT = 1000                # px
# ================================
# Streamlit page used by the root `streamlitapp.py` ASGI entry point (Streamlit Community Cloud).
# Workflow: the Streamlit ASGI app serves this script and mounts the backend -> this script embeds the
# full StarNest page (served by the FastAPI app mounted at PAGE_URL) in a full-width frame.
import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(page_title=PAGE_TITLE, layout="wide", initial_sidebar_state="collapsed")
st.markdown("<style>header,footer{visibility:hidden;height:0} .block-container{padding:0!important;max-width:100%!important}</style>", unsafe_allow_html=True)
components.iframe(PAGE_URL, height=FRAME_HEIGHT, scrolling=True)
