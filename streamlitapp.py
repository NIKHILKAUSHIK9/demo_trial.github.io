"""Single Streamlit Community Cloud entry point for StarNest / OCMass.

Select this file as the app's Main file path. It hosts the Streamlit shell and
mounts the existing FastAPI backend on the same origin, so the embedded index
page can call /starnest/api/* without a second deployment or localhost URL.
"""
from pathlib import Path

import streamlit as st
from starlette.routing import Mount

from main import app as backend_app

ROOT = Path(__file__).resolve().parent

# Streamlit's ASGI-compatible App (available in Streamlit >= 1.53) lets the
# Streamlit page and FastAPI backend share one host/port.
app = st.App(
    str(ROOT / "streamlit_page.py"),
    routes=[Mount("/starnest", app=backend_app, name="starnest-backend-and-page")],
)
