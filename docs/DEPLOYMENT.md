# StarNest / OCMass deployment guide

## Recommended: Streamlit Community Cloud with the Python backend enabled

This repository now has **`streamlitapp.py` at the project root**. Choose that file as the Streamlit Community Cloud **Main file path**. It creates a Streamlit ASGI app and mounts the FastAPI backend from `main.py` under `/starnest` on the same origin. The embedded `index.html` therefore calls `/starnest/api/*` without a second service or a `localhost` URL.

### GitHub repository layout

The files inside the `ocmass_webapp/` folder must be the repository root. When you unzip the supplied archive, open `ocmass_webapp/` and push its **contents**, not the outer folder as a nested directory.

At minimum, the root should contain:

```text
streamlitapp.py
streamlit_page.py
main.py
requirements.txt
index.html
assets/
images/
ocmass_gmm.py
ocmass_extinction_isochrone.py
ocmass_imf.py
```

### Publish steps

1. Push the project-root files to your GitHub repository's `main` branch.
2. Open Streamlit Community Cloud and create a new app (or edit the existing app's settings).
3. Select your repository and branch `main`.
4. Set **Main file path** to exactly `streamlitapp.py`.
5. Select Python 3.11 or 3.12 and deploy. The root `requirements.txt` pins Streamlit to `1.66.0` because the entry point uses `st.App` and custom ASGI routes.
6. After deployment, verify the page's engine badge says **Python backend**. In a new browser tab, open `<your-app-url>/starnest/api/health`; it should return JSON with `"status":"ok"` and `"engine":"python"`.
7. Upload a small Gaia CSV and test GMM before trying the IRSA/PARSEC network workflow.

`st.App` is Streamlit's ASGI integration with custom-route support. The dependency is pinned to `1.66.0` to keep the entry point and routing behavior reproducible. I could not validate the deployment against your actual Community Cloud instance from this archive review, so perform the health check above after publishing.

### Important limits

- The Python backend can run GMM, kinematics, IMF, and server-side data processing. The IRSA Dust and PARSEC CMD calls still depend on the external services being reachable from Streamlit's cloud environment. If those requests time out or are blocked, use the fixed E(B-V) mode as a fallback or deploy the full backend on Render / Hugging Face Spaces as described below.
- Uploaded catalogues and analysis sessions are stored in memory and will be lost when the app restarts. Do not treat this as persistent storage.
- Use only public data and team photos you have permission to publish. Check the licence and replace template repository references before release.

## Local development

From the repository root:

```bash
python -m venv .venv
# Windows PowerShell:
.\.venv\Scripts\Activate.ps1
# macOS/Linux:
# source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run streamlitapp.py
```

For testing the FastAPI backend directly, use a separate terminal:

```bash
uvicorn main:app --host 127.0.0.1 --port 8000
```

This direct FastAPI mode serves `index.html` at `http://127.0.0.1:8000/` and the API at `http://127.0.0.1:8000/api/health`.

## Alternative: Render (full backend)

1. Sign in to Render using GitHub.
2. Create a new Blueprint from the repository; it reads `render.yaml`. Or create a Python Web Service with:
   - Build command: `pip install -r requirements.txt`
   - Start command: `uvicorn main:app --host 0.0.0.0 --port $PORT`
   - Health check path: `/api/health`
3. Deploy, then check `/api/health` at the generated service URL.

Render's free instances may sleep when idle. The service uses in-memory sessions, so a restart clears work in progress.

## Alternative: Hugging Face Spaces (Docker)

1. Create a new Space with the **Docker** SDK.
2. Add this header at the top of the Space repository's `README.md`:

   ```yaml
   ---
   title: StarNest OCMass
   sdk: docker
   app_port: 7860
   ---
   ```
3. Push the project files, including `Dockerfile`, and deploy.

## Browser-only demo (not the recommended full-backend entry point)

`streamlit_app/streamlit_app.py` is the older lightweight wrapper. It deliberately sets `window.STARNEST_NO_BACKEND=true` and runs the browser engine only. It **does not run the Python backend**, so do not select this file when you need IRSA/PARSEC server-side processing or the Python backend badge.

GitHub Pages is also browser-only and cannot host this Python backend.
