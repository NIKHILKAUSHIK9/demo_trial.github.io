# Publishing the web app

## Streamlit Community Cloud: what the platform needs

Streamlit runs the main file and, if it finds an ASGI object called `app` (the FastAPI app), serves that object itself. Two things followed:
1. The platform health probe `/_stcore/health` returned 404 -> "responding 404 to health checks". `main.py` now answers it.
2. A bare FastAPI page passes the probe but is not a Streamlit app, so the platform page kept loading ("buffering"). **Fix (in `main.py` since this version):** when started by Streamlit, `app` becomes `st.App(streamlit_page.py, routes=[Mount('/starnest', FastAPI app)])`. Streamlit serves a genuine Streamlit page (`streamlit_page.py`), which embeds the complete StarNest page and its Python API (mounted at `/starnest`) in a full-width frame. Under `uvicorn main:app` (Render, Hugging Face, local) nothing changes.

Tested locally with `streamlit run main.py` (Streamlit 1.66): `/_stcore/health`, `/`, `/starnest/api/health` return 200; the badge reads "Engine: Python backend"; GMM on UPK 545 gives 34 members, K = 4. **Not tested on Community Cloud itself** (not reachable from the build environment); if it still fails there, use route B or C.

## Choose a route

| Route | What works | What does not | Effort | Cost |
|---|---|---|---|---|
| **A0. streamlit.io with `main.py`** | **Everything** (Streamlit page embeds the full app) | Untested on Community Cloud; possible IRSA/CMD limits from the cloud address | Update repo, reboot | free |
| **A. Streamlit wrapper** (`streamlit_app/streamlit_app.py`) | Page, GMM, kinematics, IMF, team photos (browser engine) | IRSA extinction and PARSEC isochrone | 5 min, on the platform you already use | free |
| **B. Render.com** (`render.yaml`) | **Everything**, including IRSA + PARSEC (the server makes those requests, so no CORS) | Free tier sleeps after about 15 min idle (first load ~30-60 s); 512 MB RAM | 10 min | free tier |
| **C. Hugging Face Spaces** (`Dockerfile`) | Everything | Same cold start behaviour; needs the README header below | 10 min | free tier |
| **D. GitHub Pages** | Page, GMM, kinematics, IMF | IRSA + PARSEC | 3 min | free |

Recommendation: **A0** to stay on streamlit.io (fallback **A**); **B** for the real research tool; **D** as a quick public demo.

## Route A0: streamlit.io with `main.py` as the main file (full app)

1. Push the repository with this `main.py`.
2. Keep the app's main file path as `main.py` (the file in the repository root).
3. **Reboot app** in the manage menu. Files needed in the repository root: `main.py`, `streamlit_page.py`, `requirements.txt` (includes `streamlit`), `assets/`, `images/`, `ocmass_*.py`, `index.html`. The app shows the StarNest page with "Engine: Python backend".
4. If the first load seems stuck, press Reboot app once and wait ~1 min (dependencies install). If the repository has an extra folder level (for example `StarNest-OCMass/main.py`), the main file path must include it.

Limits: Community Cloud memory is small; IRSA/CMD may behave differently from a cloud address. Sessions are in memory.

## Route A: keep streamlit.io (browser-only wrapper)

1. Push the repository (already done).
2. share.streamlit.io → your app → **Settings** (or create a new app).
3. **Main file path:** `streamlit_app/streamlit_app.py`. Branch: `main`.
4. Python 3.11 or newer. The dependency file `streamlit_app/requirements.txt` (only `streamlit`) is found next to the main file.
5. Save / **Reboot app**. The page appears inside the Streamlit page; the badge reads "Engine: in-browser".
6. Remove the old app entry if it still points to another main file.

Limits: the page runs inside an iframe, so only the browser engine is available.

## Route B: Render (full app)

1. https://render.com → sign in with GitHub.
2. **New → Blueprint** (reads `render.yaml`) or **New → Web Service** with:
   - Build command: `pip install -r requirements.txt`
   - Start command: `uvicorn main:app --host 0.0.0.0 --port $PORT`
   - Health check path: `/api/health`
   - Environment variable `PYTHON_VERSION` = `3.11.9`
3. Deploy; open the `onrender.com` URL. The badge reads "Engine: Python backend (PARSEC + IRSA)".
4. If IRSA or CMD time out from the cloud IP, press Apply again (results are cached); use E(B-V) mode `fixed` as a fallback.

Public hosting notes:
- Each visitor gets a private in-memory session (cookie `sn_sid`); sessions expire after 6 h idle and at most 20 are kept. Uploaded files are never written to disk.
- Uploads are limited to 50 MB (`MAX_UPLOAD_MB` in `main.py`).
- IRSA and CMD are shared academic services: the backend caches results and rate-limits requests. Do not add automatic bulk queries.
- Free instances restart and forget all sessions.

## Route C: Hugging Face Spaces

1. huggingface.co → **New Space** → SDK **Docker**.
2. Add this header at the very top of `README.md` in the Space repository:
   ```
   ---
   title: StarNest OCMass
   sdk: docker
   app_port: 7860
   ---
   ```
3. Push the repository to the Space. The `Dockerfile` starts the server on port 7860.

## Route D: GitHub Pages

Repository → Settings → Pages → branch `main`, folder `/ (root)`. Same limits as route A.

## Before going public

| Check | Why |
|---|---|
| Team photos added with consent | They become public |
| License confirmed with supervisors | Research code policy |
| `<your-user>/<your-repo>` replaced in README and CITATION.cff | Correct links |
| Visit the published URL in a private window | Confirms a new visitor gets the page and a clean session |
