# Publishing the web app

## Why Streamlit Community Cloud shows "responding 404 to health checks"

Streamlit Community Cloud only runs **Streamlit apps** (a Python file that calls `import streamlit`). StarNest is a plain web page plus a **FastAPI** server,
so no Streamlit server ever starts on the platform and its health check (`/_stcore/health`) gets a 404. No setting in the dashboard can fix this: the entry file must be a Streamlit app.

## Choose a route

| Route | What works | What does not | Effort | Cost |
|---|---|---|---|---|
| **A. Streamlit wrapper** (`streamlit_app/streamlit_app.py`) | Page, GMM, kinematics, IMF, team photos (browser engine) | IRSA extinction and PARSEC isochrone | 5 min, on the platform you already use | free |
| **B. Render.com** (`render.yaml`) | **Everything**, including IRSA + PARSEC (the server makes those requests, so no CORS) | Free tier sleeps after about 15 min idle (first load ~30-60 s); 512 MB RAM | 10 min | free tier |
| **C. Hugging Face Spaces** (`Dockerfile`) | Everything | Same cold start behaviour; needs the README header below | 10 min | free tier |
| **D. GitHub Pages** | Page, GMM, kinematics, IMF | IRSA + PARSEC | 3 min | free |

Recommendation: **B** for the real research tool, **A or D** as a quick public demo.

## Route A: keep streamlit.io (browser-only mode)

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
