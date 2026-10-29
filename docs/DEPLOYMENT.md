# Publishing the web app

## Why Streamlit Community Cloud shows "responding 404 to health checks"

Recent Streamlit versions detect an ASGI app (such as the FastAPI `app` in `main.py`) in the main file and serve it themselves. The platform then probes
`/_stcore/health`; the FastAPI app had no such route and answered 404. Reproduced locally with `streamlit run main.py` (Streamlit 1.66).
**Fix (in `main.py` since this version): the route `/_stcore/health` now returns `ok`.** With the main file set to `main.py` the **full app** (including IRSA + PARSEC) is served on streamlit.io.

## Choose a route

| Route | What works | What does not | Effort | Cost |
|---|---|---|---|---|
| **A0. streamlit.io with `main.py`** | **Everything** | Possible IRSA/CMD limits from the cloud address | Update `main.py`, reboot | free |
| **A. Streamlit wrapper** (`streamlit_app/streamlit_app.py`) | Page, GMM, kinematics, IMF, team photos (browser engine) | IRSA extinction and PARSEC isochrone | 5 min, on the platform you already use | free |
| **B. Render.com** (`render.yaml`) | **Everything**, including IRSA + PARSEC (the server makes those requests, so no CORS) | Free tier sleeps after about 15 min idle (first load ~30-60 s); 512 MB RAM | 10 min | free tier |
| **C. Hugging Face Spaces** (`Dockerfile`) | Everything | Same cold start behaviour; needs the README header below | 10 min | free tier |
| **D. GitHub Pages** | Page, GMM, kinematics, IMF | IRSA + PARSEC | 3 min | free |

Recommendation: **A0** if you want to stay on streamlit.io; **B** for the real research tool, **A or D** as a quick public demo.

## Route A0: streamlit.io with `main.py` as the main file (full app)

1. Push the repository with this `main.py`.
2. Keep the app's main file path as `main.py` (the file in the repository root).
3. **Reboot app** in the manage menu. The health check now gets 200 and `https://<your-app>.streamlit.app` shows StarNest with "Engine: Python backend".
4. If the repository has an extra folder level (for example `StarNest-OCMass/main.py`), the main file path must include it.

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
