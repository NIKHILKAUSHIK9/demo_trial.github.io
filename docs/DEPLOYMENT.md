# Publishing the web app

## Streamlit Community Cloud: one file runs page + backend

Main file: **`streamlit_app/streamlit_app.py`** (branch `main`). It chooses the engine itself:
1. `STARNEST_BACKEND_URL` secret (or `BACKEND_URL` in the file) set -> the page calls that external backend (Render / Hugging Face).
2. Otherwise **bridge mode (default): Streamlit runs the backend in-process.** Community Cloud exposes only one port, so the page cannot call a FastAPI server. The page is a Streamlit component; each API call goes page -> Streamlit script -> FastAPI app of `main.py` (in memory) -> page. Uploads larger than 2 MB are sent in parts (the browser-to-Streamlit message limit is ~25 MB). Per-visitor state travels in the `X-Session-Id` header.
3. If `main.py` cannot be imported (missing package), the page falls back to the browser-only engine.

**Do not use `main.py` as the Streamlit main file**: it is a plain FastAPI app; Community Cloud shows it as an endless "buffering" page.

Dependencies come from `streamlit_app/requirements.txt` (next to the main file; it lists the backend packages too).

Tested locally with `streamlit run streamlit_app/streamlit_app.py` (no secret set): badge "Engine: Python backend", a 29 MB CSV upload, GMM on UPK 545 (34 members, K = 4), file download. **Not tested on Community Cloud itself** (not reachable from the build environment), so IRSA/CMD from Cloud addresses, memory (about 1 GB) and the first-load time (backend import ~10-20 s) are unverified. If it fails there, use route E or B.

## Choose a route

| Route | What works | What does not | Effort | Cost |
|---|---|---|---|---|
| **A. Streamlit only** (`streamlit_app/streamlit_app.py`) | Page + in-process Python backend (bridge) incl. IRSA + PARSEC | Untested on Cloud; shared 1 GB RAM; slow first load | 5 min | free |
| **E. Streamlit front end + Render/HF backend** | Everything, more robust | Two deployments; backend cold start | B or C, then one secret | free |
| **B. Render** (`render.yaml`), whole site | **Everything** | Free tier sleeps after ~15 min idle (first load 30-60 s) | 10 min | free tier |
| **C. Hugging Face Space** (`Dockerfile`), whole site | Everything | Same cold start; README header below | 10 min | free tier |
| **D. GitHub Pages** | Page, GMM, kinematics, IMF (browser engine) | IRSA + PARSEC | 3 min | free |

## Route A: Streamlit only

1. share.streamlit.io -> Create app -> your repo, branch `main`, main file `streamlit_app/streamlit_app.py`, Python 3.11+.
2. No secrets needed. Deploy / Reboot.
3. Wait for the badge "Engine: Python backend (PARSEC + IRSA)" (first load up to ~20 s). "Engine: in-browser" means the backend import failed: read **Manage app -> logs** for `StarNest backend unavailable`.

## Route E: Streamlit front end + backend elsewhere

1. Deploy the backend with route B or C; `https://<backend>/api/health` must return JSON.
2. Same Streamlit app as route A, plus Settings -> Secrets: `STARNEST_BACKEND_URL = "https://<backend>"` (no trailing slash, no `/api`).
3. Reboot. Backend cold start 30-60 s on the first call.

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
