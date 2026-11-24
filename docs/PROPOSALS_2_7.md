# Phase 2.7: proposals (NOT implemented; need your decision)

Confidence tags: [Certain] = verified in this repository, [Likely] = strong inference, [Guessing] = needs checking.

## A. The in-browser JavaScript engine: delete it, or replace it with Pyodide?

### Why it is a problem now
- [Certain] The JS engine is not numerically equivalent to the Python back end. It uses the binned slope fit (no MLE), has no measured completeness limit, no approx_ML gate, no cluster-component rule, and none of the Phase 1/2 analyses (reliability, likelihood fit, velocities, binaries, dynamics). Panels built in Phase 2 already refuse to run without Python.
- [Certain] When the back end is unreachable the page silently falls back to that engine. A user can therefore publish a mass-function slope from a method the paper does not describe.
- [Certain] Every scientific change has to be written twice (Python and JS), and only the Python copy is tested.

### Options
| | A1. Delete the JS engine | A2. Pyodide (run the real Python in the browser) | A3. Keep, but label |
|---|---|---|---|
| What the user gets offline / on a static host | Page shows "start the backend" | Full pipeline, same code | Approximate numbers |
| Numerical parity | n/a (one engine) | Same source files; [Likely] differences only at floating-point level (BLAS) | None |
| Cost | Delete ~2k lines of `app.js`; update README and tests | Spike 1-2 days; then ship `starnest/` as a wheel loaded by Pyodide | Banner only |
| Stage 2 / 3 online services | Backend does them | [Likely] IRSA and CMD do not send CORS headers, so a browser cannot call them; Pyodide mode would need the local grid (2.2) and a manual E(B-V) (or a small proxy) | n/a |
| Load time / size | smallest | [Likely] 30-60 MB first load (numpy, pandas, scipy, scikit-learn, matplotlib, astropy are available as Pyodide packages) | none |
| galpy (orbit panel) | backend only | [Guessing] no Pyodide build; orbit panel stays backend-only | n/a |

### Recommendation
1. Do **not** delete anything yet: the local PARSEC grid is not built (needs your internet run) and Pyodide is unmeasured.
2. Run a Pyodide spike gated on two numbers: first load under about 15 s on a normal connection, and UPK 545 results identical to the Python run to 1e-6.
3. If the spike passes, replace the JS engine with Pyodide (A2) and delete the JS maths. If it fails, take A1 (delete) and rely on the Streamlit bridge or the external back end (section B).
4. Until then, A3 is cheap and removes the silent-fallback risk: a red "Browser engine: approximate, not the published method" banner on every result and no export of IMF slopes from it. I can do A3 in Phase 3 if you want it.

Risk of A1: users who open `index.html` from a file with no Python lose the page's functionality (they keep the CLI). Risk of A2: scope creep (Pyodide packaging, CORS, worker threads to keep the page responsive).

## B. Streamlit with an external back end: job + poll design

### Why
- [Likely] Streamlit Community Cloud runs the in-process bridge on a small shared machine. The long steps (binaries Monte Carlo, bootstraps, likelihood fit) hold the script thread for minutes, and a sleeping or restarted app loses every session's state.
- [Certain] `main.py` sessions live in process memory, so any restart wipes them. This is already documented in the README.

### Design
```
Browser (Streamlit page)                       External back end (uvicorn starnest.api:app, one process)
  POST /api/jobs {endpoint, body}  ---------->  validate with the same Pydantic model, create job, 202 {job_id}
  GET  /api/jobs/{id}   (every 2 s) -------->  {state: queued|running|done|failed|cancelled, elapsed_s, result?, error?}
  DELETE /api/jobs/{id}                        cooperative cancel (checked between bootstrap / MC batches)
```
- Worker: `ThreadPoolExecutor(max_workers=2)`; the worker sets the session id context variable before calling the existing route function, so **no scientific code changes**. Per-session limit of 1 running job, global queue limit, 10-minute timeout, result kept 1 h.
- Session state: stays in `_SESSIONS` (now lock-protected, 2.1); add a per-session lock so two jobs of one visitor cannot interleave stage 2 and stage 3.
- Capability flag: `/api/health` returns `{"jobs": true}`. `app.js` wraps `apiCall`: for the long endpoints, if `jobs` is advertised, submit and poll with a progress line; otherwise call directly as today. Short endpoints are unchanged.
- Streamlit: the page talks to the back end directly from the browser (`API_BASE`), so Streamlit only hosts the HTML. Needs a CORS allow-list for `https://starnest.streamlit.app`; no cookies (the existing `X-Session-Id` header already supports this).
- Security: upload size cap (exists), per-IP rate limit, no secrets in the page, sessions isolated by an unguessable id, jobs never write outside a per-session temp directory.

### Open decisions for you
1. Host. [Guessing] Hugging Face Spaces (Docker, free CPU tier, more RAM than Streamlit Cloud) or a university/institute server; Render's free tier sleeps and has small RAM. Check current limits before choosing.
2. Persistence: in-memory is enough for a demo; if sessions must survive restarts, add SQLite or Redis (extra work, not in this proposal).
3. Who may use it: public with rate limits, or behind a token.

### Effort and risk
About 2-3 days: ~150 lines for the job layer, ~60 in `app.js`, tests for queueing, cancel and timeout. Main risk: a single-process back end is a single point of failure and a shared-CPU contention point; limits above bound it but do not remove it.
