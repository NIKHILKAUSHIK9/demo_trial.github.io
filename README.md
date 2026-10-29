# StarNest / OCMass: photometry and kinematics of young open clusters (Gaia DR3)

A web interface plus a Python backend that takes a Gaia DR3 cone-search CSV of an open cluster and runs the full analysis chain:
**membership (GMM) → extinction + isochrone (IRSA Dust, PARSEC CMD 3.8) → kinematics → mass function (IMF)**.
Target cluster used during development: UPK 545.

> Project documentation with every term, parameter and formula: [`docs/StarNest_Project_Documentation.docx`](docs/StarNest_Project_Documentation.docx)

## Quick start

```bash
git clone https://github.com/<your-user>/<your-repo>.git
cd <your-repo>
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
uvicorn main:app --port 8000
```
Open **http://localhost:8000** in a browser.

| Mode | How | What works |
|---|---|---|
| Full (recommended) | `uvicorn main:app --port 8000`, open `http://localhost:8000` | Everything, incl. IRSA E(B-V) and the PARSEC isochrone from CMD 3.8 |
| Browser only | open `index.html` (or GitHub Pages) | GMM, kinematics, IMF, isochrone plot with a **typed** E(B-V) but **no** PARSEC (IRSA and CMD refuse direct browser requests: CORS) |

**Publishing online:** for Streamlit Community Cloud with the Python backend, choose the root `streamlitapp.py` as the Main file path. See [`docs/DEPLOYMENT.md`](docs/DEPLOYMENT.md) for GitHub setup, health checks, Render / Hugging Face alternatives, and the browser-only demo distinction.

The page loads Tailwind, GSAP, Plotly, Font Awesome and the Inter font from CDNs, so it needs internet access.

## Input file

A VizieR Gaia DR3 CSV with at least these columns: `Plx, e_Plx, pmRA, pmDE, e_pmRA, e_pmDE, Gmag, BP-RP, RUWE, gofAL, Nper`;
the isochrone stage also needs `RA_ICRS, DE_ICRS, BPmag, RPmag`.

## Workflow

| # | Page | What it does | Code |
|---|---|---|---|
| 1 | Membership | Quality cuts, distance window, GMM in (pmRA, pmDE, Plx), cluster component = smallest det(Σ), BIC choice of K, star classification | `ocmass_gmm.py`, `main.py` (`/run-gmm`) |
| 2 | Isochrone | IRSA Dust E(B-V) → A_V, A_G, A_BP, A_RP → CMD 3.8 PARSEC file (Gaia EDR3 = DR3, Vega) → isochrone over the members, fit statistics, downloadable `.dat` | `ocmass_extinction_isochrone.py`, `main.py` (`/compute-cmd`) |
| 3 | Kinematics | 2-D proper-motion vectors and a 3-D map with parallax-error diagnostics | `main.py` (`/compute-kinematics`), `assets/js/app.js` |
| 4 | IMF | Log-mass bins, Salpeter / Kroupa / Chabrier comparison, reduced χ² | `ocmass_imf.py`, `main.py` (`/compute-imf`) |

## Repository layout

```
index.html                  page markup
assets/css/styles.css       styles
assets/js/app.js            UI logic + in-browser GMM / IMF engine
assets/js/tailwind-config.js
images/team/                team photos (see the README in that folder)
main.py                     FastAPI backend (serves the page and the /api/* endpoints)
ocmass_gmm.py               stage 1: quality cuts + GMM membership
ocmass_extinction_isochrone.py  stage 2: IRSA + CMD 3.8 + fit
ocmass_imf.py               stage 4: mass function
tools/netcheck.py           tests the network path to IRSA and CMD
tools/cmd_probe.py          prints the live CMD form fields (debugging)
streamlitapp.py              Streamlit ASGI entry point: mounts the full FastAPI backend
streamlit_page.py            Streamlit page shell that embeds the StarNest HTML at /starnest/
streamlit_app/              Legacy lightweight wrapper (browser-only mode; not full backend)
Dockerfile, render.yaml     deployment of the full app
docs/                       project documentation, GitHub guide, deployment guide, development notes
```
Every Python module has a `GLOBAL CONFIGURATION` block at the top: change thresholds and paths there, not inside functions.
The Python modules contain no UI code and can be imported on their own, e.g. `import ocmass_gmm as g; g.run_gmm(...)`.

## Network notes (backend)

IRSA and CMD are shared public services: results are cached, requests are rate-limited, nothing is fetched on a slider move.
Environment switches: `OCMASS_NO_PROXY=1`, `OCMASS_PROXY=http://host:port`, `OCMASS_IPV4_ONLY=1`. Run `python tools/netcheck.py` if a request times out.

## Known limitations

- Star classification thresholds (white dwarf / giant / pre-main-sequence) are uncalibrated heuristics.
- IRSA's "SFD mean" is a regional mean (≥ 2°) of the total line-of-sight dust: it over-estimates E(B-V) for nearby clusters.
- The isochrone stage compares a **chosen** logAge and [M/H]; it does not search for the best age.
- Backend sessions live in memory (one per visitor, expire after 6 h idle, at most 20); a restart forgets them.
- Browser-only mode cannot reach IRSA or CMD.

## Citing

See `CITATION.cff`. Data: Gaia DR3 (ESA/Gaia/DPAC); dust: Schlafly & Finkbeiner (2011) via IRSA; isochrones: PARSEC (Bressan et al. 2012) through the CMD service of INAF-OAPd.

## License

MIT, see `LICENSE`.
