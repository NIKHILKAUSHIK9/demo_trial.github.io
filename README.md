# StarNest

**Photometry and kinematics of open star clusters with Gaia DR3**: from a raw Gaia cone-search CSV to cluster members, extinction, an isochrone fit, 3-D / proper-motion kinematics and a mass-function comparison, in one web app and as plain Python modules.

| | |
|---|---|
| Live app | <https://starnest.streamlit.app/> |
| Example cluster | UPK 545 (Gaia DR3 cone search, 71 830 sources) |
| Language / stack | Python 3 (NumPy, pandas, SciPy, scikit-learn, Astropy/astroquery, FastAPI) + HTML/JS front end (Plotly) |
| Licence | MIT (see `LICENSE`) |

---

## Contents

1. [What the app does](#1-what-the-app-does)
2. [Workflow](#2-workflow)
3. [Parameters by stage](#3-parameters-by-stage)
4. [External websites and services](#4-external-websites-and-services)
5. [Repository layout](#5-repository-layout)
6. [Using the modules from Python](#6-using-the-modules-from-python)
7. [Limitations and assumptions](#7-limitations-and-assumptions)
8. [Team](#8-team)
9. [How to cite](#9-how-to-cite)
10. [References](#10-references)
11. [Run the app offline on your own machine](#11-run-the-app-offline-on-your-own-machine)

---

## 1. What the app does

| Stage | Question answered | Module |
|---|---|---|
| 1. Membership | Which stars belong to the cluster? | `gmm.py` |
| 2. Extinction | How much dust lies in front of the cluster? | `extinction_isochrone.py` |
| 3. Isochrone | Which PARSEC isochrone matches the colour-magnitude diagram (CMD)? | `extinction_isochrone.py` |
| 4. Kinematics | What is the 3-D shape of the cluster and how do the stars move on the sky? | `kinematics.py` |
| 5. Mass function | Does the stellar mass distribution follow Salpeter, Kroupa or Chabrier? | `imf.py` |

`main.py` is the FastAPI back end that connects the modules to the web page (`index.html`, `assets/`). Each module has a **GLOBAL CONFIGURATION** block at the top, so thresholds and defaults can be changed in one place.

## 2. Workflow

```
Gaia DR3 CSV (VizieR)
   |
   v
[1] Quality cuts -> distance window -> standardise (pmRA, pmDE, Plx) -> Gaussian Mixture Model
   |                                                                  (K by BIC, tightest component = cluster)
   v  members (PMemb >= cutoff)
[2] IRSA Dust: E(B-V) -> A_V -> A_G, A_BP, A_RP -> de-reddened M_G, (BP-RP)_0
   |
   v
[3] CMD 3.8 (PARSEC, Gaia EDR3 Vega) isochrone at logAge, [M/H], A_V -> overlay + fit quality
   |        -> star types and masses (isochrone interpolation)
   v
[4] Kinematics: 3-D map (distance filter / shrinkage) + residual proper-motion vectors
   |
   v
[5] Mass function: log-mass bins -> Poisson errors -> Salpeter / Kroupa / Chabrier -> reduced chi^2
```

Extinction is applied **once**: either to the isochrone (`reddened`, stars plotted as observed) or to the stars (`intrinsic`, isochrone at A_V = 0). The two are never mixed.

Typical use in the web app: upload the CSV on the **GMM** tab, run the analysis, then open **Isochrone**, **Kinematics** and **IMF**; every tab reads the member table of the previous step, and each stage offers CSV downloads.

---

## 3. Parameters by stage

### 3.0 Input data

The input is a Gaia DR3 CSV exported from VizieR (catalogue `I/355/gaiadr3`) with these columns (names as VizieR writes them):

| Column | Meaning | Unit | Used in |
|---|---|---|---|
| `RA_ICRS`, `DE_ICRS` | Sky position (ICRS) | deg | extinction, kinematics, sky plots |
| `Plx`, `e_Plx` | Parallax and its error | mas | distance, quality cut, GMM, kinematics |
| `pmRA`, `pmDE`, `e_pmRA`, `e_pmDE` | Proper motion (RA includes cos Dec) and errors | mas/yr | quality cut, GMM, kinematics |
| `Gmag`, `BPmag`, `RPmag`, `BP-RP` | Gaia magnitudes and colour (Vega) | mag | CMD, extinction |
| `RUWE`, `gofAL`, `Nper` | Astrometric quality flags | - | quality cut |

### 3.1 Stage 1: quality cuts and GMM membership (`gmm.py`)

**Quality cuts** keep stars with a well-behaved astrometric solution. A star must pass all of them.

| Parameter | Default | What it is | How it is used |
|---|---|---|---|
| `RUWE_MAX` | 1.6 | Renormalised unit weight error; about 1 for a good single-star fit [3] | keep `RUWE < 1.6` |
| `E_PMRA_MAX`, `E_PMDE_MAX` | 0.1 mas/yr | Proper-motion uncertainty | keep `e_pmRA < 0.1` and `e_pmDE < 0.1` |
| `PLX_SNR_MIN` | 5 | Parallax signal-to-noise | keep `Plx / e_Plx > 5` |
| `GOFAL_MAX` | 3 | Goodness of fit of the along-scan astrometry | keep `gofAL < 3` |
| `NPER_MIN` | 10 | Number of visibility periods used | keep `Nper > 10` |

**Distance window and GMM.** A Gaussian Mixture Model [14] is fitted in the standardised space (pmRA, pmDE, Plx). The component with the smallest covariance determinant (the tightest blob in proper motion and parallax) is the cluster; every star gets `PMemb`, the posterior probability of belonging to it.

| Parameter | Default | What it is | How it is used |
|---|---|---|---|
| `DIST_LIT_PC` (`Dist_Lit`) | 341 pc | Literature distance of the cluster (UPK 545 in the example; **set your own**) | centre of the distance window |
| `HALF_WIDTH_PC` (`Half_Width`) | 195.5 pc | Half-width of the window | keep `abs(1000/(Plx - zp) - Dist_Lit) < Half_Width` |
| `PLX_ZP` (`Plx ZP`) | 0.0 mas | Parallax zero-point offset; the Gaia DR3 global value is about -0.017 mas [4] | distance = `1000 / (Plx - zp)` |
| `CUTOFF` (`PMemb Cutoff`) | 0.5 | Membership probability limit | member: `PMemb >= cutoff`; non-member: `PMemb <= 1 - cutoff` |
| `FEATURES` | pmRA, pmDE, Plx | Dimensions of the mixture | standardised as `(x - median) / std` |
| `K_DEFAULT` (`GMM Components K`) | `None` = auto | Number of Gaussians | auto: fit K = 2 to 6 and keep the lowest BIC; or force one value |
| `K_RANGE` | (2, 6) | Range searched in auto mode | - |
| BIC | formula | Bayesian Information Criterion [14]: `-2 ln L + p ln N`, with `p = 10K - 1` free parameters for full covariances in 3-D | lower is better; penalises extra components |
| `OPTIMIZER` (`EM optimizer`) | `converged` | EM stopping rule. `colab`: scikit-learn defaults (tol 1e-3, 100 iterations), reproduces the original Colab script. `converged`: `GMM_TOL` 1e-6, `GMM_MAX_ITER` 1000 | `converged` is recommended for papers |
| `N_INIT`, `RANDOM_STATE` | 5, 42 | Random restarts and seed | reproducible fit |
| `MIN_COMPONENT_STARS` | 5 | A fitted component with fewer stars is a collapsed (singular) solution | that K is excluded from the BIC choice |
| `MIN_STARS_PER_COMPONENT` | 10 | Minimum window stars per component | error if `N < 10 x K_max` |
| MSS | output | Mean Separation Score: mean over features of `(std_nonmember - std_member) / max(std)` | 1 = compact members on a broad field, 0 = no separation |

Stage 1 applies **no** extinction: the member table carries only observed-plane columns (`M_G_obs`, `bp_rp_obs`) and a provisional `Star_Type` (`Star_Type_plane` = observed). `M_G_corr`, `bp_rp_corr` and `MC_Mass` appear after stage 2, together with `E_BV_used` (the E(B-V) actually applied) and `Mass_quality` (`ok` / `approx`). An explicit `ebv_ui` in the stage-1 request is still honoured and stamped in `E_BV_used`.

### New analysis endpoints (additive; existing routes are unchanged)

| Endpoint | Module | What it returns |
|---|---|---|
| `POST /api/gmm-diagnostics` | `gmm.py` | HDBSCAN cross-check (Jaccard overlap with the GMM members) and per-star PMemb scatter from `N_RESAMPLE` draws within the Gaia errors |
| `POST /api/compute-binaries` | `binaries.py` | `MC_Mass`, `MC_CompMass`, `is_binary`, `q` per star (nearest synthetic system in G, BP, RP, `N_MC` repeats) and the cluster mass `M_T`. `is_binary` is a probability under the assumed `F_BIN`, not a detection |
| `POST /api/compute-velocities` | `kinematics.py` | `v = 4.74 mu / plx`, perspective-corrected (Gaia RV), median `v_out` / `v_az` with bootstrap errors, ML `sigma_1D` with an outlier mixture chosen by delta-BIC, virial sigma |
| `POST /api/compute-dynamics` | `dynamics.py` | Mass-segregation ratio sweep, tidal radius, relaxation time, optional galpy orbit and tail check |
| `POST /api/fit-isochrone-likelihood` | `extinction_isochrone.py` | Likelihood fit of A_V, DM, equal-mass-binary fraction and outlier fraction; `dust_source` = `sfd` / `bayestar` / `free` |
| `GET /api/hunt-reffert?name=` | `gmm.py` | Optional literature prefill; never applied automatically |

`/api/run-gmm` also accepts optional `cluster_rule` (`smallest_det` or `closest_to_literature`, with `lit_pm` and `lit_plx`).
Everything that needs an isochrone (binaries, dynamics, isochrone fit) requires the "Apply Extinction & Fetch PARSEC" step first.

### 3.2 Stage 2: extinction (`extinction_isochrone.py`)

Dust between us and the cluster dims and reddens the stars. The colour excess E(B-V) is read from the SFD dust map [7] through IRSA (see section 4).

| Parameter | Default | What it is | How it is used |
|---|---|---|---|
| `EBV_MODE` | `center` | `center`: one IRSA query at the cluster centre. `per_star`: one query per 0.05 deg sky cell, then the mean. `fixed`: your own number | gives the mean E(B-V) |
| `EBV_FIXED` | none | E(B-V) you supply (mag) | used when mode is `fixed` |
| `EBV_KEY` | `ext SFD mean` | IRSA column: the SFD value averaged over a region of a few degrees around the position (the total dust along the line of sight, an **upper limit** for a nearby cluster) | source of E(B-V) |
| `AV_PER_EBV` | 2.742 | `A_V / E(B-V)` for SFD maps [8] (3.1 would over-correct SFD by about 13 %) | `A_V = 2.742 x E(B-V)` |
| `AG_AV`, `ABP_AV`, `ARP_AV` | 0.789, 1.002, 0.589 | Gaia-band to V-band extinction ratios, fixed (values of the original analysis script, not recomputed) | `A_G = 0.789 A_V`, `A_BP = 1.002 A_V`, `A_RP = 0.589 A_V` |
| `USE_CLUSTER_DM` | True | One distance modulus from the median member parallax | `DM = 5 log10(d_pc) - 5` |
| `IRSA_ROUND_DEG`, `IRSA_MAX_QUERIES` | 0.05 deg, 150 | Sky-cell size and query budget for `per_star` | stars in one cell share one IRSA call |
| `BPRP0_RANGE` | (-0.8, 5.0) | Plausible intrinsic BP-RP | outside it the star is flagged `bprp_corr_ok = False` (over-correction) |

Derived columns:

| Column | Definition |
|---|---|
| `M_G_corr` | `G - DM - A_G` (absolute, de-reddened G) |
| `bp_rp_corr` | `(BP - A_BP) - (RP - A_RP)` (intrinsic colour) |
| `M_G_obs`, `bp_rp_obs` | `G - DM` and observed `BP - RP` (not de-reddened) |

### 3.3 Stage 3: PARSEC isochrone and fit (`extinction_isochrone.py`)

The isochrone is requested from the CMD web service (PARSEC [10], [11]) in the Gaia EDR3 Vega system, which is the correct system for DR3 photometry (G, BP, RP are identical in EDR3 and DR3) [5].

| Parameter | Default / range | What it is | How it is used |
|---|---|---|---|
| `logAge` | 7.5 in the app; accepted 6.6 to 10.13 | log10 of the cluster age in years | isochrone age |
| `[M/H]` | 0.0; sanity range -2.5 to 0.8 | Metallicity | isochrone composition |
| `iso_extinction` | `reddened` | `reddened`: mean A_V is sent to CMD and stars are plotted as observed. `intrinsic`: A_V = 0 is sent and the stars are de-reddened | extinction is applied once |
| `PHOTSYS`, `PHOTSYS_VERSION` | `gaiaEDR3`, `YBCnewVega` | Photometric system and bolometric-correction tables (YBC [12], Vega spectrum [13]) | magnitudes G, G_BP, G_RP |
| `CMD_URL` | `.../cgi-bin/cmd_3.8` | CMD 3.8 CGI. The generic `/cgi-bin/cmd` now serves CMD 3.9, which this client does not use | request target |
| `CMD_BASE_PARAMS` | form defaults | Every field of the CMD 3.8 form is sent (PARSEC v1.2S tracks, Cardelli extinction curve [9], Kroupa initial mass function file, ...) so nothing is left to server guessing | request body |
| `ISO_MAX_LABEL` | 7 | PARSEC evolutionary labels drawn (0 pre-MS, 1 MS, 2 sub-giant, 3 RGB, 4 core-He, 5 to 7 early AGB) | polyline used for the fit |
| `FIT_NEAR_MAG` | 0.2 mag | A star "fits" when its distance to the isochrone polyline in the CMD plane is below this | fit-quality percentage |
| `CHECK_DG_RANGE`, `CHECK_DC_RANGE` | 0.60 to 1.00; 0.25 to 0.65 | Expected shift of the reddened isochrone per unit A_V in G and in BP-RP | self-check that the server really applied A_V |
| `CMD_TIMEOUT_S`, `ISO_RETRIES`, `ISO_WAIT_S`, `CMD_MIN_INTERVAL_S` | 180 s, 4, 5 s, 2 s | Network politeness for the shared server | retry only on connection errors |

**Star types and masses** (`classify_star`, `derive_columns` in `main.py`): rough labels from the de-reddened CMD, used for colours and for excluding evolved stars from the mass function.

| Type | Rule (`M_G` = `M_G_corr`, `c` = `bp_rp_corr`) |
|---|---|
| White Dwarf | `M_G > 10` and `c < 0.6` |
| Red Giant / Supergiant | `M_G < 3` and `c > 0.8` |
| Pre-Main Sequence | `4 < M_G < 3c + 2` |
| Main Sequence | everything else |

Mass (`MC_Mass`): the isochrone's main-sequence mass is interpolated at each star's `M_G`; outside the isochrone range a piecewise mass-luminosity relation is the fallback (flagged in `MC_Mass_src`). These thresholds are heuristic and not calibrated.

### 3.4 Stage 4: kinematics (`kinematics.py`)

| Parameter | Default | What it is | How it is used |
|---|---|---|---|
| `KIN3D_MODE` | `precision` | `precision`: keep stars with `e_Plx/Plx <= KIN3D_MAX_FRAC`. `shrink`: empirical-Bayes parallax pulled toward the cluster value. `raw`: `d = 1000/Plx` for all stars | controls the line-of-sight smearing of the 3-D map |
| `KIN3D_MAX_FRAC` | 0.05 | Largest allowed parallax fractional error | `precision` filter |
| `KIN3D_DEPTH_PC` | 5.0 pc | Assumed cluster depth that sets the prior width `tau = plx_c^2 x depth / 1000` (mas) | `shrink` weights: `1/e_Plx^2` (data) and `1/tau^2` (prior) |
| `KIN3D_FRAME` | `galactic` | `galactic` rotates ICRS (x, y, z) to Galactic axes with the Hipparcos matrix [16]; `icrs` keeps equatorial axes | 3-D coordinates |
| `KIN3D_CENTRE` | True | Subtract the median position | cluster at the origin |
| `MIN_STARS_3D` | 3 | Stars that must survive the distance filter | otherwise an error |

Outputs: x, y, z, distance and `e_Plx/Plx` per star; the residual proper motion `(pmRA - median, pmDE - median)` per star (drawn as unit arrows); and shape diagnostics:

| Diagnostic | Meaning |
|---|---|
| `elongation` | `sqrt(largest / mean of the other two eigenvalues)` of the position covariance; about 1 means round |
| `angle_to_los_deg` | Angle between the long axis and the line of sight (a needle along the line of sight signals parallax-error smearing) |
| `expected_err_sd_pc`, `implied_depth_pc` | Spread expected from parallax errors alone, and the remaining true depth `sqrt(observed^2 - expected^2)` |

### 3.5 Stage 5: mass function (`imf.py`)

| Parameter | Default | What it is | How it is used |
|---|---|---|---|
| `mass_min`, `mass_max` (`M_MIN_FLOOR`, `M_MAX_CAP`) | 0.08, 50 Msun | Mass range accepted | masses outside are dropped |
| `include_companions` | off | Add companion masses (`MC_CompMass`) to the sample. The box is unchecked by default: no module produces `MC_CompMass` yet | needs your own `MC_CompMass` column |
| `include_approx_ml` (`INCLUDE_APPROX_ML`) | off | Masses flagged `MC_Mass_src = approx_ML` (crude luminosity fallback, `Mass_quality = approx`) are excluded from the fit unless you opt in | quality gate |
| `method` (`METHOD`) | `mle` | `mle`: truncated power-law maximum likelihood on [m_min, m_max], Gamma = 1 - alpha, sigma from `N_BOOT` = 1000 seeded bootstrap resamples. `binned`: the old Poisson-weighted line fit | slope and error |
| `N_BINS` (`bins`) | auto | Number of log-mass bins; auto = `min(15, max(8, ceil(log2 N + 1)))` | equal-width bins in log10 M |
| `COMPLETE_FRAC` (`gmm.py`) | 0.9 | Retained fraction (quality-cut survivors / raw stars in the distance window) that still counts as complete. `G_lim` = faintest G reached before the fraction first drops below it; converted to a mass with `m_min = interp(G_lim - DM - A_G, M_G, Mass)` on the A_V = 0 isochrone | **measured** completeness limit (needs stage 1 + a fetched isochrone) |
| `FIT_MIN_LOGM` (`fit_min_logM`) | -0.7 (about 0.2 Msun) | Fallback completeness limit, used only when no measured limit exists. Sending `fit_min_logM` in the request (or moving the slider) overrides the measured value | bins below the limit are shown (open markers) but not fitted. **Do not tune it to favour a model** (circular) |
| `MIN_MASSES` | 20 | Minimum stars | otherwise an error |
| `dN/dlogM` | formula | Stars per dex of mass, `count / bin_width`; error `sqrt(count) / bin_width` (Poisson) | observed points |
| Models | fixed | Salpeter [15]: slope `1 - alpha` with `alpha = 2.35`. Kroupa [15]: broken power law (alpha 0.3 / 1.3 / 2.3, breaks 0.08 and 0.5 Msun). Chabrier [15]: log-normal below 1 Msun (`m_c = 0.08`, `sigma = 0.69` dex), power law above | compared with the data |
| Scaling | formula | Each model is scaled so its integral over the complete range equals the observed count there | no free parameter |
| Reduced chi^2 / Cash | formula | Model counts are integrated over each bin. `sum((obs - model)^2 / obs) / (n_bins - 1)` when every complete bin holds >= 5 stars, otherwise the Cash statistic `2 sum(mu - n + n ln(n / mu))` / (n_bins - 1) | compare models **relatively** |
| Slope fit | formula | `mle` (default): see `method` above; the binned fit is always reported alongside (binned R^2). The Salpeter slope is -1.35 | reported with error |
| `MF_LABEL_INITIAL`, `MF_LABEL_PRESENT`, `IMF_LOGAGE_SPLIT` | IMF / Present-Day Mass Function, 7.5 | Heading and plot title follow the age: log(age/yr) <= 7.5 -> "Initial Mass Function (IMF)", above -> "Present-Day Mass Function" (stellar evolution and dynamics have modified the observed mass function) | label |
| `COMPLETENESS_G50`, `COMPLETENESS_W` | 17.2, 0.65 | Illustrative logistic completeness curve; **an assumption**, shown (flagged) only when no measured curve exists | display only |

### 3.6 Web-service and session settings (`starnest/api.py`, entry point `main.py`)

| Parameter | Default | Meaning |
|---|---|---|
| `MAX_UPLOAD_MB` | 50 | Largest accepted upload |
| `SESSION_MAX`, `SESSION_TTL_S` | 20, 6 h | Visitor sessions kept in memory, and idle time before one is dropped; every visitor has a private state |
| `STARNEST_NO_PROXY`, `STARNEST_PROXY`, `STARNEST_IPV4_ONLY` | environment variables | Network switches for IRSA and CMD calls: ignore the system proxy, force a proxy URL, use IPv4 only |

---

## 4. External websites and services

| Service | Address | Used for | Needed when |
|---|---|---|---|
| Gaia DR3 via VizieR [1], [2], [6] | https://vizier.cds.unistra.fr (catalogue `I/355/gaiadr3`) | Source of the input CSV (cone search around the cluster) | before the app, by you |
| IRSA Dust (SFD) [7], [8] | https://irsa.ipac.caltech.edu/applications/DUST/ (`nph-dust`) | E(B-V) at the cluster position | stage 2 (internet) |
| CMD 3.8, PARSEC isochrones [10], [11] | https://stev.oapd.inaf.it/cgi-bin/cmd_3.8 | Isochrone in Gaia EDR3 photometry | stage 3 (internet) |
| Streamlit Community Cloud | https://streamlit.io | Hosting of the public app | online use |
| Front-end libraries (CDN) | Tailwind CSS (`cdn.tailwindcss.com`), Plotly.js 2.27 (`cdn.plot.ly`), GSAP 3.12 and Font Awesome 6.5 (`cdnjs.cloudflare.com`), Google Fonts (Inter) | Page style, plots, animation, icons | whenever the page loads (internet) |

IRSA and CMD are shared public services: the app caches results, spaces its requests, never retries a request that failed because of bad parameters, and sends nothing while a slider moves.

---

## 5. Repository layout

```
starnest/                       Python package (pyproject.toml): every module below lives here
  api.py                          FastAPI back end: sessions (locked), upload, routes, star types, masses
  schemas.py                      Pydantic request models (HTTP 400 naming the bad fields)
  gmm.py                          Stage 1: quality cuts + GMM membership (+ reliability diagnostics)
  extinction_isochrone.py         Stages 2-3: IRSA extinction + PARSEC isochrone + fit quality + likelihood fit
  parsec_grid.py                  Local PARSEC grid: read, bilinear interpolation, A_V shift
  build_grid.py                   One-off download of that grid from CMD 3.8 (python build_grid.py)
  kinematics.py                   Stage 4: 3-D map, proper-motion vectors, velocities
  imf.py                          Stage 5: mass function and model comparison
  binaries.py, dynamics.py        Binary-aware masses; mass segregation, tidal radius, relaxation time, orbit
  cli.py, report.py               Batch runner (python -m starnest run config.yaml); figure export and one-page report
main.py, gmm.py, imf.py, ...    Thin compatibility shims: `import gmm` / `uvicorn main:app` still work (same module objects)
index.html                      Web page (markup)
assets/css/styles.css           Page style
assets/js/app.js                Page logic, plots, in-browser fallback engine
assets/js/panels.js             Python-backend-only panels (reliability, cluster rule, likelihood fit, velocities, binaries, dynamics)
assets/js/footer.js             Two-mode footer + citation modal (BibTeX constants)
assets/js/tailwind-config.js    Tailwind theme
images/team/                    Team photos (file name = member name)
streamlit_app/                  Streamlit Community Cloud entry point + bridge
tests/                          pytest suite (runs offline)
examples/                       Batch config, literature table template and verify.py (section 12.4)
.github/workflows/ci.yml        pytest + ruff on every push
pyproject.toml                  Package metadata, extras (api, orbit, dust, dev), ruff / pytest settings
requirements.txt                Packages for the back end
CITATION.cff, .zenodo.json      Citation metadata (GitHub "Cite this repository") and Zenodo archive metadata
LICENSE                         MIT licence
```

All modules are free of user-interface code and can be imported on their own.

## 6. Using the modules from Python

```python
import gmm as GMM, kinematics as KIN, imf as IMF

raw   = GMM.load_raw("UPK545.csv")
clean = GMM.apply_quality_cuts(raw)                                   # section 3.1 cuts
res   = GMM.run_gmm(clean, distance_lit_pc=341.0, half_width_pc=195.5, optimizer="converged")
print(res.k, res.summary["n_members"])                                # chosen K, members

kin = KIN.compute_kinematics(res.members, mode="precision", max_frac=0.05)
KIN.save_kinematics(kin, "UPK545", out_dir="starnest_output")           # CSV files (index=False)

imf = IMF.compute_imf(masses)                                         # masses: array in Msun
print(imf.chi2)                                                       # reduced chi^2 per model
```

Stage 2-3 functions (`apply_extinction`, `get_parsec`, `fit_quality`) are in `extinction_isochrone.py` and need internet access.

## 7. Limitations and assumptions

- The membership is a three-dimensional Gaussian mixture; with a wide distance window the member list can contain field stars. Always inspect the distance distribution of the members.
- `ext SFD mean` is the total dust along the line of sight, so for a nearby cluster the extinction is an upper limit.
- The Gaia extinction ratios (0.789, 1.002, 0.589) are fixed and independent of stellar colour.
- Star-type rules are heuristics. The mass-function completeness limit is measured from the fraction of stars that survive the quality cuts versus G (`COMPLETE_FRAC`); `FIT_MIN_LOGM` (-0.7) is only the fallback when no measured limit exists. That fraction counts losses to the quality cuts, not Gaia's intrinsic incompleteness (crowding, scanning law).
- The mass function endpoint (`/compute-imf`) uses `MC_Mass` from `derive_columns`: a single-star main-sequence interpolation on the isochrone, so unresolved binaries are not modelled there and the masses depend on the age you set. The binary-aware masses and total cluster mass (`binaries.py`) are a separate endpoint and are not fed back into the mass-function slope.
- `binaries.py` assumes a fixed synthetic binary fraction `F_BIN` = 0.5 with companions drawn from U(0, M_primary), injects no photometric errors into the synthetic cluster, and uses the assumed `F_BIN` (not the fraction measured from the observed stars) in `M_T = M_obs + M_unseen (1 + F_bin q) + M_WD`. The unseen mass is a single power law extrapolated to 0.09 Msun, so `M_T` carries a large systematic error.
- Results depend on the Gaia zero-point (`PLX_ZP`); the default is 0.

## 8. Team

| Member | Role |
|---|---|
| Neelam Panwar | Research supervisor, Scientist "D", ARIES Nainital |
| Richa Kundu | Research supervisor, Assistant Professor, Miranda House (University of Delhi) |
| Nikhil Kaushik | Researcher, M.Sc. Physics, IIT Roorkee |
| Pryanshu Basur | Researcher, M.Sc. Physics, IIT Mandi |
| Monu | Researcher, M.Sc. Physics, Dept. of Physics and Astrophysics (University of Delhi) |
| Shruti Prasad | Researcher, Int. PhD Physics, IISc Bengaluru |

## 9. How to cite

**Software.** Use `CITATION.cff` (GitHub shows a "Cite this repository" button). A Zenodo DOI is added to it after the first archived release (`.zenodo.json` holds the Zenodo metadata).

**Paper.** Our photometric and kinematic study of six open clusters (Haffner 14, Melotte 111, NGC 2910, Ruprecht 147, UBC 52, UPK 545) is in preparation and has no journal reference or DOI yet. It is listed under `references` in `CITATION.cff` as an unpublished working paper and is not the preferred citation. When it is published, move it to `preferred-citation`.

**Data and models.** Cite the data and model papers of the stages you use (section 10): Gaia [1-5], PARSEC [10-12], the dust map [7-8] and the cluster reference COCD [20].

### Relation to Almeida et al. (2023)

The binary-aware mass method (`binaries.py`) is **inspired by** Almeida, Monteiro & Dias (2023, MNRAS 525, 2315; [13]) and shares its core idea: a synthetic cluster (10 000 stars, Chabrier IMF, 50 % binaries, companion mass uniform in 0 to M_primary) is matched star by star to the observed stars by the nearest Euclidean distance in (G, BP, RP), repeated many times, with the median as the mass and the standard deviation as the error. **StarNest is not a strict replication.** What differs, as implemented:

| Step | Almeida et al. 2023 | StarNest |
|---|---|---|
| Cluster parameters | Age, distance, extinction and [Fe/H] taken from the Dias et al. 2021 catalogue, updated with Gaia eDR3 | Membership by the GMM of this app; E(B-V) from the IRSA dust map; age from the user (or the likelihood fit) |
| Synthetic cluster | Gaia photometric errors are added to the synthetic stars | No photometric errors are added |
| Binary fraction in `M_T` | Measured from the observed stars | The assumed input `F_BIN` (0.5) |
| Present-day mass function | Two-part segmented line fitted to a Sturges-binned histogram | Truncated power-law maximum likelihood above a measured completeness limit (one slope, bootstrap error); the binned fit is kept for comparison |
| Cluster mass | Two variants (integrated mass function; detailed, with separate single / primary / secondary mass functions); flat 20 % error | One formula with a single power-law extrapolation; no separate primary / secondary mass functions; error is not propagated |
| White-dwarf mass | Cummings et al. 2018 [21] relation for progenitors of 0.85 to 7.50 Msun | Cummings et al. 2018 [21] piecewise relation coded for 0.87 to 8.2 Msun |
| Mass segregation | Ratio of mean radii of stars above and below 1 Msun, plus a Kolmogorov-Smirnov test | Minimum-spanning-tree ratio Lambda_MSR (Allison et al. 2009) for the N most massive stars |

Results from StarNest are therefore not expected to reproduce numbers published with the OCMass app, and should not be presented as such.

## 10. References

1. Gaia Collaboration, Prusti T. et al., 2016, A&A 595, A1 (the Gaia mission).
2. Gaia Collaboration, Vallenari A. et al., 2023, A&A 674, A1 (Gaia DR3).
3. Lindegren L. et al., 2021, A&A 649, A2 (Gaia EDR3 astrometric solution; RUWE).
4. Lindegren L. et al., 2021, A&A 649, A4 (parallax bias and zero-point).
5. Riello M. et al., 2021, A&A 649, A3 (Gaia EDR3 photometry).
6. Ochsenbein F., Bauer P., Marcout J., 2000, A&AS 143, 23 (VizieR).
7. Schlegel D. J., Finkbeiner D. P., Davis M., 1998, ApJ 500, 525 (SFD dust map).
8. Schlafly E. F., Finkbeiner D. P., 2011, ApJ 737, 103 (recalibration, `A_V = 2.742 E(B-V)` for SFD).
9. Cardelli J. A., Clayton G. C., Mathis J. S., 1989, ApJ 345, 245 (extinction law).
10. Bressan A. et al., 2012, MNRAS 427, 127 (PARSEC).
11. Chen Y. et al., 2014, MNRAS 444, 2525; Tang J. et al., 2014, MNRAS 445, 4287; Chen Y. et al., 2015, MNRAS 452, 1068; Marigo P. et al., 2017, ApJ 835, 77 (PARSEC extensions).
12. Chen Y. et al., 2019, A&A 632, A105 (YBC bolometric corrections).
13. Almeida A., Monteiro H., Dias W. S., 2023, MNRAS 525, 2315 (arXiv:2307.15182), "Revisiting the mass of open clusters with Gaia data" (the OCMass method and web app, https://ocmass.streamlit.app). StarNest is inspired by, not a replication of, this method (section 9).
14. Dempster A. P., Laird N. M., Rubin D. B., 1977, JRSS B 39, 1 (EM algorithm); Schwarz G., 1978, Ann. Statist. 6, 461 (BIC); Pedregosa F. et al., 2011, JMLR 12, 2825 (scikit-learn).
15. Salpeter E. E., 1955, ApJ 121, 161; Kroupa P., 2001, MNRAS 322, 231; Chabrier G., 2003, PASP 115, 763 (initial mass functions).
16. ESA, 1997, The Hipparcos and Tycho Catalogues, ESA SP-1200 (ICRS to Galactic rotation matrix).
17. Bohlin R. C., Hubeny I., Rauch T., 2020, AJ 160, 21 (CALSPEC / Vega reference spectrum).
18. Astropy Collaboration, 2022, ApJ 935, 167; Ginsburg A. et al., 2019, AJ 157, 98 (astroquery).
19. Sim G. et al., 2019, J. Korean Astron. Soc. 52, 145 (the UPK cluster catalogue; source of the example cluster).
20. Kharchenko N. V., Piskunov A. E., Roser S., Schilbach E., Scholz R.-D., 2005, A&A 438, 1163 (COCD, the Catalogue of Open Cluster Data; the "Cluster parameters: COCD" label in the app footer). Note that COCD (2005) predates the UPK and UBC cluster discoveries, so it cannot be the source of parameters for those clusters.
21. Cummings J. D., Kalirai J. S., Tremblay P.-E., Ramirez-Ruiz E., Choi J., 2018, ApJ 866, 21 (initial-final mass relation).

---

## 11. Run the app offline on your own machine

Use this guide if the hosted app is unavailable or too slow. "Offline" here means **on your own computer**: the page still downloads its display libraries (Tailwind, Plotly, GSAP, fonts) from the internet, and stages 2-3 contact IRSA and CMD, so an internet connection is needed for those. GMM, kinematics and the mass function run entirely on your machine.

### 11.1 Requirements

| Item | Needed |
|---|---|
| Python | 3.11 or newer (developed and tested with 3.13) |
| Packages | listed in `requirements.txt` |
| Browser | Any current Chrome, Edge, Firefox or Safari |
| Data | your Gaia CSV (section 3.0) |

### 11.2 Install (once)

**Linux / macOS**

```bash
git clone https://github.com/NIKHILKAUSHIK9/StarNest.git
cd StarNest
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

**Windows (PowerShell)**

```powershell
git clone https://github.com/NIKHILKAUSHIK9/StarNest.git
cd StarNest
py -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

No git? Download the repository as a ZIP from GitHub, extract it, and open a terminal in the extracted folder.

### 11.3 Start the app

```bash
uvicorn main:app --port 8000
```

Then open **http://localhost:8000** in your browser. The badge at the bottom of the page must read **"Engine: Python backend (PARSEC + IRSA)"**. Stop the app with `Ctrl + C`.

### 11.4 Workflow in the browser

| Step | Tab | Action | Result |
|---|---|---|---|
| 1 | GMM | Upload the CSV, set the quality cuts, `Dist_Lit` and `Half_Width`, press **Run GMM Analysis** | members, K, BIC table, sky / proper-motion / parallax plots |
| 2 | GMM | **Download Filtered Members CSV** (optional) | table with `PMemb`, star type, `M_G_corr`, `bp_rp_corr` |
| 3 | Isochrone | Enter `logAge`, `[M/H]`, E(B-V) mode and extinction handling, run | IRSA E(B-V), PARSEC isochrone, CMD overlay, fit quality |
| 4 | Kinematics | Pick the 3-D mode and frame | 3-D map and proper-motion vectors |
| 5 | IMF | Set mass range, bins and `fit_min_logM` | mass function with Salpeter, Kroupa, Chabrier and reduced chi^2 |

The workflow tracker on the Isochrone tab shows each network step (IRSA, PARSEC) in green or red with the server's message.

### 11.5 Other ways to start

| Goal | Command | Notes |
|---|---|---|
| Different port | `uvicorn main:app --port 8010` | open `http://localhost:8010` |
| Auto-reload while editing code | `uvicorn main:app --reload` | development only |
| Streamlit front end on your machine | `pip install -r streamlit_app/requirements.txt` then `streamlit run streamlit_app/streamlit_app.py` | same page, backend runs inside Streamlit |
| Page only, no server | open `index.html` in the browser | GMM, kinematics and IMF run in the browser; IRSA and PARSEC are unavailable (badge: "Engine: in-browser") |
| Use the back end with a page opened from a file | start `uvicorn main:app --port 8000`, then open `index.html` | the page finds `http://localhost:8000/api` by itself |

### 11.6 Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError` | Packages not installed or the virtual environment is not active | activate `.venv`, run `pip install -r requirements.txt` |
| `Address already in use` | Port 8000 is taken | use `--port 8010` |
| Badge says "in-browser" | Back end not reachable | check the terminal for errors; open `http://localhost:8000/api/health` (should return JSON) |
| Page looks unstyled or plots are missing | No internet for the CDN libraries | connect to the internet and reload |
| IRSA or PARSEC step red, "connect timeout" | Network, proxy or IPv6 problem | set an environment variable before starting: `STARNEST_NO_PROXY=1` (ignore the system proxy), `STARNEST_PROXY=http://host:port` (force a proxy) or `STARNEST_IPV4_ONLY=1`; on Windows PowerShell use `$env:STARNEST_NO_PROXY="1"` |
| PARSEC error "did not return an output file" | CMD server busy or rejected the parameters | read the message in the tracker; check `logAge` (6.6 to 10.13) and `[M/H]`; retry after a minute |
| "Missing required columns" | The CSV lacks columns from section 3.0 | export all of them from VizieR |
| Upload refused | File larger than 50 MB | raise `MAX_UPLOAD_MB` in `main.py` |

### 11.7 Output files

Downloads from the page are saved by your browser. The Python helpers (`save_members`, `save_kinematics`, `save_imf`, `save_outputs`) write to `./starnest_output/` (set by `OUT_DIR` in each module), always without the row index.

## 12. Local PARSEC grid, batch runs and tests (Phase 2)

### 12.1 Local PARSEC grid (no CMD request at run time)

```
python build_grid.py            # once, needs internet: 11 CMD requests (+3 checks), resumable with --resume
```

Writes `data/parsec_grid.parquet`: logAge 6.6-10.1 (step 0.05) x [M/H] -0.5..+0.5 (step 0.1), A_V = 0, Gaia EDR3 (= DR3) Vega system. `get_parsec()` then blends the 4 surrounding isochrones bilinearly (per evolutionary phase, equal normalised arc length in the HR diagram) and applies A_V analytically with the A_band/A_V coefficients **measured from CMD at build time** (stored in the file). The reddening self-check is kept, but for a grid isochrone it only tests that the shift lands in the expected range: it is not independent. The build also compares the interpolation with three off-grid CMD isochrones and stores the errors in the file (`interp_check`).

| Setting | Values | Effect |
|---|---|---|
| `STARNEST_PARSEC_SOURCE` | `auto` (default), `grid`, `cmd` | `auto`: grid if it covers (logAge, [M/H]), else CMD; `grid`: never contact CMD; `cmd`: ignore the grid |
| `STARNEST_PARSEC_GRID` | path | grid file (default `data/parsec_grid.parquet`) |
| `STARNEST_ALLOW_MOCK_GRID` | `1` | tests only: a grid whose metadata says `source: MOCK` is refused otherwise |

`.gitignore` keeps `data/` out of the repository except `data/parsec_grid.parquet`, so you can commit the finished grid (needed if the Streamlit app should not call CMD).

### 12.2 Batch runs

```
pip install -e ".[api]"
python -m starnest run config.yaml [--out DIR] [--only NAME ...]
```

`config.yaml` has `output_dir`, `steps`, `defaults` (one dict per step: gmm, cmd, imf, binaries, velocities, dynamics, isofit) and `clusters` (name, csv, optional per-step overrides). See `examples/config.example.yaml`. Per cluster it writes the member / extinction / isochrone / IMF / binary CSVs (`index=False`), PDF + PNG figures and `<name>_report.pdf` (one A4 page, warnings included). Per run it writes `summary.csv` (N, d, E(B-V), logAge, Gamma, sigma_Gamma, MSR, v_out, sigma_1D plus status, isochrone source, warnings) and `run_manifest.json` (parameters, package versions, SHA-256 of every input and of the grid). The CLI calls the same FastAPI app in-process, so its numbers are the web page's numbers. Exit code 0 = all ok, 1 = a cluster or step failed (others still written), 2 = config error.

### 12.3 Tests and style

```
pip install -e ".[api]" pytest ruff
pytest -q          # offline; galpy tests skip when galpy is absent
ruff check .
```

GitHub Actions (`.github/workflows/ci.yml`) runs both on every push.

### 12.4 Examples and cross-verification against the literature

`examples/` holds the batch config template, a literature-table template and `verify.py`. No expected outputs are committed: age, E(B-V), masses and the mass-function slope depend on the isochrone, so they are only meaningful when produced with the real PARSEC grid (section 12.1) on your machine. See `examples/README.md` for the run order.

```
python -m starnest run my_config.yaml --out starnest_output
python examples/verify.py starnest_output/summary.csv my_literature_table.csv
```

`verify.py` compares your run with a literature table **you** supply (it never uses StarNest's own earlier output as the reference), using tolerances fixed in the `GLOBAL CONFIGURATION` block at the top of the script; it exits 1 if any checked value is outside its tolerance or missing.
