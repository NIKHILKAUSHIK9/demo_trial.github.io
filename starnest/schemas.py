"""
schemas.py -- Pydantic request models for the StarNest API (Phase 2.1).

TERMINOLOGY
  Request model     : typed description of one JSON body. Pydantic checks types / ranges BEFORE the route runs.
  extra='allow'     : unknown keys are kept (the in-browser client sends a few legacy keys such as dist_mod, av) instead of rejected.
  exclude_unset     : routes read the body with `req.model_dump(exclude_unset=True)`, so a key the client did NOT send stays absent and the
                      route's `p.get(key, default)` logic keeps its exact pre-2.1 meaning (an explicit JSON null stays present as None).
  Literal           : closed set of allowed strings (e.g. dust_source).

WORKFLOW
  1) FastAPI parses the JSON body into the model; 2) a type / range violation raises RequestValidationError;
  3) main.py turns it into HTTP 400 whose `detail` names every offending field ("field: reason; ...") -- the front end already shows
  string details verbatim; 4) the route receives a plain dict identical in shape to the old `p: dict`.
"""
from __future__ import annotations

from typing import Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field

# ===================== GLOBAL CONFIGURATION =====================
GMM_OPTIMIZERS = ('colab', 'converged')                 # gmm.OPTIMIZER
CLUSTER_RULES = ('smallest_det', 'closest_to_literature')   # gmm.CLUSTER_RULE
EBV_MODES = ('center', 'fixed', 'per_star')             # /compute-cmd E(B-V) source
ISO_EXTINCTION = ('reddened', 'intrinsic')              # where extinction is applied once
KIN3D_MODES = ('precision', 'shrink', 'raw')            # kinematics.KIN3D_MODE
KIN3D_FRAMES = ('galactic', 'icrs')                     # kinematics.KIN3D_FRAME
IMF_METHODS = ('mle', 'binned')                         # imf.METHOD
PLX_MODES = ('star', 'cluster')                         # kinematics.VEL_PLX_MODE
DUST_SOURCES = ('sfd', 'bayestar', 'free')              # extinction_isochrone.DUST_SOURCE
# =================================================================


class _Req(BaseModel):
    """Base: tolerate unknown keys; strings like "0.5" are coerced to numbers (lax mode), nulls are allowed on every optional field."""
    model_config = ConfigDict(extra='allow')


class RunGmmReq(_Req):
    filename: str = Field(min_length=1)
    ruwe_max: Optional[float] = Field(None, gt=0)
    epmra_max: Optional[float] = Field(None, gt=0)
    epmde_max: Optional[float] = Field(None, gt=0)
    plx_err_min: Optional[float] = Field(None, ge=0)
    gofal_max: Optional[float] = Field(None, gt=0)
    nper_min: Optional[float] = Field(None, ge=0)
    distance_lit: Optional[float] = Field(None, gt=0)
    half_width: Optional[float] = Field(None, gt=0)
    cutoff: Optional[float] = Field(None, gt=0, lt=1)
    k: Optional[int] = Field(None, ge=2, le=12)
    plx_zp: Optional[float] = Field(None, ge=-1, le=1)
    gmm_opt: Optional[Literal['colab', 'converged']] = None
    cluster_rule: Optional[Literal['smallest_det', 'closest_to_literature']] = None
    lit_pm: Optional[list[float]] = Field(None, min_length=2, max_length=2)
    lit_plx: Optional[float] = Field(None, gt=0)
    ebv_ui: Optional[float] = None
    complete_frac: Optional[float] = Field(None, gt=0, lt=1)
    complete_prefilter: Optional[bool] = None


class ComputeCmdReq(_Req):
    members_filename: str = Field(min_length=1)
    age: float = Field(ge=0, le=11)                      # log10(age / yr); the PARSEC floor message is raised by the route
    metallicity: Optional[float] = Field(None, ge=-2.5, le=1.0)
    ebv_mode: Optional[Literal['center', 'fixed', 'per_star']] = None
    ebv: Optional[float] = Field(None, ge=0, le=20)
    iso_extinction: Optional[Literal['reddened', 'intrinsic']] = None
    allow_network: Optional[bool] = None
    plx_zp: Optional[float] = Field(None, ge=-1, le=1)
    use_cluster_dm: Optional[bool] = None


class KinematicsReq(_Req):
    filename: str = Field(min_length=1)
    kin3d_mode: Optional[Literal['precision', 'shrink', 'raw']] = None
    kin3d_max_frac: Optional[float] = Field(None, ge=0)
    kin3d_depth: Optional[float] = Field(None, ge=0)
    kin3d_frame: Optional[Literal['galactic', 'icrs']] = None
    kin3d_centre: Optional[bool] = None


class ImfReq(_Req):
    filename: str = Field(min_length=1)
    mass_min: Optional[float] = Field(None, ge=0)
    mass_max: Optional[float] = Field(None, gt=0)
    include_approx_ml: Optional[bool] = None
    include_companions: Optional[bool] = None
    fit_min_logM: Optional[float] = Field(None, ge=-3, le=3)
    method: Optional[Literal['mle', 'binned']] = None
    bins: Optional[int] = Field(None, ge=2, le=200)
    log_age: Optional[float] = Field(None, ge=0, le=11)


class GmmDiagReq(_Req):
    n_resample: Optional[int] = Field(None, ge=2, le=5000)
    seed: Optional[int] = None
    min_cluster_size: Optional[int] = Field(None, ge=2)
    cluster_rule: Optional[Literal['smallest_det', 'closest_to_literature']] = None


class BinariesReq(_Req):
    filename: str = Field(min_length=1)
    f_bin: Optional[float] = Field(None, ge=0, le=1)
    n_synth: Optional[int] = Field(None, ge=100, le=1_000_000)
    n_mc: Optional[int] = Field(None, ge=1, le=2000)
    seed: Optional[int] = None


class VelocitiesReq(_Req):
    filename: str = Field(min_length=1)
    mass: Optional[float] = Field(None, ge=0)
    eta: Optional[float] = Field(None, gt=0)
    plx_mode: Optional[Literal['star', 'cluster']] = None
    epsi_max: Optional[float] = Field(None, gt=0)
    perspective: Optional[bool] = None
    n_boot: Optional[int] = Field(None, ge=0, le=20000)


class DynamicsReq(_Req):
    filename: str = Field(min_length=1)
    mass: Optional[float] = Field(None, ge=0)
    n_massive: Optional[Union[int, list[int]]] = None      # one subset size or a list (dynamics.N_MASSIVE)
    n_random: Optional[int] = Field(None, ge=1, le=100000)
    sigma_1d: Optional[float] = Field(None, ge=0)
    n_total: Optional[int] = Field(None, ge=1)
    orbit: Optional[bool] = None


class IsoFitReq(_Req):
    filename: str = Field(min_length=1)
    dust_source: Optional[Literal['sfd', 'bayestar', 'free']] = None
    fit_binaries: Optional[bool] = None
