"""Entry point kept for `uvicorn main:app`, the Streamlit bridge and the tests: it IS starnest.api (same module object)."""
import sys
import starnest.api as _api
sys.modules[__name__] = _api
