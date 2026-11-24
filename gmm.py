"""Compatibility shim: `import gmm` still works and returns the packaged module starnest.gmm (same object)."""
import sys
import starnest.gmm as _m
sys.modules[__name__] = _m
