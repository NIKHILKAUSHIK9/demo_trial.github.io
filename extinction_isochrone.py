"""Compatibility shim: `import extinction_isochrone` still works and returns the packaged module starnest.extinction_isochrone (same object)."""
import sys
import starnest.extinction_isochrone as _m
sys.modules[__name__] = _m
