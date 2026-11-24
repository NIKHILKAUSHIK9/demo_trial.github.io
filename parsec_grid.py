"""Compatibility shim: `import parsec_grid` still works and returns the packaged module starnest.parsec_grid (same object)."""
import sys
import starnest.parsec_grid as _m
sys.modules[__name__] = _m
