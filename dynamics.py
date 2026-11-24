"""Compatibility shim: `import dynamics` still works and returns the packaged module starnest.dynamics (same object)."""
import sys
import starnest.dynamics as _m
sys.modules[__name__] = _m
