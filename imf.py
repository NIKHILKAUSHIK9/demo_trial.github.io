"""Compatibility shim: `import imf` still works and returns the packaged module starnest.imf (same object)."""
import sys
import starnest.imf as _m
sys.modules[__name__] = _m
