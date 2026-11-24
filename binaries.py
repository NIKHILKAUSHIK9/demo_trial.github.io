"""Compatibility shim: `import binaries` still works and returns the packaged module starnest.binaries (same object)."""
import sys
import starnest.binaries as _m
sys.modules[__name__] = _m
