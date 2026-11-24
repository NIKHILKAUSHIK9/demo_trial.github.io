"""Compatibility shim: `import schemas` still works and returns the packaged module starnest.schemas (same object)."""
import sys
import starnest.schemas as _m
sys.modules[__name__] = _m
