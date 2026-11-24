"""Compatibility shim: `import kinematics` still works and returns the packaged module starnest.kinematics (same object)."""
import sys
import starnest.kinematics as _m
sys.modules[__name__] = _m
