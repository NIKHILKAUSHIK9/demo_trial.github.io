"""Compatibility shim: `python build_grid.py ...` runs starnest.build_grid; `import build_grid` returns that module."""
import sys

if __name__ == '__main__':
    from starnest.build_grid import main
    main()
else:
    import starnest.build_grid as _m
    sys.modules[__name__] = _m
