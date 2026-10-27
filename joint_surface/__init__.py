"""Continuous-delta feature-only selection, provisional farm branch.

Numerical dependencies load only when a numerical export is requested. Campaign
planning and browser-viewer helpers can therefore run without NumPy or SciPy.
"""
from importlib import import_module

__all__ = ['DEFAULT_CONFIG', 'FAMILIES', 'fit_axis', 'make_frame', 'predict_axis',
           'rotate', 'surface_derivative', 'surface_values']


def __getattr__(name):
    if name not in __all__:
        raise AttributeError("module %r has no attribute %r" % (__name__, name))
    value = getattr(import_module('.engine', __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
