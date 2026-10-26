"""Provisional geometric alternatives to the campaign X/Y band selectors.

These numerical engines depend only on NumPy and SciPy, not on ROOT.
"""
from .x_discovery import fit_x
from .x_geometry import global_centers, global_coordinates
from .occupancy import detect_x, predict_x
from .y_geometry import fit_y, predict_y

__all__ = ['fit_x', 'detect_x', 'predict_x', 'global_centers',
           'global_coordinates', 'fit_y', 'predict_y']
