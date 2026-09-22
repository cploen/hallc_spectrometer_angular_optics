"""Separate source-derived constructions; never used to define reference rays."""
from geometry import Ray, hms_default_mispointing
import math


def fit_target_ray(angle_deg, zfoil_cm, xhole_cm, yhole_cm, reacty_cm, xbpm_cm, L_cm):
    """Current spectrometer_config.h HMS branch, preserved including mismatches.

    xbpm is the RAW exported raster branch. Its mapping to lab X is not assumed
    here. reactx is unused by this fit branch, which is part of the audit.
    """
    a = math.radians(angle_deg)
    c, s = math.cos(a), math.sin(a)
    mx, my, _ = hms_default_mispointing(angle_deg)
    b = -xbpm_cm
    z0 = zfoil_cm*c
    yc = zfoil_cm*s+b*c-my
    q = (yhole_cm-yc)/(L_cm-z0)
    p = xhole_cm/(L_cm-z0)
    return Ray(-reacty_cm-mx-p*z0,
               zfoil_cm*(s-q*c)+b*(c+q*s)-my, p, q)


def extcor_xtar_update(reacty_cm, reactz_cm, pointing_x_cm, xptar, theta_rad):
    """ExtTarCor expression recorded in prior pinned audit; isolated for testing.

    No optics polynomial, iteration, beam provider, or actual replay is inferred.
    The full updater must call the real matrix before repeating this expression.
    """
    return -reacty_cm-pointing_x_cm-xptar*reactz_cm*math.cos(theta_rad)


def hcana_geometry_step(lab_vertex_cm, pointing_x_cm, theta_rad,
                        p_old, p_new, ytar_matrix_m, q_new):
    """Trace the HCANA geometric interface with explicit synthetic matrix outputs.

    Source: sieve_afterburner.py:47-57 (local mirror), and the prior pinned
    upstream audit. No matrix is evaluated. p_old constructs xt; p_new projects
    it. Supplying reference yt/100 and q tests units/projection only in y.
    """
    xt=extcor_xtar_update(lab_vertex_cm[1],lab_vertex_cm[2],pointing_x_cm,p_old,theta_rad)
    return Ray(xt,100*ytar_matrix_m,p_new,q_new)
