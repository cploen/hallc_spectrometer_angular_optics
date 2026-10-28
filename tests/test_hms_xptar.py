"""Vertex-to-hole closure and rejection of obsolete HMS fitting values."""
import math
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from spectrometer_config import hms_xptar_truth


class HMSRayTests(unittest.TestCase):
    def test_vertex_and_sieve_closure(self):
        for angle in (12.49, -12.49, 55.):
            degree=min(abs(angle),50.)
            xmis=.1*(2.37-.086*degree+.0012*degree**2)
            for zfoil in (-8.,-3.,0.,3.,8.):
                zvertex=zfoil*math.cos(math.radians(angle))
                for xs in (-10.16,0.,10.16):
                    reacty=np.array([-.1,0.,.08])
                    xptarT=hms_xptar_truth(xs,zfoil,angle,reacty)
                    xvertex=-reacty-xmis
                    xtarT=xvertex-xptarT*zvertex
                    np.testing.assert_allclose(xtarT+168*xptarT,xs,atol=1e-13)
                    np.testing.assert_allclose(xtarT+zvertex*xptarT,xvertex,atol=1e-13)

    def test_central_hole_residual_change(self):
        reacty=np.array([-.1,0.,.08])
        corrected=hms_xptar_truth(0.,0.,12.49,reacty)
        np.testing.assert_allclose(1000*np.diff(corrected)/np.diff(reacty),1000/168)


if __name__=='__main__':
    unittest.main()
