
import numpy as np
import unittest

from openquake.globalm.grid.refgrid.triangular_grid import (ll2xyz, xyz2ll,
                                                            get_midp)


class TestCoordinateConversion(unittest.TestCase):

    def test_fwd_inv_00(self):
        ilo = 30.
        ila = 0.
        x, y, z = ll2xyz(ilo, ila)
        lo, la = xyz2ll(x, y, z)
        self.assertAlmostEqual(ilo, lo)
        self.assertAlmostEqual(ila, la)

    def test_fwd_inv_01(self):
        ilo = 30.
        ila = 60.
        x, y, z = ll2xyz(ilo, ila)
        lo, la = xyz2ll(x, y, z)
        self.assertAlmostEqual(ilo, lo)
        self.assertAlmostEqual(ila, la)

    def test_fwd_inv_02(self):
        ilo = -30.
        ila = 60.
        x, y, z = ll2xyz(ilo, ila)
        lo, la = xyz2ll(x, y, z)
        self.assertAlmostEqual(ilo, lo)
        self.assertAlmostEqual(ila, la)

    def test_fwd_01(self):
        ilo = 30.
        ila = 60.
        x, y, z = ll2xyz(ilo, ila)
        ho = np.cos(np.deg2rad(ila))
        self.assertAlmostEqual(x, np.cos(np.deg2rad(ilo))*ho)
        self.assertAlmostEqual(y, np.sin(np.deg2rad(ilo))*ho)
        self.assertAlmostEqual(z, np.sin(np.deg2rad(ila)))

    def test_fwd_02(self):
        ilo = 0.
        ila = 60.
        x, y, z = ll2xyz(ilo, ila)
        self.assertAlmostEqual(x, 0.5)
        self.assertAlmostEqual(y, 0)
        self.assertAlmostEqual(z, np.sin(np.deg2rad(ila)))

    def test_fwd_03(self):
        ilo = 0.
        ila = -60.
        x, y, z = ll2xyz(ilo, ila)
        self.assertAlmostEqual(x, 0.5)
        self.assertAlmostEqual(y, 0)
        self.assertAlmostEqual(z, np.sin(np.deg2rad(ila)))

    def test_fwd_04(self):
        ilo = -45.
        ila = 0.
        x, y, z = ll2xyz(ilo, ila)
        self.assertAlmostEqual(x, 1/np.sqrt(2))
        self.assertAlmostEqual(y, -1/np.sqrt(2))
        self.assertAlmostEqual(z, 0)

    def test_fwd_05(self):
        ilo = 144.0
        ila = 26.56505
        x, y, z = ll2xyz(ilo, ila)
        ho = np.cos(np.deg2rad(ila))
        ve = np.sin(np.deg2rad(ila))
        self.assertAlmostEqual(x, -ho * np.cos(np.deg2rad(180-ilo)))
        self.assertAlmostEqual(y, ho*np.sin(np.deg2rad(180-ilo)))
        self.assertAlmostEqual(z, ve)


class TestMidPointCalculation(unittest.TestCase):

    def test_mid_point_01(self):
        lo1 = 90.0
        la1 = 0.0
        lo2 = 90.0
        la2 = 90.0
        mlo, mla = get_midp(lo1, la1, lo2, la2)
        self.assertAlmostEqual(90, mlo)
        self.assertAlmostEqual(45, mla)

    def test_mid_point_02(self):
        lo1 = 0.0
        la1 = 90.0
        lo2 = 0.0
        la2 = 26.56505
        mlo, mla = get_midp(lo1, la1, lo2, la2)
        self.assertAlmostEqual(0, mlo)
        self.assertAlmostEqual(58.282525, mla)

    def test_mid_point_03(self):
        lo1 = 0.0
        la1 = 90.0
        lo2 = 144.0
        la2 = 26.56505
        mlo, mla = get_midp(lo1, la1, lo2, la2)
        self.assertAlmostEqual(144, mlo)
        self.assertAlmostEqual(58.282525, mla)
