import os
import unittest

import numpy as np
from osgeo import gdal

from openquake.sep.utils import (
    make_2d_array_strides,
    rolling_array_operation,
    rolling_raster_operation,
    relief,
    make_local_relief_raster,
)


BASE_DATA_PATH = os.path.join(os.path.dirname(__file__), "data")
test_dem = os.path.join(BASE_DATA_PATH, "dem_small.tif")
test_relief = os.path.join(BASE_DATA_PATH, "relief_out.tif")


class test_array_funcs_super_simple(unittest.TestCase):
    def setUp(self):
        self.array = np.arange(25).reshape((5, 5))

    def test_make_2d_array_strides(self):
        strides = make_2d_array_strides(self.array, 1)

        strides_result = np.array(
            [
                [
                    [
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        0.0,
                        1.0,
                        np.nan,
                        5.0,
                        6.0,
                    ],
                    [np.nan, np.nan, np.nan, 0.0, 1.0, 2.0, 5.0, 6.0, 7.0],
                    [np.nan, np.nan, np.nan, 1.0, 2.0, 3.0, 6.0, 7.0, 8.0],
                    [np.nan, np.nan, np.nan, 2.0, 3.0, 4.0, 7.0, 8.0, 9.0],
                    [
                        np.nan,
                        np.nan,
                        np.nan,
                        3.0,
                        4.0,
                        np.nan,
                        8.0,
                        9.0,
                        np.nan,
                    ],
                ],
                [
                    [np.nan, 0.0, 1.0, np.nan, 5.0, 6.0, np.nan, 10.0, 11.0],
                    [0.0, 1.0, 2.0, 5.0, 6.0, 7.0, 10.0, 11.0, 12.0],
                    [1.0, 2.0, 3.0, 6.0, 7.0, 8.0, 11.0, 12.0, 13.0],
                    [2.0, 3.0, 4.0, 7.0, 8.0, 9.0, 12.0, 13.0, 14.0],
                    [3.0, 4.0, np.nan, 8.0, 9.0, np.nan, 13.0, 14.0, np.nan],
                ],
                [
                    [np.nan, 5.0, 6.0, np.nan, 10.0, 11.0, np.nan, 15.0, 16.0],
                    [5.0, 6.0, 7.0, 10.0, 11.0, 12.0, 15.0, 16.0, 17.0],
                    [6.0, 7.0, 8.0, 11.0, 12.0, 13.0, 16.0, 17.0, 18.0],
                    [7.0, 8.0, 9.0, 12.0, 13.0, 14.0, 17.0, 18.0, 19.0],
                    [8.0, 9.0, np.nan, 13.0, 14.0, np.nan, 18.0, 19.0, np.nan],
                ],
                [
                    [
                        np.nan,
                        10.0,
                        11.0,
                        np.nan,
                        15.0,
                        16.0,
                        np.nan,
                        20.0,
                        21.0,
                    ],
                    [10.0, 11.0, 12.0, 15.0, 16.0, 17.0, 20.0, 21.0, 22.0],
                    [11.0, 12.0, 13.0, 16.0, 17.0, 18.0, 21.0, 22.0, 23.0],
                    [12.0, 13.0, 14.0, 17.0, 18.0, 19.0, 22.0, 23.0, 24.0],
                    [
                        13.0,
                        14.0,
                        np.nan,
                        18.0,
                        19.0,
                        np.nan,
                        23.0,
                        24.0,
                        np.nan,
                    ],
                ],
                [
                    [
                        np.nan,
                        15.0,
                        16.0,
                        np.nan,
                        20.0,
                        21.0,
                        np.nan,
                        np.nan,
                        np.nan,
                    ],
                    [
                        15.0,
                        16.0,
                        17.0,
                        20.0,
                        21.0,
                        22.0,
                        np.nan,
                        np.nan,
                        np.nan,
                    ],
                    [
                        16.0,
                        17.0,
                        18.0,
                        21.0,
                        22.0,
                        23.0,
                        np.nan,
                        np.nan,
                        np.nan,
                    ],
                    [
                        17.0,
                        18.0,
                        19.0,
                        22.0,
                        23.0,
                        24.0,
                        np.nan,
                        np.nan,
                        np.nan,
                    ],
                    [
                        18.0,
                        19.0,
                        np.nan,
                        23.0,
                        24.0,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                    ],
                ],
            ]
        )

        np.testing.assert_array_almost_equal(strides, strides_result)

    def test_relief_w_nan(self):
        strides = make_2d_array_strides(self.array, 1)
        relief_result = np.array(
            [
                [np.nan, np.nan, np.nan, np.nan, np.nan],
                [np.nan, 12.0, 12.0, 12.0, np.nan],
                [np.nan, 12.0, 12.0, 12.0, np.nan],
                [np.nan, 12.0, 12.0, 12.0, np.nan],
                [np.nan, np.nan, np.nan, np.nan, np.nan],
            ]
        )

        np.testing.assert_array_almost_equal(relief_result, relief(strides))

    def test_relief_wo_nan(self):
        strides = make_2d_array_strides(self.array, 1)
        strides = np.ma.masked_array(strides, mask=np.isnan(strides))

        relief_result = np.array(
            [
                [6.0, 7.0, 7.0, 7.0, 6.0],
                [11.0, 12.0, 12.0, 12.0, 11.0],
                [11.0, 12.0, 12.0, 12.0, 11.0],
                [11.0, 12.0, 12.0, 12.0, 11.0],
                [6.0, 7.0, 7.0, 7.0, 6.0],
            ]
        )

        np.testing.assert_array_almost_equal(relief_result, relief(strides))

    def test_rolling_array_operation(self):
        relief_result = np.array(
            [
                [6.0, 7.0, 7.0, 7.0, 6.0],
                [11.0, 12.0, 12.0, 12.0, 11.0],
                [11.0, 12.0, 12.0, 12.0, 11.0],
                [11.0, 12.0, 12.0, 12.0, 11.0],
                [6.0, 7.0, 7.0, 7.0, 6.0],
            ]
        )

        relief_res = rolling_array_operation(
            self.array, relief, window_size=3, trim=False
        )

        np.testing.assert_array_almost_equal(relief_res, relief_result)


class test_make_local_relief_raster(unittest.TestCase):
    def setUp(self):
        self.test_relief_raster = gdal.Open(test_relief)
        self.lrr = make_local_relief_raster(
            test_dem, 5, outfile=None, write=False, trim=False
        )

    def test_make_local_relief_raster_geo_transform(self):

        np.testing.assert_allclose(
            self.lrr.GetGeoTransform(),
            self.test_relief_raster.GetGeoTransform(),
        )

    def test_make_2d_local_relief_raster_array_vals(self):
        np.testing.assert_array_almost_equal(
            self.lrr.GetRasterBand(1).ReadAsArray(),
            self.test_relief_raster.GetRasterBand(1).ReadAsArray(),
        )

    def tearDown(self):
        os.remove("tmp.tiff")
