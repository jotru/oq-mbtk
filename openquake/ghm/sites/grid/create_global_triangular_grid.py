#!/usr/bin/python

"""
"""

import sys
from openquake.globalm.grid.refgrid.triangular_grid import \
    RegularTriangularGrid

DST = {0:0000, 1: 3200, 2: 1600, 3:880, 4: 440, 5: 220, 6: 110, 7: 55,
       8: 27, 9: 13, 10: 6}
REF_GRID_PATH = './out/'


def main(argv):
    """
    """
    num_splits = int(argv[0])

    ref_grid_hdf5_name = 'trigrd_split_{:d}_spacing_{:d}.hdf5'.format(
        num_splits, DST[num_splits])

    tri_grd = RegularTriangularGrid(1.0)
    tri_grd.create_initial_icosahedron()
    for i in range(0, num_splits):
        tri_grd.split_facets()

    tri_grd.add_facet_centroids(num_splits)
    tri_grd.serialize_to_hdf5(REF_GRID_PATH, ref_grid_hdf5_name)
    print("done")


if __name__ == '__main__':
    main(sys.argv[1:])
