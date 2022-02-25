
from pyproj import Geod
from openquake.globalm.grid.refgrid.triangular_grid import (
    RegularTriangularGrid, dst)
ELLIPSOID = "sphere"


def compute_distances(tgrd, level):
    """
    This computes the distances between the nodes of the facets

    :param tgrd:
        An instance of :class:`RegularTriangularGrid`
    :param int level:
        The splitting level
    """

    geod = Geod(ellps=ELLIPSOID)
    for pts in tgrd.iter_facets(level):
        #_, _, d1 = geod.inv(pts[0][0], pts[0][1], pts[1][0], pts[1][1])
        #_, _, d2 = geod.inv(pts[0][0], pts[0][1], pts[2][0], pts[2][1])
        #_, _, d3 = geod.inv(pts[1][0], pts[1][1], pts[2][0], pts[2][1])
        #for i in range(3):
        #    print(pts[i][0], pts[i][1])
        d1 = dst(pts[0][0], pts[0][1], pts[1][0], pts[1][1])
        d2 = dst(pts[0][0], pts[0][1], pts[2][0], pts[2][1])
        d3 = dst(pts[1][0], pts[1][1], pts[2][0], pts[2][1])
        print(d1, d2, d3)
        assert abs(d1-d2) < d1*0.01
        assert abs(d1-d3) < d1*0.01
        assert abs(d2-d3) < d1*0.01


def main():
    tri_grd = RegularTriangularGrid(1.0)
    tri_grd.create_initial_icosahedron()
    compute_distances(tri_grd, 0)

    tri_grd.split_facets()
    print('-------- split 1')
    compute_distances(tri_grd, 1)

    tri_grd.split_facets()
    print('-------- split 2')
    compute_distances(tri_grd, 2)



if __name__ == '__main__':
    main()

