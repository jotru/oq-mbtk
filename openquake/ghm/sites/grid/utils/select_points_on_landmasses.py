

from rtree import index

from shapely.geometry import Polygon, Point
from shapely.prepared import prep

from openquake.globalm.grid.refgrid.landmasses_poly import LandmassesPolygons


class SelectPointOnLandmasses():
    """
    :param pnt_lst:
    :param shapefile:
    """

    def __init__(self, pnt_lst, shapefile):
        self.lm_poly = LandmassesPolygons(shapefile)
        self.nodes = pnt_lst
        self.idx = index.Index()
        self._create_index()

    def _create_index(self):
        for i, pnt in enumerate(self.nodes):
            self.idx.add(i, (pnt[0], pnt[1], pnt[0], pnt[1]))

    def iter_nodes(self):
        """
        """
        idx = self.idx
        for poly in self.lm_poly.iter_polygons(True):
            #
            # loading the vertexes of the polygon
            vtxs = []
            for p in range(poly.GetGeometryRef(0).GetPointCount()):
                lon, lat, _ = poly.GetGeometryRef(0).GetPoint(p)
                vtxs.append([lon, lat])
            #
            # create Shapely polygon
            poly = Polygon(vtxs)
            poly_pre = prep(poly)
            bounds = poly.bounds

            # Find index of points in the bounding box
            hits = list(idx.intersection(bounds, objects=False))
            pnt_in_bbox = [Point(self.nodes[j]) for j in hits]

            if len(pnt_in_bbox):
                selected_points = filter(poly_pre.contains,
                                         pnt_in_bbox)
                for p in selected_points:
                    yield p.x, p.y
