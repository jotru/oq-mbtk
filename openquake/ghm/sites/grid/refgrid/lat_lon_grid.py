"""
Module
"""

from openquake.globalm.grid.refgrid.base import ReferenceGrid
from openquake.globalm.grid.refgrid.landmasses_poly import LandmassesPolygons
from shapely.prepared import prep
from rtree import index
from shapely.geometry import Polygon, Point


class RegularLatLonGrid(ReferenceGrid):
    """
    Creates a grid of points equally spaced in units of longitude, latitude.
    This is clearly an irregularly spaced grid in terms of distances.

    :parameter grid_spacing:
        As per abstract class definition
    :parameter bounds:
        List [lon_min, lat_min, lon_max, lat_max]
    """

    def __init__(self, grid_spacing, bounds=[-180.0, -90.0, 180.0, 90.0]):
        """
        This object is described by
        """
        super(RegularLatLonGrid, self).__init__(grid_spacing)
        self.coords = []
        self.index = index.Index()
        """
        self.bounds = bounds
        # Rounding bounds
        self._round_bounds()
        """

    def iter_nodes(self, within_landmasses_flag):
        """
        """
        self.create_spatial_index()

        if within_landmasses_flag:
            name = "/Users/marcop/Documents/work/201201_gar_smoothed_psha/" + \
                "dat/gshhs/GSHHS_shp/h/GSHHS_h_L1.shp"
            lml = LandmassesPolygons(name)

            for i, poly in enumerate(lml.iter_polygons(True)):

                # Create the spatial index
                vtxs = []
                # Loop over the polygons in the shapefile
                for p in xrange(poly.GetGeometryRef(0).GetPointCount()):
                    lon, lat, _ = \
                            poly.GetGeometryRef(0).GetPoint(p)
                    vtxs.append([lon, lat])
                poly = Polygon(vtxs)
                poly_pre = prep(poly)
                bounds = poly.bounds

                # Find index of points in the bounding box
                hits = list(self.index.intersection(bounds, objects=False))
                idxs_pnt_in_bbox = [self.points[i] for i in hits]

                if len(idxs_pnt_in_bbox):
                    selected_points = filter(poly_pre.contains,
                                             idxs_pnt_in_bbox)
                    for p in selected_points:
                        yield p.x, p.y
        else:
            for idx in list(self.index.intersection(self.index.bounds)):
                yield self.coords[idx][0], self.coords[idx][1]

    def create_spatial_index(self):
        """
        Create a rtree spatial index for the nodes of the grid
        """
        # Create the initial grid
        list_of_nodes = []
        list_of_points = []
        cnt = 0
        lon = -180.0+self.grid_spacing/2
        while lon <= 180:
            lat = -90.0+self.grid_spacing/2
            while lat <= 90.0:
                self.index.insert(cnt, (lon, lat, lon, lat))
                cnt += 1
                list_of_nodes.append([lon, lat])
                list_of_points.append(Point([lon, lat]))
                lat += self.grid_spacing
            lon += self.grid_spacing
        self.coords = list_of_nodes
        self.points = list_of_points
