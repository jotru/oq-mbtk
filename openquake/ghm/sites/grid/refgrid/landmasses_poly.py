import osgeo.ogr


class LandmassesPolygons:

    def __init__(self, shapefile_name):
        """
        Loads continents from a shapefile
        """
        self.shapefile_name = shapefile_name

    def iter_polygons(self):
        """
        Polygon iterator
        """
        # Open shapeData
        shapeData = osgeo.ogr.Open(self.shapefile_name)
        if not shapeData:
            raise Exception("Can't open %s" % self.shapefile_name)
        # Get the first layer
        layer = shapeData.GetLayer()
        # Get layer projection and print info
        _ = layer.GetSpatialRef().ExportToProj4()
        # For each polygon
        for index in range(layer.GetFeatureCount()):
            feature = layer.GetFeature(index)
            geometry = feature.GetGeometryRef()
            # Make sure that it is a polygon
            if geometry.GetGeometryType() != osgeo.ogr.wkbPolygon:
                raise Exception('This module can only load polygons')
            yield geometry
            feature.Destroy()
        # Cleanup
        shapeData.Destroy()
