#!/usr/bin/python

"""
"""

import re
import os
import sys
import h5py
import shapely
import numpy as np

from osgeo import ogr
from shapely.wkt import loads
from rtree.index import Index

from openquake.baselib import sap

from pyproj import Proj


def read_polygon_shapefile(fname_shp, field=None, values=None):
    """
    Read the polygons in a shapefile and exports into WKT the ones whose field
    value is within the options provided by the user with the `field` and
    `values` parameters.

    :param str fname_shp:
        Name of the shapefile
    :param str field:
        Field name in the shapefile
    :param list values:
        Values admitted for the in the `field` of the attribute table
    :return:
        A list with polygons on the WKT format
    """
    #
    #
    polys = []
    #
    # Open shapeData
    shapedata = ogr.Open(fname_shp)
    if not shapedata:
        raise Exception("Can't open {:s}".format(fname_shp))
        #
        # Get the first layer
    layer = shapedata.GetLayer()
    #
    # Get layer projection and print info
    _ = layer.GetSpatialRef().ExportToProj4()
    #
    # For each polygon
    for index in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(index)
        geometry = feature.GetGeometryRef()
        #
        # Make sure that it is a polygon
        if geometry.GetGeometryName() == 'POLYGON':
            #
            # TODO add check on selection
            if field is None:
                polys.append(geometry.ExportToWkt())
        elif geometry.GetGeometryName() == 'MULTIPOLYGON':
            #
            # TODO add check on selection
            for geom in geometry:
                if field is None:
                    polys.append(geom.ExportToWkt())
        else:
            print('geometry:', geometry.GetGeometryName())
            raise Exception('This module can only load polygons')
        #
        # Destroying
        feature.Destroy()
    # Cleanup
    shapedata.Destroy()
    return polys


def select_nodes(fname_hdf5, polys):
    """
    Select the points in the mesh included in the polygons of the `polys` list

    :param fname_hdf5:
        The name of the hdf5 containing the grid
    :param list polys:
    """
    #
    # get the spatial index filename
    fname = re.split('\.', os.path.basename(fname_hdf5))[0]
    fpath = os.path.dirname(fname_hdf5)
    #
    # mesh
    f = h5py.File(fname_hdf5, 'r')
    mesh = f['centroids'][:]
    #
    # create the spatial index
    filename_rtree = os.path.join(fpath, fname)
    sidx = Index(filename_rtree)
    #
    # process polygons
    out = []
    for i, poly in enumerate(polys):
        if i > 0:
            exit()
        #
        # find polygon bounding box
        poa = loads(poly)
        lomi, lami, loma, lama = poa.bounds
        #
        # check if the polygon crosses the international date line
        # TODO proper selection
        if (loma - lomi) > 50.:
            idx1 = list(sidx.intersection((loma, lami, 180., lama)))
            idx2 = list(sidx.intersection((-180., lami, lomi, lama)))
            idx = idx1 + idx2
            midlo = 180
        else:
            idx = list(sidx.intersection((lomi, lami, loma, lama)))
            midlo = (lomi + loma) / 2
        #
        # create projection
        p = Proj('+proj=lcc +lon_0={:f}'.format(midlo))
        #
        # updating the list of points
        if len(idx):
            #
            # create a list of pre-selected points
            tmp = []
            for i in idx:
                tmp.append((mesh[i, 0], mesh[i, 1]))
            tmp = np.array(tmp)
            px, py = p(tmp[:, 0], tmp[:, 1])
            #
            # creating the polygon
            lopoly, lapoly = poa.exterior.coords.xy
            xpoly, ypoly = p(lopoly, lapoly)
            #
            # ppoly = shapely.geometry.Polygon([[x, y] for x, y in zip(xpoly,
            #                                                          ypoly)])
            # Create ring
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in zip(xpoly, ypoly):
                ring.AddPoint(x, y)
            # Create polygon
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            #
            #
            if len(xpoly) > 1000:
                fou = open('poly.txt', 'w')
                llo, lla = p(xpoly, ypoly, inverse=True)
                for x, y in zip(llo, lla):
                    fou.write('{:f},{:f}\n'.format(x, y))
                fou.close()
            #
            # selecting the points inside each polygon
            cnt = 0
            for i, (x, y) in enumerate(zip(px, py)):
                # pnt = shapely.geometry.Point(x, y)
                # chk = pnt.within(ppoly)
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(x, y)
                chk = point.Intersects(poly)
                if chk:
                    out.append((tmp[i, 0], tmp[i, 1]))
                    cnt += 1
            #
            # info
            tmps = 'Number of vertexes composing the polygon '
            tmps += ' {:6d} '.format(len(xpoly))
            tmps += '- {:6d} preselected nodes '.format(len(px))
            tmps += '- {:6d} selected nodes '.format(cnt)
            print(tmps)
    #
    # info
    print('Found {:d} points'.format(len(out)))
    #
    # saving csv file
    fou = open('out.csv', 'w')
    for cc in out:
        fou.write('{:.4f},{:.4f}\n'.format(cc[0], cc[1]))
    fou.close()
    return out


def get_pnts(shapefile, hdf5_fname, field_name=None, field_values=None):
    """
    """
    polys = read_polygon_shapefile(shapefile, field_name, field_values)
    select_nodes(hdf5_fname, polys)


def main(argv):
    """
    """
    p = sap.Script(get_pnts)
    p.arg(name='shapefile', help='Shapefile filename')
    p.arg(name='hdf5_fname', help='Name of the .hdf5 file with points')
    p.opt(name='field_name', help='Name of the field in the attribute table')
    p.opt(name='field_values', help='Values admitted for the selection')

    if len(argv) < 1:
        print(p.help())
    else:
        p.callfunc()

if __name__ == '__main__':
    main(sys.argv[1:])
