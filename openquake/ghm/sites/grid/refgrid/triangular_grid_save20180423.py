"""
Module
"""

import h5py
import pyproj
import pickle
import numpy as np

from os import path
from pyproj import Proj
from rtree import index
from shapely.prepared import prep
from shapely.geometry import Polygon, Point
from openquake.globalm.grid.refgrid.landmasses_poly import LandmassesPolygons
from openquake.globalm.grid.refgrid.triangular_grid_aux import (FacetList,
                                                                VertexList)
from openquake.globalm.grid.refgrid.lat_lon_grid import ReferenceGrid
from openquake.globalm.grid.utils.polygonsphere import CoordsConverter
from openquake.globalm.grid.utils.utilities import (resample_grand_arc,
                                                    normalize_lon_lat,
                                                    need_to_normalize)

INFO = False
OPTINS = 3
ELLIPSOID = "sphere"
LANDMASSES = "/Users/marcop/Documents/work/201201_gar_smoothed_psha/" + \
             "dat/gshhs/GSHHS_shp/h/GSHHS_h_L1.shp"


class RegularTriangularGrid(ReferenceGrid):
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
        sons
            It's a hash map that given a facet index contains a list of indexes
            specifying the facets included the next level
        level
            List
        index
            Hash map with a list of the spatial indexes available. The hash key
            specifies the type of index (e.g. 'xy' indicates that the index is
            on the xy plane)
        """
        super(RegularTriangularGrid, self).__init__(grid_spacing)
        self.coords = VertexList()
        self.facets = FacetList()
        self.sons = {}
        self.level = []
        self.index = {}
        self.nodes_interdistance = 3526.82

    def serialize_to_cpickle_file(self,
                                  cpickle_file_path,
                                  cpickle_file_name,
                                  override=False):
        """Serialize a triangular grid to cPickle file"""
        if ((path.exists(cpickle_file_path) and
             path.isfile(cpickle_file_name)) and not override):
            print("File {0:s} already exists" % (cpickle_file_path +
                                                 cpickle_file_name))
        else:
            output = open(cpickle_file_path+cpickle_file_name, 'wb')
            #
            # Serializing object
            pickle.dump(self, output)
            #
            # Info
            output.close()

    def serialize_to_hdf5(self,
                          hdf5_file_path,
                          hdf5_file_name,
                          override=False):
        if ((path.exists(hdf5_file_path) and
             path.isfile(hdf5_file_name)) and not override):
            print("File {0:s} already exists" % (hdf5_file_path +
                                                 hdf5_file_name))
        else:
            f = h5py.File(hdf5_file_path+hdf5_file_name, 'w')
            #
            # nodes group
            g_nodes = f.create_group("nodes")
            #
            # nodes - coordinates dataset
            vtx = np.empty((len(self.coords.vtx), 2))
            for i, tple in enumerate(self.coords.vtx):
                vtx[i, 0] = float(tple[0])
                vtx[i, 1] = float(tple[1])
            g_nodes.create_dataset("vertexes", data=vtx)
            #
            # nodes - indexes dataset
            idxs = np.empty(len(self.coords.idx), dtype=[('key', 'S10'),
                                                         ('vtx_index', 'i4')])
            for i, key in enumerate(self.coords.idx):
                idxs['key'][i] = key
                idxs['vtx_index'][i] = self.coords.idx[key]
            g_nodes.create_dataset("indexes", data=idxs)
            #
            # facets dataset
            f.create_dataset("facets",
                             data=np.array(self.facets.node_indexes))
            #
            # centroids of each facet
            if 'centroids' in self.__dict__:
                f.create_dataset("centroids",
                                 data=np.array(self.centroids))
            f.close()

    def add_facet_centroids(self, split_level=None):
        """
        This computes for every facet the corresponding centroid
        """
        #
        #
        if split_level is None:
            split_level = self.level[-1]
        #
        # initialise the centroids list
        clist = []
        #
        # compute the centroid
        for fct in self.iter_facets(split_level):
            vtxs = np.array(fct)
            p = Proj('+proj=lcc +lon_0={:f}'.format(vtxs[0, 0]))
            x, y = p(vtxs[:, 0], vtxs[:, 1])
            #
            # mid-point on one side
            mx1 = (x[0]+x[1])/2.
            my1 = (y[0]+y[1])/2.
            #
            # centroid
            cx = mx1 + (x[2] - mx1) * 1./3
            cy = my1 + (y[2] - my1) * 1./3
            mlo, mla = p(cx, cy, inverse=True)
            #
            # old procedure
            #mx = np.mean(x)
            #my = np.mean(y)
            #mlo, mla = p(mx, my, inverse=True)
            clist.append([mlo, mla])
        self.centroids = clist
        self.centroids_split_level = split_level

    def iter_facets(self, split_level):
        """
        Iterator for the facets belonging to a specified splitting level

        :parameter split_level:
            Splitting level

        :rtype:
            list of list
        """
        facet_indexes = self.level[split_level]
        for i in facet_indexes:
            # Find the mid point coordinates - first segment
            idx = self.facets.node_indexes[i]
            pnt0 = self.coords.vtx[idx[0]]
            pnt1 = self.coords.vtx[idx[1]]
            pnt2 = self.coords.vtx[idx[2]]
            yield [pnt0, pnt1, pnt2]

    def iter_nodes(self, within_landmasses_flag=False):
        """
        Nodes iterator
        """
        if within_landmasses_flag:
            idx = self.index["ll"]
            # Name of the shapefile containing landmasses
            lml = LandmassesPolygons(LANDMASSES)
            # Looping over the polygons
            for i, poly in enumerate(lml.iter_polygons(True)):
                # Create the spatial index
                vtxs = []
                #
                # Loop over the vertexes in the shapefile
                for p in range(poly.GetGeometryRef(0).GetPointCount()):
                    lon, lat, _ = poly.GetGeometryRef(0).GetPoint(p)
                    vtxs.append([lon, lat])
                poly = Polygon(vtxs)
                poly_pre = prep(poly)
                bounds = poly.bounds
                #
                # Find index of points in the bounding box
                hits = list(idx.intersection(bounds, objects=False))
                pnt_in_bbox = [Point(self.coords.vtx[i]) for i in hits]

                if len(pnt_in_bbox):
                    selected_points = filter(poly_pre.contains,
                                             pnt_in_bbox)
                    for p in selected_points:
                        yield p.x, p.y
        else:
            for vtx in self.coords.vtx:
                yield vtx[0], vtx[1]

    def ins_vtx(self, point):
        """
        Method for inserting a vertex
        """
        idx = -1
        lon = point[0]
        lat = point[1]
        processed = {}
        if len(self.level) > 0:
            geod = pyproj.Geod(ellps=ELLIPSOID)
            # Computing 6 nodes around
            for i in range(-1, 6):
                dst = self.nodes_interdistance/4*1000
                if i > -1:
                    ang = 30 + 60*i
                    lon, lat, unused = geod.fwd(point[0], point[1], ang, dst)
                _, fct_idx = self.find_nearest_node([lon, lat])
                if fct_idx is not(None) and not(processed. has_key(fct_idx)):
                    processed[fct_idx] = 1
                    idx_vtx_facet = self.facets.node_indexes[fct_idx]
                    for vtx_idx in idx_vtx_facet:
                        vtx = self.coords.vtx[vtx_idx]
                        _, _, dist = geod.inv(vtx[0], vtx[1],
                                              point[0], point[1])
                        if dist < 10000:
                            idx = vtx_idx
                            break
        if idx < 0:
            self.coords.vtx.append([point[0], point[1]])
            idx = len(self.coords.vtx)-1
        return idx

    def create_initial_icosahedron(self):
        """
        Create the initial icosahedron
        """
        # North pole
        if OPTINS == 1 or OPTINS == 2:
            self.coords.vtx.insert(0, [0.0, 90.0])
        else:
            self.coords.insvtx([0.0, 90.0], 0)
        # Create the list of longitudes at lat +26.565 and -26.565
        longitude_array = np.linspace(0, 4, 5)*72.0
        for i, lon in enumerate(longitude_array):
            # Upper hemisphere
            lat = 26.56505
            if lon > 180:
                lon = -360 + lon

            if OPTINS == 1 or OPTINS == 2:
                self.coords.vtx.insert(i+1, [lon, lat])
            else:
                self.coords.insvtx([lon, lat], i+1)
            # Lower hemisphere
            lon_lower = lon + 72.0/2.0
            lat = -26.56505
            #
            if OPTINS == 1 or OPTINS == 2:
                self.coords.vtx.insert(i+6, [lon_lower, lat])
            else:
                self.coords.insvtx([lon_lower, lat], i+6)
        # Adding vertexes
        if OPTINS == 1 or OPTINS == 2:
            self.coords.vtx.insert(11, [0.0, -90.0])
        else:
            self.coords.insvtx([0.0, -90.0], 11)
        # Creating facets
        self.facets.node_indexes.insert(0, [0, 1, 2])
        self.facets.node_indexes.insert(1, [0, 2, 3])
        self.facets.node_indexes.insert(2, [0, 3, 4])
        self.facets.node_indexes.insert(3, [0, 4, 5])
        self.facets.node_indexes.insert(4, [0, 5, 1])
        self.facets.node_indexes.insert(5, [1, 6, 2])
        self.facets.node_indexes.insert(6, [2, 6, 7])
        self.facets.node_indexes.insert(7, [2, 7, 3])
        self.facets.node_indexes.insert(8, [3, 7, 8])
        self.facets.node_indexes.insert(9, [3, 8, 4])
        self.facets.node_indexes.insert(10, [4, 8, 9])
        self.facets.node_indexes.insert(11, [4, 9, 5])
        self.facets.node_indexes.insert(12, [5, 9, 10])
        self.facets.node_indexes.insert(13, [5, 10, 1])
        self.facets.node_indexes.insert(14, [1, 10, 6])
        self.facets.node_indexes.insert(15, [6, 11, 7])
        self.facets.node_indexes.insert(16, [7, 11, 8])
        self.facets.node_indexes.insert(17, [8, 11, 9])
        self.facets.node_indexes.insert(18, [9, 11, 10])
        self.facets.node_indexes.insert(19, [10, 11, 6])
        #
        # Creating splitting level
        self.level.insert(0, [i for i in range(20)])

    def split_facets(self):
        """
        Split each facet into 4 smaller equilateral triangular facets
        """
        # Assuming a spheric geoid with radius = 6370997m
        geod = pyproj.Geod(ellps='sphere')
        if INFO:
            print("Splitting level ", len(self.level)-1)
        #
        # angular distance
        angd = 26.56505 / 2**(len(self.level)-1)
        print(angd)
        #
        # Finding indexes of the facets to split
        facet_indexes = self.level[len(self.level)-1]
        facets_in_this_level = []
        # Now starting to split the facets in the uppermost splitting level
        first = True
        for i in facet_indexes:
            #
            # Find the mid point coordinates - first segment
            idx = self.facets.node_indexes[i]
            pnt0 = self.coords.vtx[idx[0]]
            pnt1 = self.coords.vtx[idx[1]]
            pnt2 = self.coords.vtx[idx[2]]
            #
            # Find the mid point of each side
            mid_pnt_0a = geod.npts(pnt0[0], pnt0[1], pnt1[0], pnt1[1], 1)
            mid_pnt_1a = geod.npts(pnt1[0], pnt1[1], pnt2[0], pnt2[1], 1)
            mid_pnt_2a = geod.npts(pnt2[0], pnt2[1], pnt0[0], pnt0[1], 1)

            mid_pnt_0 = [((pnt0[0]+pnt1[0])/2, (pnt0[1]+pnt1[1])/2)]
            mid_pnt_1 = [((pnt1[0]+pnt2[0])/2, (pnt1[1]+pnt2[1])/2)]
            mid_pnt_2 = [((pnt0[0]+pnt2[0])/2, (pnt0[1]+pnt2[1])/2)]

            print(mid_pnt_0, mid_pnt_0a)
            print(mid_pnt_1, mid_pnt_1a)
            print(mid_pnt_2, mid_pnt_2a)

            if 1:
                #print('>>>>', pnt0, pnt1, pnt2)
                #
                # check mid-point calculation
                _, _, d0 = geod.inv(mid_pnt_0[0][0], mid_pnt_0[0][1],
                                    pnt0[0], pnt0[1])
                _, _, d1 = geod.inv(mid_pnt_0[0][0], mid_pnt_0[0][1],
                                    pnt1[0], pnt1[1])
                assert abs(d0-d1) < d0*0.01
                _, _, d0 = geod.inv(mid_pnt_1[0][0], mid_pnt_1[0][1],
                                    pnt1[0], pnt1[1])
                _, _, d1 = geod.inv(mid_pnt_1[0][0], mid_pnt_1[0][1],
                                    pnt2[0], pnt2[1])
                assert abs(d0-d1) < d0*0.01
                _, _, d0 = geod.inv(mid_pnt_2[0][0], mid_pnt_2[0][1],
                                    pnt0[0], pnt0[1])
                _, _, d1 = geod.inv(mid_pnt_2[0][0], mid_pnt_2[0][1],
                                    pnt2[0], pnt2[1])
                assert abs(d0-d1) < d0*0.01
                #
                # check distance between vertexes
                _, _, d0 = geod.inv(pnt0[0], pnt0[1], pnt1[0], pnt1[1])
                _, _, d1 = geod.inv(pnt0[0], pnt0[1], pnt2[0], pnt2[1])
                _, _, d2 = geod.inv(pnt1[0], pnt1[1], pnt2[0], pnt2[1])
                
                try:
                    assert abs(d0-d1) < d0*0.01
                    assert abs(d0-d2) < d0*0.01
                    assert abs(d1-d2) < d1*0.01
                except:
                    print('facet index:', i)
                    print('facet node indexes:', idx)
                    print(pnt0, pnt1, pnt2)
                    print('mid point 0:', mid_pnt_0)
                    print('mid point 1:', mid_pnt_1)
                    print('mid point 2:', mid_pnt_2)
                    print('d0-d1', abs(d0-d1))
                    print('d0-d2', abs(d0-d2))
                    print('d1-d2', abs(d1-d2))
                    exit(0)

            if OPTINS == 1:
                # Update the vertex
                idx0 = self.coords.insert_vtx([mid_pnt_0[0][0],
                                               mid_pnt_0[0][1]])
                idx1 = self.coords.insert_vtx([mid_pnt_1[0][0],
                                               mid_pnt_1[0][1]])
                idx2 = self.coords.insert_vtx([mid_pnt_2[0][0],
                                               mid_pnt_2[0][1]])
            elif OPTINS == 2:
                idx0 = self.ins_vtx([mid_pnt_0[0][0], mid_pnt_0[0][1]])
                idx1 = self.ins_vtx([mid_pnt_1[0][0], mid_pnt_1[0][1]])
                idx2 = self.ins_vtx([mid_pnt_2[0][0], mid_pnt_2[0][1]])
            else:
                idx0 = self.coords.insvtx([mid_pnt_0[0][0], mid_pnt_0[0][1]])
                idx1 = self.coords.insvtx([mid_pnt_1[0][0], mid_pnt_1[0][1]])
                idx2 = self.coords.insvtx([mid_pnt_2[0][0], mid_pnt_2[0][1]])
            #
            # Update the facets array
            idx_low = len(self.facets.node_indexes)
            self.facets.node_indexes.append([idx[0], idx0, idx2])
            self.facets.node_indexes.append([idx0, idx[1], idx1])
            self.facets.node_indexes.append([idx2, idx1, idx[2]])
            self.facets.node_indexes.append([idx0, idx1, idx2])
            idx_upp = len(self.facets.node_indexes)
            lst = []
            # Updating the list of facets in this level
            for j in range(idx_low, idx_upp):
                facets_in_this_level.append(j)
                lst.append(j)
            #
            # Updating the hash map
            self.sons["%d" % (i)] = lst
            #
            # Compute approximate distance between two vertexes
            if first:
                unused, unused, distance = geod.inv(
                        pnt0[0], pnt0[1], mid_pnt_0[0][0], mid_pnt_0[0][1])
                first = False
            # self.level[level_idx] = facets_in_this_level
            # print level_idx, len(self.level)
        # Final info
        print("distance between nodes: %.2f [km]" % (distance/1000.0))
        self.level.append(facets_in_this_level)
        # self.nodes_interdistance = distance/1000.0
        return distance/1000.0

    def create_spatial_indexes(self):
        """
        """
        conv = CoordsConverter()
        idx = index.Index()
        idxll = index.Index()
        cnt = 0
        for vtx in self.coords.vtx:
            lo = vtx[0]
            la = vtx[1]
            xyz = conv.get_xyz([[lo, la]])
            # Insert point left, bottom, right, top
            idx.add(cnt, (xyz[0][0], xyz[0][1], xyz[0][0], xyz[0][1]))
            idxll.add(cnt, (lo, la, lo, la))
            cnt += 1
        self.index["xy"] = idx
        self.index["ll"] = idxll
        print("created spatial indexes", self.index.keys())

    def find_nearest_node(self, point):
        """
        :param point:
        """
        tmp_fct_idx = None
        # Map projection string
        frmt_str = "+proj=stere +lat_0=%.1f +lon_0=%.1f"
        # Finding the parent facet
        for fct_idx in self.level[0]:
            ppp = np.array(point)
            node_idx = self.facets.node_indexes[fct_idx]
            pnt = self.coords.vtx[node_idx[0]]
            # Create the map projection
            pstr = frmt_str % (pnt[0], pnt[1])
            proj = pyproj.Proj(pstr)
            # Resample the borders of the triangle so as to avoid point in
            # polygon problems due to map projection distortions
            flag = need_to_normalize(self.coords.vtx[node_idx[0]][0],
                                     self.coords.vtx[node_idx[0]][1],
                                     self.coords.vtx[node_idx[1]][0],
                                     self.coords.vtx[node_idx[1]][1],
                                     self.coords.vtx[node_idx[2]][0],
                                     self.coords.vtx[node_idx[2]][1])
            # Creating the polygon border
            lonlats = resample_grand_arc(self.coords.vtx[node_idx[0]][0],
                                         self.coords.vtx[node_idx[0]][1],
                                         self.coords.vtx[node_idx[1]][0],
                                         self.coords.vtx[node_idx[1]][1], 1e5)
            if flag:
                lonlats = normalize_lon_lat(lonlats)
            coo = lonlats
            lonlats1 = resample_grand_arc(self.coords.vtx[node_idx[1]][0],
                                          self.coords.vtx[node_idx[1]][1],
                                          self.coords.vtx[node_idx[2]][0],
                                          self.coords.vtx[node_idx[2]][1], 1e5)
            if flag:
                lonlats1 = normalize_lon_lat(lonlats1)
            coo = np.concatenate((coo, lonlats1), axis=0)
            lonlats2 = resample_grand_arc(self.coords.vtx[node_idx[2]][0],
                                          self.coords.vtx[node_idx[2]][1],
                                          self.coords.vtx[node_idx[0]][0],
                                          self.coords.vtx[node_idx[0]][1], 1e5)
            if flag:
                lonlats2 = normalize_lon_lat(lonlats2)
            coo = np.concatenate((coo, lonlats2), axis=0)
            # Projecting point and polygon
            ppoint = np.array((proj(ppp[0], ppp[1])))
            pcoo = np.transpose(np.array(proj(coo[:, 0], coo[:, 1])))
            # Correcting the point coordinates
            norm_pnt = ppp
            if flag:
                norm_pnt = normalize_lon_lat(ppp)
            ppoint = norm_pnt
            pcoo = coo

            if Polygon(pcoo).contains(Point(ppoint)):
                if INFO:
                    print("POINT    :", norm_pnt[0], norm_pnt[1])
                    print("FACET IDX:", fct_idx)
                tmp_fct_idx = fct_idx
                break

        # Now drilling
        for unused in range(0, len(self.level)-1):
            if self.sons. has_key("%d" % (fct_idx)):
                if INFO:
                    print("drilling from facet ", fct_idx)
                for tmp_fct_idx in self.sons["%d" % (fct_idx)]:
                    node_idx = self.facets.node_indexes[tmp_fct_idx]
                    pnt = self.coords.vtx[node_idx[0]]
                    pstr = frmt_str % (pnt[0], pnt[1])
                    proj = pyproj.Proj(pstr)
                    pnt0 = np.array((proj(pnt[0], pnt[1])))
                    pnt = self.coords.vtx[node_idx[1]]
                    pnt1 = np.array((proj(pnt[0], pnt[1])))
                    pnt = self.coords.vtx[node_idx[2]]
                    pnt2 = np.array((proj(pnt[0], pnt[1])))
                    ppoint = np.array((proj(norm_pnt[0], norm_pnt[1])))
                    if Polygon([pnt0, pnt1, pnt2]).contains(Point(ppoint)):
                        fct_idx = tmp_fct_idx
                        if INFO:
                            print("FACET:", fct_idx)
                            print("POINT:", norm_pnt[0], norm_pnt[1])
                            print("FACET VTX IDX:",)
                            print(node_idx[0], node_idx[1], node_idx[2])
                            print("  ", self.coords.vtx[node_idx[0]])
                            print("  ", self.coords.vtx[node_idx[1]])
                            print("  ", self.coords.vtx[node_idx[2]])
                        break
        # This is temporary
        vtx_idx = None
        return vtx_idx, tmp_fct_idx
