import pyproj


class FacetList(object):
    """
    """
    def __init__(self):
        self.node_indexes = []


class VertexList(object):
    """
    Stores information about the vertexes used to create the triangles at
    various scales. The information collected about vertexes includes their
    coordinates and a dictionary that given a label representing the
    coordinates of a point provides the corresponding index.

    :param list vtx:
    :param dict idx:
    """
    def __init__(self):
        self.vtx = []
        self.idx = {}

    def insert_vtx(self, vtx):
        """
        Adding a new vertex to the list making sure that vertexes are unique
        :param vtx:
        """
        #
        # projection
        geod = pyproj.Geod(ellps='sphere')
        #
        # search for a closeby vertex. The threshold distance is in meters.
        # TODO consider reducing the threshold distance.
        idx = -1
        for i, tvtx in enumerate(self.vtx):
            _, _, dist = geod.inv(tvtx[0], tvtx[1], vtx[0], vtx[1])
            if dist < 10000:
                idx = i
                break
        #
        # if  this condition is satisfied it means that the vertex is a new
        # one and therefore it is appended to the `vtx` list
        if idx < 0:
            self.vtx.append([vtx[0], vtx[1]])
            idx = len(self.vtx)-1
        return idx

    def insvtx(self, vtx, idx=None):
        """
        Adding a new vertex to the list making sure that vertexes are unique

        :param list vtx:
        :param int idx:
        """
        lat = round(vtx[1], 2)
        lon = round(vtx[0], 2)
        if abs(lat) < 0.01:
            lat = +0.0
        elif abs(abs(lat) - 180.0) < 0.01:
            lat = 180.0
        if abs(lon) < 0.01:
            lon = +0.0
        #
        # this is the key used in the `idx` dictionary
        key = "%+.2f_%+.2f" % (lon, lat)
        if idx is None:
            if key in self.idx:
                idx = self.idx[key]
            else:
                self.vtx.append([vtx[0], vtx[1]])
                self.idx[key] = idx = len(self.vtx)-1
        else:
            self.vtx.insert(idx, [vtx[0], vtx[1]])
            self.idx[key] = len(self.vtx)-1
        return idx
