# -*- coding: utf-8 -*-

import numpy as np
import pylab

from pyproj import Geod
from shapely.geometry import Polygon, Point

# Earth radius (from proj manual)
RADIUS = 6370997
GEOID = "WGS84"

def get_azimuths(pnt0,pnt1):
    geod = Geod("ellps="+GEOID)
    aziA, aziB, dst = geod.inv(pnt0[0], pnt0[1], pnt1[0], pnt1[1], 
                               radians=False)
    # Same latitude
    if pnt0[1] == pnt1[1]:
        aziA = 90; aziB = -90;
        if pnt0[0] > pnt1[0]: aziA = -90; aziB = 90;
    # Same longitude
    elif pnt0[0] == pnt1[0]:
        aziA = 0; aziB = 180;
        if pnt0[1] > pnt1[1]: aziA = 180; aziB = 0;
    return aziA, aziB, dst            

class CoordsConverter(object):
    
    def __init__(self):
        pass
    
    def get_xyz(self,vtx):
        """
        This converts from lon, lat coordinates to x,y,z coordinates. 
        The origin is at the earth center in [0,0,0]
        
        :parameter vtx:
        :type :
        
        Returns a point in a xyz reference (unit of measure is meters) 
        """
        out = []
        for v in vtx:
            z = np.sin(np.radians(v[1])) * RADIUS 
            t = np.cos(np.radians(v[1])) * RADIUS
            x = np.cos(np.radians(v[0])) * t
            y = np.sin(np.radians(v[0])) * t
            out.append([x, y, z])
        return out
            
    def get_ll(self,xyz):
        """
        This converts a x,y,z triple into the corresponding lon and lat 
        coordinates
        """
        out = []
        for x in xyz:
            lon = np.degrees(np.arctan(x[1]/x[0]))
            dst = np.sqrt(x[0]*x[0]+x[1]*x[1])
            lat = np.degrees(np.arctan(x[2]/dst))
            out.append([lon,lat])
        return out

class RectangleSphere(object):
    """
    This class 
    """
    
    def __init__(self,oppvtx,buffer_dst=0):
        """
        :parameter oppvtx:
            A list containing the ll_lon, ll_lat, ur_lon, ur_lat
        :parameter delta:
            The buffer distance to be used in km
        """
        # Rectangle vertexes in lon, lat
        self.vtx = [[oppvtx[0],oppvtx[1]]]
        self.vtx.append([oppvtx[2],oppvtx[1]])
        self.vtx.append([oppvtx[2],oppvtx[3]])
        self.vtx.append([oppvtx[0],oppvtx[3]])
        # Container for the xyz coordinates
        self.xyz = []
        self.vtxb = []
        self.xyzb = []
        self._convert_rectangle_into_xyz()
        self.has_buffer = False
        
    def _convert_rectangle_into_xyz(self):
        """
        Converts the coordinates of the vertexes from longitude, latitude to 
        x, y, z coordinates
        """
        con = CoordsConverter()
        for v in self.vtx:
            self.xyz.append(con.get_xyz([v])[0])
            
    def add_buffer(self,buffer_dst):
        """
        :parameter buffer_dst:
            Buffer distance [in km]
        """
        geod = Geod("ellps="+GEOID)
        vtxb = []
        # Compute vertexes of the buffer: LL
        azi02, azi20, _ = geod.inv(self.vtx[0][0], self.vtx[0][1], 
                                   self.vtx[2][0], self.vtx[2][1], 
                                   radians=False)
        azi13, azi31, _ = geod.inv(self.vtx[1][0], self.vtx[1][1], 
                                   self.vtx[3][0], self.vtx[3][1], 
                                   radians=False)  
        azi01, azi10, _ = geod.inv(self.vtx[0][0], self.vtx[0][1], 
                                   self.vtx[1][0], self.vtx[1][1], 
                                   radians=False)
        azi12, azi21, _ = get_azimuths([self.vtx[1][0], self.vtx[1][1]], 
                                       [self.vtx[2][0], self.vtx[2][1]]) 
        azi23, azi32, _ = geod.inv(self.vtx[2][0], self.vtx[2][1], 
                                   self.vtx[3][0], self.vtx[3][1], 
                                   radians=False)   
        azi30, azi03, _ = geod.inv(self.vtx[3][0], self.vtx[3][1], 
                                   self.vtx[0][0], self.vtx[0][1], 
                                   radians=False)
        # Compute the vertexes of the buffer: LL Vertex
        lon = [0,0,0]; lat = [0,0,0]
        # Define insert sequence
        seq = [0, 1, 2]
        # Compute coordinates around the LL vertex
        lon[0], lat[0], _ = geod.fwd(self.vtx[0][0], self.vtx[0][1], 
                               azi10, buffer_dst*1e3,
                               radians=False)
        lon[1], lat[1], _ = geod.fwd(self.vtx[0][0], self.vtx[0][1], 
                               azi20, buffer_dst*1e3,
                               radians=False)
        lon[2], lat[2], _ = geod.fwd(self.vtx[0][0], self.vtx[0][1], 
                               azi30, buffer_dst*1e3,
                               radians=False)
        for i in seq: vtxb.append([lon[i],lat[i]])
        # Compute vertexes of the buffer: LR
        lon[0], lat[0], _ = geod.fwd(self.vtx[1][0], self.vtx[1][1], 
                               azi21, buffer_dst*1e3,
                               radians=False)        
        lon[1], lat[1], _ = geod.fwd(self.vtx[1][0], self.vtx[1][1], 
                               azi31, buffer_dst*1e3,
                               radians=False)        
        lon[2], lat[2], _ = geod.fwd(self.vtx[1][0], self.vtx[1][1], 
                               azi01, buffer_dst*1e3,
                               radians=False)
        for i in seq: vtxb.append([lon[i],lat[i]])        
        # Compute vertexes of the buffer: UR
        lon[0], lat[0], _ = geod.fwd(self.vtx[2][0], self.vtx[2][1], 
                               azi32, buffer_dst*1e3,
                               radians=False)           
        lon[1], lat[1], _ = geod.fwd(self.vtx[2][0], self.vtx[2][1], 
                               azi02, buffer_dst*1e3,
                               radians=False)           
        lon[2], lat[2], _ = geod.fwd(self.vtx[2][0], self.vtx[2][1], 
                               azi12, buffer_dst*1e3,
                               radians=False)           
        for i in seq: vtxb.append([lon[i],lat[i]])
        # Compute vertexes of the buffer: UL
        lon[0], lat[0], _ = geod.fwd(self.vtx[3][0], self.vtx[3][1], 
                               azi03, buffer_dst*1e3,
                               radians=False)           
        lon[1], lat[1], _ = geod.fwd(self.vtx[3][0], self.vtx[3][1], 
                               azi13, buffer_dst*1e3,
                               radians=False)           
        lon[2], lat[2], _ = geod.fwd(self.vtx[3][0], self.vtx[3][1], 
                               azi23, buffer_dst*1e3,
                               radians=False)  
        for i in seq: vtxb.append([lon[i],lat[i]])  
        self.vtxb = vtxb
        # Adding intermediate points
        num = 20
        lonlat = geod.npts(self.vtxb[2][0], self.vtxb[2][1], 
                           self.vtxb[3][0], self.vtxb[3][1], num)
        for i,v in enumerate(lonlat):
            self.vtxb.insert(i+3,v)
        lonlat = geod.npts(self.vtxb[5+num][0], self.vtxb[5+num][1], 
                           self.vtxb[6+num][0], self.vtxb[6+num][1], num)
        for i,v in enumerate(lonlat):
            self.vtxb.insert(i+6+num,v)
        lonlat = geod.npts(self.vtxb[8+2*num][0], self.vtxb[8+2*num][1], 
                           self.vtxb[9+2*num][0], self.vtxb[9+2*num][1], num)
        for i,v in enumerate(lonlat):
            self.vtxb.insert(i+9+2*num,v)
        
        lonlat = geod.npts(self.vtxb[len(self.vtxb)-1][0], 
                            self.vtxb[len(self.vtxb)-1][1], 
                           self.vtxb[0][0], self.vtxb[0][1], num)
        for i,v in enumerate(lonlat):
            self.vtxb.append(v)
        # Updating the class parameter
        self.vtxb = vtxb
        self.has_buffer = True
        # Converting coordinates into xyz
        cnv = CoordsConverter()
        self.xyzb = cnv.get_xyz(vtxb)

    def get_bb(self):
        """
        """
        if self.has_buffer: 
            minx = +1e100; maxx = -1e100
            miny = +1e100; maxy = -1e100
            minz = +1e100; maxz = -1e100
            for x in self.xyzb:
                minx = min(minx,x[0])
                maxx = max(maxx,x[0])
                miny = min(miny,x[1])
                maxy = max(maxy,x[1])
                minz = min(minz,x[2])
                maxz = max(maxz,x[2])
        else:
            raise StandardError('rectangle hasn\'t a buffer')
        return [minx,maxx,miny,maxy,minz,maxz]

    def get_points_inside(self,hashkey,idx,xyzp):
        """
        :parameter hash key:
            It's a string specifying the type of spatial index 
                - xy
                - xz
                - yz
                - ll
        :parameter idx: 
            It's an instance of the rtree class
        :parameter xyzp:
            A list of points
        """
        lst = []
        # Finding the points included in the bounding box
        if hashkey is 'xy':
            minx = +1e100; maxx = -1e100
            miny = +1e100; maxy = -1e100
            for x in self.xyzb:
                minx = min(minx,x[0])
                maxx = max(maxx,x[0])
                miny = min(miny,x[1])
                maxy = max(maxy,x[1])
            lst = list(idx.intersection((minx, miny, maxx, maxy)))
            # Now use shapely to find the points inside the polygon
            polygon = Polygon(self.xyzb)
            outlist = []
            for idx in lst:
                #print type(xyzp)
                #print type(xyzp[idx])
                #print xyzp[idx][0]
                #print type(xyzp[idx][0])
                if Point(list(xyzp[idx])).within(polygon):
                    outlist.append(idx)
        elif hashkey is 'xz':
            xzp = []
            for x in xyzp:
                xzp.append([x[0],x[2]])
            xzb = []
            for x in self.xyzb:
                xzb.append([x[0],x[2]]) 
            
            minx = +1e100; maxx = -1e100
            minz = +1e100; maxz = -1e100
            for x in self.xyzb:
                minx = min(minx,x[0])
                maxx = max(maxx,x[0])
                minz = min(minz,x[2])
                maxz = max(maxz,x[2])
            lst = list(idx.intersection((minx, minz, maxx, maxz)))
            # Now use shapely to find the points inside the polygon
            polygon = Polygon(xzb)
            outlist = []
            for idx in lst:
                if Point(xzp[idx]).within(polygon):
                    outlist.append(idx)

        elif hashkey is 'yz':
            yzp = []
            for x in xyzp:
                yzp.append([x[1],x[2]])
            yzb = []
            for x in self.xyzb:
                yzb.append([x[1],x[2]]) 
            
            miny = +1e100; maxy = -1e100
            minz = +1e100; maxz = -1e100
            for x in self.xyzb:
                miny = min(miny,x[1])
                maxy = max(maxy,x[1])
                minz = min(minz,x[2])
                maxz = max(maxz,x[2])
            lst = list(idx.intersection((miny, minz, maxy, maxz)))
            # Now use shapely to find the points inside the polygon
            polygon = Polygon(yzb)
            outlist = []
            for idx in lst:
                if Point(yzp[idx]).within(polygon):
                    outlist.append(idx)  
          
        else:
            raise IOError('hash key not recognized')
        
#        # Plotting
#        x = []; y = []
#        for v in xyzp:
#            x.append(v[0])
#            y.append(v[1])
#        pylab.plot(x,y,'g+')
#        x = []; y = []
#        for v in self.xyzb:
#            x.append(v[0])
#            y.append(v[1])
#        pylab.plot(x,y,'mx--')
#        x = []; y = []
#        for v in self.xyz:
#            x.append(v[0])
#            y.append(v[1])
#        pylab.plot(x,y,'ro')
#        pylab.ylim([-RADIUS,RADIUS])
#        pylab.xlim([-RADIUS,RADIUS])
#        pylab.show()
                
        return outlist, lst
        
