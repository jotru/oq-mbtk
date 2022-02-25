import sys
import pickle

from osgeo import ogr
from osgeo import osr
from rtree import index


class OutGrd(object):
    """
    A container for the result of a smoothing or box counting procedure

    :parameter grd:
        Each element contains longitude and latitude, number of occurrences,
        bGR and sigma bGR, aGR, sigma aGR, area
    :type grd: list of list
        TODO check the significance of aGR as obtained by the weichert
        procedure prepared by Graeme and Giuseppe
    :parameter grd_spacing:
        The distance between two nearby nodes. Can be either in decimal degrees
        as well as km.
    :type grd_spacing: float
    """

    def __init__(self, grd_spacing):
        self.grd = []
        self.grd_spacing = grd_spacing
        self.index = None

    def create_spatial_index(self):
        """
        Create a rtree spatial index for the nodes of the grid
        """
        self.index = index.Index()
        for i, node in enumerate(self.grd):
            self.index.insert(i, (node[0], node[1], node[0], node[1]))

    def serialize_to_cpickle_file(self, cpickle_file_path, cpickle_file_name):
        """ Serialize an output grid to a cPickle file"""
        # Remove spatial index since it causes ploblems
        self.index = []
        # Serialize object
        output = open(cpickle_file_path+cpickle_file_name, 'wb')
        # Serializing object
        pickle.dump(self, output)
        # Info
        output.close()

    def deserialize_cpickle_file(self, cpickle_file_path, cpickle_file_name):
        """ Serialize an output grid to a cPickle file"""
        input_pickle = open(cpickle_file_path+cpickle_file_name, 'rb')
        aa = pickle.load(input_pickle)
        self.grd = aa.grd
        self.grd_spacing = aa.grd_spacing
        input_pickle.close()

    def create_nrml_as_area_sources(self,
                                    nrml_file_path,
                                    nrml_file_name,
                                    tect_reg):
        """
        This creates a nrml file

        :parameter nrml_file_path:
            Path to the folder where to create the nrml file
        :parameter nrml_file_name:
            Name of the nrml file name
        :parameter tect_reg:
            A tectonic regionalisation object
        :parameter tect_reg_mapping:
            This is simply a hash map where keys are the tectonic regions
            supported by tect_reg and the value is the corresponding tectonic
            region in the classical OQ classification
        """
        strike = 0.0
        rake = 0.0
        dip = 90.0

        hypocentral_depth = 10.0
        magnitude = [9.0]
        depth = [0.0]
        tr_list = set()

        tect_reg_mapping = tect_reg.get_tectonic_region_mapping()

        try:
            nrml = open(nrml_file_path + nrml_file_name, 'w')
            header = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            header += "<nrml xmlns=\"http://openquake.org/xmlns/nrml/0.3\" \n"
            header += "      xmlns:gml=\"http://www.opengis.net/gml\" \n"
            header += "      gml:id=\"n1\">\n"
            header += "<sourceModel gml:id=\"sm01\">\n"
            nrml.write(header)
            for i, pnt in enumerate(self.grd):
                all_tr_list = tect_reg.get_tectonic_region([pnt[0],pnt[1]])
                for tr in all_tr_list:
                    if (len(all_tr_list)) > 1:
                        print(tr, i)
                    tmp_str = ""
                    if (len(all_tr_list)) > 1:
                        tmp_str = tr
                    tstr = "    <areaSource gml:id=\"ID%05d%s\">\n" % (i,tmp_str)
                    tstr += "        <gml:name>"+'{0:d}'.format(i).zfill(5)
                    tstr += "</gml:name>\n"
                    tstr += "        <tectonicRegion>"
                    tstr += tect_reg_mapping[tr]
                    tr_list.add(tr)
                    tstr += "</tectonicRegion>\n"
                    tstr += "        <areaBoundary>\n"
                    tstr += "             <gml:Polygon>\n"
                    tstr += "                 <gml:exterior>\n"
                    tstr += "                     <gml:LinearRing>\n"
                    tstr += "                         <gml:posList>\n"
                    tstr += "%.2f %.2f \n" % (pnt[0]-self.grd_spacing/2,
                                              pnt[1]-self.grd_spacing/2)
                    tstr += "%.2f %.2f \n" % (pnt[0]+self.grd_spacing/2,
                                              pnt[1]-self.grd_spacing/2)
                    tstr += "%.2f %.2f \n" % (pnt[0]+self.grd_spacing/2,
                                              pnt[1]+self.grd_spacing/2)
                    tstr += "%.2f %.2f \n" % (pnt[0]-self.grd_spacing/2,
                                              pnt[1]+self.grd_spacing/2)
                    tstr += "                         </gml:posList>\n"
                    tstr += "                     </gml:LinearRing>\n"
                    tstr += "                 </gml:exterior>\n"
                    tstr += "             </gml:Polygon>\n"
                    tstr += "         </areaBoundary>\n"
                    tstr += "        <ruptureRateModel>\n"
                    tstr += "            <truncatedGutenbergRichter>\n"

                    aGR = pnt[5]
                    bGR = pnt[3]
                    aGR_C = -5.30
                    bGR_C = 1.0
                    aGR_I = -7.30
                    bGR_I = 1.0

                    if (tr == 'C' or tr == 'O' or tr == 'R' or tr == 'S') and \
                       (aGR < aGR_C):
                        aGR = aGR_C
                        bGR = bGR_C
                    elif tr == 'I' and aGR < aGR_I:
                        aGR = aGR_I
                        bGR = bGR_I

                    tstr += "                <aValueCumulative>"
                    tstr += "{0:.2f}".format(aGR)
                    tstr += "</aValueCumulative>\n"
                    tstr += "                <bValue>"
                    tstr += "{0:.2f}".format(bGR)
                    tstr += "</bValue>\n"
                    tstr += "                <minMagnitude>5.0</minMagnitude>\n"
                    tstr += "                <maxMagnitude>8.5</maxMagnitude>\n"
                    tstr += "            </truncatedGutenbergRichter>\n"
                    tstr += "            <strike>%.2f</strike>\n" % (strike)
                    tstr += "            <dip>%.2f</dip>\n" % (dip)
                    tstr += "            <rake>%.2f</rake>\n" % (rake)
                    tstr += "        </ruptureRateModel>\n"
                    tstr += "        <ruptureDepthDistribution>\n"
                    tstr += "            <magnitude>"
                    for mag in magnitude:
                        tstr += "%.2f " % (mag)
                    tstr += "</magnitude>\n"
                    tstr += "            <depth>"
                    for dep in depth:
                        tstr += "%.2f " % (dep)
                    tstr += "</depth>\n"
                    tstr += "        </ruptureDepthDistribution>\n"
                    tstr += "        <hypocentralDepth>"
                    tstr += "%.2f" % (hypocentral_depth)
                    tstr += "</hypocentralDepth>\n"
                    tstr += "    </areaSource>\n"
                    nrml.write(tstr)
            tstr = "</sourceModel>"
            tstr += "</nrml>"
            nrml.write(tstr)
            nrml.close()

        except IOError as err:
            print("I/O error({0}): {1}".format(err.errno, err.strerror))
        except ValueError:
            print("Could not convert data to an integer.")
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

        tect_reg_list = []
        if 'I' in tr_list or 'O' in tr_list or 'R' in tr_list:
            tect_reg_list.append(tect_reg_mapping['I'])
        if 'S' in tr_list:
            tect_reg_list.append(tect_reg_mapping['S'])
        if 'C' in tr_list:
            tect_reg_list.append(tect_reg_mapping['C'])

        return tect_reg_list

    def create_shapefile(self, shp_file_path, shp_file_name):
        """
        This create a shapefile of points. Each point has the following
        attributes: NumOcc, number of occurrences, aGR, Gutenberg-Richter
        """
        spatialReference = osr.SpatialReference()
        spatialReference.SetWellKnownGeogCS('WGS84')
        driverName = "ESRI Shapefile"
        drv = ogr.GetDriverByName(driverName)
        if drv is None:
            print("%s driver not available.\n" % driverName)
            sys.exit(1)
        ds = drv.CreateDataSource(shp_file_path)
        if ds is None:
            print("Creation of output file failed.\n")
            sys.exit(1)
        lyr = ds.CreateLayer(
                shp_file_name,
                spatialReference,
                ogr.wkbPoint)
        if lyr is None:
            print("Layer creation failed.\n")
            sys.exit(1)
        field_defn = ogr.FieldDefn("NumOcc", ogr.OFTReal)
        field_defn.SetWidth(32)
        if lyr.CreateField(field_defn) != 0:
            print("Creating NumOcc field failed.\n")
            sys.exit(1)
        field_defn = ogr.FieldDefn("aGR", ogr.OFTReal)
        field_defn.SetWidth(32)
        if lyr.CreateField(field_defn) != 0:
            print("Creating aGR field failed.\n")
            sys.exit(1)
        field_defn = ogr.FieldDefn("bGR", ogr.OFTReal)
        field_defn.SetWidth(32)
        if lyr.CreateField(field_defn) != 0:
            print("Creating bGR field failed.\n")
            sys.exit(1)
        field_defn = ogr.FieldDefn("aGRSigma", ogr.OFTReal)
        field_defn.SetWidth(32)
        if lyr.CreateField(field_defn) != 0:
            print("Creating aGRSigma field failed.\n")
            sys.exit(1)
        field_defn = ogr.FieldDefn("bGRSigma", ogr.OFTReal)
        field_defn.SetWidth(32)
        if lyr.CreateField(field_defn) != 0:
            print("Creating bGRSigma field failed.\n")
            sys.exit(1)
        # Adding points
        for i, pnt in enumerate(self.grd):
            feat = ogr.Feature(lyr.GetLayerDefn())
            feat.SetField("NumOcc", pnt[2])
            feat.SetField("bGR", pnt[3])
            feat.SetField("bGRSigma", pnt[4])
            feat.SetField("aGR", pnt[5])
            feat.SetField("aGRSigma", pnt[6])

            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(pnt[0], pnt[1])
            feat.SetFID(i)
            feat.SetGeometry(point)
            if lyr.CreateFeature(feat) != 0:
                print("Failed to create feature in shapefile.\n")
                sys.exit(1)
            point.Destroy()
            feat.Destroy()

        ds = None
