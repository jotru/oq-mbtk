
import os
from osgeo import ogr


def create_filtered_shapefile(value, filter_field, in_shapefile,
                              out_shapefile='tmp.shp', name=None):
    """
    :param value:
        A list values
    :param str filter_field:
        The name of the field to be used for the selection
    :param in_shapefile:
        The name of the input shapefile
    :param out_shapefile:
        The name of the output shapefile
    """
    input_layer = read_shapefile(in_shapefile)
    #
    # Filter by our query
    # mylayer.SetAttributeFilter("C_info = 'g19' or C_info = 'g22' or C_info = 'g23'")
    if len(value) == 1:
        query_str = '"{}" = "{}"'.format(filter_field, value[0])
        input_layer.SetAttributeFilter(query_str)
    else:
        raise ValueError("")
    input_layer.SetAttributeFilter(query_str)
    #
    # Copy Filtered Layer and Output File
    driver = ogr.GetDriverByName('ESRI Shapefile')
    out_ds = driver.CreateDataSource(out_shapefile)
    if name == None:
        name = str(value[0])
    out_layer = out_ds.CopyLayer(input_layer, name)
    #
    # cleanup
    del input_layer, out_layer, out_ds
    #
    # results

def read_shapefile(filename):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(filename, 0) # 0 means read-only. 1 means writeable.
    #
    # Check to see if shapefile is found.
    if dataSource is None:
        print('Could not open %s' % (filename))
    else:
        print('Opened %s' % (filename))
        layer = dataSource.GetLayer()
        featureCount = layer.GetFeatureCount()
        print("Number of features: %d" % (featureCount))
    return layer


def main():

    in_shapefile = "/Users/mpagani/NC/Hazard_Charles/Data/Administrative/world_country_admin_boundary_shapefile_with_fips_codes.shp"
    #layer = read_shapefile(in_shapefile)
    #
    # select polygons
    create_filtered_shapefile(['US'], 'FIPS_CNTRY', in_shapefile,
                              out_shapefile='./us.shp')

    #
    #
    """
    bufferDistance = 500
    poly = pt.Buffer(bufferDistance)
    print "%s buffered by %d is %s" % (pt.ExportToWkt(), bufferDistance, poly.ExportToWkt())
    """



if __name__== "__main__":
  main()


