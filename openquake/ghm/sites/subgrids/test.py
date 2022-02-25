
import os
import h5py
import numpy as np
import geopandas as gpd

from rtree import index
from shapely.geometry import Point
from oqmbt.tools.geo import get_idx_points_inside_polygon

# limits are defined as: lomin, lomax, lamin, lamax
LIMITS = {'cac18': [{'CO': [-82, -81, 12, 13]}

DATA = {'alk10': ['US'],
        'arb18': ['SA', 'QA', 'YM', 'MU', 'AE', 'BA'],
        'aus18': ['AU'],
        'cac18': ['GT', 'BH', 'ES', 'NU', 'CS', 'PM', 'CU', 'JM', 'HA', 'DR',
                  'HO', 'RQ', 'GP', 'DO', 'MB', 'ST', 'VC', 'TD', 'GJ', 'BF',
                  'TK', 'SC', 'AV', 'AC', 'CJ', 'VQ', 'VI'],
        'can15': ['CA'],
        'chn15': ['CH'],
        'emm16': ['TU', 'IZ', 'SY', 'GG', 'AM', 'AJ', 'IR', 'AF', 'PK', 'CY',
                  'LE', 'JO', 'IS', 'KU', 'GZ', 'WE'],
        'esh15': ['IC', 'EI', 'UK', 'PO', 'SP', 'IT', 'FR', 'GM', 'BE', 'NL',
                  'SZ', 'AU', 'NO', 'SW', 'FI', 'DA', 'PL', 'LH', 'LG', 'EN',
                  'BK', 'HR', 'SR', 'AL', 'MK', 'GR', 'BU', 'RO', 'MD', 'UP',
                  'BO', 'EZ', 'LO', 'HU', 'SI', 'LU', 'AL'],
        'idi12': ['IN', 'BT', 'BG', 'CE', 'NP'],
        'ind17': ['ID'],
        'jap14': ['JA'],
        'mex18': ['MX'],
        'nea17': ['RS', 'MG'],
        'nzl10': ['NZ'],
        'pai18': ['BP', 'NC', 'NH', 'FG'],
        'png16': ['PP'],
        'phi17': ['RP'],
        'saf18': ['SF', 'LT', 'WZ'],
        'sar18': ['VE', 'CO', 'GY', 'NS', 'FG', 'BR', 'EC', 'PE', 'CI', 'AR',
                  'UY', 'PA', 'BL', 'FK'],
        'sea18': ['BM', 'TH', 'LA', 'CB', 'VM', 'MY'],
        'ssa16': ['ET', 'SO', 'SU', 'KE', 'TZ', 'UG', 'MI', 'MZ', 'ZI', 'ZA',
                  'ER', 'DJ', 'CG', 'RW', 'BY'],
        'tem15': ['TW'],
        'usa14': ['US'],
        }

def create_query(inpt, field, labels):
    sel = None
    for lab in labels:
        if sel is None:
            sel = inpt[field] == lab
        else:
            sel = sel | (inpt[field] == lab)
    return inpt.loc[sel]


def main():
    #
    # shapefile
    path = '/Users/mpagani/NC/Hazard_Charles/Data/Administrative/'
    fname = 'world_country_admin_boundary_shapefile_with_fips_codes.shp'
    in_file = os.path.join(path, fname)
    #
    # spatial index
    path = '/Users/mpagani/Documents/Projects/2017_global_model/py/openquake/globalm/grid/out/'
    fname_sidx = 'trigrd_split_9_spacing_13'
    sidx_file = os.path.join(path, fname_sidx)
    sidx = index.Rtree(sidx_file)
    #
    # grid points
    fname_hdf5 = 'trigrd_split_9_spacing_13.hdf5'
    hdf5_file = os.path.join(path, fname_hdf5)
    grdf = h5py.File(hdf5_file,'r')
    grdp = grdf['centroids'][:]
    grdf.close()
    #
    # set model
    #model = 'sar18'
    #model = 'esh15'
    model = 'cac18'
    #model = 'mex18'
    #model = 'usa14'
    #model = 'saf18'
    #model = 'sea18'
    #model = 'png16'
    #
    # read polygon file and set MODEL attribute
    inpt = gpd.read_file(in_file)
    inpt['MODEL'] = model
    #
    # select polygons
    selection = create_query(inpt, 'FIPS_CNTRY', DATA[model])
    #
    # merge polygons into a single one
    one_polygon = selection.dissolve(by='MODEL')
    #
    # bounding box
    bb = one_polygon.bounds.values[0]
    dlt = 0.5
    pnt_idxs = [n for n in sidx.intersection((bb[0]-dlt, bb[1]-dlt,
                                              bb[2]+dlt, bb[3]+dlt))]

    idx_all_sel = set()
    for iii, row in one_polygon.iterrows():
        for pol in row.geometry:
            pcoo = []
            for pt in list(pol.exterior.coords):
                pcoo.append(pt)
            pcoo = np.array(pcoo)
            sel_idx = get_idx_points_inside_polygon(grdp[pnt_idxs, 0],
                                                    grdp[pnt_idxs, 1],
                                                    pcoo[:, 0], pcoo[:, 1],
                                                    pnt_idxs,
                                                    buff_distance=100000.)
            idx_all_sel = set(sel_idx) | idx_all_sel
            print('num idxs in common:', len(list(set(sel_idx) & idx_all_sel)))
    #
    # add a buffer
    one_polygon.buffer(0.3)
    #
    # output file
    out_file = './shp/{:s}.shp'.format(model)
    one_polygon.to_file(out_file)
    #
    # output file
    grd_file = './grd/sites_{:s}.csv'.format(model)
    fout = open(grd_file, 'w')
    for i in idx_all_sel:
        fout.write('{:.6f},{:.6f}\n'.format(grdp[i, 0], grdp[i, 1]))
    fout.close()


if __name__== "__main__":
    main()
