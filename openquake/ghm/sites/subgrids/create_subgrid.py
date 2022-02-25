
import os
import re
import h5py
import toml
import numpy as np
import geopandas as gpd

from pathlib import Path
from rtree import index
from openquake.baselib import sap
from openquake.ghm.maps.utils import create_query, explode

from shapely.geometry import Polygon, MultiPolygon
from openquake.mbt.tools.geo import get_idx_points_inside_polygon

# limits are defined as: lomin, lomax, lamin, lamax
SUBSETS = {
           'als07': {'US': ['-180 50 -127 50 -127 73 -180 73',
                            '170 47 180 47 180 57 170 57']},
           'cca18': {'CO': ['-83 11 -80 11 -80 14 -83 14']},
           'ids18': {'MY': ['108 0 120 0 120 9 108 9']},
           'esm15': {'RS': ['25 50 25 60 16 60 16 50']},
           'nea18': {'RS': ['75 0 180 0 180 88 75 88',
                            '-180 0 -160 0 -160 88 -180 88']},
           'haw96': {'US': ['-162 18 -153 18 -153 24 -162 24']},
           'nwa18': {'RS': ['75 0 25 0 25 88 75 88']},
           'pac18': {'NZ': ['-180 -20 -160 -20 -160 -10 -180 -10'],
                     'US': ['-180 -20 -160 -20 -160 -10 -180 -10']},
           'naf18': {'AG': ['-21 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20 25 20 30 20 39 20 39 39 -21 39'],
                     'CD': ['-21 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20 25 20 30 20 39 20 39 39 -21 39'],
                     'LY': ['-21 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20 25 20 30 20 39 20 39 39 -21 39'],
                     'ML': ['-21 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20 25 20 30 20 39 20 39 39 -21 39'],
                     'MR': ['-21 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20 25 20 30 20 39 20 39 39 -21 39'],
                     'NG': ['-21 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20 25 20 30 20 39 20 39 39 -21 39'],
                     'SU': ['-21 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20 25 20 30 20 39 20 39 39 -21 39']},
           'sea18': {'MY': ['98 0 106 0 106 8 98 8']},
           'wmy18': {'MY': ['98 0 106 0 106 8 98 8']},
           'ssa16': {'SU': ['25 20 25 5 18 -5.5 20.5 -28.5 50 -28.5 50 20 40 20 35 20 30 20'],
                     'CT': ['25 20 25 5 18 -5.5 20.5 -28.5 50 -28.5 50 20 40 20 35 20 30 20'],
                     'AO': ['25 20 25 5 18 -5.5 20.5 -28.5 50 -28.5 50 20 40 20 35 20 30 20'],
                     'CF': ['25 20 25 5 18 -5.5 20.5 -28.5 50 -28.5 50 20 40 20 35 20 30 20'],
                     'CG': ['25 20 25 5 18 -5.5 20.5 -28.5 50 -28.5 50 20 40 20 35 20 30 20'],
                     'WA': ['25 20 25 5 18 -5.5 20.5 -28.5 50 -28.5 50 20 40 20 35 20 30 20']},
           'usa14': {'US': ['-128 22 -90 22 -75 22 -60 22 -60 53 -128 53']},
           'waf18': {'AG': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'AO': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'BC': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'CD': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'CG': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'CT': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'LY': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'ML': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'MR': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'NG': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'SF': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'SU': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20'],
                     'WA': ['25 20 25 5 18 -5.5 25.0 -31.0 9 -31 -29 20 -20 20 -15 20 -10 20 -5 20 0 20 5 20 10 20 15 20 20 20']}
          }

DATA = {'als07': ['US'],
        'arb18': ['SA', 'QA', 'YM', 'MU', 'AE', 'BA'],
        'aus18': ['AS'],
        'bru18': ['BX'],
        'cca18': ['GT', 'BH', 'ES', 'NU', 'CS', 'PM', 'CU', 'JM', 'HA', 'DR',
                  'HO', 'RQ', 'GP', 'DO', 'MB', 'ST', 'VC', 'TD', 'GJ', 'BF',
                  'TK', 'SC', 'AV', 'AC', 'CJ', 'VQ', 'VI', 'CO', 'BB'],
        'cea15': ['KZ', 'TX', 'UZ', 'KG', 'TI'],
        'can15': ['CA'],
        'chn15': ['CH'],
        'emm16': ['TU', 'IZ', 'SY', 'GG', 'AM', 'AJ', 'IR', 'AF', 'PK', 'CY',
                  'LE', 'JO', 'IS', 'KU', 'GZ', 'WE'],
        'esm15': ['IC', 'EI', 'UK', 'PO', 'SP', 'IT', 'FR', 'GM', 'BE', 'NL',
                  'SZ', 'AU', 'NO', 'SW', 'FI', 'DA', 'PL', 'LH', 'LG', 'EN',
                  'BK', 'HR', 'SR', 'AL', 'MK', 'GR', 'BU', 'RO', 'MD', 'UP',
                  'BO', 'EZ', 'LO', 'HU', 'SI', 'LU', 'AL', 'JE', 'MW', 'MT',
                  'RS'],
        'gre19': ['GL'],
        'haw96': ['US'],
        'ind12': ['IN', 'BT', 'BG', 'CE', 'NP', 'MV'],
        'ids18': ['ID', 'MY', 'TM', 'BX'],
        'jpn14': ['JA'],
        'kor18': ['KN', 'KS'],
        'mex18': ['MX'],
        'nea18': ['RS', 'MG'],
        'nwa18': ['RS'],
        'naf18': ['MO', 'CD', 'TS', 'MR', 'AG', 'GI', 'NG', 'GZ',
                  'SU', 'WI', 'EG', 'LY', 'ML'],
        'nzl10': ['NZ'],
        'pac18': ['AQ', 'CW', 'BP', 'NC', 'NH', 'FJ', 'NF', 'WS', 'NZ', 'TN',
                  'NE', 'KR', 'WF', 'US', 'NE', 'TV', 'TL', 'NR'],
        'png16': ['PP'],
        'phi17': ['RP'],
        'saf18': ['SF', 'LT', 'WZ'],
        'sar18': ['VE', 'CO', 'GY', 'NS', 'FG', 'BR', 'EC', 'PE', 'CI', 'AR',
                  'UY', 'PA', 'BL', 'FK', 'NT'],
        'sea18': ['BM', 'TH', 'LA', 'CB', 'VM', 'MY'],
        'wmy18': ['MY'],
        'ssa16': ['ET', 'SO', 'SU', 'KE', 'TZ', 'UG', 'MI', 'MZ', 'ZI', 'ZA',
                  'ER', 'DJ', 'CG', 'RW', 'BY', 'SU', 'CT', 'AO', 'WA', 'BC',
                  'MA', 'CN', 'RE', 'MP'],
        'tem15': ['TW'],
        'usa14': ['US'],
        'ucerf': ['CA'],
        'waf18': ['TP', 'IV', 'CD', 'MR', 'CF', 'CG', 'AG', 'GA', 'GB', 'GH',
                  'CM', 'AO', 'UV', 'SF', 'SG', 'CT', 'NG', 'SH', 'NI',
                  'SL', 'BC', 'WA', 'GV', 'LI', 'SU', 'LY', 'BN', 'EK',
                  'ML', 'TO', 'PU']}


def get_poly_from_str(tstr):
    """
    :param str tstr:
        A string with a sequence of lon, lat tuples
    """
    li = re.split('\\s+', tstr)
    coo = []
    for i in range(0, len(li), 2):
        coo.append([float(li[i]), float(li[i+1])])
    coo = np.array(coo)
    return coo


"""
def create_query(inpt, field, labels):
    sel = None
    for lab in labels:
        if sel is None:
            sel = inpt[field] == lab
        else:
            sel = sel | (inpt[field] == lab)
    return inpt.loc[sel]
"""


def create(model: str, path: str, path_shp: str, splitlevel: int):
    """
    Create the grids

    :param model:
        The model key
    :param path:
        The path to the folder with the spatial indexes of the global grid
    :param path_shp:
        The path to the folder with the reference shapefile
    """

    # Output folder
    outfolder = './tmp'

    # Path to the folder with the grid spatial indexes
    # path = '/Users/mpagani/Documents/Projects/2017_global_model/py/openquake/globalm/grid/out/'

    # Model key
    # model = m

    # Shapefile
    # path_shp = '/Users/mpagani/Repos/venv/src/oq-mbtk/openquake/ghm/data/gis/'
    #if model == 'ucerf':
    #    fname = 'cb_2017_us_state_500k.shp'
    #else:
    #    fname = 'world_country_admin_boundary_with_fips_codes_mosaic_eu_russia.shp'
    #in_file = os.path.join(path_shp, fname)

    in_file = path_shp

    # Set splitting level and create the names of input files
    splitlevel = int(splitlevel)
    if splitlevel == 9:
        fname_hdf5 = 'trigrd_split_9_spacing_13.hdf5'
        fname_sidx = 'trigrd_split_9_spacing_13'
    elif splitlevel == 6:
        fname_hdf5 = 'trigrd_split_6_spacing_110.hdf5'
        fname_sidx = 'trigrd_split_6_spacing_110'
    elif splitlevel == 5:
        fname_hdf5 = 'trigrd_split_5_spacing_220.hdf5'
        fname_sidx = 'trigrd_split_5_spacing_220'
    else:
        raise ValueError('Split level not recognized')

    sidx_file = os.path.join(path, fname_sidx)
    sidx = index.Rtree(sidx_file)

    # Grid points
    hdf5_file = os.path.join(path, fname_hdf5)
    grdf = h5py.File(hdf5_file, 'r')
    grdp = grdf['centroids'][:]
    grdf.close()

    # Read polygon file and set MODEL attribute
    tmpdf = gpd.read_file(in_file)
    inpt = explode(tmpdf)
    inpt['MODEL'] = model

    # Select polygons the countries composing the given model
    selection = create_query(inpt, 'FIPS_CNTRY_x', DATA[model])

    # Merge polygons into a single one
    one_polygon = selection.dissolve(by='MODEL')

    # Bounding box - the delta is used to create the buffer around
    bb = one_polygon.bounds.values[0]
    dlt = 1.2
    pnt_idxs = [n for n in sidx.intersection((bb[0]-dlt, bb[1]-dlt,
                                              bb[2]+dlt, bb[3]+dlt))]
    idx_all_sel = set()
    for iii, row in selection.iterrows():

        key = row['FIPS_CNTRY_x']
        print('Processing: {:s}'.format(key))

        for pol in [row.geometry]:
            pcoo = []
            for pt in list(pol.exterior.coords):
                pcoo.append(pt)
            pcoo = np.array(pcoo)

            # Select the points within one nation - with buffer
            sel_idx = get_idx_points_inside_polygon(grdp[pnt_idxs, 0],
                                                    grdp[pnt_idxs, 1],
                                                    pcoo[:, 0], pcoo[:, 1],
                                                    pnt_idxs,
                                                    buff_distance=0.001)

            # Select the points
            if model in SUBSETS and key in SUBSETS[model]:
                tmp_idx = []
                for k, tstr in enumerate(SUBSETS[model][key]):
                    tcoo = get_poly_from_str(tstr)
                    sss_idx = get_idx_points_inside_polygon(grdp[sel_idx, 0],
                                                            grdp[sel_idx, 1],
                                                            tcoo[:, 0],
                                                            tcoo[:, 1],
                                                            sel_idx,
                                                            buff_distance=0.)
                    tmp_idx += sss_idx
                sel_idx = tmp_idx

            # Update the list of points selected i.e. the union of what already
            # selected and new points
            idx_all_sel = set(sel_idx) | idx_all_sel
            num1 = len(list(set(sel_idx) & idx_all_sel))
            num2 = len(sel_idx)
            #print('Number of points in common: {:d}/{:d}'.format(num1, num2))
            print('Number of selected points: {:d}'.format(len(idx_all_sel))

    # Add a buffer to the polygon
    one_polygon.buffer(0.3)

    # Output folders
    out_shp = os.path.join(outfolder, 'shp')
    Path(out_shp).mkdir(parents=True, exist_ok=True)
    out_grd = os.path.join(outfolder, 'grd')
    Path(out_grd).mkdir(parents=True, exist_ok=True)

    # Output shapefile
    out_file = os.path.join(out_shp, '{:s}.shp'.format(model))
    one_polygon.to_file(out_file)

    # Output file
    grd_file = 'sites_{:s}_s{:d}.csv'.format(model, splitlevel)
    grd_file = os.path.join(out_grd, grd_file)
    fout = open(grd_file, 'w')
    for i in idx_all_sel:
        fout.write('{:.6f},{:.6f}\n'.format(grdp[i, 0], grdp[i, 1]))
    fout.close()


def main(model_key, split_level, config_file):

    conf = toml.load(config_file)
    path_grd = conf['path_grd']
    path_shp = conf['path_shp']
    create(model_key, path_grd, path_shp, split_level)


main.model_key = "The key identifying a model e.g. esm15"
main.split_level = "An integer specifying the resolution of the grid"
main.config_file = "The name of a .toml file with main settings"

if __name__ == "__main__":
    sap.run(main)
