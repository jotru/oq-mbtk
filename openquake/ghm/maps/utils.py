"""
Module :module:`~openquake.ghm.utils`
"""

import re
import geopandas as gpd

from shapely.geometry import Polygon, MultiPolygon


def explode_old(indf):
    """
    Implements what's suggested here: http://goo.gl/nrRpdV

    :param indf:
        A geodataframe instance
    :returns:
        A geodataframe instance
    """
    outdf = gpd.GeoDataFrame(columns=indf.columns)
    for idx, row in indf.iterrows():
        if type(row.geometry) == Polygon:
            outdf = outdf.append(row, ignore_index=True)
        if type(row.geometry) == MultiPolygon:
            multdf = gpd.GeoDataFrame(columns=indf.columns)
            recs = len(row.geometry)
            multdf = multdf.append([row]*recs, ignore_index=True)
            for geom in range(recs):
                multdf.loc[geom, 'geometry'] = row.geometry[geom]
            outdf = outdf.append(multdf, ignore_index=True)
    return outdf


def explode(gdf):
    """
    Explode multi-polygons into polygons. Implements what's suggested
    here: http://goo.gl/nrRpdV

    :param indf:
        A geodataframe instance
    :returns:
        A geodataframe instance
    """
    gs = gdf.explode(index_parts=True)
    gdf2 = gs.reset_index().rename(columns={0: 'geometry'})
    gdf_out = gdf2.merge(gdf.drop('geometry', axis=1), left_on='level_0',
                         right_index=True)
    gdf_out = gdf_out.set_index(['level_0', 'level_1']).set_geometry('geometry')
    gdf_out.crs = gdf.crs
    return gdf_out


def read_hazard_map_csv(fname):
    """
    Read the content of a .csv file with the mean hazard map computed for a
    given hazard model.

    :param str fname:
        The name of the file containing the results
    :return:
        A dictionary with key the sting in the header and values the floats
        in the csv file
    """
    data = {}
    for line in open(fname):
        if re.search('^#', line):
            pass
        else:
            if re.search('lon', line):
                labels = re.split('\,', line)
            else:
                aa = re.split('\,', line)
                for l, d in zip(labels, aa):
                    if l in data:
                        data[l].append(float(d))
                    else:
                        data[l] = [float(d)]
    return data


def create_query(inpt, field, labels, nott=False):
    """
    This function returns a dataframe containing the rows with values in the
    `field` column that match (or do not match) the values in labels.

    :param inpt:
        A :class:`pandas.DataFrame` instance
    :param field:
        A string with the name of the column
    :param labels:
        A list with the values to match
    :param nott:
        A boolean indicating the the values to keep are the ones matching the
        values in `labels` [False] or not matching [True]
    :returns:
        A :class:`pandas.DataFrame` instance with a subset of the row in the
        input DataFrame.
    """
    sel = None
    if nott:
        sel = inpt.copy()
        sel['idx'] = range(0, len(sel))
        sel.set_index(['idx'], inplace=True)

    for lab in labels:
        if sel is None:
            sel = inpt[field] == lab
        else:
            if nott:
                sel = sel.drop(sel[sel[field] == lab].index, axis=0)
            else:
                sel = sel | (inpt[field] == lab)
    if nott:
        return sel
    return inpt.loc[sel]
