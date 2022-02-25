import math
import numpy as np

from pyproj import Geod

INFO = False


def normalize_points(pnt):
    """
    Normalize the point coordinates. Normalize in this context means that
    each point coordinates are constrained within the domains [-180:lon:180]
    and [-90:lat:90]

    :parameter pnt: The point coordinates eventually to be normalized.
    :type pnt: list with 2 elements

    :rtype:
        list, list
    """
    normalisation_map = [0, 0, 0, 0]
    delta = 1e-6
    # Normalizing longitude
    if pnt[0] > 180+delta:
        if INFO:
            print("A from", pnt[0], end='')
        pnt[0] = -360+pnt[0]
        normalisation_map[0] = 1
        if INFO:
            print(" to", pnt[0], end='')
    elif pnt[0] < -180-delta:
        if INFO:
            print("B from", pnt[0], end='')
        pnt[0] = 360+pnt[0]
        normalisation_map[1] = 1
        if INFO:
            print(" to", pnt[0], end='')
    # Normalizing latitude
    if pnt[1] > 90.0+delta:
        if INFO:
            print("C from", pnt[1], end='')
        pnt[1] = 180-pnt[1]
        normalisation_map[2] = 1
        if INFO:
            print(" to", pnt[1], end='')
    elif pnt[1] < -90-delta:
        if INFO:
            print("D from", pnt[1], end='')
        pnt[1] = -180-pnt[1]
        normalisation_map[3] = 1
        if INFO:
            print(" to", pnt[1], end='')

    return pnt, normalisation_map


def need_to_normalize(lon0, lat0, lon1, lat1, lon2, lat2):
    """
    """
    flag = False
    if (((math.copysign(1, lon0) != math.copysign(1, lon1)) and
            max(lon0, lon1) > 90.0) or
        ((math.copysign(1, lon1) != math.copysign(1, lon2)) and
            max(lon1, lon2) > 90.0) or
        ((math.copysign(1, lon2) != math.copysign(1, lon0)) and
            max(lon2, lon0) > 90.0)):
        flag = True
    return flag


def normalize_lon_lat(coo):
    """
    """
    if len(coo.shape) < 2:
        coo[0] += 360
    elif isinstance(coo, np.ndarray):
        coo[:, 0][coo[:, 0] < 0.0] = coo[:, 0][coo[:, 0] < 0.0]+360
    return coo


def resample_grand_arc(lon1, lat1, lon2, lat2, nodes_spacing):
    """
    :parameter nodes_spacing:
        Spacing between nodes created by the resampling procedure

    :raises ValueError:
        If ``nodes_spacing`` is too small i.e. the number of resampled points
        is too large
    """
    geod = Geod(ellps='sphere')
    flag = False
    if ((math.copysign(1, lon1) != math.copysign(1, lon2)) and
            max(lon1, lon2) > 90.0):
        flag = True
    _, _, dst = geod.inv(lon1, lat1, lon2, lat2)
    npts = int(round(dst / nodes_spacing))
    if npts > 200:
        raise ValueError('number of points used to resample is too large')
    coo = np.array(geod.npts(lon1, lat1, lon2, lat2, npts))
    if flag:
        coo = normalize_lon_lat(coo)
    return coo
