#!/usr/bin/python

"""
"""

import os
import re
import sys
import h5py

from rtree.index import Index


def generator_function(mesh):
    """
    Generator function for quick loading the spatial index
    """
    for i in range(0, mesh.shape[0]):
        yield (i, (mesh[i, 0], mesh[i, 1], mesh[i, 0], mesh[i, 1]), None)


def create_mesh_spatial_index(filename):
    """
    :param str filename:
        The name of the .hdf5 file containing the regular mesh
    """
    #
    # create the spatial index filename
    fname = re.split('\.', os.path.basename(filename))[0]
    fpath = os.path.dirname(filename)
    #
    # mesh
    f = h5py.File(filename, 'r')
    mesh = f['centroids'][:]
    #
    # create the spatial index
    filename_rtree = os.path.join(fpath, fname)
    print('Processing:', filename_rtree)
    # p = rtree.index.Property()
    # p.filename = filename_rtree
    sidx = Index(filename_rtree, generator_function(mesh))


def main(argv):
    """
    argv[0] - name of the .hdf5 file
    """
    create_mesh_spatial_index(argv[0])


if __name__ == '__main__':
    main(sys.argv[1:])
