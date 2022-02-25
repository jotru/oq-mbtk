
import os
import re
import h5py
import struct
import numpy as np


def get_centroids(filename, outfilename):
    """
    :param filename:
    :param outfilename:
    """

    fi = h5py.File(filename, 'r')
    cens = fi['centroids'][:]
    #
    #
    fo = open(outfilename, 'wb')
    #
    #
    for i in range(cens.shape[0]):
        binary = struct.pack('ff', cens[i, 0], cens[i, 1])
        fo.write(binary)
    #
    #
    fo.close()
    fi.close()


def get_triangles(filename, outfilename):
    """
    :param filename:
    :param outfilename:
    """

    fi = h5py.File(filename, 'r')
    nodes = fi['nodes']
    vtx = nodes['vertexes'][:]
    facets = fi['facets'][:]

    fo = open(outfilename, 'wb')
    #
    #
    binary = struct.pack('ff', np.nan, np.nan)
    fo.write(binary)
    #
    #
    for i in range(facets.shape[0]):
        for j in range(3):
            idx = facets[i, j]
            #
            #
            binary = struct.pack('ff', vtx[idx, 0], vtx[idx, 1])
            fo.write(binary)
        #
        #
        binary = struct.pack('ff', np.nan, np.nan)
        fo.write(binary)
    fo.close()
    fi.close()


def main():
    filename = './../out/trigrd_split_4_spacing_432.hdf5'
    filename = './../out/trigrd_split_6_spacing_108.hdf5'
    filename = './../out/trigrd_split_5_spacing_216.hdf5'
    outpath = './../out/'

    tfname = re.split('\.', os.path.basename(filename))[0] + '_tri.bin'
    cfname = re.split('\.', os.path.basename(filename))[0] + '_cen.bin'
    get_triangles(filename, os.path.join(outpath, tfname))
    get_centroids(filename, os.path.join(outpath, cfname))


if __name__ == "__main__":
    main()
