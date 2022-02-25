"""
This module defines a base class for reference grid to use for calculations
"""

import abc
import numpy as np

class ReferenceGrid(object):
    """
    A base class for the calculation reference grid.

    :param grid_spacing:
        Spacing between nodes [decimal degrees or km]
    :parameter bounds:
        Lower left and upper right nodes represented in terms of a list 
        [lon_min, lat_min, lon_max, lat_max] 
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, grid_spacing, bounds=[-180.0,-90.0,180.0,90.0]):
        if not grid_spacing > 0.0:
            raise ValueError('grid spacing must be positive')
        self.grid_spacing = grid_spacing
        self.bounds = bounds
        self._round_bounds()

    @abc.abstractmethod
    def iter_nodes(self,within_landmasses_flag):
        """
        Get a generator object that yields calculation grid nodes 

        :param within_landmasses_flag:
            Boolean variable. If true the method return only nodes 
            within landmasses are provided 
        :returns:
            A list of nodes in [x,y] format
        """
    
    def _round_bounds(self):
        """
        This method makes the bounds a multiple of the grid spacing
        """
        arr = np.array(self.bounds)
        arr[0:1] = np.floor(arr[0:1]/self.grid_spacing)*self.grid_spacing
        arr[2:3] = np.ceil(arr[2:3]/self.grid_spacing)*self.grid_spacing
        self.bounds = list(arr)   

