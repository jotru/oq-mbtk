#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2015-2018 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.
#
# Authors: Julio Garcia, Richard Styron, Valerio Poggi
# Last modify: 18/07/2018

# -----------------------------------------------------------------------------

import sys
import ast
import json
import configparser
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import openquake.mbt.tools.fault_modeler.fault_modeling_utils as fmu
from openquake.hazardlib.sourcewriter import write_source_model
from openquake.baselib import sap

# -----------------------------------------------------------------------------

# Parameters required from the fault modeler
key_list = ['source_id',
            'name',
            'average_dip',
            'average_rake',
            'net_slip_rate',
            'vert_slip_rate',
            'strike_slip_rate',
            'shortening_rate',
            'dip_dir',
            'dip_slip_rate',
            'slip_type']

"""
key_list += ['mfd',
             'M_max',
             'M_min',
             'fault_trace',
             'rupture_aspect_ratio',
             'rupture_mesh_spacing',
             'magnitude_scaling_relationship',
             'tectonic_region_type',
             'temporal_occurrence_model',
             'upper_seismogenic_depth',
             'lower_seismogenic_depth'] 
"""

# Conversion map for geojson input parameters
param_map = {'source_id': 'ogc_fid',
             'name': 'ns_name',
             'average_dip': 'ns_average_dip',
             'average_rake': 'ns_average_rake',
             'net_slip_rate': 'ns_net_slip_rate',
             'vert_slip_rate': 'ns_vert_slip_rate',
             'strike_slip_rate': 'ns_strike_slip_rate',
             'shortening_rate': 'ns_shortening_rate',
             'dip_dir': 'ns_dip_dir',
             'dip_slip_rate': 'ns_dip_slip_rate'}

# Default values
defaults = {}

# -----------------------------------------------------------------------------

def build_fault_model(cfg_file=None,
                      geojson_file=None,
                      xml_output=None,
                      project_name=None,
                      black_list=None,
                      select_list=None,
                      param_map=param_map,
                      defaults=defaults,
                      **kwargs):
    """
    Note: we have to set the priority when using both ini file and
    optional arguments
    """

    if cfg_file is not None:
        cfg_dict = read_config_file(cfg_file)

        if 'config' in cfg_dict:
            if 'geojson_file' in cfg_dict['config']:
                geojson_file = cfg_dict['config']['geojson_file']
            if 'xml_output' in cfg_dict['config']:
                xml_output = cfg_dict['config']['xml_output']

        if 'param_map' in cfg_dict:
            param_map.update(cfg_dict['param_map'])

        if 'defaults' in cfg_dict:
            defaults.update(cfg_dict['defaults'])

    for key in kwargs:
        defaults[key] = kwargs[key]

    # Import the fault database from geojson
    if geojson_file is not None:
        fault_db = fault_database()
        fault_db.import_from_geojson(geojson_file,
                                     black_list=black_list,
                                     select_list=select_list,
                                     param_map=param_map)
    else:
        print('Geojson file no specified')
        return

    # Create the fault source model in xml_format
    build_model_from_db(fault_db, xml_output, project_name, defaults=defaults)

# -----------------------------------------------------------------------------

def read_config_file(cfg_file):
    """
    """

    cfg = configparser.ConfigParser(dict_type=dict)
    cfg.read(cfg_file)

    cfg_dict = {}

    for key in ['config', 'param_map', 'defaults']:

        if cfg.has_section(key):
            cfg_dict[key] = {k:v for k, v in cfg.items(key)}

    return cfg_dict

# -----------------------------------------------------------------------------

def build_model_from_db(fault_db,
                        xml_output=None,
                        project_name=None,
                        defaults=defaults,
                        **kwargs):
    """
    """

    for key in kwargs:
        defaults[key] = kwargs[key]

    srcl = []

    for fl in fault_db.db:

        sfs_dict = fmu.construct_sfs_dict(fl, defaults=defaults)
        sfs = fmu.make_fault_source(sfs_dict)
        srcl.append(sfs)

    if xml_output is not None:
        # Write the final fault model
        write_source_model(xml_output, srcl, project_name)
    else:
        return srcl

# -----------------------------------------------------------------------------

class fault_database():
    """
    """

    def __init__(self, geojson_file=None):
        """
        """

        # Initialise an empty fault list
        self.db = []

        if geojson_file:
            self.import_from_geojson(geojson_file)

    def import_from_geojson(self, geojson_file,
                                  black_list=None,
                                  select_list=None,
                                  param_map=param_map,
                                  defaults=defaults):
        """
        """

        # Import database
        with open(geojson_file, 'r') as f:
            data = json.load(f)

            # Loop over faults
            for feature in data['features']:

                # Get fault properties
                fault = {}
                for key in key_list:

                    pm_key = param_map[key] if key in param_map else key
                    default = defaults[key] if key in defaults else None

                    if pm_key in feature['properties'].keys():
                        fault[key] = feature['properties'][pm_key]
                    else:
                        fault[key] = default

                # Process only faults in the selection list
                if select_list is not None:
                    if fault['source_id'] not in select_list:
                        continue

                # Skip further processing for blacklisted faults
                if black_list is not None:
                    if fault['source_id'] in black_list:
                        continue

                # Get fault geometry
                fault['trace_coordinates'] = feature['geometry']['coordinates']

                self.db.append(fault)

            # Temporary adjustment to make the code running....
            self.add_property('name', 'XXX')


    def add_property(self, property, value=None, id=None):
        """
        """

        for fault in self.db:
            if fault['source_id'] is id or id is None:
                fault[property] = value

# -----------------------------------------------------------------------------

def main(argv):
    """
    """

    p = sap.Script(build_fault_model)
    p.arg(name='geojson_file',
          help='Fault database in geojson format')
    p.arg(name='xml_output',
          help='Output xml containing the fault model')
    p.opt(name='project_name',
          help='Name of the current project', abbrev='-h')
    p.opt(name='black_list',
          help='List of fault IDs NOT to be processed [id1,id2]', 
          type=ast.literal_eval, abbrev='-h')
    p.opt(name='select_list',
          help='List of selected fault IDs to be processed [id1,id2]', 
          type=ast.literal_eval, abbrev='-h')

    if len(argv) < 1:
        print(p.help())
    else:
        p.callfunc()

if __name__ == "__main__":
    main(sys.argv[1:])
