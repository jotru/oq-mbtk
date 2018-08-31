import unittest
import os
import json

import openquake.mbt.tools.fault_modeler.fault_source_modeler as fsm

# -----------------------------------------------------------------------------

class TestDatabaseIO(unittest.TestCase):

    geojson_file = os.path.join('Data', 'ne_asia_faults_rates.geojson')

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

    def test_fault_database(self):

        # Import the database
        fault_db = fsm.FaultDatabase()
        fault_db.import_from_geojson(self.geojson_file,
                                     param_map=self.param_map,
                                     select_list=[1, 2])

        # Adding a key/value to all faults
        fault_db.add_property('lower_seismogenic_depth', value=25)

        # Adding key/value to given faults
        fault_db.add_property('M_max', value=7., id=1)
        fault_db.add_property('M_max', value=7.5, id=2)

        fault_db.remove_property('name')

        # Target and reference files
        tar_file = os.path.join('Data', 'out_db.tar.geojson')
        ref_file = os.path.join('Data', 'out_db.ref.geojson')

        # Export the augmented database
        fault_db.export_to_geojson(tar_file)

        with open(tar_file, 'r') as f:
            tar = json.load(f)

        with open(ref_file, 'r') as f:
            ref = json.load(f)

        for fr, ft in zip(ref['features'], tar['features']):
            self.assertTrue(fr == ft)

    def test_build_model_from_db(self):

        # Import the database
        fault_db = fsm.FaultDatabase()
        fault_db.import_from_geojson(self.geojson_file,
                                     param_map=self.param_map,
                                     select_list=[1, 2])

        # Create and export the model
        fsm.build_model_from_db(fault_db,
                                xml_output='Data/fault_model.xml')

    def test_build_source_model_single_args(self):

        fsm.build_fault_model(geojson_file=self.geojson_file,
                              xml_output='Data/fault_model.xml',
                              black_list=[1,2,3],
                              param_map=self.param_map,
                              M_max=8.2,
                              lower_seismogenic_depth=30.)

    def test_build_source_model_dictionary(self):

        fsm.build_fault_model(geojson_file=self.geojson_file,
                              xml_output='Data/fault_model.xml',
                              param_map=self.param_map,
                              defaults={'upper_seismogenic_depth':10.,
                                        'lower_seismogenic_depth':30.})

    def test_build_source_model_config_file(self):

        fsm.build_fault_model(cfg_file='Data/config.ini')


if __name__ == "__main__":
    unittest.main()
