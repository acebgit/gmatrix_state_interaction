__author__ = 'Antonio Cebreiro'

import sys
from projection_method.parsers_OPO.gmatrix_classes import OutToJsonConverter, GTensorConfig, GTensorPipeline

# Get the .out file from command line arguments
if len(sys.argv) != 2:
    print("Usage: python new_example.py <qchem_output_file.out>")
    sys.exit(1)
qchemout = sys.argv[1] 

# Obtain the json file from the .out file
converter = OutToJsonConverter(qchemout)
converter.read_file()
converter.parse()
converter.save_json()

# Create a configuration object
config = GTensorConfig()
config.state_selection = 0  # 0: use "state_ras" ; 1: use all states_selected
config.initial_states = [1, 5, 6]  # select only states 1,2,3,4

# Create the pipeline using this configuration
pipeline = GTensorPipeline(qchemout, config=config)
pipeline.load_data()                    # Step 1: load the JSON
pipeline.select_states()                # Step 2: select the states
pipeline.run_gshift_calculation()       # Step 3: calculate g-shift


