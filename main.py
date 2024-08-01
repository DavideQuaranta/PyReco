import sys

import pandas as pd
import numpy as np
import awkward as ak

sys.path.append('PyReco')
from PyReco.chamber import Chamber, create_chambers
from PyReco.gas import Gas
from PyReco.utils import read_configuration_file, print_event
from PyReco.rootio import open_tree


# Read configuration file with chambers settings
config = read_configuration_file("config.yaml")

# Open ROOT TTree containing data
run_number = sys.argv[1]
raw_dir = config["paths"]["raw_dir"]
filename = "".join([raw_dir, "run", run_number, "_xtalk1_xtalk2", ".root"])
treename = "apv_raw"
tree = open_tree(filename, treename)

# Define Gas
gas = Gas(config["gas"]["name"] , config["gas"]["name"], config["gas"]["drift_speed"])
print(f"- {gas}")

# Create chambers
chambers = create_chambers(config)


# Define the chunk size
chunk_size = 1000  # Adjust based on your memory constraints

# Number of events to process
nentries = int(tree.num_entries*config['general']['max_events'])

# Iterate over the tree in chunks
for data in tree.iterate(entry_start=0, entry_stop=nentries, step_size=chunk_size):
    for event in data:

        print(event['evt'])
        
        if config['general']['verbose'] == True : print_event(event)

        for chamb in chambers:
            for layer in chamb.strips.keys():
                for readout in chamb.strips[layer].keys():
                    chamber_cut = event.mmChamber == chamb.name
                    layer_cut = event.mmLayer == int(layer)
                    readout_cut = event.mmReadout == int(readout)
                    chamb.strips[layer][readout] = np.array(event.mmStrip[chamber_cut & layer_cut & readout_cut])



    


