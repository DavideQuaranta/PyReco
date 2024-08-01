import sys

import pandas as pd

sys.path.append('PyReco')
from PyReco.chamber import Chamber, create_chambers
from PyReco.gas import Gas
from PyReco.utils import read_configuration_file
from PyReco.rootio import open_tree


# Read configuration file with chambers settings
config = read_configuration_file("config.yaml")

# Open ROOT TTree containing data
run_number = sys.argv[1]
raw_dir = config["paths"]["raw_dir"]
filename = "".join([raw_dir, "run", run_number, "_xtalk1_xtalk2", ".root"])
treename = "apv_raw"
tree = open_tree(filename, treename)

# Create chambers
gas = Gas(config["gas"]["name"] , config["gas"]["name"], config["gas"]["drift_speed"])
print(f"- {gas}")

chambers = create_chambers(config)


# print(chambers)

# # Define the chunk size
# chunk_size = 1  # Adjust based on your memory constraints

# # Iterate over the tree in chunks
# for data in tree.iterate(entry_start=0, entry_stop=10, step_size=chunk_size):
#     print(f"\n#### EVENT {data['evt']} ####")


