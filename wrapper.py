# launcher.py
import subprocess
from tqdm import tqdm
import os
import yaml
import re

# Read the configuration file
with open('config.yaml', 'r') as file:
    config = yaml.safe_load(file)

# List of scripts and their respective arguments
script_name = 'PyReco.py'
# args = ['2136', '2124', '2130', '2131', '2132', '2133', '2134', '2135', '2137', '2121']
# args = ['2130', '2131', '2132', '2133', '2134', '2135', '2136, '2137', '2121', '2124]
# args = ['2135', '2121', '2130', '2131', '2132', '2133', '2134', '2136', '2137', '2124']

# get the runs in the raw dir to process
files = os.listdir(config['paths']['raw_dir'])
run_numbers = [re.findall(r'\d\d\d\d', file)[0] for file in files]

# Function to run the script with a single argument
def run_script(script_path, arg):
    command = ['python', script_path, arg]
    result = subprocess.run(command, capture_output=True, text=True)
    print(f"Output of {script_path} with argument '{arg}':")
    print(result.stdout)
    print(f"Errors of {script_path} with argument '{arg}':")
    print(result.stderr)

# Run the script with each argument in sequence
for run in run_numbers:
    print(f'\nRunning {script_name} for run {run}...\n')
    run_script(script_name, run)
