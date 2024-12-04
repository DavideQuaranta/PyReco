# PyReco

This repository contains the code i used to analyze the data of May 2024 test beam. It is structured like the following:

- **PyReco.py** : it is the main program. To run it navigate to the Pyreco directory and type:
    ```sh
    python3 PyReco.py 2124
    ```
    for example if you want to process run 2124.
- **PyMegas.py** : it is the library where all the function and the classes used in PyReco.py are defined.
- **config.yaml** : this is a configuration file in YAML format where alll the parameters of the chambers are included, togheter with some analysis parameters to change them easily from the same place. Feel free to add in there whatever you need.
- **wrapper.py** : a utility script to run PyReco.py on multiple runs at once.

Togheter with the main code you find also a series of folders, in which the results of PyReco.py will be saved. Keep in mind, however that the directory where you save your results can be changed easily from the config.yaml file.

- **data** : ideally this is where you put your ROOT file. The subdirectory `data\raw` will contain the raw apv files, while `data\timefit` will contain the results of the macros contained in the `macros` folder, needed to perform the time fit on the chamber strips (see the README.md file in that folder).
- **macros** : contains the ROOT macros to perform the Fermi-Dirac fit on the strips.
- **output** : this is where the results of Pyreco.py are saved, `output\plots` contains the debug plots that the program produces at the end of its running; `output\histos` contains the produced in histograms in `.pickle` format, for them to be used in further analysis. In the output folder the `.txt` files containing the results on drift speed, efficiencies and alignement are also saved.
- **notebooks** : contains some jupyter notebooks that i used for specific tasks.


## Installation

- Install the requirements by running: 
    ```sh
    pip install -r requirements.txt
    
- Clone this repository by clicking on the green button "Code" and copying the HTTPS link, then in your terminal:
    ```sh
    # Clone the repository
    git clone https://github.com/username/repository-name.git PyReco

    # Navigate into the repository
    cd PyReco

    # Create a new branch
    git checkout -b my-new-branch

    # Make your changes to the files

    # When you want to save a version with your changes
    git add .

    # Commit the changes
    git commit -m "Describe your changes"

    # Push your branch to GitHub
    git push -u origin my-new-branch
    ```



