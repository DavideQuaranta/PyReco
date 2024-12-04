import os
import sys
import pickle

import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

from scipy.optimize import curve_fit
from scipy import stats
from scipy.stats import norm
from LatexConstants import *


run_number = sys.argv[1]
title = "Centroid"

histo_dir = "".join(["../reco/histos/", run_number, "/"])
histo_filenames = [name for name in os.listdir(histo_dir) if name.startswith(title + "P")]

chambers = ["TMM", "ExMe"]

data = {}

for i, histo in enumerate(histo_filenames):
    with open(histo_dir + histo, "rb") as file:
        data[chambers[i]] = np.array(pickle.load(file))


# nbins = 100
# n, bins, patches = plt.hist(data["TMM"], bins = nbins, density = True)

mean = np.mean(data["ExMe"])
std = np.std(data["ExMe"])
cut = (data["ExMe"] < mean+0.5*std) & (data["ExMe"] > mean-0.3*std)
sliced_data = data["ExMe"][cut]

result = stats.fit(norm, sliced_data, bounds=[(0, 500), (0, 60)])

result.plot()
print(result.params)
plt.show()
