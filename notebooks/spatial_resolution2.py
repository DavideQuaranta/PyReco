import os
import sys
import pickle

import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

from scipy.optimize import curve_fit
from scipy import stats
from scipy.stats import norm
from uncertainties import ufloat

from LatexConstants import *

def fit_histogram(func, n, bins, p0):
    y = n
    x = bins[0:-1] + 0.5*(bins[1]-bins[0])
    par, cov = curve_fit(func, x, y, p0=p0)
    return par, np.sqrt(np.diag(cov))

def gaussian(x, N, mu, sigma):
    return N*np.exp(-(x-mu)**2/(2*sigma**2))

def bigaus(x, N, mu, sigma, N1, mu1, sigma1):
    return gaussian(x, N, mu, sigma) + gaussian(x, N1, mu1, sigma1)

def print_fit_results(par, err):
    for i, p in enumerate(par):

        print(f"par[{i}] = {ufloat(p, err[i]):.2uL}")


if __name__ == "__main__":
    
    run_number = sys.argv[1]
    title = "Centroid"
    histo_dir = "".join(["../reco/histos/", run_number, "/"])
    histo_filenames = [name for name in os.listdir(histo_dir) if name.startswith(title)]
    histo_filenames=histo_filenames[3:]
    histo_filenames = np.roll(histo_filenames, 1)
    print(histo_filenames)

    names = ["TMM (x)", "TMM (y)", "ExMe (x)"]
    colors = ["royalblue", "royalblue", "darkorange"]
    fitcolor = ["red", "red", "blue"]
    # names = ["ExMe"]

    data = {}

    for name, histo in zip(names, histo_filenames):
        with open(histo_dir + histo, "rb") as file:
            data[name] = np.array(pickle.load(file))
            print(len(data[name]))


    plt.style.use(hep.style.ATLAS)
    plt.figure(figsize=(2*textwidth, 0.5*textheight))
    
    for i,name in enumerate(names):
        nbins = 100
        mean = np.mean(data[name])
        median = np.median(data[name])
        std = 0.5*np.std(data[name])
        lims  = [mean-2*std, mean*2+std]

        plt.subplot(1, len(histo_filenames), i+1)
        n, bins, patches = plt.hist(data[name], bins=nbins, range = lims, histtype='step', color = colors[i], label = f'{name}')
        par, errors = fit_histogram(gaussian, n, bins, p0 = [100,mean,std])
        x = np.linspace(lims[0],lims[1], 500)
        plt.plot(x, gaussian(x, *par), color = fitcolor[i], lw = 2)
        print_fit_results(par, errors)

        # if abs(par[5]) < abs(par[2]):
        #     errors = np.roll(errors, 3)
        #     par = np.roll(par, 3)

        # plt.plot(x, gaussian(x, *par[0:3]), color = 'royalblue', label = 'core', lw = 1)
        # plt.plot(x, gaussian(x, *par[3:]), color = 'darkorange', label = 'tail', lw = 1)


        plt.xlabel(r"Centroid  [mm]", fontsize = fontsize)
        plt.xlim(par[1]-8*par[2], par[1]+8*par[2])
        plt.ylim(0, 1.3*par[0])
        y = plt.gca().get_ylim()[1]*(1-0.1)

        plt.text(mean + 0.4*std, y, s = r'$\mu$ = ' + f"{par[1]:.2f}" + r' $\pm$ ' f'{errors[1]:.2f}', fontsize = 0.6*fontsize)
        plt.text(mean + 0.4*std, y*(1-0.06), s = r'$\sigma$ = ' + f"{abs(par[2]):.2f}" + r' $\pm$ ' f'{errors[2]:.2f}', fontsize = 0.6*fontsize)
        # plt.text(mean + 0.4*std, y, s = r'$\mu_{core}$ = ' + f"{par[1]:.2f}" + r' $\pm$ ' f'{errors[1]:.2f}', fontsize = 0.5*fontsize)
        # plt.text(mean + 0.4*std, y*(1-0.05), s = r'$\sigma_{core}$ = ' + f"{abs(par[2]):.2f}" + r' $\pm$ ' f'{errors[2]:.2f}', fontsize = 0.5*fontsize)
        # plt.text(mean + 0.4*std, y*(1-2*0.05), s = r'$\mu_{tail}$ = ' + f"{par[4]:.2f}" + r' $\pm$ ' f'{errors[4]:.2f}', fontsize = 0.5*fontsize)
        # plt.text(mean + 0.4*std, y*(1-3*0.05), s = r'$\sigma_{tail}$ = ' + f"{abs(par[5]):.2f}" + r' $\pm$ ' f'{errors[5]:.2f}', fontsize = 0.5*fontsize)


        if i == 0:
            plt.ylabel("Entries", fontsize = fontsize)
        plt.tick_params(labelsize = 0.9*fontsize)
        # plt.tight_layout()
        plt.legend(loc = 'upper left',fontsize = 0.7*fontsize)
        plt.savefig("centroid" + f"_run{run_number}.pdf", dpi = 600)
        plt.savefig("centroid" + f"_run{run_number}.png", dpi = 600)
    plt.show()





