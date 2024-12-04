import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import pickle
from scipy.optimize import curve_fit
from scipy.stats import norm
from LatexConstants import *

 
def gaussian(x, mu, sigma, N):
    return N*np.exp(-(x-mu)**2/(2*sigma**2))

def bigaus(x, N, mu, sigma, N1, mu1, sigma1):
    # N = 1/np.sqrt(2*np.pi*sigma**2)
    return N*np.exp(-(x-mu)**2/(2*sigma**2)) + N1*np.exp(-(x-mu1)**2/(2*sigma1**2))

def normal(x, N, mu, sigma):
    # N = 1/np.sqrt(2*np.pi*sigma**2)
    return N*np.exp(-(x-mu)**2/(2*sigma**2))

 
# open file with data
with open("../reco/histos/2124/Z ResiduesP2_M01.pkl", "rb") as file:
    res = pickle.load(file)


plt.style.use(hep.style.ATLAS)
fig, ax = plt.subplots(figsize=(textwidth, 0.5*textheight), dpi = 100) #plt.figure(figsize=(6,4), dpi = 100)


limits = (-10,10)

n, bins, patches = ax.hist(res, bins = 200, histtype='step', color = 'blue', range = limits, lw = 2)
w = 0.5*(bins[1]-bins[0])
a = bins[0]
b = bins[-2]
bin_centers = np.linspace(a+w, b+w, len(bins)-1)
par, cov = curve_fit(bigaus, bin_centers, n, p0 = (1000, 0, 1, 100, 0, 0.1))

# plt.errorbar(bin_centers, n, yerr = np.sqrt(n), fmt = '.', color = 'k')

if par[2] < par[5]:
    par_core = par[0:3] 
    par_tail = par[3:]
    err_tail = [np.sqrt(cov[i][i]) for i in range(3, 6)]
    err_core = [np.sqrt(cov[i][i]) for i in range(0, 3)]
else:
    par_tail = par[0:3] 
    par_core = par[3:]
    err_core = [np.sqrt(cov[i][i]) for i in range(3, 6)]
    err_tail = [np.sqrt(cov[i][i]) for i in range(0, 3)]

x = np.linspace(limits[0],limits[1], 1000)
ax.plot(x, bigaus(x, *par), color = 'red', lw = 1, label = 'Total')
ax.plot(x, normal(x, *par_core), lw = 1, label = 'Core')
ax.plot(x, normal(x, *par_tail), lw = 1, label = 'Tail')
ax.set_xlim(limits)
ax.annotate(text = r' $\mu_{tail}$' + f' = {abs(par_tail[1]):.3f} $\pm$ {err_tail[1]:.3f}', xy = (-8.9, 2.2*10**4), fontsize = 0.8*fontsize)
ax.annotate(text = r'$\mu_{core}$' + f' = {abs(par_core[1]):.3f} $\pm$ {err_core[1]:.3f}', xy = (-8.7, 2.05*10**4), fontsize = 0.8*fontsize)
ax.annotate(text = r' $\sigma_{tail}$' + f' = {abs(par_tail[2]):.3f} $\pm$ {err_tail[2]:.3f}', xy = (-8.9, 1.9*10**4), fontsize = 0.8*fontsize)
ax.annotate(text = r'$\sigma_{core}$' + f' = {abs(par_core[2]):.3f} $\pm$ {err_core[2]:.3f}', xy = (-8.7, 1.75*10**4), fontsize = 0.8*fontsize)
ax.set_xlabel('Z residuals [mm]', fontsize = fontsize)
ax.set_ylabel('Entries', fontsize = fontsize)
ax.legend(fontsize = 0.8*fontsize)
plt.savefig('ZResolution_2124_timefit.png', dpi = 600)
plt.savefig('ZResolution_2124_timefit.pdf', dpi = 600)
# plt.show()

