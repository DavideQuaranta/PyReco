# import libraries
import numpy as np
import pickle
import sys
import mplhep as hep
import matplotlib.pyplot as plt


from scipy.optimize import curve_fit
from LatexConstants import *


def normal(x, N, mu, sigma):
    # N = 1/np.sqrt(2*np.pi*sigma**2)
    return N*np.exp(-(x-mu)**2/(2*sigma**2))


# read file
run = sys.argv[1]
title = 'xhalf'
chamb = 'P0_RM1'
filename = '../reco_spatres/histos/'+run+'/'+title+chamb+'.pkl'

with open(filename, "rb") as file:
    xhalf = pickle.load(file)

# plot constants
nbins = 100
# hrange = (-0.5*np.std(xhalf), 0.5*np.std(xhalf))
hrange = (-20, 20)


plt.style.use(hep.style.ATLAS)
fig, ax = plt.subplots(figsize = (textwidth, 0.5*textheight))
n, bins, patches = ax.hist(xhalf - np.mean(xhalf), bins = nbins, range = hrange, histtype = 'step', color = 'blue')
bin_center = bins[0:-1] + 0.5*(bins[1]-bins[0])

par, cov = curve_fit(normal, bin_center, n)
ax.plot(bin_center, normal(bin_center, *par), color = 'red', lw = 1)
ax.annotate(text = rf' $\sigma$ = {abs(par[2]):.2f} $\pm$ {np.sqrt(cov[2][2]):.2f}', 
            xy = (hrange[0]+0.6*(hrange[1]-hrange[0]), np.max(n)), fontsize = fontsize)


plt.title(chamb)
ax.set_xlabel(r'$x_{half}$ [mm]', fontsize = fontsize)
ax.set_ylabel('Entries', fontsize = fontsize)
plt.savefig('xhalf_TMM_2017.pdf', dpi  = 600)
plt.savefig('xhalf_TMM_2017.png', dpi  = 600)
plt.tight_layout()
plt.show()


# read file
run = sys.argv[1]
title = 'xhalf'
chamb = 'P2_M01'
filename = '../reco_spatres/histos/'+run+'/'+title+chamb+'.pkl'

with open(filename, "rb") as file:
    xhalf = pickle.load(file)

# plot constants
nbins = 100
# hrange = (-0.5*np.std(xhalf), 0.5*np.std(xhalf))
hrange = (-20, 20)


plt.style.use(hep.style.ATLAS)
fig, ax = plt.subplots(figsize = (textwidth, 0.5*textheight))
n, bins, patches = ax.hist(xhalf - np.mean(xhalf), bins = nbins, range = hrange, histtype = 'step', color = 'blue')
bin_center = bins[0:-1] + 0.5*(bins[1]-bins[0])

par, cov = curve_fit(normal, bin_center, n)
ax.plot(bin_center, normal(bin_center, *par), color = 'red', lw = 1)
ax.annotate(text = rf' $\sigma$ = {abs(par[2]):.2f} $\pm$ {np.sqrt(cov[2][2]):.2f}', 
            xy = (hrange[0]+0.6*(hrange[1]-hrange[0]), np.max(n)), fontsize = fontsize)


plt.title(chamb)
ax.set_xlabel(r'$x_{half}$ [mm]', fontsize = fontsize)
ax.set_ylabel('Entries', fontsize = fontsize)
plt.savefig('xhalf_ExMe_2017.pdf', dpi  = 600)
plt.savefig('xhalf_ExMe_2017.png', dpi  = 600)
plt.tight_layout()
plt.show()

# read file
run = sys.argv[1]
title = 'xhalf_res'
chamb = 'P2_M01'
filename = '../reco_spatres/histos/'+run+'/'+title+chamb+'.pkl'

with open(filename, "rb") as file:
    xhalf = pickle.load(file)

# plot constants
nbins = 200
# hrange = (-0.5*np.std(xhalf), 0.5*np.std(xhalf))
hrange = (-40, 40)

plt.style.use(hep.style.ATLAS)
fig, ax = plt.subplots(figsize = (textwidth, 0.5*textheight))
n, bins, patches = ax.hist(xhalf - np.mean(xhalf), bins = nbins, range = hrange, histtype = 'step', color = 'blue')
bin_center = bins[0:-1] + 0.5*(bins[1]-bins[0])

par, cov = curve_fit(normal, bin_center, n)
ax.plot(bin_center, normal(bin_center, *par), color = 'red', lw = 1)
ax.annotate(text = rf' $\sigma$ = {abs(par[2]):.2f} $\pm$ {np.sqrt(cov[2][2]):.2f}', 
            xy = (hrange[0]+0.6*(hrange[1]-hrange[0]), np.max(n)), fontsize = 0.9*fontsize)

ax.set_xlabel(r'$x_{half}$ residues [mm]', fontsize = 0.9*fontsize)
ax.set_ylabel('Entries', fontsize = 0.9*fontsize)
plt.savefig('xhalf_residues_2017.pdf', dpi  = 600)
plt.savefig('xhalf_residues_2017.png', dpi  = 600)
plt.tight_layout()
plt.show()

