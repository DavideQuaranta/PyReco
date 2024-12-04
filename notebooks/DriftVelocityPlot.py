import numpy as np
import pandas as pd
import mplhep as hep
import matplotlib.pyplot as plt

from LatexConstants import *
from scipy.optimize import curve_fit
from scipy.interpolate import make_interp_spline


# open file with drift speeds
filename = '../reco/drift_speed.txt'
speed = pd.read_csv(filename, header=None, sep=' ', names=['run', 'speed', 'error'])

# Read Garfield data
garfield = pd.read_csv('../TestBeamMMTPC/ExMe_DriftTime_Garfield.txt', header=None, names=['E', 'v'], sep=' ')

# Drift voltage scan runs
gap = 50.128
err_gap = 1 #mm
real_gap = gap
print(real_gap)
runs = ['2125', '2126', '2124', '2127', '2128', '2129']
V = [2500, 2750, 3000, 3250, 3500, 3750]
E = [(v / (0.1 * real_gap)) / 1000 for v in V]


# match units used on x-axis by garfield 
drift_scan = np.array([100 * (real_gap/gap) * speed['speed'][speed['run'] == int(run)].values[0] for run in runs])
drift_scan_err = np.array([100 * (real_gap/gap) * speed['error'][speed['run'] == int(run)].values[0] for run in runs])

# try to include +/- 1mm error on the drift gap thickness
rel_t_err = drift_scan_err/drift_scan
rel_gap_err = np.repeat(err_gap/gap, repeats=len(drift_scan))

new_error = drift_scan * np.sqrt(rel_gap_err**2 + rel_t_err**2)


print(drift_scan)
print(new_error)

plt.style.use(hep.style.ATLAS)

fig, (ax1, ax2) = plt.subplots(2, 1, height_ratios=[6,2], sharex=True)
# plt.figure(figsize=(0.5*textwidth, 0.5*0.5*textheight), dpi = 50)

ax1.grid(zorder = 0, alpha = 0.4)
ax1.errorbar(E, drift_scan, yerr = new_error, color = 'k', lw = 1, capsize = 3,markersize = 10,fmt = '.', label = 'ExMe', zorder =2, alpha = 0.7)

garfield = pd.read_csv('../TestBeamMMTPC/ExMe_DriftTime_Garfield.txt', header = None, names = ['E', 'v'], sep  =' ')
ax1.plot(garfield['E'], garfield['v'], marker = '.', markersize = 5,zorder = 1, label = 'Garfield', lw = 1, alpha = 0.7)
ax1.set_ylabel(r'Drift Speed [cm/$\mu$s]', fontsize = fontsize)
ax1.legend(loc = 'upper right', fontsize = fontsize)
ax1.tick_params(top = False, labeltop = False,  which = 'both', labelsize = fontsize)


from scipy.interpolate import make_interp_spline
x = garfield['E'][0:15]
y = garfield['v'][0:15]
bspl = make_interp_spline(x, y, k=3)

x1 = np.linspace(0.4,0.8, 100)
ax2.grid(zorder = 0, alpha = 0.4)
ax2.errorbar(E, 100*(drift_scan-bspl(E))/bspl(E),yerr = [100*new_error[i]/bspl(E[i]) for i in range(len(drift_scan_err))], color = 'k', lw = 1, capsize = 3,markersize = 10,fmt = '.', zorder = 2)
ax2.set_xlabel(r'E[V/cm]$\times 10^3$', fontsize = 15)
ax2.set_ylabel(r'Residues $\%$', fontsize = 15)
ax2.set_ylim(-10, 10)
ax2.tick_params(labelsize = fontsize)
ax2.hlines(0, 0, 1.6, color = 'red', linestyle = 'solid', alpha = 0.8, lw = 1, zorder = 1)

print(f'squared res = {np.sum((drift_scan-bspl(E))**2/new_error**2)}')

def E2kV(x):
    return x*1000*(0.1*50.128) / 1000


def kV2E(x):
    return (x/(0.1*50.128))/1000

print(E2kV(E))

# ax = plt.gca()
secax = ax1.secondary_xaxis('top', functions=(E2kV, kV2E))
secax.set_xlabel(r'Voltage [kV]', fontsize = fontsize)
secax.tick_params(labelsize = fontsize)


# plt.tick_params(axis='both', which='major', labelsize=fontsize)
# ax = plt.gca()
ax1.xaxis.offsetText.set_fontsize(fontsize)
ax1.set_ylim(5, 12)
ax1.set_xlim(0.21,1.5)
ax1.annotate(text = r'ArCF4Iso (88:10:2)', xy = (1.1, 1))


plt.tight_layout()
# plt.savefig('DriftSpeed_w_Garfield.pdf', dpi = 600)
# plt.savefig('DriftSpeed_w_Garfield.png', dpi = 600)

plt.show()
