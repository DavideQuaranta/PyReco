import yaml
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from sklearn.cluster import DBSCAN
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
from iminuit import Minuit
from iminuit.cost import LeastSquares
from numba import njit
import re
import time

# Read the configuration file
with open('config.yaml', 'r') as file:
    config = yaml.safe_load(file)

# Finds the right TTree in the TFile
def get_highest_apv_raw(filenames):
    
    # Regular expression to match filenames of the pattern apv_raw{i}
    pattern = re.compile(r'apv_raw;(\d+)')
    max_number = -1
    highest_filename = None
    
    for filename in filenames:
        match = pattern.match(filename)
        if match:
            number = int(match.group(1))
            if number > max_number:
                max_number = number
                highest_filename = filename
    
    return highest_filename


class Gas:
    def __init__(self, name, composition, drift_speed):
        self.name = name
        self.composition = composition
        self.drift_speed = drift_speed
    
    def GetName(self):
        return self.name
    
    def GetComposition(self):
        return self.composition
    
    def GetDriftSpeed(self):
        return self.drift_speed
    
    def Print(self):
        print(f"Gas : {self.name} {self.composition}")
        print(f" - drift velocity : {self.drift_speed} {config['units']['speed']}")


class Chamber:
    def __init__(self, name, gas):
        self.name = name
        self.gas = gas
        self.__nstrips = 0
        self.__pitch = 0
        self.__MinStripCharge = 0
        self.__MaxStripCharge = 0
        self.__MinStripsInCluster = 0
        self.__MaxStripsInCluster = 0
        self.__MinClusterCharge = 0
        self.__MaxNConsecutiveHolesInCluster = 0
        self.__MaskedStrips = []

    def SetNStrips(self, nstrips):
        self.__nstrips = nstrips

    def SetPitch(self, pitch):
        self.__pitch = pitch

    def SetMinStripCharge(self, cmin):
        self.__MinStripCharge = cmin

    def SetMaxStripCharge(self, cmax):
        self.__MaxStripCharge = cmax
    
    def SetMinStripsInCluster(self, smin):
        self.__MinStripsInCluster = smin

    def SetMaxStripsInCluster(self, smax):
        self.__MaxStripsInCluster = smax

    def SetMinClusterCharge(self, cmin):
        self.__MinClusterCharge = cmin
    
    def SetMaxNConsecutiveHolesInCluster(self, hmax):
        self.__MaxNConsecutiveHolesInCluster = hmax

    def SetMaskedStrips(self, masked):
        self.__MaskedStrips = masked

    def GetMaskedStrips(self):
        return self.__MaskedStrips

    def GetName(self):
        return self.name
    
    def GetGas(self):
        return self.gas
    
    def GetNStrips(self):
        return self.__nstrips

    def GetPitch(self):
        return self.__pitch

    def GetMinStripCharge(self):
        return self.__MinStripCharge
    
    def GetMaxStripCharge(self):
        return self.__MaxStripCharge
    
    def GetMinStripsInCluster(self):
        return self.__MinStripsInCluster 

    def GetMaxStripsInCluster(self):
        return self.__MaxStripsInCluster

    def GetMinClusterCharge(self):
        return self.__MinClusterCharge 
    
    def GetMaxNConsecutiveHolesInCluster(self):
        return self.__MaxNConsecutiveHolesInCluster
    
    def Print(self):
        print(f'Chamber : {self.name}')
        print(f' - Gas : {self.gas.GetName()}')
        print(f' - N° strips : {self.__nstrips}')
        print(f' - Pitch : {self.__pitch}')
    

def CreateChambers(config, nchambers, verbose = False):

    gas = Gas(config['gas']['name'], config['gas']['composition'], config['gas']['drift_speed'])
    
    if verbose == True:
        gas.Print()
    
    names = list(config.keys())[0:nchambers]
    chambers = []
    
    for name in names:
        chamb = Chamber(name, gas)
        chamb.SetNStrips(config[name]['nstrips'])
        chamb.SetPitch(config[name]['pitch'])
        chamb.SetMinStripCharge(config[name]['MinStripCharge'])
        chamb.SetMaxStripCharge(config[name]['MaxStripCharge'])
        chamb.SetMinStripsInCluster(config[name]['MinStripsInCluster'])
        chamb.SetMinClusterCharge(config[name]['MinClusterCharge'])
        chamb.SetMaxStripsInCluster(config[name]['MaxStripsInCluster'])
        chamb.SetMaxNConsecutiveHolesInCluster(config[name]['MaxNConsecutiveHolesInCluster'])
        chamb.SetMaskedStrips([i for i in range(config[name]['MaskedStrips'][0], config[name]['MaskedStrips'][1])])
        chambers.append(chamb)
        if verbose == True:
            chamb.Print()
    
    return chambers

class Strip:
    def __init__(self, id, charge, time = 0, time_error = 0, chamber = None):
        self.id = id
        self.charge = charge
        self.chamber = chamber
        self.time = time
        self.time_error = time_error
    
    def isGoodStrip(self):
        # print(self.charge)
        if self.chamber.GetMinStripCharge() < self.charge:
            if self.chamber.GetMaxStripCharge() > self.charge:
                if self.id not in self.chamber.GetMaskedStrips():
                    if (self.time > 100) & (self.time < 400): # added to fit only central region
                        if self.time_error < 10: # the value is assigned after wathing the distribution of fit times error
                            return True
        else:
            return False
    
    def Print(self):
        print(f"Strip : ID = {self.id}; Charge  = {self.charge}; Chamber = {self.chamber.GetName()}")
    
    def GetID(self):
        return self.id
    
    def GetCharge(self):
        return self.charge
    
    def ParentChamber(self):
        return self.chamber
    
    def SetParentChamber(self, chamb):
        self.chamber = chamb
    
    def GetTime(self):
        return self.time
    def SetTime(self, t):
        self.time = t
    




class Cluster:
    def __init__(self, chamber):
        self.chamber = chamber
        self.__strips = []
    
    def AddStrip(self, strip):
        self.__strips.append(strip)

    def __AddAllStrips(self, strips):
        self.__strips = strips
    
    def __add__(self, other):
        if isinstance(other, Cluster):
            if self.chamber == other.chamber:
                newCluster = Cluster(other.chamber);
                newCluster.__AddAllStrips(self.__strips + other.__strips)
                return newCluster
        return NotImplemented
    
    def isGoodCluster(self):
        min = self.chamber.GetMinStripsInCluster()
        if self.Size() >= min:
            if self.Size() <=self.chamber.GetMaxStripsInCluster():
                if self.GetTotalCharge() >= self.chamber.GetMinClusterCharge():
                    return True
        else:
            return False
        
    def Size(self):
        return len(self.__strips)

    def Center(self):
        sum = 0
        for i in range(self.Size()):
            sum+=self.__strips[i].GetID()
        return sum/self.Size()
    
    # def GetTotalCharge(self):
    #     qtot = 0
    #     for i in range(self.Size()):
    #         qtot+= self.__strips[i].GetCharge()
    #     return qtot  

    def GetTotalCharge(self):
        return np.sum([strip.GetCharge() for strip in self.__strips])

    def GetEarliestStripTime(self):
        return np.min([strip.GetTime() for strip in self.__strips])
    
    def GetEarliestStrip(self):
        tmin = 900
        minstrip = None
        for strip in self.__strips:
            if strip.GetTime() < tmin:
                minstrip = strip
        return minstrip
    
    def GetDeltaT(self):
        t_min = 1000
        t_max = 0.1
        for i in range(self.Size()):
            t =  self.__strips[i].GetTime()
            if t < t_min:
                t_min = t
            if t > t_max:
                t_max = t
        return t_max - t_min + 0.00000000000001 # pecionata vera

    # def Centroid(self):
    #     qtot = 0
    #     centroid = 0
    #     pitch = self.ParentChamber().GetPitch()
    #     for i in range(self.Size()):
    #         qtot+= self.__strips[i].GetCharge()
    #         centroid += self.__strips[i].GetCharge() * (self.__strips[i].GetID()*pitch)
    #         return centroid/qtot  

    def Centroid(self):
        pitch = self.ParentChamber().GetPitch()
        qtot = np.sum([strip.GetCharge() for strip in self.__strips])
        centroid = np.sum([strip.GetCharge() * (strip.GetID()*pitch) for strip in self.__strips])
        return centroid/qtot if qtot != 0 else 0  # Avoid division by zero   
    
    def GetAllStrips(self):
        return self.__strips
    
    def ParentChamber(self):
        return self.chamber
    
    def Print(self):
        print(f"Cluster size = {len(self.__strips)}")
        for i in range(len(self.__strips)):
            self.__strips[i].Print()



def neighbor_distance(strip1, strip2):
    return abs(strip1.GetID()-strip2.GetID())


def CreateClusters(strips):

    if len(strips) == 0: 
        return {}

    # Calculate pairwise distances between strips using the custom distance metric
    distances = np.array([[neighbor_distance(strip1, strip2) for strip2 in strips] for strip1 in strips])

    # Use DBSCAN for clustering
    chamber = strips[0].ParentChamber()
    dbscan = DBSCAN(eps=chamber.GetMaxNConsecutiveHolesInCluster()+1, min_samples=1, metric='precomputed')
    labels = dbscan.fit_predict(distances)

    # Output clusters
    clusters = {}
    for strip, label in zip(strips, labels):
        if label not in clusters:
            clusters[label] = Cluster(strip.ParentChamber())
        clusters[label].AddStrip(strip)
    
    goodclusters = [clusters[i] for i in range(len(clusters)) if clusters[i].isGoodCluster()]

    return goodclusters

# def FindHighestCluster(clusters):
#     charge = 0
#     highest_cluster = 0
#     for i in range(len(clusters)):
#         if clusters[i].GetTotalCharge() > charge:
#             highest_cluster = clusters[i]
#     return highest_cluster

def FindHighestCluster(clusters):
    highest_cluster = max(clusters, key=lambda cluster: cluster.GetTotalCharge())
    return highest_cluster

@njit
def line(x, m, q):
    return m*x + q


'''The Fermi Dirac function can be used to fit the rise 
    of the APV charge signal and extract the time of the hit
    - max : is the value of the function at infinity
    - thalf: is the value for which the function value is half the maximum
    - slop: defines the steepness of the rise
    - shift: represents the charge baseline'''

@njit
def FermiDirac(x, max, thalf, slope, shift):
    return (max/(1+np.exp(-(x-thalf)/slope))) + shift

@njit
def fermi_dirac_jacobian(x, max, thalf, slope, shift):
    exp_term = np.exp((x - thalf) / slope)
    df_dmu = -max * exp_term / (slope * (1 + exp_term) ** 2)
    df_dkT = max * (x - thalf) * exp_term / (slope ** 2 * (1 + exp_term) ** 2)
    df_dA = 1 / (1 + exp_term)
    df_dB = np.ones_like(x)
    return np.column_stack((df_dmu, df_dkT, df_dA, df_dB))



def TimeFit(apv_times, apv_signal, q_error, maxindex, plot_apv_signal = False):

    ls = LeastSquares(apv_times[0:maxindex], apv_signal[0:maxindex], q_error[0:maxindex], FermiDirac)
    m = Minuit(ls, apv_signal[maxindex], maxindex, 50, apv_signal[0])
    m.limits = [(0, 3000), (0, 675), (1, 100), (-50, 50)]
    m.migrad()

    if m.fmin.is_valid:  # if fit is successful use fitted times
        t = m.values[1]
        err_t = m.errors[1]
        if plot_apv_signal==True:
            PlotAPVSignal(apv_times, apv_signal, t, err_t, m.values)

    else:                # otherwise compute time with maximum bin time
        t = 25*maxindex + 0.5*25
        err_t = 25/np.sqrt(12)   # assuming uniform probability in the bin
        if plot_apv_signal==True:
            PlotAPVSignal(apv_times, apv_signal, t, err_t, m.values)
    
    return t




def PlotAPVSignal(apv_times, apv_signal, t, err_t, fit_params = None):
    plt.plot(apv_times,apv_signal, marker = 's', alpha = 0.6)
    plt.xlabel(f"Time [{config['units']['time']}]")
    plt.ylabel(f"Charge [{config['units']['charge']}]")
    plt.vlines(t,0, 1000, color = 'darkgrey', lw = 2)
    plt.axvspan(t-err_t, t+err_t, alpha = 0.5)
    if fit_params != None:
        plt.plot(apv_times, FermiDirac(apv_times, *fit_params), color = 'red', lw = 2, alpha = 0.7)    
    plt.show()




def PlotTrackFromClusters(clusters, show_tracks = False, show_fit = False):

    scaling_factor = 0.3 #only for plotting purposes, scales the dimension of markers
    chamber = clusters[0].ParentChamber()

    pitch = chamber.GetPitch()
    drift_velocity = chamber.GetGas().GetDriftSpeed()

    pvalue = []
    thetas = []
    residues = []
        
    for i in range(len(clusters)):
        ids = []
        z = []
        charges = []

  
        ids = [strip.GetID() for strip in clusters[i].GetAllStrips()]
        z = [strip.GetTime()*drift_velocity for strip in clusters[i].GetAllStrips()]
        charges = [strip.GetCharge() for strip in clusters[i].GetAllStrips()]
        err_z = [25*drift_velocity for i in range(len(ids))] #the factor 2 is a test


        x = np.array(ids)*pitch
        z = np.array(z)
        err_z_arr = np.array(err_z)

        # cut = (z > 20) & (z < 50) 

        par, cov, fit_info, msg, ier = curve_fit(line, x, z, sigma = err_z_arr, absolute_sigma = True, full_output=True) 
        residues += list(fit_info['fvec']*err_z)
        thetas += [90 + (np.arctan(par[0]))*(180/np.pi)]
        chi_sq = np.sum(np.array(fit_info['fvec'])**2)
        pvalue.append(1-chi2.cdf(chi_sq, len(ids)-2))

        # ls = LeastSquares(x, z, err_z_arr, line)
        # m = Minuit(ls, 1, 1)
        # m.migrad()

        if show_tracks == True:
            plt.scatter(np.array(ids)*pitch, z, marker = 'o', edgecolor = 'royalblue', alpha = 0.5, color = 'lightskyblue', linewidth = 2,s = scaling_factor*np.array(charges), label = 'Charge Deposit')
            # plt.plot(np.array([max(ids), min(ids)])*pitch, line(np.array([max(ids), min(ids)])*pitch, *par), color = 'midnightblue', label = 'Best Fit')
            if show_fit == True:
                plt.plot(np.array([max(x), min(x)]), line(np.array([max(x), min(x)]), *par), color = 'midnightblue', label = 'Best Fit')
            plt.xlabel('x [mm]')
            plt.ylabel('z [mm]')
            # plt.legend()
            
    
    plt.show()

    return pvalue, thetas, residues


def FitTrackFromCluster(cluster, show_tracks = False):

    chamber = cluster.ParentChamber()
    pitch = chamber.GetPitch()
    drift_velocity = chamber.GetGas().GetDriftSpeed()

    strips = cluster.GetAllStrips()

    z = [strip.GetTime()*drift_velocity for strip in strips]
    x = [strip.GetID()*pitch for strip in strips]
    charges = [strip.GetCharge() for strip in strips]
    err_z = [25*drift_velocity for i in range(len(x))] #the factor 2 is a test

    # par, cov, fit_info, msg, ier = curve_fit(line, x, z, sigma = err_z, absolute_sigma = True, full_output=True) 
    par, cov, fit_info, msg, ier = curve_fit(line, x, z, absolute_sigma=True, sigma = err_z, full_output=True) 
    # residues = list(fit_info['fvec']*err_z)
    residues = [z[i]-line(x[i], *par) for i in range(len(x))]
    xhalf = (0.5*config['P2_M01']['gap'] - par[1])/par[0]
    chi_sq = np.sum(np.array(fit_info['fvec'])**2)
    pvalue = 1-chi2.cdf(chi_sq, len(x)-2)

    if show_tracks == True:
        PlotTrack(x, z, charges, par)
        plt.plot(xhalf, 0.5*config['P2_M01']['gap'])
        plt.vlines(xhalf, 0, 0.5*config['P2_M01']['gap'])
        plt.hlines(0.5*config['P2_M01']['gap'], xhalf-10, xhalf+10)

    # print(f"theta from pymegas : {90 + (np.arctan(par[0]))*(180/np.pi)}")

    return par, xhalf, pvalue, residues, z, charges



def PlotTrack(x, z, charges, pars = []):
    scaling_factor = 0.3 #only for plotting purposes, scales the dimension of markers
    plt.scatter(x, z, marker = 'o', edgecolor = 'royalblue', alpha = 0.5, color = 'lightskyblue', linewidth = 2,s = scaling_factor*np.array(charges), label = 'Charge Deposit')
    if len(pars) != 0:
        plt.plot(np.array([max(x), min(x)]), line(np.array([max(x), min(x)]), *pars), color = 'midnightblue', label = 'Best Fit')
    plt.xlabel('x [mm]')
    plt.ylabel('z [mm]')

def PlotCluster(cluster):
    chamber = cluster.ParentChamber()
    pitch = chamber.GetPitch()
    drift_velocity = chamber.GetGas().GetDriftSpeed()

    strips = cluster.GetAllStrips()

    x = [strip.GetID()*pitch for strip in strips]
    z = [strip.GetTime()*drift_velocity for strip in strips]
    charges = [strip.GetCharge() for strip in strips]
    err_z = [2*25*drift_velocity for i in range(len(x))] #the factor 2 is a test

    PlotTrack(x, z, charges)


# pars are the fit parameters for the whole track
# we use it as a seed for the fit to decrease the number of iterations
def HitEfficiency(cluster, par):
    chamber = cluster.ParentChamber()
    drift_velocity = chamber.GetGas().GetDriftSpeed()

    strips = cluster.GetAllStrips()

    ids = np.array([strip.GetID() for strip in strips])
    z = np.array([strip.GetTime()*drift_velocity for strip in strips])
    err_z = np.array([2*25*drift_velocity for i in range(len(ids))]) #the factor 2 is a test
    first_strip = min(ids)
    last_strip = max(ids) + 1

    cluster_strips_id =  []
    tot = []
    hit = []

    for i in range(first_strip, last_strip):
        ignore_strip = i
        mask = (ids != ignore_strip)
        par, cov = curve_fit(line, ids[mask], z[mask], absolute_sigma= True, sigma = err_z[mask], p0 = (par[0]/0.4, par[1]))
        fitx_end = -par[1]/par[0]
        fitx_start = (config['P2_M01']['gap']-par[1])/par[0]
        if i<fitx_end and i > fitx_start:
            cluster_strips_id.append(i)
            if i in ids:
                hit.append(1)
            else:
                hit.append(0)
    
    return cluster_strips_id, hit


def DoubleFermiDiracFit(times, show_fit = False):
    
    nbins = 500
    hist_range = (-100, 675)
    counts, bin_edges = np.histogram(times, bins=nbins, range=hist_range, density = False)

    bin_width = bin_edges[1]-bin_edges[0]
    bin_centers = bin_edges[0:-1] + 0.5*bin_width

    # this is needed otherwise we get nan for fit parameter errors.
    # Why? I have absolutely no idea.
    x = [bin_centers[i] for i in range(len(bin_centers)) if counts[i] != 0]
    counts = [counts[i] for i in range(len(counts)) if counts[i] != 0]
    err_y = [0.1*count for count in counts]

    # plt.scatter(x, counts)


    #for the rising part fit we take the first 50% of the data
    fraction = 0.5
    end_idx = np.argmax(x[0:int(fraction*len(x))])
    ls = LeastSquares(x[0:end_idx], counts[0:end_idx], err_y[0:end_idx], FermiDirac)
    m = Minuit(ls, max = counts[end_idx], thalf= x[end_idx],slope = 10, shift = counts[0])
    m.fixed["shift"] = True
    m.migrad()
    m.hesse()

    start_idx = np.argmax(x[int(len(x)*(1-fraction)):])
    ls1 = LeastSquares(x[start_idx:], counts[start_idx:], err_y[start_idx:], FermiDirac)
    m1 = Minuit(ls1, max = -m.values[0], thalf= 500,slope = 10, shift = 2*m.values[0])
    m1.fixed["max"] = True
    m1.migrad()
    m1.hesse()

    if show_fit == True:
        x_low = np.linspace(-100, 250, 100)
        x_high = np.linspace(250, 775, 200)
        plt.hist(times, bins=nbins, range=hist_range, density = False, histtype = 'step')
        # plt.scatter(x, counts)
        plt.plot(x_low, FermiDirac(x_low, *m.values), color = 'red', lw = 2)
        plt.plot(x_high, FermiDirac(x_high, *m1.values), color = 'red', lw = 2)
        plt.xlabel('Time [ns]')
        plt.ylabel('Entries (Normalized)')
        # plt.show()
    
    return m.values, m.errors, m1.values, m1.errors














