import yaml
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from sklearn.cluster import DBSCAN
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
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
        chamb.SetMinStripsInCluster(config[name]['MinStripsInCluster'])
        chamb.SetMinClusterCharge(config[name]['MinClusterCharge'])
        chamb.SetMaxStripsInCluster(config[name]['MaxStripsInCluster'])
        chamb.SetMaxNConsecutiveHolesInCluster(config[name]['MaxNConsecutiveHolesInCluster'])
        chamb.SetMaskedStrips(config[name]['MaskedStrips'])
        chambers.append(chamb)
        if verbose == True:
            chamb.Print()
    
    return chambers

class Strip:
    def __init__(self, id, charge, chamber = None, time = 0):
        self.id = id
        self.charge = charge
        self.chamber = chamber
        self.time = time
    
    def isGoodStrip(self):
        if self.charge >= self.chamber.GetMinStripCharge():
            # if self.id not in self.chamber.GetMaskedStrips():
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
    
    def GetDeltaT(self):
        t_min = 1000
        t_max = 0.1
        for i in range(self.Size()):
            t =  self.__strips[i].GetTime()
            if t < t_min:
                t_min = t
            if t > t_max:
                t_max = t
        return t_max - t_min + 0.00000000000001

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

# def Moyal(x, a, m, s, base):


def TimeFit(apv_times, apv_signal, maxindex, plot_apv_signal = False):
    # bin = config['apv']['bin_width']
    # nbins = config['apv']['nbins']
    # maxindex = np.argmax(apv_signal)
    # apv_times = np.arange(0+(0.5*bin), bin*nbins+(0.5*bin), bin)
    q_error = np.repeat(config['apv']['charge_error'], len(apv_signal))

    '''If the number of points to fit is lower than the 
    number of parameters we cannot do the fit, in this
    case we return the time of the maximum charge bin'''

    t = 0
    err_t = 0

    bounds = ((1, 1, 1, -10),(3000,1000,100,200))
    if maxindex > 3:
        try:
            par,cov = curve_fit(FermiDirac, apv_times[0:maxindex], apv_signal[0:maxindex], p0=(apv_signal[maxindex], maxindex, 50, 0), bounds=bounds, absolute_sigma=True, sigma = q_error[0:maxindex], jac = fermi_dirac_jacobian)
            t = par[1]
            err_t = np.sqrt(cov[1][1])
            if plot_apv_signal==True:
                plt.plot(apv_times,apv_signal, marker = 's', alpha = 0.6)
                plt.xlabel(f"Time [{config['units']['time']}]")
                plt.ylabel(f"Charge [{config['units']['charge']}]")
                plt.vlines(t,0, 1000, color = 'darkgrey', lw = 2)
                plt.axvspan(t-err_t, t+err_t, alpha = 0.5)
                plt.plot(apv_times, FermiDirac(apv_times, *par), color = 'red', lw = 2, alpha = 0.7)    
                plt.show()   
            return t
        except (RuntimeError, RuntimeWarning, OptimizeWarning) as e:
            t = 25*maxindex + 0.5*25
            err_t = 25/2
            if plot_apv_signal==True:
                plt.plot(apv_times,apv_signal, marker = 's')
                plt.xlabel(f"Time [{config['units']['time']}]")
                plt.ylabel(f"Charge [{config['units']['charge']}]")
                plt.vlines(t,0, 1000)
                plt.axvspan(t-err_t, t+err_t, alpha = 0.5)
                plt.show()       
            return t
    else:
        t = 25*maxindex + 0.5*25
        err_t = 25/2
        if plot_apv_signal==True:
            plt.plot(apv_times,apv_signal, marker = 's')
            plt.xlabel(f"Time [{config['units']['time']}]")
            plt.ylabel(f"Charge [{config['units']['charge']}]")
            plt.vlines(t,0, 1000)
            plt.axvspan(t-err_t, t+err_t, alpha = 0.5)
            plt.show()       
        return t






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
        charges = [strip.GetCharge()*scaling_factor for strip in clusters[i].GetAllStrips()]
        err_z = [2*25*drift_velocity for i in range(len(ids))] #the factor 2 is a test


        x = np.array(ids)*pitch
        z = np.array(z)
        err_z_arr = np.array(err_z)

        # cut = (z > 20) & (z < 50) 


        par, cov, fit_info, msg, ier = curve_fit(line, x, z, sigma = err_z_arr, absolute_sigma = True, full_output=True) 
        residues += list(fit_info['fvec']*err_z)
        thetas += [90 + (np.arctan(par[0]))*(180/np.pi)]
        chi_sq = np.sum(np.array(fit_info['fvec'])**2)
        pvalue.append(1-chi2.cdf(chi_sq, len(ids)-2))

        if show_tracks == True:
            plt.scatter(np.array(ids)*pitch, z, marker = 'o', edgecolor = 'royalblue', alpha = 0.5, color = 'lightskyblue', linewidth = 2,s = charges, label = 'Charge Deposit')
            # plt.plot(np.array([max(ids), min(ids)])*pitch, line(np.array([max(ids), min(ids)])*pitch, *par), color = 'midnightblue', label = 'Best Fit')
            if show_fit == True:
                plt.plot(np.array([max(x), min(x)]), line(np.array([max(x), min(x)]), *par), color = 'midnightblue', label = 'Best Fit')
            plt.xlabel('x [mm]')
            plt.ylabel('z [mm]')
            # plt.legend()

        # if len(z[cut]) >= 2 :
        #     par, cov, fit_info, msg, ier = curve_fit(line, x[cut], z[cut], sigma = err_z_arr[cut], absolute_sigma = True, full_output=True) 
        #     residues += list(fit_info['fvec']*err_z_arr[cut])
        #     thetas += [90 + (np.arctan(par[0]))*(180/np.pi)]
        #     chi_sq = np.sum(np.array(fit_info['fvec'])**2)
        #     pvalue.append(1-chi2.cdf(chi_sq, len(ids)-2))

        #     if show_tracks == True:
        #         plt.scatter(np.array(ids)*pitch, z, marker = 'o', edgecolor = 'royalblue', alpha = 0.5, color = 'lightskyblue', linewidth = 2,s = charges, label = 'Charge Deposit')
        #         # plt.plot(np.array([max(ids), min(ids)])*pitch, line(np.array([max(ids), min(ids)])*pitch, *par), color = 'midnightblue', label = 'Best Fit')
        #         plt.plot(np.array([max(x[cut]), min(x[cut])]), line(np.array([max(x[cut]), min(x[cut])]), *par), color = 'midnightblue', label = 'Best Fit')
        #         plt.xlabel('x [mm]')
        #         plt.ylabel('z [mm]')
        #         plt.legend()
            
    
    plt.show()

    return pvalue, thetas, residues 







