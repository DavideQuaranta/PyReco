import os
import sys
import yaml
import pickle
import numpy as np
import PyMegas as mm
import mplhep as hep
import uproot as root
import matplotlib.pyplot as plt

# number of chambers in the setup
NCHAMBERS = 2

# Read the configuration file
with open('config.yaml', 'r') as file:
    config = yaml.safe_load(file)


# Open the ROOT file of the run
run_number = sys.argv[1]
filename = config['paths']['raw_dir'] + "run" + run_number + "_xtalk1_xtalk2.root"
raw = root.open(filename)
print("\n------------- Test Beam Reconstruction -------------")
print(f"- opening run {run_number} from {filename}...")


# Check if object exist in ROOT file
print('- Available objects in your TFile:\n')
print(raw.keys())
treename = input("Insert name of object to read : ")
while treename not in raw.keys():
    print('ERROR : the object you are searching for does not exist, please try again!')
    treename = input("Insert name of object to read : ")
print(f"- creating tree from {treename} and reading branches...")
raw_tree = raw[treename]


# Reading branches of ROOT TTree
branch = raw_tree.keys()
evt = raw_tree[branch[0]].array()
error = raw_tree[branch[1]].array()
daqTimeSec = raw_tree[branch[2]].array()
daqTimeMicroSec = raw_tree[branch[3]].array()
srsTimeStamp = raw_tree[branch[4]].array()
srsTrigger = raw_tree[branch[5]].array()
srsFec = raw_tree[branch[6]].array()
srsChip = raw_tree[branch[7]].array()
srsChan = raw_tree[branch[8]].array()
mmChamber = raw_tree[branch[9]].array()
mmLayer = raw_tree[branch[10]].array()
mmReadout = raw_tree[branch[11]].array()
mmStrip = raw_tree[branch[12]].array()
raw_q = raw_tree[branch[13]].array()

# Create chambers and initialize their parameters
TMM, ExMe = mm.CreateChambers(config, NCHAMBERS)

# Here define the histograms you want to fill
# and initialize a dictionary that will contain them
histo_dict = {}
chamber_keys = [TMM.GetName(), ExMe.GetName()]
histo_keys =  [None] * NCHAMBERS
histo_keys[0] = ['Strips', 'Charge', 'Time', 'Centroid', 'Centroid Residues', 'Cluster Size', 'NClusters', 'Cluster Charge']
histo_keys[1] = ['Strips', 'Charge', 'Time', 'Centroid', 'Centroid Residues', 'Cluster Size', 'NClusters', 'Cluster Charge', 'Theta', 'Z residues', 'Drift Speed']

for i, chamb_name in enumerate(chamber_keys):
    histo_dict[chamb_name] = {}
    for key in histo_keys[i]:
        histo_dict[chamb_name][key] = []



# These variables are needed to keep count of
# how many times we have good events in TMM and ExMe
# to compute efficiency 
TMM_hit = 0
ExMe_hit = 0

nentries = len(evt)

# for loop on events
for i in range(nentries):

    if i%100 == 0:
        print(f'Reading event {i}/{nentries}')
        print(f'Progress {100*(i/nentries):.2f} %')

    strips = mmStrip[i]
    charges = [max(charges) for charges in raw_q[i]]
    chamber = mmChamber[i]

    TMM_strips = []
    ExMe_strips = []

    t_min_ExMe = 200
    t_max_ExMe = 500


    for j in range(len(strips)):
        strip = mm.Strip(strips[j], charges[j])
        if chamber[j] == "P0_RM1":
            strip.SetParentChamber(TMM)
            if mmReadout[i][j] == 89: # only read y coordinate      
                if strip.isGoodStrip():
                    strip.SetTime(25*np.argmax(raw_q[i][j])+0.5*25)
                    TMM_strips.append(strip)
        elif chamber[j] == "P2_M01":
            strip.SetParentChamber(ExMe)
            if mmReadout[i][j] == 89: # only read y coordinate      
                if strip.isGoodStrip():
                    if config['general']['do_time_fit'] == True:
                        strip.SetTime(mm.TimeFit(raw_q[i][j], config['general']['plot_apv_signals']))
                    else:
                        strip.SetTime(25*np.argmax(raw_q[i][j])+0.5*25)
                    
                    ExMe_strips.append(strip)
                    
                    if strip.GetTime() > t_max_ExMe:
                        t_max_ExMe = strip.GetTime()
                    if strip.GetTime() < t_min_ExMe:
                        t_min_ExMe = strip.GetTime()



    histo_dict['P0_RM1']['Strips'] += [TMM_strips[k].GetID() for k in range(len(TMM_strips))]
    histo_dict['P2_M01']['Strips'] += [ExMe_strips[k].GetID() for k in range(len(ExMe_strips))]

    histo_dict['P0_RM1']['Charge'] += [TMM_strips[k].GetCharge() for k in range(len(TMM_strips))]
    histo_dict['P2_M01']['Charge'] += [ExMe_strips[k].GetCharge() for k in range(len(ExMe_strips))]

    histo_dict['P0_RM1']['Time'] += [TMM_strips[k].GetTime() for k in range(len(TMM_strips))]
    histo_dict['P2_M01']['Time'] += [ExMe_strips[k].GetTime() for k in range(len(ExMe_strips))]

    histo_dict['P2_M01']['Drift Speed'].append(config['P2_M01']['gap']/(t_max_ExMe-t_min_ExMe))

    TMMClust = mm.CreateClusters(TMM_strips)
    ExMeClust = mm.CreateClusters(ExMe_strips)

    # print(f'Event {i}')

    TMM_fitinfo = 0
    ExMe_fitinfo = 0
    residues = []
    theta = []

    if len(TMMClust) == 1 and len(ExMeClust) >= 1: 
        TMM_fitinfo, _ , _ = mm.PlotTrackFromClusters(TMMClust, show_tracks = config['general']['show_tracks'])
        ExMe_fitinfo, angle, res = mm.PlotTrackFromClusters(ExMeClust, show_tracks = config['general']['show_tracks'])
        residues+=res
        theta+=angle
        # plt.show()
        for l in range(len(TMMClust)):
            if abs(TMMClust[l].Centroid() - config['P0_RM1']['offset']) < 3*6.60: #6.60 is the std of the tmm centroid computed on all events
                TMM_hit += 1                
                AlreadyFoundGoodClusterExMe = 0
                for m in range(len(ExMeClust)):
                    # print(f'fit info : {ExMe_fitinfo[m]}')
                    # print(f'distance : {abs(TMMClust[l].Centroid() - config['P0_RM1']['offset'] + config['P2_M01']['offset'] - ExMeClust[m].Centroid())}')
                    if ExMe_fitinfo[m] > 0.05:
                        if abs(TMMClust[l].Centroid() - config['P0_RM1']['offset'] + config['P2_M01']['offset'] - ExMeClust[m].Centroid()) < config['P2_M01']['EFF_WINDOW']:
                            if AlreadyFoundGoodClusterExMe == 0:
                                ExMe_hit += 1
                                AlreadyFoundGoodClusterExMe = 1
                                histo_dict['P0_RM1']['Centroid Residues'] += [TMMClust[l].Centroid() - config['P0_RM1']['offset'] + config['P2_M01']['offset'] - ExMeClust[m].Centroid()]
                                histo_dict['P2_M01']['Centroid Residues'] += [TMMClust[l].Centroid() - config['P0_RM1']['offset'] + config['P2_M01']['offset'] - ExMeClust[m].Centroid()]


    # print(f'TMM_hits : {TMM_hit} ExMe_hits : {ExMe_hit}')


    # histo_dict['P0_RM1']['Clusters'].append(TMMClust)
    # histo_dict['P2_M01']['Clusters'].append(ExMeClust)

    histo_dict['P0_RM1']['NClusters'] += [len(TMMClust)]
    histo_dict['P2_M01']['NClusters'] += [len(ExMeClust)]

    histo_dict['P0_RM1']['Cluster Charge'] += [TMMClust[k].GetTotalCharge() for k in range(len(TMMClust))]
    histo_dict['P2_M01']['Cluster Charge'] += [ExMeClust[k].GetTotalCharge() for k in range(len(ExMeClust))]

    histo_dict['P0_RM1']['Cluster Size'] += [TMMClust[k].Size() for k in range(len(TMMClust))]
    histo_dict['P2_M01']['Cluster Size'] += [ExMeClust[k].Size() for k in range(len(ExMeClust))]

    TMMCentroid = [TMMClust[k].Centroid() for k in range(len(TMMClust))]
    ExMeCentroid = [ExMeClust[k].Centroid() for k in range(len(ExMeClust))]
    histo_dict['P0_RM1']['Centroid'] += TMMCentroid
    histo_dict['P2_M01']['Centroid'] += ExMeCentroid

    histo_dict['P2_M01']['Z residues'] += residues
    histo_dict['P2_M01']['Theta'] += theta



print(f'Clusters created for run {run_number}...')

efficiency = ExMe_hit/TMM_hit
efficiency_error = np.sqrt(efficiency*(1-efficiency))/np.sqrt(TMM_hit)

print(f'Efficiency = {efficiency} +/- {efficiency_error}')



plots_dir = config['paths']['plots_dir'] + run_number + '/'

if not os.path.exists(plots_dir):
    # Create the directory
    os.makedirs(plots_dir)
    print(f"Saving plots in {plots_dir}...")

histo_dir = config['paths']['histo_dir'] + run_number + '/'

if not os.path.exists(histo_dir):
    # Create the directory
    os.makedirs(histo_dir)
    print(f"Saving histograms in {histo_dir}...")

# # plt.style.use(hep.style.ROOT) # to plot with a nicer style
# plt.style.use('default') # to plot with matplotlib default style


xlabels = ['Strip ID', 'Charge [ADC Counts]', 'Time [ns]', 'Centroid [mm]', 'Centroid Residues [mm]', 'Cluster Size', 'NClusters', 'Cluster Charge [ADC Counts]']

for i,key in enumerate(histo_keys[0]):
    plt.hist(histo_dict[chamber_keys[0]][key], bins = 100, histtype= 'step', label = chamber_keys[0])
    plt.hist(histo_dict[chamber_keys[1]][key], bins = 100, histtype= 'step', label = chamber_keys[1])
    plt.xlabel(xlabels[i])
    plt.ylabel("Entries")
    plt.legend()
    with open(histo_dir+key+chamber_keys[0]+'.pkl', 'wb') as file:
        pickle.dump(histo_dict[chamber_keys[0]][key], file)
    with open(histo_dir+key+chamber_keys[1]+'.pkl', 'wb') as file:
        pickle.dump(histo_dict[chamber_keys[1]][key], file)
    plt.savefig(plots_dir+key+'.pdf')
    plt.show()

xlabels = ['Theta [deg]', '[mm]', 'Drift Speed [mm/ns]']
for i, key in enumerate(histo_keys[1][-3:-1]):
    plt.hist(histo_dict[chamber_keys[1]][key], bins = 30, histtype= 'step', label = chamber_keys[1])
    plt.xlabel(xlabels[i])
    plt.ylabel("Entries")
    plt.legend()
    with open(histo_dir+key+chamber_keys[1]+'.pkl', 'wb') as file:
        pickle.dump(histo_dict[chamber_keys[1]][key], file)
    plt.savefig(plots_dir+key+'.pdf')
    plt.show()






TMM_mean = np.array(histo_dict[chamber_keys[0]]['Centroid']).mean()
TMM_std = np.array(histo_dict[chamber_keys[0]]['Centroid']).std()

print(f'TMM centroid : {TMM_mean:.2f} +/- {TMM_std:.2f} mm')
print('Done!')

exit(1)


