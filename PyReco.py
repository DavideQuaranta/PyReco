import os                         # to launch system commands
import re                         # to search for patterns in string
import sys                        # to deal with command-line argument
import yaml                       # to read YAML files
import shutil
import pickle                    
import numpy as np
import uproot as root
from tqdm import tqdm
import matplotlib.pyplot as plt
from PyMegas import *

print("\n------------- Test Beam Reconstruction -------------")

NCHAMBERS = 2   # number of chambers in the setup

CONFIG_FILENAME = 'config.yaml'     # name of the configuration file to use
print(f"- reading configuration file {CONFIG_FILENAME} ...")
# Read the configuration file
with open(CONFIG_FILENAME, 'r') as file:
    config = yaml.safe_load(file)
    print(f"- configuration file {CONFIG_FILENAME} opened with success ...")


# Open the ROOT file of the run
run_number = sys.argv[1]
filename = config['paths']['raw_dir'] + "run" + run_number + "_xtalk1_xtalk2.root"
raw = root.open(filename)
print(f"- opening run {run_number} from {filename}...")


# Check if object exist in ROOT file
print('- available objects in your TFile:')
print(f'- {raw.keys()}')
# treename = input("Insert name of object to read : ")
# Get the filename with the highest number
treename = get_highest_apv_raw(raw.keys())  
print(f"- creating tree from {treename} and reading branches...\n")
raw_tree = raw[treename]


# Reading branches of ROOT TTree
branch = raw_tree.keys()
evt = raw_tree[branch[0]].array()
mmChamber = raw_tree[branch[9]].array()
mmReadout = raw_tree[branch[11]].array()
mmStrip = raw_tree[branch[12]].array()
raw_q = raw_tree[branch[13]].array()

# Reading ROOT TTree containing strip fitted times
timefit_filename = config['paths']['timefit_dir']+ 'run' + run_number + '_timefit.root'
timefit_file = root.open(timefit_filename)
print(f"- opening timefit file from {filename} and reading branches...\n")
time_tree = timefit_file['output_tree']
mmStrip_time = time_tree['fitted_time'].array()
mmStrip_time_err = time_tree['fitted_time_err'].array()

# create directories for output
plots_dir = config['paths']['plots_dir'] + run_number + '/'

if not os.path.exists(plots_dir):
    # Create the directory
    os.makedirs(plots_dir)
    print(f"Plots will be saved in {plots_dir}...")

histo_dir = config['paths']['histo_dir'] + run_number + '/'

if not os.path.exists(histo_dir):
    # Create the directory
    os.makedirs(histo_dir)
    print(f"Histograms will be saved in {histo_dir}...")

results_dir = config['paths']['results_dir'] + run_number + '/'

if not os.path.exists(results_dir):
    # Create the directory
    os.makedirs(results_dir)
    print(f"Results will be saved in {results_dir}...")



# Create chambers and initialize their parameters
TMM, ExMe = CreateChambers(config, NCHAMBERS)

# Here define the histograms you want to fill
# and initialize a dictionary that will contain them
histo_dict = {}
chamber_keys = [TMM.GetName(), ExMe.GetName()]
histo_keys =  [None] * NCHAMBERS
histo_keys[0] = ['Strips', 'Charge', 'Time', 'Centroid', 'Centroid Residues',
                 'Cluster Size', 'NClusters', 'Cluster Charge']
histo_keys[0] += [name + '_1st_coord' for name in histo_keys[0]]
histo_keys[1] = ['Strips', 'Charge', 'Time', 'Centroid', 'Centroid Residues', 
                 'Cluster Size', 'NClusters', 'Cluster Charge', 'Theta',
                 'Z residues', 'Drift Speed', 'Time Residues', 
                 'GoodClustCharges', 'GoodClustZ', 'xhalf', 'xhalf_res']


for i, chamb_name in enumerate(chamber_keys):
    histo_dict[chamb_name] = {}
    for key in histo_keys[i]:
        histo_dict[chamb_name][key] = []

# we need to compute the alignement to calculate the efficiency
# but we should run the program twice...
# since the damn thing has to run i save the clusters in lists
# for the chamber so i can retrieve them after the alignement
# i have to use lists for now because numpy array are of fixed size once created
TMM_allclust = []
TMM_allclust_1st_coord = []
ExMe_allclust = []


# These variables are needed to keep count of
# how many times we have good events in TMM and ExMe
# to compute efficiency
TMM_hit = 0
ExMe_hit = 0

# number of entries to process
nentries = int(config['general']['NMAXEVENTS']*len(evt))

# these we compute once and for all
bin = config['apv']['bin_width']
nbins = config['apv']['nbins']
apv_times = np.arange(0+(0.5*bin), bin*nbins+(0.5*bin), bin)

# array for quantities needed to compute strip efficiency
strip_efficiency = {}
for i in range(ExMe.GetNStrips()):
    strip_efficiency[i] = [0,0]



for i in tqdm(range(nentries)):

    if config['general']['verbose'] == True:
        if i%100 == 0:
            print(f'Reading event {i}')
            print(f'Progress {100*(i/nentries):.2f} %')

    charges = np.max(raw_q[i], axis = 1)           # for each strip take the highest bin as the charge
    max_indices = np.argmax(raw_q[i], axis = 1)    # the bin of the maximum

    # maybe creating these outside the for and only emptying them before each event is faster
    TMM_strips = np.zeros(config['P0_RM1']['nstrips'], dtype=Strip)              
    TMM_strips_1st_coord = np.zeros(config['P0_RM1']['nstrips'], dtype=Strip)
    ExMe_strips = np.zeros(config['P2_M01']['nstrips'], dtype=Strip)
    
    # keep track of how many good strips
    current_size_TMM = 0
    current_size_TMM_1st = 0
    current_size_ExMe = 0
    
    # quantities of current event
    event_strips = np.array(mmStrip[i])
    event_strip_times = np.array(mmStrip_time[i])
    event_strip_times_error = np.array(mmStrip_time_err[i])
    event_chamber = np.array(mmChamber[i])
    event_readout = np.array(mmReadout[i])

    for j in range(len(event_strips)):
        strip_id = event_strips[j]
        strip_time = event_strip_times[j]
        strip_time_error = event_strip_times_error[j]
        charge = charges[j]
        max_index = max_indices[j]
        chamber_name = event_chamber[j]
        readout = event_readout[j]

        strip = Strip(strip_id, charge, strip_time, strip_time_error)

        if chamber_name == "P0_RM1":
            strip.SetParentChamber(TMM)
            if readout == 89: # only read y coordinate      
                if strip.isGoodStrip():
                    # TMM_strips.append(strip)
                    TMM_strips[current_size_TMM] = strip
                    current_size_TMM += 1
                    histo_dict['P0_RM1']['Strips'] += [strip_id]
                    histo_dict['P0_RM1']['Charge'] += [charge]
                    histo_dict['P0_RM1']['Time'] += [strip_time]
            if readout == 88: # read x coordinate      
                if strip.isGoodStrip():
                    # TMM_strips_1st_coord.append(strip)
                    TMM_strips_1st_coord[current_size_TMM_1st] = strip
                    current_size_TMM_1st += 1
                    histo_dict['P0_RM1']['Strips_1st_coord'] += [strip_id]
                    histo_dict['P0_RM1']['Charge_1st_coord'] += [charge]
                    histo_dict['P0_RM1']['Time_1st_coord'] += [strip_time]
        elif chamber_name == "P2_M01": # change to P1_M01 for november runs!!!!
            strip.SetParentChamber(ExMe)
            if readout == 89: # only read y coordinate  change to 88 for november runs!!!!
                if strip.isGoodStrip():
                    # ExMe_strips.append(strip)
                    ExMe_strips[current_size_ExMe] = strip
                    current_size_ExMe += 1
                    histo_dict['P2_M01']['Strips'] += [strip_id]
                    histo_dict['P2_M01']['Charge'] += [charge]
                    histo_dict['P2_M01']['Time'] += [strip_time]

    # TMM_strips = np.array(TMM_strips, dtype=Strip)
    # TMM_strips_1st_coord = np.array(TMM_strips_1st_coord, dtype=Strip)
    # ExMe_strips = np.array(ExMe_strips, dtype=Strip)


    # histo_dict['P0_RM1']['Strips'] += [TMM_strips[k].id for k in range(len(TMM_strips))]
    # histo_dict['P0_RM1']['Strips_1st_coord'] += [TMM_strips_1st_coord[k].GetID() for k in range(len(TMM_strips_1st_coord))]
    # histo_dict['P2_M01']['Strips'] += [ExMe_strips[k].id for k in range(len(ExMe_strips))]

    # histo_dict['P0_RM1']['Charge'] += [TMM_strips[k].charge for k in range(len(TMM_strips))]
    # histo_dict['P0_RM1']['Charge_1st_coord'] += [TMM_strips_1st_coord[k].GetCharge() for k in range(len(TMM_strips_1st_coord))]
    # histo_dict['P2_M01']['Charge'] += [ExMe_strips[k].charge for k in range(len(ExMe_strips))]

    # histo_dict['P0_RM1']['Time'] += [TMM_strips[k].time for k in range(len(TMM_strips))]
    # histo_dict['P0_RM1']['Time_1st_coord'] += [TMM_strips_1st_coord[k].GetTime() for k in range(len(TMM_strips_1st_coord))]
    # histo_dict['P2_M01']['Time'] += [ExMe_strips[k].time for k in range(len(ExMe_strips))]

    TMM_allclust.append(CreateClusters(TMM_strips[0:current_size_TMM]))
    TMM_allclust_1st_coord.append(CreateClusters(TMM_strips_1st_coord[0:current_size_TMM_1st]))
    ExMe_allclust.append(CreateClusters(ExMe_strips[0:current_size_ExMe]))

#end of for on all events

# now we should be able to compute the average position and
# t0 correction for each chamber
# actually we also compute another time useful for the computation of drift speed
# now this part actually seems to work

print('- computing time alignement for chambers...')

ExMe_timefit = DoubleFermiDiracFit(histo_dict['P2_M01']['Time'], show_fit=True)
plt.savefig(plots_dir + 'TimeFit.pdf')
plt.close()
# plt.show()


# TMM_1st_coord_timefit = DoubleFermiDiracFit(histo_dict['P0_RM1']['Time_1st_coord'], show_fit=True)
# TMM_timefit = DoubleFermiDiracFit(histo_dict['P0_RM1']['Time'], show_fit=True)

T0_ExMe = ExMe_timefit[0][1]
T0_ExMe_err = ExMe_timefit[1][1]

T1_ExMe = ExMe_timefit[2][1]
T1_ExMe_err = ExMe_timefit[3][1]

# T0_TMM = TMM_timefit[0][1]
# T0_TMM_err = TMM_timefit[1][1]

# T0_TMM_1st_coord = TMM_1st_coord_timefit[0][1]
# T0_TMM_1st_coord_err = TMM_1st_coord_timefit[1][1]

print(f'- ExMe T0 = {T0_ExMe} +/- {T0_ExMe_err} [ns]')
# print(f'- TMM T0 = {T0_TMM} +/- {T0_TMM_err} [ns]')
# print(f'- TMM_1st_coord T0 = {T0_TMM_1st_coord} +/- {T0_TMM_1st_coord_err} [ns]')

# we can now compute the drift velocity
delta_t = T1_ExMe-T0_ExMe
delta_t_err = np.sqrt(T0_ExMe_err**2 + T1_ExMe_err**2)
drift_speed = config['P2_M01']['gap']/delta_t
drift_speed_err = drift_speed*(delta_t_err/delta_t)

print(f'- drift speed (ExMe) : {drift_speed} +/- {drift_speed_err} [mm/ns]')

speed_filename = results_dir+'drift_speed.txt'
with open(speed_filename, 'a') as speed_file:
    print(f'- writing drift speed to {speed_filename}...')
    speed_file.write(f'{run_number} {drift_speed} {drift_speed_err}\n')

# now we compute the spatial alignement for the chambers
# we do this by finding the mean of the number of on-strips
# this part seems to work 
print(f'- computing spatial alignement for the chambers...')
ExMe_offset = np.mean(histo_dict['P2_M01']['Strips']) * config['P2_M01']['pitch'] # 216.82976111166826 
ExMe_offset_std = np.std(histo_dict['P2_M01']['Strips']) * config['P2_M01']['pitch']
TMM_offset = np.mean(histo_dict['P0_RM1']['Strips']) * config['P0_RM1']['pitch'] # 37.07723433715333 
TMM_offset_std = np.std(histo_dict['P0_RM1']['Strips']) * config['P0_RM1']['pitch']
TMM_1st_coord_offset = np.mean(histo_dict['P0_RM1']['Strips_1st_coord']) * config['P0_RM1']['pitch']
TMM_1st_coord_offset_std = np.std(histo_dict['P0_RM1']['Strips_1st_coord']) * config['P0_RM1']['pitch']

print(f'- ExMe offset = {ExMe_offset} +/- {ExMe_offset_std} [mm]')
print(f'- TMM offset = {TMM_offset} +/- {TMM_offset_std} [mm]')
print(f'- TMM_1st_coord offset = {TMM_1st_coord_offset} +/- {TMM_1st_coord_offset_std} [mm]')

align_filename = results_dir+'alignement.txt'
with open(align_filename, 'a') as align_file:
    print(f'- writing alignement to {align_filename}...')
    align_file.write(f'{run_number} ExMe {ExMe_offset} {ExMe_offset_std}\n')
    align_file.write(f'{run_number} TMM {TMM_offset} {TMM_offset_std}\n')
    align_file.write(f'{run_number} TMM_1st_coord {TMM_1st_coord_offset} {TMM_1st_coord_offset_std}\n')



# now we loop on events again for track fitting
print("- looping on events for track fitting...")
ExMe_n_clusters = len(ExMe_allclust)
for i, ExMeClust in enumerate(ExMe_allclust):
    
    if i%100 == 0 : print(f'   - Event {i} / {ExMe_n_clusters}')
    
    # check if current event has clusters
    if len(ExMeClust) > 0:
        # if it has fit all clusters found
        par, xhalf, pvalue, residues, z, charges = zip(*[FitTrackFromCluster(clust, show_tracks=config['general']['show_tracks']) for clust in ExMeClust])

    # show the fit of the tracks in real time for debug
    if config['general']['show_tracks'] == True:
        plt.show()
        for cluster in TMM_allclust[i]:
            PlotCluster(cluster)
        plt.show()


    # Check if we have at least one cluster in reference chamber
    if len(TMM_allclust[i]) == 1: 
        # Check if the reference chamber clusters are in the region where there was beam
        for l in range(len(TMM_allclust[i])):
            if abs(TMM_allclust[i][l].Centroid() - TMM_offset) < 3*TMM_offset_std: #3*std of the tmm centroid computed on all events
                TMM_hit += 1                
                AlreadyFoundGoodClusterExMe = 0     # to check if i already found a compatible cluster in exme
                    
                if len(ExMeClust) == 0:
                    continue
                    
                # If there are clusters in the tested chamber, fit them
                if len(ExMeClust) >= 1:
                    for m in range(len(ExMeClust)):
                        theta = 90 + (np.arctan(par[m][0]))*(180/np.pi)
                            
                        if config['general']['verbose'] == True:
                            print(f'fit info : {pvalue}')
                            print(f"distance : {abs(TMM_allclust[i][l].Centroid() - TMM_offset + ExMe_offset - ExMeClust[m].Centroid())}")
                            

                        if pvalue[m] > config['P2_M01']['min_pvalue_for_good_track']: 
                            # check if cluster centroid is near a TMM cluster
                            if abs(TMM_allclust[i][l].Centroid() - TMM_offset + ExMe_offset - ExMeClust[m].Centroid()) < config['P2_M01']['EFF_WINDOW']:
                                if AlreadyFoundGoodClusterExMe == 0:
                                    ExMe_hit += 1
                                    AlreadyFoundGoodClusterExMe = 1
                                    strips_in_cluster, hit_strips = HitEfficiency(ExMeClust[m], par[m])
                                    for k in range(len(strips_in_cluster)):
                                        strip_efficiency[strips_in_cluster[k]][0] += 1
                                        strip_efficiency[strips_in_cluster[k]][1] += hit_strips[k] 
                                    histo_dict['P0_RM1']['Centroid Residues'] += [TMM_allclust[i][l].Centroid() - TMM_offset + ExMe_offset - ExMeClust[m].Centroid()]
                                    histo_dict['P2_M01']['Centroid Residues'] += [TMM_allclust[i][l].Centroid() - TMM_offset + ExMe_offset - ExMeClust[m].Centroid()]
                                    histo_dict['P2_M01']['Time Residues'] += [TMM_allclust[i][l].GetEarliestStripTime() - ExMeClust[m].GetEarliestStripTime()]
                                    histo_dict['P2_M01']['Theta'] += [theta]
                                    histo_dict['P2_M01']['xhalf'] += [xhalf[m]]
                                    histo_dict['P2_M01']['xhalf_res'] += [xhalf[m] - TMM_allclust[i][l].Centroid()]
                                    histo_dict['P2_M01']['GoodClustCharges'] += charges[m]
                                    histo_dict['P2_M01']['Z residues'] += residues[m]
                                    histo_dict['P2_M01']['GoodClustZ'] += z[m]
                    if config['general']['show_tracks'] == True:
                        plt.show()

    if config['general']['verbose']:
        print(f'TMM_Clust = {len(TMM_allclust[i])} | ExMe_Clust = {len(ExMeClust)}')
        print(f'TMM_hits : {TMM_hit} ExMe_hits : {ExMe_hit}')




    # histo_dict['P0_RM1']['Clusters'].append(TMMClust)
    # histo_dict['P2_M01']['Clusters'].append(ExMeClust)

    histo_dict['P0_RM1']['NClusters'] += [len(TMM_allclust[i])]
    histo_dict['P0_RM1']['NClusters_1st_coord'] += [len(TMM_allclust_1st_coord[i])]
    histo_dict['P2_M01']['NClusters'] += [len(ExMeClust)]

    histo_dict['P0_RM1']['Cluster Charge'] += [TMM_allclust[i][k].GetTotalCharge() for k in range(len(TMM_allclust[i]))]
    histo_dict['P0_RM1']['Cluster Charge_1st_coord'] += [TMM_allclust_1st_coord[i][k].GetTotalCharge() for k in range(len(TMM_allclust_1st_coord[i]))]
    histo_dict['P2_M01']['Cluster Charge'] += [ExMeClust[k].GetTotalCharge() for k in range(len(ExMeClust))]

    histo_dict['P0_RM1']['Cluster Size'] += [TMM_allclust[i][k].Size() for k in range(len(TMM_allclust[i]))]
    histo_dict['P0_RM1']['Cluster Size_1st_coord'] += [TMM_allclust_1st_coord[i][k].GetTotalCharge() for k in range(len(TMM_allclust_1st_coord[i]))]
    histo_dict['P2_M01']['Cluster Size'] += [ExMeClust[k].Size() for k in range(len(ExMeClust))]

    TMMCentroid = [TMM_allclust[i][k].Centroid() for k in range(len(TMM_allclust[i]))]
    TMMCentroid_1st_coord = [TMM_allclust_1st_coord[i][k].Centroid() for k in range(len(TMM_allclust_1st_coord[i]))]
    ExMeCentroid = [ExMeClust[k].Centroid() for k in range(len(ExMeClust))]
    histo_dict['P0_RM1']['Centroid'] += TMMCentroid
    histo_dict['P0_RM1']['Centroid_1st_coord'] += TMMCentroid_1st_coord
    histo_dict['P2_M01']['Centroid'] += ExMeCentroid




print(f'\nClusters created for run {run_number}...')

if TMM_hit == 0:
    print('TMM_hit = 0 : Could not compute efficiency!')
elif TMM_hit >= 1:
    efficiency = ExMe_hit/TMM_hit
    efficiency_error = np.sqrt(efficiency*(1-efficiency))/np.sqrt(TMM_hit)
    print(f'Efficiency = {efficiency:.3f} +/- {efficiency_error:.3f}')
    with open(results_dir+'efficiency.txt', 'a') as file:
        print(f'- writing efficiency to efficiency.txt')
        file.write(f'{run_number} {efficiency} {efficiency_error}\n')
    




strip_efficiencies = [strip_efficiency[i][1]/strip_efficiency[i][0] for i in strip_efficiency.keys() if strip_efficiency[i][0] != 0]
errors = [np.sqrt((strip_efficiency[i][1]/strip_efficiency[i][0])*(1-(strip_efficiency[i][1]/strip_efficiency[i][0]))/strip_efficiency[i][0]) for i in list(strip_efficiency.keys()) if strip_efficiency[i][0] != 0]
strips = [i for i in strip_efficiency.keys() if strip_efficiency[i][0] != 0]
plt.errorbar(strips, strip_efficiencies, yerr= errors, fmt = 'o')
plt.savefig(plots_dir+'HitEfficiency.pdf')
plt.close()
hitefficiency_data = {
    'efficiencies': strip_efficiencies, 
    'errors': errors,
    'strips': strips
}
with open(histo_dir+'HitEfficiency.pkl', 'wb') as file:
        pickle.dump(hitefficiency_data, file)


# labels to show on histogram x-axis
xlabels = ['Strip ID', 'Charge [ADC Counts]', 'Time [ns]', 'Centroid [mm]', 'Centroid Residues [mm]', 'Cluster Size', 'NClusters', 'Cluster Charge [ADC Counts]']
xlabels += [name + '_1st_coord' for name in xlabels]

for i,key in enumerate(histo_keys[0]):

    # plotting histograms
    plt.hist(histo_dict[chamber_keys[0]][key], bins = 100, histtype= 'step', label = chamber_keys[0])
    # for keys common to the two chambers superimpose them
    if key in histo_keys[1]:
        plt.hist(histo_dict[chamber_keys[1]][key], bins = 100, histtype= 'step', label = chamber_keys[1])
    plt.xlabel(xlabels[i])
    plt.ylabel("Entries")
    plt.legend()
    
    # also save histograms in pickle files
    with open(histo_dir+key+'_'+chamber_keys[0]+'.pkl', 'wb') as file:
        pickle.dump(histo_dict[chamber_keys[0]][key], file)
    
    if key in histo_keys[1]:
        with open(histo_dir+key+'_'+chamber_keys[1]+'.pkl', 'wb') as file:
            pickle.dump(histo_dict[chamber_keys[1]][key], file)
    
    plt.savefig(plots_dir+key+'.pdf')
    # plt.show()
    plt.close()
    

xlabels = ['Theta [deg]', '[mm]', 'Drift Speed [mm/ns]', '[ns]', '[ADC Counts]', '[mm]', '[mm]',' [mm]']
for i, key in enumerate(histo_keys[1][8:16]):
    # print(xlabels[i])
    plt.hist(histo_dict[chamber_keys[1]][key], bins = 30, histtype= 'step', label = chamber_keys[1])
    plt.xlabel(xlabels[i])
    plt.ylabel("Entries")
    plt.legend()
    with open(histo_dir+key+'_'+chamber_keys[1]+'.pkl', 'wb') as file:
        pickle.dump(histo_dict[chamber_keys[1]][key], file)
    plt.savefig(plots_dir+key+'.pdf')
    # plt.show()
    plt.close()






TMM_mean = np.array(histo_dict[chamber_keys[0]]['Centroid']).mean()
TMM_std = np.array(histo_dict[chamber_keys[0]]['Centroid']).std()

print(f'TMM centroid : {TMM_mean:.2f} +/- {TMM_std:.2f} mm')






# Copy the file 
shutil.copy(CONFIG_FILENAME, results_dir)

print(f"- {CONFIG_FILENAME} copied from {os.getcwd()} to {results_dir}")

print('Done!')

exit(1)


