### This is the main script to calculate the effective correlation times ###
### originally developed by Hanne Antila, Samuli Ollila, Markus Miettinen, Hector ... ##
### initial adaptation to the NMRLipids Databank by Batuhan Kav ###

# intended use #

#teff_calculator_databank(topology, trajectory, dt, lipid_type, residue_ids, op_list, write_directory, def_file, ID = 0)

###############
###############
## IMPORTANT REMARKS ##
# topology = topology file for MDAnalysis (psf, tpr, top, pdb, gro, ...)
# trajectory = trajectory file for MDAnalysis
# dt = saving frequency in ps
# lipid_type = name of the lipid molecules, i.e. POPC. this is used later to save files with proper naming
# residue_ids = zero-based list of the lipid_type
# op_list = a numpy array that contains the label, atom1, atom2 for a given order parameter
# ex: np.array(['beta','C32','H2A'])
# write_directory = directory where the output will be written
# def_file = Order parameter definition file
# ID = simulation ID. could be removed in the future

import MDAnalysis as mda
from OrderParameter import *
import warnings
from corrtimesDatabank import *
import time

def teff_calculator_databank(topology, trajectory, dt, lipid_type, residue_ids, def_file, write_directory = './', ID = 0):

    print('Starting the effective correlation times')
    start_time = time.time() # to get the performance now, could be removed

    u = mda.Universe(topology, trajectory, dt=dt, in_memory = False)
   # u = mda.Universe(topology, trajectory, in_memory = True)
    nframes = u.trajectory.n_frames
    # note that we only calculate the correlations upto lag_time = n_frames/2
    # this is hardcoded now, could be changed in the future
    nframes_half = np.int(np.floor(nframes/2))

    # correlation function will be calculated for nlipids number of lipids
    nlipids = len(residue_ids)

    # initialize the correlation times array
    correlation_times = np.zeros((nframes_half, nlipids))
    
    # getting the time stamps for each frame
    #HANNE: should dt be passed as input or read from the traj? 
    time_axis = np.array([i * u.trajectory.dt for i in range(0, nframes_half)])
    # calculating the order parameters for normalization #
    # Eq. 7 from https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.0c01299 #
    print('Calculating order parameters')
    OrdParam = find_OP(def_file, topology, trajectory)






    outfile = open(write_directory + '/times_' + str(ID) + '_' + lipid_type + '_R1.txt', 'a')
    outfile_teffs = open(write_directory + '/times_' + str(ID) + '_' + lipid_type + '_teffs.txt', 'a')

    for i, op in enumerate(OrdParam.values()):

                print(op.name)

                resops = op.get_op_res
                (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP
                outfile_teffs.flush()
                outfile.flush()  
                convs = []
                Teffs = []
                Teffs_area = []
                R1s = []

                for j in range(1, nlipids + 1):
                    
                    atom1 = u.universe.select_atoms('resid ' + str(j) + ' and name ' + op.atAname)
                    #print(atom1)
                    atom2 = u.universe.select_atoms('resid ' + str(j) + ' and name ' + op.atBname)
                    #print(atom2)
                    correlation_times[:, j-1] = calc_tau(len(u.trajectory),u.trajectory, atom1.ix[0], atom2.ix[0])
 
                    Teff, tau_eff_area, R1, conv = calc_corrtime_noread(correlation_times[:,j-1], time_axis, resops[j - 1])
                    #the teffs from here are obsolete, can be saved to compare the area and the decomposition estimate for teff for individual lipids
                    #Teffs.append(Teff)
                    #Teffs_area.append(tau_eff_area)
                    R1s.append(R1)
                    convs.append(conv)
                np.savetxt(op.name+'_ftions.dat', np.hstack((time_axis.reshape(nframes_half,1),correlation_times)) )			
                line2=str(op.name)+" "+str(op.avg)+" "+str(op.stem)+" "+str(np.mean(R1s))+" "+str(np.std(R1s, ddof=1)/np.sqrt(len(R1s)-1))+" "+str(sum(convs)/float(len(convs)))+'\n'
                outfile.write(line2)

                #analysis of Teff error starts here
                means=np.mean(correlation_times,axis=1)
                stems=np.std(correlation_times,axis=1,ddof=1)/np.sqrt(int(nlipids)-1)
                teff,teff_min,teff_max=calc_corrtime_withee(time_axis,means, stems, op.avg, op.stem)
                line3=str(op.name)+" "+str(op.avg)+" "+str(teff)+" "+str(teff_min)+" "+str(teff_max)+'\n'
                outfile_teffs.write(line3)



    print('Finished with the effective correlation time analysis')
    print('It took ' + str(time.time() - start_time) + ' seconds')
    print('If you wish to publish these results, please cite https://dx.doi.org/10.1021/acs.jcim.0c01299')
#    np.savetxt(write_directory + '/correlation_function.csv',correlation_times,delimiter=',')
     #Hanne: this returns the correlation time functions for the last op only?		
    return correlation_times, OrdParam

######## demonstration ########
topology = './step2_drude.psf'
trajectory = './wrapped_full_remove_first_5000frames.xtc'
lipid_type = 'POPC'
residue_ids = np.arange(128) + 1
#op_list = np.array([['gamma', 'C15','H15C'],['beta1','C12','H12A'],['beta2','C12','H12B'],['alpha1','C11','H11A'],['alpha2','C11','H11B'],['g1_1','C3','HX'],['g1_2','C3', 'HY'],['g2', 'C2','HS'],['g3_1','C1','HA'],['g3_2','C1','HB']])
def_file = "./Slipids_DOPC_reduced.def"
write_directory = './'
dt=10
cor = teff_calculator_databank(topology, trajectory,dt, lipid_type, residue_ids, def_file, write_directory, ID = 0)

