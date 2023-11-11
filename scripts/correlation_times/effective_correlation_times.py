### This is the main script to calculate the effective correlation times ###
### originally developed by Hanne Antila ###
### further improved by Batuhan Kav ###


import MDAnalysis as mda
from OrderParameter import *
import warnings
from corrtimesDatabank import *
import time
import argparse

def teff_calculator_databank(topology, trajectory, dt, lipid_type, residue_ids, def_file, write_directory = './', ID = 0):

    print('Starting the effective correlation times')
    start_time = time.time() 

    u = mda.Universe(topology, trajectory, dt=dt, in_memory = False)
    nframes = u.trajectory.n_frames
    # note that we only calculate the correlations upto lag_time = n_frames/2
    # this is hardcoded now, could be changed in the future
    nframes_half = np.int(np.floor(nframes/2))

    # correlation function will be calculated for nlipids number of lipids
    nlipids = len(residue_ids)

    # initialize the correlation times array
    correlation_times = np.zeros((nframes_half, nlipids))
    
    # getting the time stamps for each frame
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
    print(f'It took {(time.time() - start_time)/60} minutes')
    print('If you wish to publish these results, please cite https://dx.doi.org/10.1021/acs.jcim.0c01299')

    return correlation_times, OrdParam

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate effective correlation times as described in https://dx.doi.org/10.1021/acs.jcim.0c01299')

    # Required string arguments
    parser.add_argument('--topology', type=str, required=True, help='Path to the topology file, must be compatible with MDAnalysis file formats')
    parser.add_argument('--trajectory', type=str, required=True, help='Path to the trajectory file, must be compatible with MDAnalysis file formats')
    parser.add_argument('--lipid_type', type=str, required=True, help='Type of lipid as defind in the definition file')
    parser.add_argument('--def_file', type=str, required=True, help='Path to the definition file that contains the atom and residue names for the order parameters')
    parser.add_argument('--write_directory', type=str, required=True, help='Directory to write output files')
    parser.add_argument('--number_of_lipid_residues', type = int, required = True, help = 'Number of lipid_type residues in the trajectory')

    # Required float argument
    parser.add_argument('--dt', type=float, required=True, help='Time step for the analysis')

    args = parser.parse_args()

    topology = args.topology
    trajectory = args.trajectory
    lipid_type = args.lipid_type
    def_file = args.def_file
    write_directory = args.write_directory
    dt = args.dt
    number_of_lipid_residues = args.number_of_lipid_residues
    residue_ids = np.arange(number_of_lipid_residues) + 1 # we're keeping this variable for historical reasons. can be changed in the future

    teff_calculator_databank(topology, trajectory,dt, lipid_type, residue_ids, def_file, write_directory, ID = 0)

