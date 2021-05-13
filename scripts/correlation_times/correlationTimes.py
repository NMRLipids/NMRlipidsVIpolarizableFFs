#!/usr/bin/env python

"""
    calculation of effective correlation times from a MD trajectory

    meant for use with NMRLipids projects

    based on https://github.com/hsantila/Corrtimes

    Initial development by H. Antila, O. Ollila, H. Martinez-Seara, M. Miettinen



    ------------------------------------------------


"""

import numpy as np
import MDAnalysis as mda
import time
from OrderParameter import *
from numba import jit, njit
import numpy as np
from scipy import optimize


def calculate_teff(def_file, top, traj, residue_ids, output_file, addMemory = True):

    print('Starting the effective correlation times')
    start_time = time.time()  # to get the performance now, could be removed

    mol = mda.Universe(top, traj)
    nframes = mol.trajectory.n_frames

    if addMemory:
        print("transfering trajectory to memory")
        mol.transfer_to_memory()
    else:
        print("trajectory won't be transfered to memory. this calculation will probably be painfully slow")


    # note that we only calculate the correlations up to lag_time = n_frames/2
    # this is hardcoded now, could be changed in the future
    nframes_half = np.int(np.floor(nframes / 2))

    # correlation function will be calculated for nlipids number of lipids
    nlipids = len(residue_ids)

    # initialize the correlation times array
    correlations = np.zeros((nframes_half, nlipids))

    # getting the time stamps for each frame
    time_axis = np.array([i * mol.trajectory.dt for i in range(0, nframes_half)])

    OrdParam = find_OP(def_file, top, traj)

    correlation_times = np.zeros((len(OrdParam),6))

    w = 0
    for i, op in enumerate(OrdParam.values()):

        resops = op.get_op_res
        (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP

        lipid_counter = 0

        for lipid_id in residue_ids:

            atom1 = mol.universe.select_atoms('resid ' + str(lipid_id) + ' and name ' + op.atAname)
            atom2 = mol.universe.select_atoms('resid ' + str(lipid_id) + ' and name ' + op.atBname)

            print("calculating correlation times for lipid " + str(lipid_counter + 1) + "/" + str(len(residue_ids)) \
                   + " and pair " + op.atAname + "," + op.atBname)

            vector = prepare_corr_vector(mol, atom1.ids[0], atom2.ids[0])

            correlations[:, lipid_counter] = calc_corr(vector)

            lipid_counter = lipid_counter + 1


        correlation_times[w,:] = calc_teff_r1(time_axis, correlations, resops, op.avg, op.std, op.stem)
        w = w + 1

    outfile = open(output_file, 'w')

    line = "Atom1" + "  " + "Atom2" + "  " + "  " + "<R1>" + "  " + "std_err(R1)" + "  " + "convs" + "  " \
           + "Teff" + "  " + "Teff_min" + "  " + "Teff_max" + '\n'

    outfile.write(line)

#    np.mean(R1s), np.std(R1s, ddof=1) / np.sqrt(len(R1s) - 1), sum(convs) / float(len(convs)), teff, teff_min, teff_max
    w = 0
    for i, op in enumerate(OrdParam.values()):

        line = op.atAname + "  " + op.atBname + "  " + str(correlation_times[w,0]) + "  " + str(correlation_times[w,1]) \
                + "  " + str(correlation_times[w, 2]) + "  " + str(correlation_times[w,3]) + "  " + str(correlation_times[w,4]) \
                + "  " + str(correlation_times[w, 5]) + '\n'

        outfile.write(line)

        w = w + 1
    print('Finished with the effective correlation time analysis')
    print('It took ' + str(time.time() - start_time) + ' seconds')
    print('If you wish to publish these results, please cite https://dx.doi.org/10.1021/acs.jcim.0c01299')

#    print('timer for prepare corr ' + str(timer_prepare_corr_vector))
#    print('timer for calc corr ' + str(timer_calc_corr))
    return

def prepare_corr_vector(mol, atom1, atom2):

    vector = np.zeros((len(mol.trajectory), 3))

    w = 0
    for ts in mol.trajectory:

        vector[w,:] = mol.coord[atom1] - mol.coord[atom2]
        vector[w,:] = vector[w,:] / np.linalg.norm(vector[w,:])

        w = w + 1

    return vector

@njit()
def calc_corr(vector):

    nframes = len(vector)
    nframes_half = np.int(len(vector)/2)

    correlation = np.zeros((nframes_half))

    for tau in range(0, nframes_half):

        w = 0
        tmp = 0

        for t in range(0, nframes-tau):

            theta = np.dot(vector[t], vector[t+tau])

            tmp = tmp + p2(theta)

            w = w+1

        correlation[tau] = tmp/w

    return correlation

def calc_teff_r1(time_axis, correlation, resops, avg, std, stem):

    # op (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP

    convs = []
    Teffs = []
    Teffs_area = []
    R1s = []

    nlipids = correlation.shape[1]

    for j in range(1, nlipids + 1):
        if j == 1:
            out, times = [correlation[:, j - 1], time_axis]
            fvals = np.asmatrix(out).T
        else:
            out, times = [correlation[:, j - 1], time_axis]
            fvals = np.concatenate((fvals, np.asmatrix(out).T), axis=1)
        Teff, tau_eff_area, R1, conv = calc_corrtime_noread(out, times, resops[j - 1])
        Teffs.append(Teff)
        Teffs_area.append(tau_eff_area)
        R1s.append(R1)
        convs.append(conv)

    # analysis of Teff error starts here
    means = np.mean(fvals, axis=1)
    stems = np.std(fvals, axis=1, ddof=1) / np.sqrt(int(nlipids) - 1)
    teff, teff_min, teff_max = calc_corrtime_withee(times, means, stems, avg, stem)

    return np.mean(R1s), np.std(R1s, ddof=1)/np.sqrt(len(R1s)-1),sum(convs) / float(len(convs)), teff, teff_min, teff_max

@njit()
def p2(x):

    return 1/2 * (3*x*x-1) # 2nd order Legendre polynomial

def read_data(datafile):
    # for reading the correlation function data
    opf = open(datafile, 'r')
    lines = opf.readlines()
    data_times = []
    data_F = []
    for line in lines:
        if '#' in line:
            continue
        if '&' in line:
            continue
        if 'label' in line:
            continue
        parts = line.split()
        data_F.append(float(parts[1]))
        data_times.append(float(parts[0]))

    # data_out=np.empty((n, m))
    data_Fout = np.array(data_F)
    times_out = np.array(data_times)
    return data_Fout, times_out


def calc_corrtime_noread(corrF, Stimes, OP):
    # flag for convergence
    conv = 0
    # read in simulation data from gromacs g_rotacf
    # corrF, Stimes=read_data(corrfile)

    # normalized correlation fuction
    NcorrF = (corrF - OP ** 2) / (1 - OP ** 2);

    # Create correlation times from 1ps to 1micros
    Ctimes = 10 ** np.arange(0, 6 + 0.1, 0.1)

    # First, no forcing the plateou
    # create exponential functions and put them into a matrix
    n = len(Stimes)
    m = len(Ctimes)
    Cexp_mat = np.zeros((n, m))

    for i in range(0, n):
        for j in range(0, m):
            Cexp_mat[i, j] = np.exp(-Stimes[i] / Ctimes[j])

    Coeffs, res = optimize.nnls(Cexp_mat, NcorrF)

    # Effective correlation time from components, in units of sec

    Teff = sum(Coeffs * Ctimes * 0.001 * 10 ** (-9))

    # calculate t_eff from area
    dt = Stimes[2] - Stimes[1]
    pos = np.argmax(NcorrF < 0);

    if pos > 0:
        tau_eff_area = sum(NcorrF[0:pos]) * dt * 0.001 * 10 ** (-9);
        conv = 1
    else:
        tau_eff_area = sum(NcorrF) * dt * 0.001 * 10 ** (-9);
        conv = 0

    # Constants for calculating R1

    wc = 2 * np.pi * 125.76 * 10 ** 6;
    wh = wc / 0.25;

    # changin the unit of time permanently
    Ctimes = Ctimes * 0.001 * 10 ** (-9);

    J0 = 0
    J1 = 0
    J2 = 0
    Jw1 = 0

    for i in range(0, m):
        w = wh - wc

        J0 = J0 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = wc
        J1 = J1 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = wc + wh
        J2 = J2 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

    # R1=(2.1*10**9)*(J0+3*J1+6*J2)
    # note! R1's are additive. Nh from the Ferreira2015 paper correctly omitted here
    R1 = (22000 * 2 * np.pi) ** 2 / 20.0 * (1 - OP ** 2) * (J0 + 3 * J1 + 6 * J2)

    return Teff, tau_eff_area, R1, conv


def calc_corrtime_withee(times, mean, std, OP, deltaOP):
    conv = 0
    tau_eff_mean = -1
    tau_eff_min = -1
    tau_eff_max = -1
    # normalized correlation fuction
    # think about the factorial here
    NcorrF = (mean - OP ** 2) / (1 - OP ** 2);
    delta = std / (1 - OP ** 2) + deltaOP * np.absolute(2 * OP * (mean - 1.0)) / (1 - OP ** 2) ** 2
    NcorrF_max = NcorrF + delta
    NcorrF_min = NcorrF - delta

    # Create correlation times from 1ps to 1micros
    Ctimes = 10 ** np.arange(0, 6 + 0.1, 0.1)

    # First, no forcing the plateou
    # create exponential functions and put them into a matrix
    n = len(times)
    m = len(Ctimes)
    Cexp_mat = np.zeros((n, m))

    # calculate t_eff from area
    dt = times[2] - times[1]
    pos = np.argmax(NcorrF < 0);
    pos2 = np.argmax(NcorrF_min < 0);

    T = 10 ** -6
    if pos > 0:
        tau_eff_mean = sum(NcorrF[0:pos]) * dt * 0.001 * 10 ** (-9);
        tau_eff_min = sum(NcorrF_min[0:pos2]) * dt * 0.001 * 10 ** (-9);

        eps = NcorrF_max[pos]
        tend = times[pos] * 0.001 * 10 ** (-9)
        mu = T * (1 - np.exp(-(T - tend) / T))
        tau_eff_max = sum(NcorrF_max[0:pos]) * dt * 0.001 * 10 ** (-9) + mu * eps;

    elif pos2 > 0:
        tau_eff_min = sum(NcorrF_min[0:pos2]) * dt * 0.001 * 10 ** (-9);

        eps = NcorrF_max[-1]
        tend = times[-1] * 0.001 * 10 ** (-9)
        mu = T * (1 - np.exp(-(T - tend) / T))
        tau_eff_max = sum(NcorrF_max) * dt * 0.001 * 10 ** (-9) + mu * eps;
    else:
        tau_eff_min = sum(NcorrF_min) * dt * 0.001 * 10 ** (-9);

        eps = NcorrF_max[-1]
        tend = times[-1] * 0.001 * 10 ** (-9)
        mu = T * (1 - np.exp(-(T - tend) / T))
        tau_eff_max = sum(NcorrF_max) * dt * 0.001 * 10 ** (-9) + mu * eps;

    return float(tau_eff_mean), float(tau_eff_min), float(tau_eff_max)



### demonstration ###

topology = './tmp/0/traj_nowater_test.tpr'
trajectory = './tmp/0/traj_nowater.xtc'
residue_ids = np.arange(128)
def_file = "./POPC_lipid14_short.def"
output_file = './test/all_corr_times.txt'

## these are required for the initial compilation with Numba
calc_corr(np.zeros((10,3)))
p2(10)

## call to the function ##
calculate_teff(def_file, topology, trajectory, residue_ids, output_file)