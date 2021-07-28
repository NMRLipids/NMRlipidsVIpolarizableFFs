### for completeness, this is the corrtimes.py file for the databank port ###
### i'll do my updates/changes on this file instead of the corrtimes.py ###

import sys
import numpy as np
from scipy import optimize
from numba import jit


# script for calculating effective correlation time and R1 by H. Antila, with help from S. Ollila and T. Ferreira
# usage: python corrtimes.py rotac.dat opvalue
# files needed: rotac.dat = the rotation correlation ftio that corresponds to the opvalue. Expects 2 column format "time" "f_value"


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

def calc_tau(nframes, trajectory, atom1, atom2):

    
    nframes_half = np.int(np.floor(nframes/2));

    norm_coor = np.zeros((nframes, 3))
    for i in range(0,nframes):

        norm_coor[i,:] = (trajectory[i][atom1] - trajectory[i][atom2]) / np.linalg.norm(trajectory[i][atom1] - trajectory[i][atom2])

    #res = calc_tau_fast(nframes, atom1_coor, atom2_coor)
    res = calc_tau_faster(nframes, norm_coor)
    #print(res)
    return res

@jit(nopython=True)
def calc_tau_fast(nframes,atom1_coor,atom2_coor):

    nframes_half = np.int(np.floor(nframes/2))
    corr = np.zeros((nframes_half));
    for tau in range(0, nframes_half):

        w = 0;
        tmp = 0;

        for t in range(0, nframes - tau):

            vec1 = (atom1_coor[t] - atom2_coor[t]) / np.linalg.norm(atom1_coor[t] - atom2_coor[t]);
            vec2 = (atom1_coor[t + tau] - atom2_coor[t + tau]) / np.linalg.norm(atom1_coor[t + tau] - atom2_coor[t + tau]);

            theta = np.dot(vec1, vec2)

            tmp = tmp + p2(theta)

            w = w + 1

        corr[tau] = tmp/w # correlation function at a given lag time tau for a given atom pair atom1 atom2

    print(corr)
    return corr # correlation function for all possible tau for a given atom pair atom1/atom2


def calc_tau_faster(nframes, atom_coor):

    nframes_half = np.int(np.floor(nframes/2))
    corr = np.zeros((nframes_half))

    #coor = (atom1_coor - atom2_coor) / np.linalg.norm(atom1_coor - atom2_coor, axis = 1)[:,None]


    for tau in range(0, nframes_half):

        res = np.einsum('ij,ij->i', atom_coor[0:nframes-tau,:], atom_coor[tau:nframes,:])

        corr[tau] = np.sum( p2(res)/(nframes-tau) )

    return corr

@jit(nopython=True)
def p2(x):

    return 1/2 * (3*x*x-1)
