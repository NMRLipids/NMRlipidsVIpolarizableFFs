#!/usr/bin/env python
"""
 calculation of order parameters of lipid bilayers
 from a MD trajectory
 
 meant for use with NMRlipids projects
 ------------------------------------------------------------
 Made by Joe,  Last edit 2017/02/02
------------------------------------------------------------
 input: Order parameter definitions
        gro and xtc file (or equivalents)
 output: order parameters (2 textfiles)
--------------------------------------------------------
"""

# coding: utf-8

import MDAnalysis as mda
import numpy as np
import math
import os, sys
import warnings
from optparse import OptionParser
from collections import OrderedDict

#k_b = 0.0083144621  #kJ/Mol*K
#f_conc=55430  # factor for calculating concentrations of salts from numbers of ions/waters; in mM/L

bond_len_max=1.5  # in A, max distance between atoms for reasonable OP calculation
bond_len_max_sq=bond_len_max**2

#%%
class OrderParameter:
    """
    Class for storing&manipulating
    order parameter (OP) related metadata (definition, name, ...)
    and OP trajectories
    and methods to evaluate OPs.
    """
    def __init__(self, name, resname, atom_A_name, atom_B_name, *args):
        """
        it doesn't matter which comes first,
        atom A or B, for OP calculation.
        """
        self.name = name             # name of the order parameter, a label
        self.resname = resname       # name of residue atoms are in
        self.atAname = atom_A_name
        self.atBname = atom_B_name
        for field in self.__dict__:
            if not isinstance(field, str):
                warnings.warn("provided name >> {} << is not a string! \n \
                Unexpected behaviour might occur.")#.format(field)
            else:
                if not field.strip():
                    raise RuntimeError("provided name >> {} << is empty! \n \
                    Cannot use empty names for atoms and OP definitions.")#.format(field)
        # extra optional arguments allow setting avg,std values -- suitable for reading-in results of this script
        if len(args) == 0:
            self.avg = None
            self.std = None
            self.stem = None
        elif len(args) == 2:
            self.avg = args[0]
            self.std = args[1]
            self.stem = None
        else:
            warnings.warn("Number of optional positional arguments is {len}, not 2 or 0. Args: {args}\nWrong file format?")
        self.traj = []  # for storing OPs
        

    def calc_OP(self, atoms):
        """
        calculates Order Parameter according to equation
        S = 1/2 * (3*cos(theta)^2 -1)
        """
        vec = atoms[1].position - atoms[0].position
        d2 = np.square(vec).sum()
	
        if d2>bond_len_max_sq:
            at1=atoms[0].name
            at2=atoms[1].name
            resnr=atoms[0].resid
            d=math.sqrt(d2)
            warnings.warn("Atomic distance for atoms \
            {at1} and {at2} in residue no. {resnr} is suspiciously \
            long: {d}!\nPBC removed???")
        cos2 = vec[2]**2/d2
        S = 0.5*(3.0*cos2-1.0)
        return S


    @property
    def get_avg_std_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        # convert to numpy array
        return (np.mean(self.traj), np.std(self.traj))
    @property
    def get_avg_std_stem_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        std=np.std(self.traj)
        # convert to numpy array
        return (np.mean(self.traj),std,std/np.sqrt(len(self.traj)-1))
    @property
    def get_op_res(self):
        """
	Provides average and stddev of all OPs in self.traj
        """

	# convert to numpy array
        return self.traj
     


def read_trajs_calc_OPs(ordPars, top, trajs):
    """
    procedure that
    creates MDAnalysis Universe with top,
    reads in trajectories trajs and then
    goes through every frame and
    evaluates each Order Parameter "S" from the list of OPs ordPars.
    ordPars : list of OrderParameter class
       each item in this list describes an Order parameter to be calculated in the trajectory
    top : str
        filename of a top file (e.g. conf.gro)
    trajs : list of strings
        filenames of trajectories
    """
    # read-in topology and trajectory
    mol = mda.Universe(top, trajs)

    # make atom selections for each OP and store it as its attribute for later use in trajectory
    for op in ordPars.values():
        # selection = pairs of atoms, split-by residues
        selection = mol.select_atoms("resname {rnm} and name {atA} {atB}".format(
                                    rnm=op.resname, atA=op.atAname, atB=op.atBname)
                                    ).atoms.split("residue")
        for res in selection:
            # check if we have only 2 atoms (A & B) selected
            if res.n_atoms != 2:
                print(res.resnames, res.resids)
                for atom in res.atoms:
                    print(atom.name, atom.id)
                atA=op.atAnam
                atB=op.atBname
                nat=res.n_atoms
                warning.warn("Selection >> name {atA} {atB} << \
                contains {nat} atoms, but should contain exactly 2!")
        op.selection = selection

    # go through trajectory frame-by-frame
    Nres=len(op.selection)
    Nframes=len(mol.trajectory)	
    for op in ordPars.values():
        op.traj=[0]*Nres		
    for frame in mol.trajectory:
	
        for op in ordPars.values():
            for i in range(0,Nres):
                residue=op.selection[i]	
                S = op.calc_OP(residue)
                op.traj[i]=op.traj[i]+S/Nframes
            #op.traj.append(np.mean(tmp))
		
        #print "--", mol.atoms[0].position
 #   for op in ordPars.values():
 #	op.traj=op.traj/Nframes   
	  	

def parse_op_input(fname):
    """
    parses input file with Order Parameter definitions
    file format is as follows:
    OP_name    resname    atom1    atom2
    (flexible cols)
    fname : string
        input file name
    returns : dictionary 
        with OrderParameters class instances
    """
    ordPars = OrderedDict()
    try:
        with open(fname,"r") as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    items = line.split()
                	
                    ordPars[items[0]] = OrderParameter(*items)

    except:
        inpf=opts.inp_fname
        raise RuntimeError("Couldn't read input file >> {inpf} <<")
    return ordPars



#%%

def find_OP(inp_fname, top_fname, traj_fname):
    ordPars = parse_op_input(inp_fname)
		
 #   for file_name in os.listdir(os.getcwd()):
 #       if file_name.startswith(traj_fname):
 #           trajs.append(file_name)
            
    read_trajs_calc_OPs(ordPars, top_fname, traj_fname)
    
    return ordPars
