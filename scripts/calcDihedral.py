#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np

import yaml
import json
import matplotlib.pyplot as plt
import mdtraj
import urllib.request
import seaborn as sns

from OrderParameter import *

# Download link
def download_link(doi, file):
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        return 'https://zenodo.org/record/' + zenodo_entry_number + '/files/' + file
    else:
        print ("DOI provided: {0}".format(doi))
        print ("Repository not validated. Please upload the data for example to zenodo.org")
        return ""
    
# read mapping file
def read_mapping_file(mapping_file, atom1):
    with open(mapping_file, 'rt') as mapping_file:
            for line in mapping_file:
                if atom1 in line:
                    m_atom1 = line.split()[1]
    return m_atom1

def read_mapping_filePAIR(mapping_file, atom1, atom2):
    with open(mapping_file, 'rt') as mapping_file:
            print(mapping_file)
            for line in mapping_file:
                if atom1 in line:
                    m_atom1 = line.split()[1]
#                    print(m_atom1)
                if atom2 in line: 
                    m_atom2 = line.split()[1]
#                    print(m_atom2)
    return m_atom1, m_atom2

def make_positive_angles(x):
    for i in range(len(x)):
        if x[i] < 0:
            x[i] = np.degrees(x[i]) + 360
        else:
            x[i] = np.degrees(x[i])
    return x


def calcDihedrals(lipids):
    colors = {'POPC' :'black','POPS':'red','POPE':'blue','POPG':'green'}
    for subdir, dirs, files in os.walk(r'../../Data/Simulations/'):
        for filename in files:
            filepath = subdir + os.sep + filename
            if filepath.endswith("README.yaml"):
                READMEfilepath = subdir + '/README.yaml'
                with open(READMEfilepath) as yaml_file:
                    readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
                    for molname in lipids:
                        doi = readme.get('DOI')
                        trj = readme.get('TRJ')
                        tpr = readme.get('PDB')
                        trj_name = subdir + '/' + readme.get('TRJ')[0][0]
                        tpr_name = subdir + '/' + readme.get('PDB')[0][0]
                        gro_name = subdir + '/conf.gro'
                        trj_url = download_link(doi, trj[0][0])
                        tpr_url = download_link(doi, tpr[0][0])
                        #Download tpr and xtc files to same directory where dictionary and data are located
                        if (not os.path.isfile(tpr_name)):
                            response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
                        if (not os.path.isfile(trj_name)):
                            response = urllib.request.urlretrieve(trj_url, trj_name)
                
#                        if sum(readme['N' + molname]) > 0:
#                            print('Analyzing '+molname+' in '+filepath)
#                            #fig= plt.figure(figsize=(12,9))
#                            if (not os.path.isfile(gro_name)):
#                                get_ipython().system('echo System | gmx trjconv -f {trj_name} -s {tpr_name}  -dump 0 -o {gro_name}')
#                        
                        xtcwhole= trj_name
#                            if (not os.path.isfile(xtcwhole)):
#                                get_ipython().system('echo System | gmx trjconv -f {trj_name} -s {tpr_name} -o {xtcwhole} -pbc mol ')
#                        
#                            try:
#                                traj = mdtraj.load(xtcwhole, top = gro_name)
#                                #print(lipid)
#                            except FileNotFoundError or OSError:
#                                continue
                        mapping_file = './mapping_files/'+readme['MAPPING_DICT'][molname] # readme.get('MAPPING')[0][0]


                        dihedrals = [['M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M', 'M_G3N6_M'],
           ['M_G3O3_M','M_G3C4_M','M_G3C5_M','M_G3C6_M'],
           ['M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M'],
           ['M_G3O1_M', 'M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M'],
           ['M_G3_M', 'M_G3O1_M', 'M_G3P2_M', 'M_G3O3_M'],
           ['M_G2_M', 'M_G3_M', 'M_G3O1_M', 'M_G3P2_M'],
           ['M_G1_M', 'M_G2_M', 'M_G3_M', 'M_G3O1_M']
           ]

                        gro_name = tpr_name 
                        traj = mdtraj.load(xtcwhole, top = gro_name)
                        for i in dihedrals:
                            DIHatoms = i
                            try:
                                atom1 = read_mapping_file(mapping_file, DIHatoms[0])
                                atom2 = read_mapping_file(mapping_file, DIHatoms[1])
                                atom3 = read_mapping_file(mapping_file, DIHatoms[2])
                                atom4 = read_mapping_file(mapping_file, DIHatoms[3])
                                print(atom1,atom2,atom3,atom4)
                            except:
                                print(atom1 + " and " + atom2 + " not found in the mapping file.")
                                print("Some atom not found in the mapping file.")
                                continue
                            index = [ [] for i in range(traj.topology.n_residues)]
                            dihRESULT = []
                            for residue in traj.topology.residues:
                             if residue.name == 'POPC':
                                atom1ind=traj.topology.select("name == " + atom1 + " and resid == " + str(residue.index))
                                atom2ind=traj.topology.select("name == " + atom2 + " and resid == " + str(residue.index))
                                atom3ind=traj.topology.select("name == " + atom3 + " and resid == " + str(residue.index))
                                atom4ind=traj.topology.select("name == " + atom4 + " and resid == " + str(residue.index))
                                if(len(atom1ind) > 0):
                                    index[residue.index].append(atom1ind[0])
                                    index[residue.index].append(atom2ind[0])
                                    index[residue.index].append(atom3ind[0])
                                    index[residue.index].append(atom4ind[0])
                                    dihRESULT.append(mdtraj.compute_dihedrals(traj,[index[residue.index]]))
                            dihRESULT = [make_positive_angles(x) for x in dihRESULT ]
                            dist = [ 0 for i in range(len(dihRESULT))]
                            distSUM = [ 0 for i in range(359)]
                            for i in range(len(dihRESULT)):
                                dist[i] =  plt.hist(dihRESULT[i], range(360),density=True);
                                distSUM = np.add(distSUM,dist[i][0])
                            
                        
                            distSUM = [x / len(dihRESULT) for x in distSUM]
                            xaxis = [ 0 for i in range(len(dist[0][1])-1)]
                            for i in range(len(dist[0][1])-1):
                                xaxis[i]=(dist[0][1][i])
                            
                            #plt.plot(xaxis,distSUM,color=colors[molname], label = readme.get('SYSTEM'))[0] 
                            #plt.legend()
                            #plt.xlabel("Angle (°)")
                            dihedralFOLDERS = subdir.replace("Simulations","dihedral")
                            get_ipython().system('mkdir -p {dihedralFOLDERS}')
                            get_ipython().system('cp {READMEfilepath} {dihedralFOLDERS}')
                            outfile=open(str(dihedralFOLDERS) + '/' + DIHatoms[0] + DIHatoms[1] + DIHatoms[2] + DIHatoms[3] +'.dat','w')
                            print(outfile)
                            for i in range(len(xaxis)):
                                outfile.write(str(xaxis[i]) + " " + str(distSUM[i])+'\n')
                            outfile.close()
                            plt.close()



dihedrals=[['M_G3_M', 'M_G3O1_M', 'M_G3P2_M', 'M_G3P2O1_M'],
           ['M_G3_M', 'M_G3O1_M', 'M_G3P2_M', 'M_G3P2O2_M'],
           ['M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M'],
           ['M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M', 'M_G3N6_M'],
           ['M_G3O1_M', 'M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M'],
           ['M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M'],
           ['M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M', 'M_G3N6_M'],
           ['M_G3O3_M','M_G3C4_M','M_G3C5_M','M_G3C6_M'],
           ['M_G3O3_M','M_G3C4_M','M_G3C5_M','M_G3C5O1_M'],
           ['M_G1C17_M', 'M_G1C16_M', 'M_G1C15_M', 'M_G1C14_M'],
           ['M_G1C16_M', 'M_G1C15_M', 'M_G1C14_M', 'M_G1C13_M'],
           ['M_G1C15_M', 'M_G1C14_M', 'M_G1C13_M', 'M_G1C12_M'],
           ['M_G1C14_M', 'M_G1C13_M', 'M_G1C12_M', 'M_G1C11_M'],
           ['M_G1C13_M', 'M_G1C12_M', 'M_G1C11_M', 'M_G1C10_M'],
           ['M_G1C12_M', 'M_G1C11_M', 'M_G1C10_M', 'M_G1C9_M'],
           ['M_G1C11_M', 'M_G1C10_M', 'M_G1C9_M', 'M_G1C8_M'],
           ['M_G1C10_M', 'M_G1C9_M', 'M_G1C8_M', 'M_G1C7_M'],
           ['M_G1C9_M', 'M_G1C8_M', 'M_G1C7_M', 'M_G1C6_M'],
           ['M_G1C8_M', 'M_G1C7_M', 'M_G1C6_M', 'M_G1C5_M'],
           ['M_G1C7_M', 'M_G1C6_M', 'M_G1C5_M', 'M_G1C4_M'],
           ['M_G1C6_M', 'M_G1C5_M', 'M_G1C4_M', 'M_G1C3_M'],
           ['M_G1C5_M', 'M_G1C4_M', 'M_G1C3_M', 'M_G1C2O1_M'],
           ['M_G1C5_M', 'M_G1C4_M', 'M_G1C3_M', 'M_G1C2_M'],
           ['M_G1C4_M', 'M_G1C3_M', 'M_G1C2_M', 'M_G1O1_M'],
           ['M_G1C3_M', 'M_G1C2_M', 'M_G1O1_M', 'M_G1_M'],
           ['M_G2C19_M', 'M_G2C18_M', 'M_G2C17_M', 'M_G2C16_M'],
           ['M_G2C18_M', 'M_G2C17_M', 'M_G2C16_M', 'M_G2C15_M'],
           ['M_G2C17_M', 'M_G2C16_M', 'M_G2C15_M', 'M_G2C14_M'],
           ['M_G2C16_M', 'M_G2C15_M', 'M_G2C14_M', 'M_G2C13_M'],
           ['M_G2C15_M', 'M_G2C14_M', 'M_G2C13_M', 'M_G2C12_M'],
           ['M_G2C14_M', 'M_G2C13_M', 'M_G2C12_M', 'M_G2C11_M'],
           ['M_G2C13_M', 'M_G2C12_M', 'M_G2C11_M', 'M_G2C10_M'],
           ['M_G2C12_M', 'M_G2C11_M', 'M_G2C10_M', 'M_G2C9_M'],
           ['M_G2C11_M', 'M_G2C10_M', 'M_G2C9_M', 'M_G2C8_M'],
           ['M_G2C10_M', 'M_G2C9_M', 'M_G2C8_M', 'M_G2C7_M'],
           ['M_G2C9_M', 'M_G2C8_M', 'M_G2C7_M', 'M_G2C6_M'],
           ['M_G2C8_M', 'M_G2C7_M', 'M_G2C6_M', 'M_G2C5_M'],
           ['M_G2C7_M', 'M_G2C6_M', 'M_G2C5_M', 'M_G2C4_M'],
           ['M_G2C6_M', 'M_G2C5_M', 'M_G2C4_M', 'M_G2C3_M'],
           ['M_G2C5_M', 'M_G2C4_M', 'M_G2C3_M', 'M_G2C2O1_M'],
           ['M_G2C5_M', 'M_G2C4_M', 'M_G2C3_M', 'M_G2C2_M'],
           ['M_G2C4_M', 'M_G2C3_M', 'M_G2C2_M', 'M_G2O1_M'],
           ['M_G2C3_M', 'M_G2C2_M', 'M_G2O1_M', 'M_G2_M'],
           ['M_G1O1_M', 'M_G1_M', 'M_G2_M', 'M_G3_M'],
           ['M_G1O1_M', 'M_G1_M', 'M_G2_M', 'M_G2O1_M'],
           ['M_G2C2_M', 'M_G2O1_M', 'M_G2_M', 'M_G1_M'],
           ['M_G1_M', 'M_G2_M', 'M_G3_M', 'M_G3O1_M'],
           ['M_G2O1_M', 'M_G2_M', 'M_G3_M', 'M_G3O1_M'],
           ['M_G2_M', 'M_G3_M', 'M_G3O1_M', 'M_G3P3_M'],
           ['M_G1C2O1_M', 'M_G1C2_M', 'M_G1O1_M', 'M_G1_M'],           
           ['M_G2C2O1_M', 'M_G2C2_M', 'M_G2O1_M', 'M_G2_M']
]

lipids = {'POPC'}
#lipids = {'POPG'}
#lipids = {'POPG','POPE','POPC'}
#for i in dihedrals:
#    DIHatoms = i
calcDihedrals(lipids)

#mapping_file = "./mapping_files/mappingPOPEcharmm.txt"
#dihedrals = parseDihedralInput(mapping_file)
#parseDihedralInput(mapping_file)

