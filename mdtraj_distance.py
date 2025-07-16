import natsort
import mdtraj
from pathlib import Path
import argparse  

import os

def quickly_list_pdb_files(path):
    file_list = []
    with os.scandir(path) as entries:
        for entry in entries:
            if entry.is_file() and entry.name.endswith(".pdb") and entry.name.startswith("S_"):
                file_list.append(entry.name)
    file_list = natsort.natsorted (file_list) 
    return file_list 


def mdtraj_distance (filename, indices): 
    structure = mdtraj.load(filename)
    #convert to angstrom.
    distance = mdtraj.compute_distances(structure, [indices]) * 10
    return distance

parser = argparse.ArgumentParser()

parser.add_argument('--directory',dest='directory') 
args = parser.parse_args()
pdb_files = quickly_list_pdb_files(args.directory)

first_structure = mdtraj.load(os.path.join(args.directory,pdb_files[0]))
indices = first_structure.topology.select('name P and (resid 64 or resid  2)')
index1=indices[0]
index2=indices[1]


for i in range(len(pdb_files)) :
    print(mdtraj_distance(os.path.abspath(os.path.join(args.directory,pdb_files[i])),indices))

