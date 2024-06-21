import MDAnalysis as mda
import tqdm
import argparse
import jax
import jax.numpy as  jnp
import sys
import pdb
import rlcompleter

pdb.Pdb.complete=rlcompleter.Completer(locals()).complete

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument ('--topology',type=str) 
parser.add_argument ('--trajectory',type=str) 
parser.add_argument ('--n_batches',default = 10000,type=int) 
parser.add_argument ('--selection',default = 'protein',type=str) 
parser.add_argument ('-o',default='bond_lengths.npy',type=str) 
parser.add_argument ('--stride',default=10,type=int) 
args = parser.parse_args()

u = mda.Universe(args.topology, args.trajectory, in_memory=True)
selected_atoms = u.select_atoms(args.selection)

stride = args.stride

frames_array = range(u.trajectory.n_frames)[::stride]
total_length = np.zeros(len(frames_array))
for i in tqdm.tqdm(range(len(frames_array))):
    u.trajectory[frames_array[i]]
    total_length[i] = np.sum([ bond.length() for bond in selected_atoms.atoms.bonds])
    
    # Calculate grid dimensions
np.save(args.o, total_length)
