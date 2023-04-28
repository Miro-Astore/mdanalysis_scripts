import numpy as np

import pandas as pd 
import MDAnalysis as mda 
from MDAnalysis.analysis import pca, align
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction
import argparse 
import sys 
import os

parser = argparse.ArgumentParser()

parser.add_argument('--topology', '-t', dest='top' , help='topology file',type=str) 
parser.add_argument('--trajectory', '-j', dest='traj', help='trajectory file',type=str) 
parser.add_argument('--out', '-o', dest='out_file',default="pca", help='Output file',type=str) 
parser.add_argument('--stride', '-dt', dest='stride' , default=1, help='Stride through trajectory skipping this many frames.',type=int) 
parser.add_argument('--n-components', '-n', dest='n_components' , default = 2, help='calculate this many Principal Components of the trajectory.',type=int) 
parser.add_argument('--selection', '-s', dest='selection_string', help='Selection string for fitting. Will be applied to both target and reference structures.',type=str, default="name CA") 
parser.add_argument('--n-frames', '-fn', dest='n_vis_frames', help='Number of frames to visualise',type=str, default=30) 
parser.add_argument('--in-mem', dest='in_mem', action="store_true", help='Do alignment processing in memory. I advise against doing this if the trajectory is large') 

args = parser.parse_args()

selection_string = args.selection_string 

universe = mda.Universe(args.top,args.traj)
#universe.trajectory = universe.trajectory[0::args.stride] 

selection_object = universe.select_atoms(selection_string)

print('Aligning Trajectory')

if args.in_mem==True:
    aligner = align.AlignTraj(universe, universe, select = args.selection_string, in_memory=True).run(start=0,stop=-1,step=args.stride)
else:
    aligner = align.AlignTraj(universe, universe, select = args.selection_string, in_memory=False).run(start=0,stop=-1,step=args.stride)

print('Calculating covariance.')

pc = pca.PCA(universe, select=args.selection_string,n_components=args.n_components,verbose=True).run(start=0,stop=-1,step=args.stride)


transformed = pc.transform(selection_object, 1)

print(pc.cov)
print(np.shape(pc.p_components))
np.save(args.out_file + '_eigenvectors.npy', pc.results.p_components)
np.save(args.out_file + '_eigenvalues.npy', np.sqrt(pc.results.variance))
print(np.shape(pc.results.variance))

for i in range(args.n_components):

    # visualise loadings
    pc_vector = pc.results.p_components[:, i] * np.sqrt(pc.results.variance[i])
    print(len(pc_vector))

    #pc1 = pc.results.p_components[:, 0] 
    #print(pc.results.variance[0])
    origin = np.zeros(len(pc_vector))

    pc_traj = np.linspace(-pc_vector,pc_vector,args.n_vis_frames)

    trans1 = transformed[:, 0]

    #projected = np.outer([dummy_coords,dummy_coords2], pc1) + pc.mean.flatten()
    projected = pc_traj + pc.mean.flatten()

    coordinates = projected.reshape(len(pc_traj), -1, 3)

    print(np.shape(coordinates))
    print(np.shape(selection_object.atoms.positions))

    proj1 = mda.Merge(selection_object)

    proj1.load_new(coordinates, order="fac")

    new_sel_object = proj1.select_atoms(selection_string)

    new_sel_object.write (args.out_file + str(i) + '.pdb')

    with mda.Writer((args.out_file + str(i) + '.xtc'), new_sel_object.n_atoms) as W:
        for ts in proj1.trajectory:
            W.write(new_sel_object)
