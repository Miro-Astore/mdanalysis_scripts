import numpy as np

import pandas as pd 
import pdb
import rlcompleter
import MDAnalysis as mda 
from MDAnalysis.analysis import pca, align
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction
import argparse 
import sys 
import os
#TODO. do more things in memory

pdb.Pdb.complete=rlcompleter.Completer(locals()).complete

parser = argparse.ArgumentParser()

parser.add_argument('--topology', '-t', dest='top' , help='topology file',type=str) 
parser.add_argument('--trajectory', '-j', dest='traj', help='trajectory file',type=str) 
parser.add_argument('--out', '-o', dest='out_file',default="pca", help='Output file',type=str) 
parser.add_argument('--stride', '-dt', dest='stride' , default=1, help='Stride through trajectory skipping this many frames.',type=int) 
parser.add_argument('--n-components', '-n', dest='n_components' , default = 2, help='calculate this many Principal Components of the trajectory.',type=int) 

parser.add_argument('--s-components', dest='symmetry_list',  help='List of selections which compose a single symmetry group. This script expects symmetric groups to have the same number of atoms when the selection string is applied. For example you might have 4 identical chains so you would pass the argument "--s-components \'segid A\' \'segid B\' \'segid C\' \'segid D\'" to this script. The user is also warned that the groups will be treated cyclically. So in the previous example A will move to B, B to C, C to D and D to A. This will keep the relative arrangement of components so long as the user names the groups in the correct order. ',nargs="+") 

parser.add_argument('--selection', '-s', dest='selection_string', help='Selection string for fitting. Will be applied to both target and reference structures.',type=str, default="name CA") 
parser.add_argument('--n-frames', '-fn', dest='n_vis_frames', help='Number of frames to visualise',type=str, default=30) 
parser.add_argument('--in-mem', dest='in_mem', action="store_true", help='Do alignment processing in memory. I advise against doing this if the trajectory is large.') 

args = parser.parse_args()

selection_string = args.selection_string 

universe = mda.Universe(args.top,args.traj)
#universe.trajectory = universe.trajectory[0::args.stride] 


selection_object = universe.select_atoms(selection_string)


def copy_coords(ag):
        return ag.positions.copy()


if args.symmetry_list != None:
    #segids_coords = [ [selection_object.select_atoms()]]
    n_atoms_set =set(([selection_object.select_atoms(current_symmetry).n_atoms for current_symmetry in args.symmetry_list]))
    if len(n_atoms_set)  > 1 :
        raise ValueError('The selections you have chosen to delineate symmetry groups results in atom groups with different numbers of atoms. Check your structure and your selection strings.')

    temp_coords1 = [[selection_object.select_atoms(current_symmetry).positions for ts in universe.trajectory[::args.stride]] for current_symmetry in args.symmetry_list]


    segids = list(set(selection_object.segids))
    segids_indices = selection_object.atoms.segindices

    names = selection_object.atoms.names
    resids = selection_object.atoms.resids
    resnames = selection_object.atoms.resnames
    selection_object.write('temp_pdb_file.pdb')

    #empty_universe = mda.Universe.empty(selection_object.n_atoms, atom_resindex=resids, n_residues=selection_object.n_residues, n_segments=selection_object.n_segments, trajectory=True)
    analysis_universe = mda.Universe('temp_pdb_file.pdb')

    #empty_universe.add_TopologyAttr('name', names)
    #empty_universe.add_TopologyAttr('resid', resids)
    #empty_universe.add_TopologyAttr('resname', resnames)
    #empty_universe.add_TopologyAttr('segid', segids)

    coordinates = np.empty((len(universe.trajectory[::args.stride]) * len(args.symmetry_list), selection_object.n_atoms, 3))

    #coordinates = np.empty((universe.trajectory.n_frames * len(args.symmetry_list), selection_object.n_atoms, 3))

    #new_universe = selection_object.load_new(coordinates, order='fac')
    #print(np.shape(np.concatenate(temp_coords1 , axis = 1)))
    coordinates [0:len(universe.trajectory[::args.stride])] = np.concatenate(temp_coords1 , axis = 1)

    for i in range(0,len(args.symmetry_list)):
        temp_coords1 = np.concatenate (([temp_coords1[-1]],temp_coords1[0:-1]), axis=0)
        coordinates[i*len(universe.trajectory[::args.stride]):(i+1)*len(universe.trajectory[::args.stride])] =  np.concatenate(temp_coords1 , axis = 1)

#    for i in range(len(args.symmetry_list)):
#
#        
#        current_symmetry = args.symmetry_list [i]
#        temp_coords1 = [selection_object.select_atoms(current_symmetry).positions for ts in universe.trajectory[::args.stride]]
#
#        if i+1 !=  len(args.symmetry_list):
#            next_symmetry = args.symmetry_list [i+1]
#        else:
#            next_symmetry = args.symmetry_list [0]
#
#        temp_coords2 = [selection_object.select_atoms(next_symmetry).positions for ts in universe.trajectory[::args.stride]]

else:
    selection_object.write('/dev/shm/temp_pdb_file.pdb')

    analysis_universe = mda.Universe('/dev/shm/temp_pdb_file.pdb')
    coordinates = [selection_object.positions for ts in universe.trajectory[::args.stride]]
     
analysis_universe = analysis_universe.load_new (coordinates)


with mda.Writer(('/dev/shm/analysis_universe_traj.xtc'), analysis_universe.atoms.n_atoms) as W:
    for ts in analysis_universe.trajectory:
        W.write(analysis_universe.atoms)

ref_universe = mda.Universe('/dev/shm/temp_pdb_file.pdb')
analysis_universe = mda.Universe('/dev/shm/temp_pdb_file.pdb')
#load like this to over write first frame
analysis_universe.load_new('/dev/shm/analysis_universe_traj.xtc')


print('Aligning Trajectory')
if args.in_mem==True:
    aligner = align.AlignTraj(analysis_universe, analysis_universe, in_memory=True).run()
else:
    print('here')
    #pdb.set_trace()
    aligner = align.AlignTraj(analysis_universe, ref_universe, filename='thing.dcd')
    aligner.run()
    print('done')

    analysis_universe.atoms.write ('/dev/shm/analysis.pdb')

    #aligner_universe=mda.Universe(aligner,format=MemoryReader)

    #aligned_coords = AnalysisFromFunction(copy_coords, aligner_universe).run().results


analysis_universe = mda.Universe('/dev/shm/analysis.pdb')
analysis_universe = analysis_universe.load_new('/dev/shm/aligned.dcd')

#https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/aligning_trajectory.html
#analysis_universe.load_new(aligned_coords['timeseries'],format=MemoryReader)
#analysis_universe.load_new(aligned_coords['timeseries'],format=MemoryReader)


print('Calculating covariance.')

pc = pca.PCA(analysis_universe, n_components=args.n_components,verbose=True).run()


#transformed = pc.transform(selection_object, 1)

np.save(args.out_file + '_eigenvectors.npy', pc.results.p_components)
np.save(args.out_file + '_eigenvalues.npy', np.sqrt(pc.results.variance))

for i in range(args.n_components):

    # visualise loadings
    pc_vector = pc.results.p_components[:, i] * np.sqrt(pc.results.variance[i])

    #pc1 = pc.results.p_components[:, 0] 
    origin = np.zeros(len(pc_vector))

    pc_traj = np.linspace(-pc_vector,pc_vector,args.n_vis_frames)

    #trans1 = transformed[:, 0]

    #projected = np.outer([dummy_coords,dummy_coords2], pc1) + pc.mean.flatten()
    projected = pc_traj + pc.mean.flatten()

    coordinates = projected.reshape(len(pc_traj), -1, 3)


    proj1 = mda.Merge(selection_object)

    proj1.load_new(coordinates, order="fac")

    new_sel_object = proj1.select_atoms(selection_string)

    new_sel_object.write (args.out_file + str(i) + '.pdb')

    with mda.Writer((args.out_file + str(i) + '.xtc'), new_sel_object.n_atoms) as W:
        for ts in proj1.trajectory:
            W.write(new_sel_object)
