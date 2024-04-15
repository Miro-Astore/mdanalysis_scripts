import numpy as np
import shutil

import MDAnalysis as mda 
from MDAnalysis.analysis import pca, align
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction
import uuid
import os
from sklearn.decomposition import PCA
import os 

import uuid

print('imported modules') 

import argparse 
#TODO. do more things in memory

import warnings
warnings.filterwarnings('ignore')

import pdb
import rlcompleter
pdb.Pdb.complete=rlcompleter.Completer(locals()).complete

parser = argparse.ArgumentParser()

parser.add_argument('--topology', '-t', dest='topology' , help='topology file',type=str) 
parser.add_argument('--trajectory', '-j', dest='traj', help='trajectory file',type=str) 
parser.add_argument('--out', '-o', dest='out_file',default="pca", help='Output file',type=str) 
parser.add_argument('--stride', '-dt', dest='stride' , default=1, help='Stride through trajectory skipping this many frames.',type=int) 
parser.add_argument('--no-std', dest='standardise_bool', action="store_false", help='Choose whether or not to standardise coordinates. Default behaviour is to standardise the cordinates.') 
parser.add_argument('--save-proj', dest='save_projection', action="store_true", help='Choose whether or not to standardise coordinates. Default behaviour is to standardise the cordinates.',default=False) 
parser.add_argument('--n-components', '-n', dest='n_components' , default = 2, help='calculate this many Principal Components of the trajectory.',type=int) 
parser.add_argument('--vis-multiply', dest='extend' , default = 1, help='multiply the eigenvectors by this magnitude for visualisation purposes.',type=float) 
parser.add_argument('--dir-root', dest='dir_root', help='Determines where the output will go. Default is to put it all into a temporary directory which is deleted after the run.') 

parser.add_argument('--ref', dest='reference',  help='Reference structure, used for alignment') 
parser.add_argument('--ref-top', dest='ref_top',  help='Reference topology, used for alignment') 
parser.add_argument('--symmetry-list', dest='symmetry_list',  help='List of selections which compose a single symmetry group. This script expects symmetric groups to have the same number of atoms when the selection string is applied. For example you might have 4 identical chains so you would pass the argument "--symmetry-list \'segid A\' \'segid B\' \'segid C\' \'segid D\'" to this script. The user is also warned that the groups will be treated cyclically. So in the previous example A will move to B, B to C, C to D and D to A. This will keep the relative arrangement of components so long as the user names the groups in the correct order. ',nargs="+") 

parser.add_argument('--selection', '-s', dest='selection_string', help='Selection string for fitting. Will be applied to both target and reference structures.',type=str, default="name CA") 
parser.add_argument('--n-frames', '-fn', dest='n_vis_frames', help='Number of frames to visualise',type=int, default=30) 
parser.add_argument('--in-mem', dest='in_mem', action="store_true", help='Do alignment processing in memory. I advise against doing this if the trajectory is large.') 

args = parser.parse_args()

print(args.standardise_bool)
selection_string = args.selection_string 

universe = mda.Universe(args.topology,args.traj)
original_trajectory_n_frames = universe.trajectory.n_frames
#universe.trajectory = universe.trajectory[0::args.stride] 


selection_object = universe.select_atoms(selection_string)

if args.dir_root == None : 
    dir_root = '/tmp/' + str(uuid.uuid4())  + '/'
else:
    dir_root = args.dir_root

print(dir_root)
os.makedirs(dir_root,exist_ok=True)

def copy_coords(ag):
        return ag.positions.copy()


if args.symmetry_list != None:
    #segids_coords = [ [selection_object.select_atoms()]]
    n_atoms_set =set(([selection_object.select_atoms(current_symmetry).n_atoms for current_symmetry in args.symmetry_list]))


    if len(n_atoms_set)  > 1 :
        raise ValueError('The selections you have chosen to delineate symmetry groups results in atom groups with different numbers of atoms. Check your structure and your selection strings.')
    if list(n_atoms_set)[0]  == 0 :
        raise ValueError('Symmetry group selection has resulted in no atoms. Check selection string.')
    if selection_object.n_atoms  == 0 :
        raise ValueError('Selection object contains no atoms. Check selection object and try again.')
    if selection_object.n_atoms % list(n_atoms_set)[0] > 1 :
        raise ValueError('Number of atoms in symmetry group does not evenly divide full selection object. Check selection strings and try again.')

    temp_coords1 = [[selection_object.select_atoms(current_symmetry).positions for ts in universe.trajectory[::args.stride]] for current_symmetry in args.symmetry_list]

    segids = list(set(selection_object.segids))
    segids_indices = selection_object.atoms.segindices

    selection_object.write(dir_root + 'temp_pdb_file.pdb')

    analysis_universe = mda.Universe(dir_root + 'temp_pdb_file.pdb')

    coordinates = np.empty((len(universe.trajectory[::args.stride]) * len(args.symmetry_list), selection_object.n_atoms, 3),dtype=np.float32)

    #coordinates = np.empty((universe.trajectory.n_frames * len(args.symmetry_list), selection_object.n_atoms, 3))

    #new_universe = selection_object.load_new(coordinates, order='fac')
    #print(np.shape(np.concatenate(temp_coords1 , axis = 1)))

    coordinates [0:len(universe.trajectory[::args.stride])] = np.concatenate(temp_coords1 , axis = 1)

    for i in range(1,len(args.symmetry_list)):
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
    selection_object.write(dir_root + 'temp_pdb_file.pdb')

    analysis_universe = mda.Universe(dir_root + 'temp_pdb_file.pdb')
    coordinates = np.array([selection_object.positions for ts in universe.trajectory[::args.stride]])
     
analysis_universe = analysis_universe.load_new (coordinates)

with mda.Writer((dir_root + 'analysis_universe_traj.dcd'), analysis_universe.atoms.n_atoms) as W:
    for ts in analysis_universe.trajectory:
        W.write(analysis_universe.atoms)

if args.reference != None: 
    if args.ref_top != None:
        ref_universe = mda.Universe(args.ref_top, args.reference)
    else:
        ref_universe = mda.Universe(args.topology, args.reference)
else:
    #ref_universe = analysis_universe.copy()
    ref_universe = mda.Universe(dir_root + 'temp_pdb_file.pdb')
    ref_universe.load_new(coordinates)



analysis_universe = mda.Universe(dir_root + 'temp_pdb_file.pdb')
#load like this to over write first frame
analysis_universe.load_new(dir_root + 'analysis_universe_traj.dcd')


#aligns by default
analysis_universe.atoms.write (dir_root + 'analysis.pdb')
print('Aligning Trajectory')
if args.in_mem==True:
    #aligner = align.AlignTraj(analysis_universe, ref_universe, in_memory=True).run(stride=args.stride)
    # we have to aligned twice for some reason.... should do the alignment not in mdanalysis i  think.
    aligner = align.AlignTraj(analysis_universe, ref_universe, filename=dir_root + 'aligned.dcd',select=args.selection_string,verbose=True).run()
    analysis_universe = analysis_universe.load_new(dir_root + 'aligned.dcd')
    aligner = align.AlignTraj(analysis_universe, analysis_universe, filename=dir_root + 'aligned2.dcd',select = args.selection_string,verbose=True).run()
    os.remove(dir_root + 'aligned.dcd')
    analysis_universe = analysis_universe.load_new(dir_root + 'aligned2.dcd')
else:
    #need to fix this
    aligner = align.AlignTraj(analysis_universe, ref_universe, filename=dir_root + 'aligned.dcd',select = args.selection_string,verbose=True).run()

    analysis_universe = analysis_universe.load_new(dir_root + 'aligned.dcd')
    #aligner_universe=mda.Universe(aligner,format=MemoryReader)
    aligner = align.AlignTraj(analysis_universe, analysis_universe, filename=dir_root + 'aligned2.dcd',select = args.selection_string,verbose=True).run()
    analysis_universe = analysis_universe.load_new(dir_root + 'aligned2.dcd')

    #aligned_coords = AnalysisFromFunction(copy_coords, aligner_universe).run().results

print('Trajectory Aligned')

analysis_universe_coordinates = np.array([analysis_universe.atoms.positions for ts in analysis_universe.trajectory],dtype=np.float32)
analysis_universe_coordinates = analysis_universe_coordinates.reshape((analysis_universe.trajectory.n_frames,analysis_universe.atoms.n_atoms*3))

store_mean =  np.mean(analysis_universe_coordinates,axis=0)

demeaned_coords = analysis_universe_coordinates - store_mean
demeaned_coords = demeaned_coords.reshape(analysis_universe.trajectory.n_frames,analysis_universe.atoms.n_atoms,3) 
dist_from_mean = np.linalg.norm(demeaned_coords,axis=2)
print(np.shape(dist_from_mean))

#np.save ('mean_molecular_structure.npy', store_mean)
#store_std =  [num for num in store_std for _ in range(3)]
store_std = (np.std(analysis_universe_coordinates,axis=0))

if args.standardise_bool == True:
    print('using standardised coordinates')
    pca_ready_coords = (analysis_universe_coordinates - store_mean ) / store_std
else:
    print('not using standardised coordinates')
    pca_ready_coords = analysis_universe_coordinates

#print(np.shape(analysis_universe_coordinates))
#print(np.shape(store_mean))

#https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/aligning_trajectory.html
#analysis_universe.load_new(aligned_coords['timeseries'],format=MemoryReader)
#analysis_universe.load_new(aligned_coords['timeseries'],format=MemoryReader)

print('Calculating covariance.')


pc = PCA(n_components = args.n_components)
pc.fit(pca_ready_coords)

np.save(args.out_file + '_eigenvectors.npy', pc.components_)
np.save(args.out_file + '_eigenvalues.npy', np.sqrt(pc.explained_variance_))
np.save(args.out_file + '_mean.npy', store_mean)
np.save(args.out_file + '_std.npy', store_std)

if args.save_projection == True:
    #np.save(args.out_file + '_projected_traj.npy', pc.transform(pca_ready_coords[:original_trajectory_n_frames,:]))
    np.save(args.out_file + '_projected_traj.npy', pc.transform(pca_ready_coords))

mean_universe = mda.Merge(selection_object)

print(np.shape(pc.mean_))
mean_coords = np.reshape(store_mean,(analysis_universe.atoms.n_atoms,3))
mean_universe.load_new(mean_coords, order="fac")
mean_universe.atoms.write(args.out_file + '_mean.pdb')

#pdb.set_trace()
print('Writing output.')
for i in range(args.n_components):

    # visualise loadings
    pc_vector = pc.components_[i, :] * np.sqrt(pc.explained_variance_[i])
    #pc_vector = pc.results.p_components[:, i] 

    #pc1 = pc.results.p_components[:, 0] 
    origin = np.zeros(len(pc_vector))

    pc_traj = np.linspace(-pc_vector,pc_vector,args.n_vis_frames) * args.extend

    #trans1 = transformed[:, 0]

    #projected = np.outer([dummy_coords,dummy_coords2], pc1) + pc.mean.flatten()
    if args.standardise_bool == False:
        projected = (pc_traj + pc.mean_.flatten())  
    else:
        projected = ((pc_traj + pc.mean_.flatten()) + store_mean) 

    #pdb.set_trace()

    coordinates = projected.reshape(len(pc_traj), -1, 3)

    proj1 = mda.Merge(selection_object)

    proj1.load_new(coordinates, order="fac")

    new_sel_object = proj1.select_atoms(selection_string)

    new_sel_object.write (args.out_file + str(i) + '.pdb')

    with mda.Writer((args.out_file + str(i) + '.xtc'), new_sel_object.n_atoms) as W:
        for ts in proj1.trajectory:
            W.write(new_sel_object)
if args.dir_root == None : 
    #shutil.rmtree(dir_root)
    pass
