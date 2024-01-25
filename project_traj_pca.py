import numpy as np

import MDAnalysis as mda 
from MDAnalysis.analysis import pca, align
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction
from sklearn.decomposition import PCA
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
parser.add_argument('--input', '-i', dest='input_file',default="pca", help='Output file',type=str) 
parser.add_argument('--out', '-o', dest='out_file',default="pc_projected_", help='Output file',type=str) 
parser.add_argument('--stride', '-dt', dest='stride' , default=1, help='Stride through trajectory skipping this many frames.',type=int) 
parser.add_argument('--no-std', dest='standardise_bool', action="store_false", help='Choose whether or not to standardise coordinates. Default behaviour is to standardise the cordinates.') 
parser.add_argument('--n-components', '-n', dest='n_components' , default = 2, help='calculate this many Principal Components of the trajectory.',type=int) 

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
#universe.trajectory = universe.trajectory[0::args.stride] 

selection_object = universe.select_atoms(selection_string)

eigenvectors = np.load(args.input_file + '_eigenvectors.npy')
eigenvalues = np.load(args.input_file + '_eigenvalues.npy')

n_loaded_components = len(eigenvalues)
if args.n_components > n_loaded_components:
    raise ValueError('the number of loaded components is smaller than the number of components you have asked to analyze. You must either run a deeper analysis or analyse fewer coordinates in this step..')
elif args.n_components < n_loaded_components :
    eigenvectors = eigenvectors [:args.n_components]
    eigenvalues = eigenvectors [:args.n_components]

store_mean = np.load(args.input_file + '_mean.npy')
store_std = np.load(args.input_file + '_std.npy')

mean_universe = mda.Merge(selection_object)

mean_coords = np.reshape(store_mean,(mean_universe.atoms.n_atoms,3))

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

    selection_object.write('proj/temp_pdb_file.pdb')

    analysis_universe = mda.Universe('proj/temp_pdb_file.pdb')

    coordinates = np.empty((len(universe.trajectory[::args.stride]) * len(args.symmetry_list), selection_object.n_atoms, 3),dtype=np.float32)

    #coordinates = np.empty((universe.trajectory.n_frames * len(args.symmetry_list), selection_object.n_atoms, 3))

    #new_universe = selection_object.load_new(coordinates, order='fac')
    #print(np.shape(np.concatenate(temp_coords1 , axis = 1)))

    coordinates [0:len(universe.trajectory[::args.stride])] = np.concatenate(temp_coords1 , axis = 1)

    for i in range(1,len(args.symmetry_list)):
        temp_coords1 = np.concatenate (([temp_coords1[-1]],temp_coords1[0:-1]), axis=0)
        coordinates[i*len(universe.trajectory[::args.stride]):(i+1)*len(universe.trajectory[::args.stride])] =  np.concatenate(temp_coords1 , axis = 1)

else:
    selection_object.write('proj/temp_pdb_file.pdb')

    analysis_universe = mda.Universe('proj/temp_pdb_file.pdb')
    coordinates = np.array([selection_object.positions for ts in universe.trajectory[::args.stride]])
     
analysis_universe = analysis_universe.load_new (coordinates)

with mda.Writer(('proj/analysis_universe_traj.dcd'), analysis_universe.atoms.n_atoms) as W:
    for ts in analysis_universe.trajectory:
        W.write(analysis_universe.atoms)

if args.reference != None: 
    print('loading reference')
    if args.ref_top != None:
        ref_universe = mda.Universe(args.ref_top, args.reference)
    else:
        ref_universe = mda.Universe(args.topology, args.reference)
else:
    #ref_universe = analysis_universe.copy()
    ref_universe = mda.Universe('proj/temp_pdb_file.pdb')
    ref_universe.load_new(mean_coords)



analysis_universe = mda.Universe('proj/temp_pdb_file.pdb')
#load like this to over write first frame
analysis_universe.load_new('proj/analysis_universe_traj.dcd')


#aligns by default
analysis_universe.atoms.write ('proj/analysis.pdb')
print('Aligning Trajectory')

if args.in_mem==True:
    print('in mem')
    #aligner = align.AlignTraj(analysis_universe, ref_universe, in_memory=True).run(stride=args.stride)
    # we have to aligned twice for some reason.... should do the alignment not in mdanalysis i  think.
    aligner = align.AlignTraj(analysis_universe, ref_universe, filename='proj/aligned.dcd',select=args.selection_string).run()
    analysis_universe = analysis_universe.load_new('proj/aligned.dcd')
    aligner = align.AlignTraj(analysis_universe, analysis_universe, filename='proj/aligned2.dcd',select = args.selection_string).run()
    os.remove('proj/aligned.dcd')
    analysis_universe = analysis_universe.load_new('proj/aligned2.dcd')
else:
    #need to fix this
    aligner = align.AlignTraj(analysis_universe, ref_universe, filename='proj/aligned.dcd',select = args.selection_string).run()

    analysis_universe = analysis_universe.load_new('proj/aligned.dcd')
    #aligner_universe=mda.Universe(aligner,format=MemoryReader)
    aligner = align.AlignTraj(analysis_universe, analysis_universe, filename='proj/aligned2.dcd',select = args.selection_string).run()
    analysis_universe = analysis_universe.load_new('proj/aligned2.dcd')

    #aligned_coords = AnalysisFromFunction(copy_coords, aligner_universe).run().results

print('Trajectory Aligned')

analysis_universe_coordinates = np.array([analysis_universe.atoms.positions for ts in analysis_universe.trajectory],dtype=np.float32)
analysis_universe_coordinates = analysis_universe_coordinates.reshape((analysis_universe.trajectory.n_frames,analysis_universe.atoms.n_atoms*3))

if args.standardise_bool == True:
    print('using standardised coordinates')
    projection_ready_coords = (analysis_universe_coordinates - store_mean ) / store_std
else:
    print('not using standardised coordinates')
    projection_ready_coords = analysis_universe_coordinates



projected_arr = np.zeros([analysis_universe.trajectory.n_frames,args.n_components])
print(np.shape(eigenvectors))
print(np.shape(projection_ready_coords))


#eigenvectors =  np.dot(np.diag(eigenvalues), eigenvectors ) 
projected_arr =  projection_ready_coords @eigenvectors.T
print(np.shape(projected_arr))

np.save(args.out_file + '.npy',projected_arr)
