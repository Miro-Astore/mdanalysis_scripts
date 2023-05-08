import MDAnalysis as mda 
import numpy as np 
import matplotlib.pyplot as plt
from MDAnalysis.analysis import  align
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--topology', '-t', dest='topology' , help='topology file',type=str) 
parser.add_argument('--trajectory', '-j', dest='traj', help='trajectory file',type=str) 
parser.add_argument('--out', '-o', dest='out_file',default="projected.npy", help='Output file',type=str) 
parser.add_argument('--stride', '-dt', dest='stride' , default=1, help='Stride through trajectory skipping this many frames.',type=int) 
parser.add_argument('--pc-eigenvectors', '-v', dest='eigvecs' , help='Eigenvectors corresponding to the PC\'s.',type=str) 
parser.add_argument('--pc-eigenvalues', '-l', dest='eigvals' , help='Eigenvalues corresponding to the PC\'s.',type=str) 
parser.add_argument('--mean', '-m', dest='mean' , help='Mean vector.',type=str) 
parser.add_argument('--align', '-a', dest='align' , action='store_true', help='Weather or not to align the trajectory before analysis.',default=False) 

parser.add_argument('--symmetry-components', dest='symmetry_list',  help='List of selections which compose a single symmetry group. This script expects symmetric groups to have the same number of atoms when the selection string is applied. For example you might have 4 identical chains so you would pass the argument "--s-components \'segid A\' \'segid B\' \'segid C\' \'segid D\'" to this script. The user is also warned that the groups will be treated cyclically. So in the previous example A will move to B, B to C, C to D and D to A. This will keep the relative arrangement of components so long as the user names the groups in the correct order. ',nargs="+") 

parser.add_argument('--selection', '-s', dest='selection_string', help='Selection string for fitting. Will be applied to both target and reference structures.',type=str, default="name CA") 

args = parser.parse_args()

universe = mda.Universe(args.topology,args.traj)
selection_subset = universe.select_atoms(args.selection_string)

pc_eigvecs =  np.load (args.eigvecs)
pc_eigvals = np.load (args.eigvals)
pc_mean = np.load (args.mean)

n_atoms_desired = np.shape(pc_eigvecs)[0]/3
n_modes = np.shape(pc_eigvecs)[1]

if n_atoms_desired != selection_subset.n_atoms :
    raise ValueError('Looks like the number of atoms you\'ve selected is not the same as was used to produce the eigenvectors you\'ve loaded. Check your selection string for both the principal component analysis and this projection script and try again.')

#aligning universe for analysis.
ref_universe = mda.Universe('/dev/shm/pca_analysis_temp.pdb')
ref_universe = ref_universe.load_new(pc_mean.reshape(selection_subset.n_atoms,3))

if args.align == True:
    aligner = align.AlignTraj(universe, ref_universe, select=args.selection_string,  filename='/dev/shm/pca_analysis.dcd',verbose=True).run()
    analysis_universe = mda.Universe(args.topology)
    analysis_universe = analysis_universe.load_new('/dev/shm/pca_analysis.dcd')
else:
    analysis_universe = universe

analysis_selection_subset = analysis_universe.select_atoms(args.selection_string)

trajectory_array = [analysis_selection_subset.positions for ts in analysis_universe.trajectory [::args.stride]]

trajectory_array = np.reshape (trajectory_array,newshape= (analysis_universe.trajectory.n_frames,analysis_selection_subset.n_atoms*3))
selection_subset.write('/dev/shm/pca_analysis_temp.pdb')

project = np.zeros ([len(pc_eigvals),universe.trajectory.n_frames])

for i in range(n_modes):
    eigenvector = pc_eigvecs[:,i]
    eigenvalue = pc_eigvals[i]
    print(pc_mean)

    project[i] = [np.dot(curr_frame-pc_mean,eigenvector) for curr_frame in trajectory_array]

np.save(args.out_file,project)

if len(pc_eigvals) == 2:
    print('plotting')
    plt.scatter(project[0],project[1],linewidths=0.2)
    plt.show()
    plt.savefig('fig.pdf')

