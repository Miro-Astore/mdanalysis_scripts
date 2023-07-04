import MDAnalysis as mda
import sys
from MDAnalysis.analysis import align
from multiprocessing import Pool, cpu_count
import itertools
import numpy as np

def compute_rmsd_pair(args):
    topology_file, trajectory_file, i, j = args
    u1 = mda.Universe(topology_file, trajectory_file)
    u2 = mda.Universe(topology_file, trajectory_file)
    
    u1.trajectory[i]
    u2.trajectory[j]
    
    rmsd = align.alignto(u1.atoms, u2.atoms, select='protein')[1]
    return (i, j, rmsd)

def compute_rmsd_matrix(topology_file, trajectory_file):
    u = mda.Universe(topology_file, trajectory_file)

    num_frames = len(u.trajectory)
    rmsd_matrix = np.zeros((num_frames, num_frames))

    # List of arguments for each compute_rmsd_pair task
    args = [(topology_file, trajectory_file, i, j) for i in range(num_frames) for j in range(i, num_frames)]

    # Compute RMSDs in parallel
    with Pool(cpu_count()) as pool:
        results = pool.map(compute_rmsd_pair, args)

    # Update RMSD matrix with results
    for i, j, rmsd in results:
        rmsd_matrix[i, j] = rmsd
        rmsd_matrix[j, i] = rmsd

    return rmsd_matrix

if __name__ == '__main__':
    # Use the function to compute the RMSD matrix
    # Replace 'topology_file.pdb' and 'trajectory_file.xtc' with your actual files
    rmsd_matrix = compute_rmsd_matrix(sys.argv[1], sys.argv[2])
    np.savetxt(rmsd_matrix,sys.argv[3])

