import natsort
import numpy as np 
import mdtraj
from pathlib import Path
import argparse
import os
import concurrent.futures
from tqdm import tqdm
import csv

def quickly_list_pdb_files(path):
    file_list = []
    with os.scandir(path) as entries:
        for entry in entries:
            if entry.is_file() and entry.name.endswith(".pdb") and entry.name.startswith("S_"):
                file_list.append(entry.name)
    file_list = natsort.natsorted(file_list)
    return file_list

def mdtraj_bs_rmsd(structure,reference_structure, bottom_stem_indices):
    curr_rmsd = mdtraj.rmsd(structure,reference_structure)
    return curr_rmsd


def mdtraj_planarity(structure, bs_base_base_pairs_inds, ts_base_base_pairs_inds, ss_base_base_pairs_inds):
    bs_points_mr = structure.xyz[0][bs_base_base_pairs_inds] - np.mean(structure.xyz[0][bs_base_base_pairs_inds],axis=0)
    ts_points_mr = structure.xyz[0][ts_base_base_pairs_inds] - np.mean(structure.xyz[0][ts_base_base_pairs_inds],axis=0)
    ss_points_mr = structure.xyz[0][ss_base_base_pairs_inds] - np.mean(structure.xyz[0][ss_base_base_pairs_inds],axis=0)

    bs_u,bs_s,bs_vT = np.linalg.svd(bs_points_mr)
    ts_u,ts_s,ts_vT = np.linalg.svd(ts_points_mr)
    ss_u,ss_s,ss_vT = np.linalg.svd(ss_points_mr)

    bs_normal_axis=bs_vT[-1]
    ts_normal_axis=ts_vT[-1]
    ss_normal_axis=ss_vT[-1]

    # check if bs major axis is pointing inward compared to z axis.
    # check that it's displacement vector matches the sign of the major axis
    #if np.dot(ts_major_axis,structure.xyz[0][top_stem_indices]) < 0 : 
    #    ts_major_axis = ts_major_axis * -1 
    #
    #if np.dot(bs_major_axis,structure.xyz[bottom_stem_indices]) < 0 : 
    #    bs_major_axis = bs_major_axis * -1 

    bs_ts_angle = np.rad2deg(np.arccos(np.dot(bs_normal_axis,ts_normal_axis))) 
    ss_bs_angle = np.rad2deg(np.arccos(np.dot(bs_normal_axis,ss_normal_axis)))
    ss_ts_angle = np.rad2deg(np.arccos(np.dot(ss_normal_axis,ts_normal_axis)))
    return (bs_ts_angle, ss_bs_angle, ss_ts_angle)

def mdtraj_angle(structure, bottom_stem_indices, top_stem_indices, side_stem_indices):
    bs_points_mr = structure.xyz[0][bottom_stem_indices] - np.mean(structure.xyz[0][bottom_stem_indices],axis=0)
    ts_points_mr = structure.xyz[0][top_stem_indices] - np.mean(structure.xyz[0][top_stem_indices],axis=0)
    ss_points_mr = structure.xyz[0][side_stem_indices] - np.mean(structure.xyz[0][side_stem_indices],axis=0)

    
    bs_u,bs_s,bs_vT = np.linalg.svd(bs_points_mr)
    ts_u,ts_s,ts_vT = np.linalg.svd(ts_points_mr)
    ss_u,ss_s,ss_vT = np.linalg.svd(ss_points_mr)

    bs_major_axis=bs_vT[0]
    ts_major_axis=ts_vT[0]
    ss_major_axis=ss_vT[0]

    # check if bs major axis is pointing inward compared to z axis.
    # check that it's displacement vector matches the sign of the major axis
    #if np.dot(ts_major_axis,structure.xyz[0][top_stem_indices]) < 0 : 
    #    ts_major_axis = ts_major_axis * -1 
    #
    #if np.dot(bs_major_axis,structure.xyz[bottom_stem_indices]) < 0 : 
    #    bs_major_axis = bs_major_axis * -1 

    bs_ts_angle = np.rad2deg(np.arccos(np.dot(bs_major_axis,ts_major_axis))) 
    ss_bs_angle = np.rad2deg(np.arccos(np.dot(bs_major_axis,ss_major_axis)))
    ss_ts_angle = np.rad2deg(np.arccos(np.dot(ss_major_axis,ts_major_axis)))
    return (bs_ts_angle, ss_bs_angle, ss_ts_angle)

def mdtraj_distance(structure, indices):
    distance = mdtraj.compute_distances(structure, [indices]) * 10  # Convert to Ångstroms
    return  distance[0][0]  # Flatten to a single value

def mdtraj_cvs (filepath, fitting_indices, bs_indices, ts_indices, ss_indices, bs_base_base_inds,  ss_base_base_inds, ts_base_base_inds, distance_indices, state1, state2, bs_indices_ref, reference_inds):
    structure = mdtraj.load(filepath)

    rmsd_bs = mdtraj.rmsd(structure,state1,atom_indices=bs_indices, ref_atom_indices=bs_indices_ref)*10
    rmsd_bs = rmsd_bs [0]
    distance = mdtraj_distance(structure, distance_indices)

    stem_angles = mdtraj_angle (structure, bs_indices, ts_indices, ss_indices)  
    base_angles = mdtraj_planarity (structure, bs_base_base_inds, ts_base_base_inds, ss_base_base_inds)  

    rmsd_state1 = mdtraj.rmsd(structure,state1,atom_indices=fitting_indices, ref_atom_indices=reference_inds)*10
    rmsd_state1 = rmsd_state1[0]
    rmsd_state2 = mdtraj.rmsd(structure,state2,atom_indices=fitting_indices, ref_atom_indices=reference_inds)*10
    rmsd_state2 = rmsd_state2[0]

    #print(os.path.basename(filepath), rmsd_bs, angle, distance, rmsd_state1, rmsd_state2)

    return  os.path.basename(filepath), rmsd_bs, stem_angles, base_angles, distance, rmsd_state1, rmsd_state2

def process_file(filepath_and_indices):
    filepath, fitting_indices, bs_indices, ts_indices, ss_indices, bs_base_base_inds, ts_base_base_inds, ss_base_base_inds, distance_indices, state1, state2, bs_indices_ref, reference_inds =  filepath_and_indices

    structure = mdtraj.load(filepath)
    #try:
    results_tuple =  mdtraj_cvs(filepath, fitting_indices, bs_indices, ts_indices, ss_indices, bs_base_base_inds, ts_base_base_inds, ss_base_base_inds, distance_indices, state1, state2, bs_indices_ref, reference_inds) 
    return results_tuple
    #except Exception as e:
    #    return os.path.basename(filepath), f"ERROR: {str(e)}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', dest='directory', required=True)
    parser.add_argument('--output', dest='output', default="cvs.csv")
    parser.add_argument('--reference', dest='reference_structure', default=None)
    args = parser.parse_args()

    pdb_files = quickly_list_pdb_files(args.directory)

    # Load first structure to define indices
    if args.reference_structure == None:
        reference_structure = mdtraj.load(os.path.join(args.directory, pdb_files[0]))
    else:
        reference_structure = mdtraj.load(os.path.join(args.reference_structure))

    state1 = mdtraj.load('/mnt/ceph/users/mastore/sla_lorena_2025/source_files/renumbered_state1_SLA.pdb')
    state2 = mdtraj.load('/mnt/ceph/users/mastore/sla_lorena_2025/source_files/renumbered_state2_SLA.pdb')

    first_structure = mdtraj.load(os.path.join(args.directory,pdb_files[0]))
    reference_inds = state1.topology.select('name P and residue 2 to 69')

    distance_indices = first_structure.topology.select('name C5 and (residue 65 or residue 34)')
    bs_indices_ref =      state1.topology.select('name P and (residue 2 to 15 or residue 56 to 66)')
    bs_indices = first_structure.topology.select('name P and (residue 2 to 15 or residue 56 to 66)')
    bs_base_base_inds = first_structure.topology.select("(residue 15 56) and not element 'H.*'")

    ts_base_base_inds = first_structure.topology.select("(residue 20 43) and not element 'H.*' ")
    ss_base_base_inds = first_structure.topology.select("(residue 44 54) and not element 'H.*' ")

    ss_indices = first_structure.topology.select('name P and (residue 44 to 54)')

    ts_indices = first_structure.topology.select('name P and (residue 20 to 43)')
    fitting_indices = first_structure.topology.select('residue 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 and name P ')

    if len(distance_indices) != 2:
        raise ValueError(f"Expected 2 indices, got {len(indices)}: {indices}")

    # Build input for multiprocessing
    filepaths = [os.path.abspath(os.path.join(args.directory, f)) for f in pdb_files]
    input_tuple = (fitting_indices,bs_indices, ts_indices, ss_indices, bs_base_base_inds, ts_base_base_inds, ss_base_base_inds, distance_indices, state1, state2, bs_indices_ref, reference_inds)

    inputs = [(fp,) + input_tuple  for fp in filepaths]

    results = []
    print(process_file(inputs[0]))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        for result in tqdm(executor.map(process_file, inputs), total=len(inputs), desc="Computing distances"):
            results.append(result)

    print(results)
    # Save to CSV
    with open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["filename", "bottom_stem_rsmd", "stem_loop_angles", "base_pair_angles", "U34-U65 C5 distance", "state 1 RMSD", "state 2 RMSD"])
        writer.writerows(results)
    print(f"\n Distances saved to {args.output}")

if __name__ == "__main__":
    main()

