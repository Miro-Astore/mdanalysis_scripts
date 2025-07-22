import natsort
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


def mdtraj_angle(structure, bottom_stem_indices,top_stem_indices):
    
    bs_points_mr = structure.xyz[bottom_stem_indices] - np.mean(structure.xyz[bottom_stem_indices],axis=0)
    ts_points_mr = structure.xyz[top_stem_indices] - np.mean(structure.xyz[top_stem_indices],axis=0)

    bs_u,bs_s,bs_vT = np.linalg.svd(bs_points_mr.ravel())
    ts_u,ts_s,ts_vT = np.linalg.svd(ts_points_mr.ravel())

    bs_major_axis=bs_vT[0]
    ts_major_axis=ts_vT[0]

    # check if bs major axis is pointing inward compared to z axis.
    # check that it's displacement vector matches the sign of the major axis
    if np.dot(ts_major_axis,structure.xyz[top_stem_indices]) < 0 : 
        ts_major_axis = ts_major_axis * -1 
    
    if np.dot(bs_major_axis,structure.xyz[bottom_stem_indices]) < 0 : 
        bs_major_axis = bs_major_axis * -1 

    return  np.arcos(np.dot(bs_major_axis,ts_major_axis))  

def mdtraj_distance(structure, indices):
    distance = mdtraj.compute_distances(structure, [indices]) * 10  # Convert to Ångstroms
    return  distance[0][0]  # Flatten to a single value

def mdtraj_cvs (filepath, fitting_indices, bs_indices, ts_indices, distance_indices, state1, state2, bs_indices_ref, reference_inds):
    structure = mdtraj.load(filename)

    rmsd_bs = mdtraj.rmsd(structure,state1,atom_indices=bs_indices, ref_atom_indices=bs_indices_ref)*10
    rmsd_bs = rmsd_bs [0]
    distance = mdtraj_distance(structure, distance_indices)

    angle = mdtraj_angle (structure, bs_indices, ts_indices)  

    rmsd_state1 = mdtraj.rmsd(structure,state1,atom_indices=fitting_indices, ref_atom_indices=reference_inds)*10
    rmsd_state1 = rmsd_state1[0]
    rmsd_state2 = mdtraj.rmsd(structure,state2,atom_indices=fitting_indices, ref_atom_indices=reference_inds)*10
    rmsd_state2 = rmsd_state2[0]

    return  os.path.baname(filepath), rmsd_bs, angle, distance, rmsd_state1, rmsd_state2

def process_file(filepath_and_indices):
    filepath, fitting_indices, bs_indices, ts_indices, distance_indices, state1, state2, bs_indices_ref, reference_inds =  filepath_and_indices

    structure = mdtraj.load(filepath)
    try:
        rmsd_bs, angle, distance, rmsd_state1, rmsd_state2 =  md_traj_cvs(filepath, fitting_indices, bs_indices, ts_indices, distance_indices, state1, state2, bs_indices_ref, reference_inds) 
        return os.path.basename(structure), rmsd_bs, angle, distance, rmsd_state1, rmsd_state2
    except Exception as e:
        return os.path.basename(filepath), f"ERROR: {str(e)}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', dest='directory', required=True)
    parser.add_argument('--output', dest='output', default="distances.csv")
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
    reference_inds = state1.topology.select('name P and resid 2 to 69')
    bs_indices_ref = state1.topology.select('name P and (resid 2 to 16 or resid 56 to 70)')

    distance_indices = first_structure.topology.select('name C5 and (resid 65 or resid 34)')
    bs_indices = first_structure.topology.select('name P and (resid 2 to 16 or resid 56 to 70)')

    ts_indices = first_structure.topology.select('name P and (resid 20 to 43)')
    fitting_indices = first_structure.topology.select('name P and resid 2 to 69')

    if len(distance_indices) != 2:
        raise ValueError(f"Expected 2 indices, got {len(indices)}: {indices}")

    # Build input for multiprocessing
    filepaths = [os.path.abspath(os.path.join(args.directory, f)) for f in pdb_files]
    input_tuple = (fitting_indices, bs_indices, ts_indices, distance_indices, state1, state2, bs_indices_ref, reference_inds)

    inputs = [(fp,) + input_tuple  for fp in filepaths]

    results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for result in tqdm(executor.map(process_file, inputs), total=len(inputs), desc="Computing distances"):
            results.append(result)

    print(results)
    # Save to CSV
    with open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["filename", "bottom_stem_rsmd", "bs_ts_angle", "U34-U65 C5 distance", "state 1 RMSD", "state 2 RMSD"])
        writer.writerows(results)
    print(f"\n Distances saved to {args.output}")

if __name__ == "__main__":
    main()

