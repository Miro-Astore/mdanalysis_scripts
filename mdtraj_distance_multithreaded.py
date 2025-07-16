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


def mdtraj_distance(filename, indices):
    structure = mdtraj.load(filename)
    distance = mdtraj.compute_distances(structure, [indices]) * 10  # Convert to Ångstroms
    return  distance[0][0]  # Flatten to a single value


def process_file(filepath_and_indices):
    filepath, indices = filepath_and_indices
    try:
        distance = mdtraj_distance(filepath, indices)
        return os.path.basename(filepath), distance
    except Exception as e:
        return os.path.basename(filepath), f"ERROR: {str(e)}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', dest='directory', required=True)
    parser.add_argument('--output', dest='output', default="distances.csv")
    args = parser.parse_args()

    pdb_files = quickly_list_pdb_files(args.directory)

    # Load first structure to define indices
    first_structure = mdtraj.load(os.path.join(args.directory, pdb_files[0]))
    indices = first_structure.topology.select('name P and (resid 64 or resid 2)')
    if len(indices) != 2:
        raise ValueError(f"Expected 2 indices, got {len(indices)}: {indices}")

    # Build input for multiprocessing
    filepaths = [os.path.abspath(os.path.join(args.directory, f)) for f in pdb_files]
    inputs = [(fp, indices) for fp in filepaths]

    results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for result in tqdm(executor.map(process_file, inputs), total=len(inputs), desc="Computing distances"):
            results.append(result)

    print(results)
    # Save to CSV
    with open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["filename", "distance_angstroms"])
        writer.writerows(results)

    print(f"\n Distances saved to {args.output}")

if __name__ == "__main__":
    main()

