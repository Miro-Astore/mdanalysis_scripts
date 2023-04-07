import sys
import argparse
import MDAnalysis as mda

# usage 
# python write_frames.py <PSFFILE> <TRAJFILE> <selection_string>
#TODO



u=mda.Universe(str(sys.argv[1]), sys.argv[2])
file_root = str(sys.argv[2])[:-4]
atoms = u.select_atoms(str(sys.argv[3]))
frame_num = 0 
for ts in u.trajectory:
    with mda.Writer(file_root + str(frame_num) + '.pdb',atoms.n_atoms) as W:
        W.write(atoms)
        frame_num = frame_num + 1
