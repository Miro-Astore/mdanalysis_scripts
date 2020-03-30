import sys
import numpy as np 
import MDAnalysis as mda

# usage 
# python grab_last_frame.py <PSFFILE> <TRAJFILE> frame_start number_of_snapshots
#last argument optional, default output last.pdb 
#TODO

u=mda.Universe(str(sys.argv[1]), sys.argv[2])
frame_start=int(sys.argv[4])
frames_left = u.trajectory.n_frames -  int(sys.argv[3]) - 1 
frame_slice = int(frames_left/int(sys.argv[4]))

frame_indexes = np.zeros(int(sys.argv[4]))

for i in range(int(sys.argv[4])):
    frame_indexes[i] = int(frame_start + frame_slice * i)

file_name_root=str(sys.argv[2])[:-4]
print(frame_indexes)

for i in range(int(sys.argv[4])):
    u.trajectory[int(frame_indexes[i])]
    atoms = u.select_atoms('protein')

    with mda.Writer(file_name_root + '_snaps_' + str(i) + '.pdb', atoms.n_atoms) as W:
        W.write(atoms)
