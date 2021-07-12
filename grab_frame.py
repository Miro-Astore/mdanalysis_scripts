import sys
import MDAnalysis as mda

# usage 
# python grab_frame.py <PSFFILE> <TRAJFILE> <frame_num> <OUTFILE.pdb>
#TODO

u=mda.Universe(str(sys.argv[1]), sys.argv[2])
atoms = u.select_atoms('all')
with mda.Writer(sys.argv[4],atoms.n_atoms) as W:
    u.trajectory[int(sys.argv[3])]
    W.write(atoms)
