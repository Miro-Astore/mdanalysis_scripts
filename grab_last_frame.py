import sys
import MDAnalysis as mda

# usage 
# python grab_last_frame.py <PSFFILE> <TRAJFILE> <OUTFILE.pdb>
#last argument optional, default output last.pdb
#TODO

u=mda.Universe(str(sys.argv[1]), sys.argv[2])
atoms = u.select_atoms('all')
with mda.Writer(sys.argv[3],atoms.n_atoms) as W:
    u.trajectory[-1]
    W.write(atoms)
