#usage:
# python mdanalysis_scripts/rmsd.py top_file traj_file ref_struct out_name

import MDAnalysis as mda
import MDAnalysis.analysis.rms
from MDAnalysis.analysis import rms, align

import numpy as np 
import sys

print ('psf ' + str(sys.argv[1]))
print ('traj ' + str(sys.argv[2]))
print ('refpdb ' + str(sys.argv[3]))
print ('sel_text ' + str(sys.argv[4]))
print ('out_name ' + str(sys.argv[5]))

select_string = str(sys.argv[4])
u = mda.Universe(sys.argv[3], sys.argv[2])
measure_sel = u.select_atoms (select_string)
ref = mda.Universe(sys.argv[3])
print(u.atoms.names)

aligner = align.AlignTraj(u, ref, select='name CA', in_memory=True, match_atoms=False).run()


R = MDAnalysis.analysis.rms.RMSF(measure_sel)

R.run()

np.savetxt(str(sys.argv[5]),R.rmsf)
