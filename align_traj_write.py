import MDAnalysis
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

#from MDAnalysis.tests.datafiles import PDB, XTC

PDB = 'bpti_prot.pdb'
XTC = './bpti-protein/sum.dcd'

ref = MDAnalysis.Universe(PDB)
u = MDAnalysis.Universe(PDB, XTC)
protein = u.select_atoms("protein")

alignment = align.AlignTraj(u, ref, filename='rmsfit.dcd')
alignment.run()

#with MDAnalysis.Writer("protein.xtc", protein.n_atoms) as W:
#
#    for ts in u.trajectory:
#        W.write(protein)


