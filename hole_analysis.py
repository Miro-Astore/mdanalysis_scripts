#####
#usage: python mdanalysis_scripts/hole_analysis.py PDB.pdb 
# will output a file with the name PDB.sph 
# TODO add processing to output .sph file 
#####

import MDAnalysis as mda
import sys
import os 
import numpy as np 
from MDAnalysis.analysis.hole import HOLEtraj
from MDAnalysis.transformations.translate import center_in_box
from MDAnalysis.transformations.rotate import rotateby
import MDAnalysis.transformations 
import MDAnalysis as mda
from scipy.spatial.transform import Rotation as R
##finding origin ,
import functools
u=mda.Universe(sys.argv[-1])
#workflow=align_to_principal_axes(ag=pore_sel)

#align_to_principal_axes = functools.partial(align_to_principal_axes,pore_sel)
#workflow=align_to_principal_axes
#align_to_principal_axes = functools.partial(align_to_principal_axes,ag=pore_sel)
#u.trajectory.add_transformations(workflow)

#u.atoms.align_principal_axis(2,[0,0,1])
#u.atoms.align_principal_axis(1,[0,1,0])

#pore_entry_sel=u.select_atoms("resid 187 249 and name CA")
#pore_start=pore_entry_sel.center_of_mass()
#print (pore_start)
#orr=u.atoms.center_of_mass()
#print (orr)

atoms = u.select_atoms('all')
#with mda.Writer('out.pdb',atoms.n_atoms) as W:
#    u.trajectory[-1]
#    W.write(atoms)

# selecting pore entry residues



##demeaning 
#for ts in u.trajectory:
##    new_ts = mda.transformations.center_in_box(pore_sel,point=[0,0,0])(ts)
#    u.ts = mda.transformations.center_in_box(pore_sel,point=[0,0,0])(ts)
 
    
#rotating

#u2=mda.Universe('.pdb')
#annoying wt sel
#pore_entry_sel=u.select_atoms("resid  186 190 363 241 and name CA")
pore_entry_sel=u.select_atoms("resid 993 352")
#pore_entry_sel=u.select_atoms("resid 334 336 352 and name CA")

pore_start=pore_entry_sel.center_of_mass() 
#pore_start=  [160.598007, 158.057007,  163.794006]

print (pore_start)
sphpdb=str(sys.argv[-1])[:-4] + '.sph'

print (sphpdb)

H = HOLEtraj(u, cpoint=pore_start, cvect = [ 0 , 0 , 1],  executable="~/md/hole2/exe/hole")  # set path to your hole binary 

#H = HOLE('./6msm_prot.pdb', executable="~/md/hole2/exe/hole")  # set path to your hole binary 
H.run()

os.rename('hole.sph',sphpdb)
#H.collect()
H.plot(linewidth=3, color="black", label=False)
