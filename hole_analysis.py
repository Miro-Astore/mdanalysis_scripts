import MDAnalysis as mda
import numpy as np 
from MDAnalysis.analysis.hole import HOLEtraj
from MDAnalysis.transformations.translate import center_in_box
from MDAnalysis.transformations.rotate import rotateby
from MDAnalysis.transformations.rotate import rotateby
import MDAnalysis.transformations 
import MDAnalysis as mda
from scipy.spatial.transform import Rotation as R

def align_to_principal_axes(ag,ts):
    """
    align principal components of the system to the xyz axes 
    ts is a timestep argument and ag is the atom group we wish to use as the    to detect the principal components
    TODO add wrapping.
    """
    #finding origin of rotation
    origin=ag.center_of_mass()
    p_axes=ag.principal_axes()

    x=[1,0,0]
    y=[0,1,0]
    z=[0,0,1]

    #only two rotations are required because if two components are aligned the 3rd one has to be since it's all orthogonal in 3 dimenions
    rot1_axis=np.cross(p_axes[2],z)
    rot2_axis=np.cross(p_axes[1],y)
    rot1_angle=np.degrees(np.arccos(np.dot(p_axes[0],x)))
    rot2_angle=np.degrees(np.arccos(np.dot(p_axes[1],y)))
    rotated=MDAnalysis.transformations.rotate.rotateby(rot1_angle, direction=rot1_axis,point=origin)(ts)
    new_ts=MDAnalysis.transformations.rotate.rotateby(rot2_angle, direction=rot2_axis,point=origin)(rotated)

    return new_ts 
##finding origin ,
import functools
u=mda.Universe('./1unp.pdb')
#pore_sel=u.select_atoms("resid 95 104 134 190 248 334 335 352 370 and name CA")
pore_sel=u.select_atoms("all")
print (pore_sel.principal_axes())
#workflow=align_to_principal_axes(ag=pore_sel)
u.atoms.align_principal_axis(0,[0,0,1])
#u.atoms.align_principal_axis(1,[0,1,0])
u.atoms.write('aligned.pdb')

#align_to_principal_axes = functools.partial(align_to_principal_axes,pore_sel)
#workflow=align_to_principal_axes
#align_to_principal_axes = functools.partial(align_to_principal_axes,ag=pore_sel)
#u.trajectory.add_transformations(workflow)

atoms = u.select_atoms('all')
with mda.Writer('out.pdb',atoms.n_atoms) as W:
    u.trajectory[-1]
    W.write(atoms)

# selecting pore entry residues


pore_entry_sel=u.select_atoms("resid 95 104 134 190 248 334 335 352 370 and name CA")

pore_start=pore_entry_sel.center_of_mass()

##demeaning 
#for ts in u.trajectory:
##    new_ts = mda.transformations.center_in_box(pore_sel,point=[0,0,0])(ts)
#    u.ts = mda.transformations.center_in_box(pore_sel,point=[0,0,0])(ts)
 
    
#rotating

H = HOLEtraj(u, cpoint=pore_start, cvect =[ 0 , 0 , 1],  executable="~/md/hole2/exe/hole")  # set path to your hole binary 
#H = HOLE('./6msm_prot.pdb', executable="~/md/hole2/exe/hole")  # set path to your hole binary 
H.run()
#H.collect()
H.plot(linewidth=3, color="black", label=False)
