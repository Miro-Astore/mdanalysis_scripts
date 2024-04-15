import MDAnalysis as mda
from MDAnalysis import transformations
import sys

#error handling
if len(sys.argv) != 4:
    raise Exception('wrong number of arguments. Need 3')
u = mda.Universe(sys.argv[1], sys.argv[2])
prot = u.select_atoms("protein")
# we load another universe to define the reference
# it uses the same input files, but this doesn't have to be always the case
ref_u = u.copy()
reference = ref_u.select_atoms("protein")
ag = u.atoms
workflow = (transformations.unwrap(ag), transformations.center_in_box(prot, center='mass'), transformations.wrap(ag, compound='fragments'))
u.trajectory.add_transformations(*workflow)

all_as=u.select_atoms('all')
with mda.Writer(sys.argv[3], all_as.n_atoms) as W:
    for ts in u.trajectory:
        W.write(all_as)
