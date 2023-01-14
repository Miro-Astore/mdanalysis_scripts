import pdb
import numpy as np
import MDAnalysis as mda

top_file = 'window_40/ionized.psf'
traj_files = ['window_40/memb_prod1.xtc','window_41/memb_prod1.xtc']

u = mda.Universe(top_file)
u = mda.Universe(top_file)

for i in traj_files:
    u.load_new(i)

chloride_index = str(int('154659' ) - 1)
chloride_sel_string = 'index ' + chloride_index
chloride_sel = u.select_atoms (chloride_sel_string)

body_indices = np.loadtxt('indices.txt',dtype=int)

temp_str = ''  
for i in range(len(body_indices)):
    temp_num = str(body_indices [i] - 1)
    temp_str = temp_str + (temp_num) + ' '

body_indices=temp_str
body_sel_string = 'index ' + str(body_indices)

body_sel = u.select_atoms(body_sel_string)

z_pos = np.zeros (len(u.trajectory))
hydration_num = np.zeros (len(u.trajectory))


for i in range(len(u.trajectory)):
    u.trajectory[i]
    body_com = body_sel.center_of_mass()[2]
    chloride_com = chloride_sel.center_of_mass()[2]

    hydration_sel  =  u.select_atoms('resname SOL and around 3.8 (index  ' + chloride_index  +  ')' )
    resids = set(hydration_sel.resids)
    hydration_num[i] = len(resids)

    z_pos [i] =  np.abs(body_com - chloride_com)

np.save('ion_z_position.npy',z_pos)
np.save('chloride_coodination_number.npy',hydration_num)
