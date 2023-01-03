from MDAnalysis.analysis.distances import dist
import MDAnalysis as mda
import MDAnalysis.analysis.rms

import numpy as np 
import sys
#####
#USAGE: python mdanalysis_scripts/distance.py TOPFILE traj sel1 sel2 out_name
#
#TODO better detection of lysine bridges, trying to fix argenine detection
#####
print ('TOP ' + str(sys.argv[1]))
print ('traj ' + str(sys.argv[2]))
print ('sel1 ' + str(sys.argv[3]))
print ('sel2 ' + str(sys.argv[4]))
print ('out_name ' + str(sys.argv[5]))
u=mda.Universe(sys.argv[1],sys.argv[2])

file_out_name = str(sys.argv[5])

sel1txt=str(sys.argv[3])
sel2txt=str(sys.argv[4])

sel1=u.select_atoms(sel1txt)
sel2=u.select_atoms(sel2txt)


dist_arr=np.zeros([u.trajectory.n_frames,2])

i=0
for ts in u.trajectory:
   sel1_com = sel1.center_of_mass()
   sel2_com = sel2.center_of_mass()

   dist_arr[i][1] =  np.linalg.norm (sel1_com-sel2_com)
   dist_arr[i][0] =  ts.time*0.001
   i=i+1

print (dist_arr)
np.savetxt(file_out_name,dist_arr)
