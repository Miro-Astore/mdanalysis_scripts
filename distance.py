import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis.analysis.distances as dist
import sys
import multiprocessing

PSF=sys.argv[1]
TRJ=sys.argv[2]
u=mda.Universe(PSF,TRJ)
sel1_str=str(sys.argv[3])
sel2_str=str(sys.argv[4])
sel1=u.select_atoms(sel1_str)
sel2=u.select_atoms(sel2_str)
times=np.zeros(len(u.trajectory))
dist_arr=np.zeros(len(u.trajectory))
for i in range(u.trajectory.n_frames):
    u.trajectory[i]
    times[i]=u.trajectory.time*0.001
    dist_arr[i]=dist.dist(sel1,sel2)[-1]

paired=[times,dist_arr]
numpy_name=str(sys.argv[5])[:-4] + '.npy'

np.save(numpy_name,paired)
plt.plot(times,dist_arr)
plt.xlabel('time (ns)')
plt.ylabel('distance $(\AA)$')
plt.savefig(str(sys.argv[5]))

