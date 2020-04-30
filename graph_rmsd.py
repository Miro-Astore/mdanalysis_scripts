import numpy as np
import os 

import MDAnalysis as mda
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec

cwd=os.getcwd()
#system_list=np.loadtxt('system_list',dtype=str)
system_list=['_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_I37R_3_sum_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_1_sum_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_1_sum_rmsd.dat']
colors=['orange','tab:blue','red']

num_systems=len(system_list)

plt.figure(figsize=(10,5))
numblocks=1
sys_count=1

for i in range(len(system_list)):
    sys_name=system_list[i].split('_')[-3].upper()
    print (sys_name)
    data = np.loadtxt(system_list[i])
    x = data[:,1]*0.001
    y = data[:,2]
    plt.plot(x,y,label=sys_name,color=colors[i])

#plt.title('RMSD of the lasoo domain in the wild type and the I37R mutant')
plt.xlabel('time (ns)')
plt.ylabel('Transmembrane RMSD ($\AA$)')
plt.legend()


#os.chdir(cwd)
plt.tight_layout()
plt.savefig('fig.pdf')
plt.show()

