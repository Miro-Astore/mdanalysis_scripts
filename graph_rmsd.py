import numpy as np
import os 

import MDAnalysis as mda
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec

cwd=os.getcwd()
#system_list=np.loadtxt('system_list',dtype=str)
#$system_list=['WT_rmsd_tm6_1.dat','WT_rmsd_tm6_2.dat','WT_rmsd_tm6_3.dat']
system_list=['R347P_rmsd_tm6_1.dat','R347P_rmsd_tm6_2.dat','R347P_rmsd_tm6_3.dat']
colors=['orange','tab:blue','red','yellow']

num_systems=len(system_list)

plt.figure(figsize=(5,5))
numblocks=1
sys_count=1

for i in range(len(system_list)):
    sys_name=i
    print (sys_name)
    data = np.loadtxt(system_list[i])
    x = data[:,1]*0.001
    y = data[:,-1]
    plt.plot(x,y,label='simulation ' + str(sys_name),color=colors[i])

plt.ylim([0,5])
plt.xlim([0,2000])


#plt.title('RMSD of the lasoo domain in the wild type and the I37R mutant')
plt.xlabel('Time (ns)',fontsize=16)
plt.ylabel('TM6 RMSD ($\AA$)',fontsize=16)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.legend(fontsize=14)


#os.chdir(cwd)
plt.tight_layout()
plt.savefig('R347P.pdf')
plt.show()

