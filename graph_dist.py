import numpy as np
import os 

import MDAnalysis as mda
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec

cwd=os.getcwd()
#system_list=np.loadtxt('system_list',dtype=str)
#system_list=['R347P_distance_1.dat','R347P_distance_2.dat','R347P_distance_3.dat']
system_list=['WT_distance_1.dat','WT_distance_2.dat','WT_distance_3.dat']
colors=['orange','tab:blue','red','yellow']

num_systems=len(system_list)

plt.figure(figsize=(5,5))
numblocks=1
sys_count=1

for i in range(len(system_list)):
    sys_name=i
    print (sys_name)
    data = np.loadtxt(system_list[i])
    print(data)
    x = data[:,0]
    y = data[:,-1]
    plt.plot(x,y,label='simulation ' + str(sys_name),color=colors[i])

plt.ylim([8,18])
plt.xlim([0,2000])


#plt.title('RMSD of the lasoo domain in the wild type and the I37R mutant')
plt.xlabel('Time (ns)',fontsize=16)
plt.ylabel('CA Distance R347 - 924 ($\AA$)',fontsize=16)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.legend(fontsize=14)


#os.chdir(cwd)
plt.tight_layout()
plt.savefig('WT_distance.pdf')
plt.show()

