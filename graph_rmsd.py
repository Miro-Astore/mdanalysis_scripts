import numpy as np
import os 

import MDAnalysis as mda
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec

cwd=os.getcwd()
#system_list=np.loadtxt('system_list',dtype=str)
#$system_list=['WT_rmsd_tm6_1.dat','WT_rmsd_tm6_2.dat','WT_rmsd_tm6_3.dat']
system_list=['WT_rmsd_icl_1.dat','WT_rmsd_icl_2.dat','WT_rmsd_icl_3.dat']
#system_list=['L1065R_rmsd_icl_1.dat','L1065R_rmsd_icl_2.dat','L1065R_rmsd_icl_3.dat','L1065R_rmsd_icl_4.dat','L1065R_rmsd_icl_5.dat']
#system_list=['L1065R_rmsd_icl_1.dat','L1065R_rmsd_icl_2.dat','L1065R_rmsd_icl_3.dat','L1065R_rmsd_icl_4.dat','L1065R_rmsd_icl_5.dat']
system_list=['L1065P_rmsd_icl_1.dat','L1065P_rmsd_icl_2.dat','L1065P_rmsd_icl_3.dat']
colors=['orange','tab:blue','red','yellow','purple','green']

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

plt.ylim([0,6])
plt.xlim([0,2000])


#plt.title('RMSD of the lasoo domain in the wild type and the I37R mutant')
plt.xlabel('Time (ns)',fontsize=16)
plt.ylabel('ICL RMSD ($\AA$)',fontsize=16)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.legend(fontsize=14)


#os.chdir(cwd)
plt.tight_layout()
plt.savefig('L1065P_icl.pdf')
plt.show()

