import numpy as np
from functools import reduce
import os 

import MDAnalysis as mda
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec

system_list=['WT_oxidised','WT_reduced','F55A','C77F','C77F_F55A','C60A_C77A','pT34_reduced','pT34_oxidised']

system_locations=['../WT_oxidised','../WT_reduced','../F55A','../C77F','../C77F_F55A','../C60A_C77A','../pT34_reduced','../pT34_oxidised']
dat_file_name = 'R86_E17_C-alpha_distance.dat'

replicates=['1','2','3']

colors=['tab:blue','red','orange']

def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

def factors(n):    
    return set(reduce(list.__add__, 
                    ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

def gridLayout(sys_count):
    factor_array = factors(sys_count)
    factor_array = list(factor_array)
    factor_array = np.sort(factor_array)
    target = np.sqrt (sys_count)
    closeness_array = np.abs(factor_array-target)
    x = factor_array[np.argmin (closeness_array)]
    y = sys_count / x
    return np.sort([int(x),int(y)])

#cwd=os.getcwd()
#system_list=np.loadtxt('system_list',dtype=str)
#num_systems=len(system_list)

numblocks=1
sys_count=len(system_list)
grid_dims = gridLayout(sys_count)
#plt.figure(figsize=(10,5))

fig,ax = plt.subplots(grid_dims[0],grid_dims[1], figsize=(10,5))
grid_x_pos = 0 
grid_y_pos = 0

graph_count=0

temp_max = 0 

for i in range(len(system_list)):
    grid_x_pos = int(np.floor((i)/grid_dims[1]))
    grid_y_pos = int(np.mod(i,grid_dims[1]))

    #plt.subplot2grid((grid_dims[0],grid_dims[1]),(grid_x_pos,grid_y_pos),colspan=1,rowspan=1)
    sys_name=system_list[i]
    for j in range(len(replicates)):
        data = np.loadtxt(system_locations[i] + '/' + replicates[j] + '/' + dat_file_name )
        data=np.array(data)
        #convert to nanoseconds
        x = np.array(data[:,0])
        y = data[:,1]
        #plt.plot(x,y,label=sys_name,color=colors[j])
        ax[grid_x_pos,grid_y_pos].plot(x,y,label=sys_name,color=colors[j])
        if np.max(y) > temp_max:
            temp_max = np.max(y)
    ax[grid_x_pos,grid_y_pos].set_title(sys_name,fontsize=18)
    #fig.xticks(fontsize=14)
    #fig.yticks(fontsize=14)

for i in range(len(system_list)):
    grid_x_pos = int(np.floor((i)/grid_dims[1]))
    grid_y_pos = int(np.mod(i,grid_dims[1]))
#    plt.subplot2grid((grid_dims[0],grid_dims[1]),(grid_x_pos,grid_y_pos),colspan=1,rowspan=1)
    sys_name=system_list[i]
    ax[grid_x_pos,grid_y_pos].set_ylim([0,temp_max*1.1])
    ax[grid_x_pos,grid_y_pos].tick_params(axis='x',labelsize=14)
    ax[grid_x_pos,grid_y_pos].tick_params(axis='y',labelsize=14)

#plt.xlabel('time (ns)',fontsize=16)
#plt.ylabel('Minimum N-O distance ($\AA$)',fontsize=16)

fig.tight_layout()
fig.savefig('salt_bridge_distance.pdf')
plt.show()

