import numpy as np 
import sys
sys.path.insert(0, '.') # added on recomendation from https://stackoverflow.com/questions/279237/import-a-module-from-a-relative-path/6098238#6098238 probably a bad fix
import matplotlib.pyplot as plt
from mdanalysis_scripts import LoadGrace
from functools import reduce
from matplotlib.ticker import MaxNLocator
import re
import glob
import itertools
import matplotlib.gridspec as gridspec

plt.figure(figsize=(5,5))
labels = ['simulation 1','simulation 2', 'simulation 3', 'simulation 4']
color = ['tab:blue','green', 'red', 'orange']
 
num_pca_components_to_graph = 2
AX=gridspec.GridSpec(1,1)
plot_array = []
for i in range(num_pca_components_to_graph):
    plot_array.append(plt.subplot2grid((1,1),(0,0)))

data = LoadGrace.LoadGrace('cov-domain-proj1to9.xvg')
num_replicates=4
data_in_each_replicate = [100000, 100001, 100001, 109997]

start_indices = np.zeros(num_replicates+1,dtype=int)
start_indices[0]=int(start_indices[0])
for i in range(num_replicates):
    start_indices[i+1]=int(start_indices[i])+int(data_in_each_replicate[i])
total_pca_components = len(data.sets())
comps =(data.sets())
t_min = 0
t_max = 0
y_min = 0
y_max = 0

for i in range(num_pca_components_to_graph):
    x_dat = data.sets()[0]
    y_dat = data.sets()[1]
    #temp_plot=plt.subplot2grid((num_pca_components_to_graph,1),(i,0))
    for j in range(num_replicates):
        #y = curr_dat[start_indices[j]:start_indices[j+1]]
        t = x_dat[start_indices[j]:start_indices[j+1]:10,1]
        y = y_dat[start_indices[j]:start_indices[j+1]:10,1]
        if np.max(t) > t_max:
            t_max = np.max(t)
        if np.min(t) < t_min:
            t_min = np.min(t)
        if np.max(y) > y_max:
            y_max = np.max(y)
        if np.min(y) < y_min:
            y_min = np.min(y)
        #t = curr_dat[start_indices[j]:start_indices[j+1],0]
        #y = curr_dat[start_indices[j]:start_indices[j+1],1]
        plot_array[i].set_xlabel('PC 1 Projection (nm)',fontsize=20)
        plot_array[i].set_ylabel('PC 2 Projection (nm)',fontsize=20)
        plot_array[i].scatter(t,y,color=color[j],label=labels[j],marker='.',s=0.2)


        plot_array[i].legend()

for i in range(num_pca_components_to_graph):
    plot_array[i].set_ylim([1.1*y_min,1.1*y_max])
    plot_array[i].set_xlim([1.1*y_min,1.1*y_max])
    plot_array[i].tick_params(axis='both',  labelsize=16)

plt.tight_layout()
#plt.show()
plt.savefig('scatter_pca_1_2.png' )
