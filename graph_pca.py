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

plt.figure(figsize=(10,5))
labels = ['simulation 1','simulation 2', 'simulation 3', 'simulation 4']
color = ['tab:blue','green', 'red', 'orange']
 
num_pca_components_to_graph = 2
AX=gridspec.GridSpec(num_pca_components_to_graph,1)
plot_array = []
for i in range(num_pca_components_to_graph):
    plot_array.append(plt.subplot2grid((num_pca_components_to_graph,1),(i,0)))

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
    curr_dat = data.sets()[i]
    #temp_plot=plt.subplot2grid((num_pca_components_to_graph,1),(i,0))
    for j in range(num_replicates):
        #y = curr_dat[start_indices[j]:start_indices[j+1]]
        t = curr_dat[start_indices[j]:start_indices[j+1]:500,0]*0.001
        y = curr_dat[start_indices[j]:start_indices[j+1]:500,1]
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
        plot_array[i].set_xlabel('Time (ns)',fontsize=14)
        plot_array[i].set_ylabel('PC' + str(i+1) + ' Projection (nm)',fontsize=14)
        plot_array[i].plot(t,y,color=color[j],label=labels[j])


        plot_array[i].legend()

for i in range(num_pca_components_to_graph):
    plot_array[i].set_ylim([1.1*y_min,1.1*y_max])
    plot_array[i].set_xlim([t_min,250+t_max])
    plot_array[i].tick_params(axis='both',  labelsize=12)

plt.tight_layout()
plt.show()
