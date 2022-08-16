import numpy as np 
import matplotlib.pyplot as plt
import LoadGrace
from functools import reduce
from matplotlib.ticker import MaxNLocator
import re
import glob
import itertools
import matplotlib.gridspec as gridspec

labels = ['simulation 1','simulation 2', 'simulation 3', 'simulation 4']
color = ['tab:blue','green', 'red', 'orange']
 
num_pca_components_to_graph = 2
AX=gridspec.GridSpec(num_pca_components_to_graph,1)
data = LoadGrace.LoadGrace('proj1.xvg')
num_replicates=4
data_in_each_replicate = [100001, 100001, 100002, 109997]

start_indices = np.zeros(num_replicates+1,dtype=int)
start_indices[0]=int(start_indices[0])
for i in range(num_replicates):
    start_indices[i+1]=int(start_indices[i])+int(data_in_each_replicate[i])
total_pca_components = len(data.sets())

for i in range(num_pca_components_to_graph):
    curr_dat = data.sets()[i]
    plt.subplot2grid((num_pca_components_to_graph,1),(i,0))
    for j in range(num_replicates):
        #y = curr_dat[start_indices[j]:start_indices[j+1]]
        t = curr_dat[start_indices[j]:start_indices[j+1]:500,0]
        y = curr_dat[start_indices[j]:start_indices[j+1]:500,1]
        print(y)
        plt.plot(t,y,color=color[j],label=labels[j])


plt.legend()
plt.show()
