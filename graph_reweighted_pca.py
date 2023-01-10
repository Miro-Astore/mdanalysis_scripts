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

fig, ax1 = plt.subplots(figsize=(5,5),ncols=1)

#labels = ['simulation 1','simulation 2', 'simulation 3', 'simulation 4']
#color = ['tab:blue','green', 'red', 'orange']
 
num_pca_components_to_graph = 2
AX=gridspec.GridSpec(1,1)
plot_array = []
for i in range(num_pca_components_to_graph):
    plot_array.append(plt.subplot2grid((1,1),(0,0)))

data = LoadGrace.LoadGrace('cov-domain-tetramer-proj1to9.xvg')

num_replicates=1

data_in_each_replicate = [19731]
total_data_points = np.sum(data_in_each_replicate)

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
        t = x_dat[start_indices[j]:start_indices[j+1],1]
        y = y_dat[start_indices[j]:start_indices[j+1],1]
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
        #plot_array[i].set_xlabel('PC 1 Projection (nm)',fontsize=20)
        #plot_array[i].set_ylabel('PC 2 Projection (nm)',fontsize=20)
        #plot_array[i].scatter(t,y,color=color[j],label=labels[j],marker='.',s=0.2)


        #plot_array[i].legend()


side_spacing = 1.00
bounds = [-70, 70, -70, 70] 
num_side = ['']*int(len(bounds)/2)


for i in range(num_pca_components_to_graph):
    #plot_array[i].set_ylim([bounds[0],bounds[1]])
    #plot_array[i].set_xlim([bounds[2],bounds[3]])
    plot_array[i].tick_params(axis='both',  labelsize=16)

j = 0

for i in range(0,len(bounds),2) :
    num_side[j]=int(np.round(np.abs(bounds[i]-bounds[i+1])/side_spacing))
    j=j+1

count_mat = np.zeros(num_side)

for i in range(total_data_points): 
    x = x_dat [i]
    y = y_dat [i]
    x = x[1]
    y = y[1]
    x_shifted = x - bounds[0]
    y_shifted = y - bounds[2]

    #choosing x bin 
    x_coord = int(x_shifted/side_spacing)
    y_coord = int(y_shifted/side_spacing)

    if  y < bounds[3] and x < bounds[1] :

        count_mat[x_coord,y_coord] = count_mat[x_coord,y_coord] + 1

count_mat = count_mat / total_data_points
free_energy_mat=np.zeros(np.shape(count_mat))

for i in range(np.shape(count_mat)[0]):
    for j in range(np.shape(count_mat)[1]):
        if count_mat[i][j] == 0 :
            free_energy_mat[i][j] = np.nan
            continue
        free_energy_mat[i][j] = (-1/ (1.987204259e-3*310)) * np.log(count_mat[i][j])


max_non_nan = np.nanmax(np.nanmax(free_energy_mat))
min_non_nan = np.nanmin(np.nanmin(free_energy_mat))



for i in range(np.shape(count_mat)[0]):
    for j in range(np.shape(count_mat)[1]):
        if np.isnan(free_energy_mat[i][j]): 
            free_energy_mat [i][j] = np.nan
        else:
            free_energy_mat [i][j] = free_energy_mat [i][j] - min_non_nan

linear = np.reshape(free_energy_mat, [1,  np.shape(free_energy_mat)[0] * np.shape(free_energy_mat)[1]])


#for i in range(num_pca_components_to_graph):
#    plot_array[i].imshow(free_energy_mat)
    

#for i in range((num_side[0])):
#    for j in range((num_side[1])):
#        temp_lower_bound1 = bounds[0] + num_side[0]*side_spacing*i
#        temp_upper_bound1 = bounds[1] + num_side[0]*side_spacing*i
#        temp_lower_bound2 = bounds[2] + num_side[1]*side_spacing*j
#        temp_upper_bound2 = bounds[3] + num_side[1]*side_spacing*j
#
#        count_mat[i][j] =  np.sum([data_point for data_point in  x_dat < temp_upper_bound1 and x_dat > temp_lower_bound1 and y_dat ]) / total_data_points


plot = ax1.imshow(free_energy_mat,origin='lower')

#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])

#fig.colorbar(plot,cax=cax)

#plt.tight_layout()
plt.show()
#plt.savefig('boltzmann_weighted_fes_pca_1_2.pdf' )
