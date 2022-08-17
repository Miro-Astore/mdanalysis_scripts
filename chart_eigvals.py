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
special_vals = [0, 1]
color_list = ['red','blue']
 
num_special_components = len(special_vals)
num_components_chart = 9
colors=['']*num_components_chart



AX=gridspec.GridSpec(1,1)

bars_pos = list(range(1,num_components_chart+1))

data = np.loadtxt('eigvals.dat')
x = data [:num_components_chart,0]
x = [int(i) for i in x]
y = data [:num_components_chart,1]

kin_var = 0 
for i in range(num_components_chart):
    if i in special_vals:
        colors [i] = color_list[0]
        kin_var = kin_var + data[i,1] 
    else:
        colors [i] = color_list[1]
plt.bar(x,y,align='center',color=colors)
plt.xticks(x,x,fontsize=16)
plt.yticks(x,x,fontsize=16)
plt.xlabel('Principal Component', fontsize=20)
plt.ylabel('Eigenvalue (nm$^2$)',fontsize=20)
plt.ylim([0,3.0])
plt.tight_layout()
print(kin_var/np.sum(data[:,1]))
plt.savefig('eigvals.pdf')
plt.show()

