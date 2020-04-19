import numpy as np
import os 

import MDAnalysis as mda
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec
def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out


plt.figure(figsize=(10,5))
numblocks=1
sys_count=1
gridspec.GridSpec(2,1)
plt.subplot2grid((2,1),(0,0),colspan=1,rowspan=1)

labels=['5','6','9']
data1 = np.loadtxt('_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_3_sum_tms_rmsd.dat')
data2 = np.loadtxt('_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_1_sum_tms_rmsd.dat')
num_cols=len(data1[0,3:-1])
print(num_cols)

plt.subplot2grid((2,1),(0,0),colspan=1,rowspan=1)
for i in range(num_cols+1):
    x = data1[:,1]*0.001
    y = data1[:,3+i]
    plt.plot(x,y,label=labels[i])

plt.subplot2grid((2,1),(1,0),colspan=1,rowspan=1)
for i in range(num_cols+1):
    x = data2[:,1]*0.001
    y = data2[:,3+i]
    print(3+i)
    plt.plot(x,y,label=labels[i])
    plt.title('R352Q')

#plt.title('TM domains')
plt.xlabel('time (ns)')
plt.ylabel('RMSD $(\AA)$')
plt.legend()
plt.tight_layout()
plt.savefig('tms.pdf')
plt.show()

