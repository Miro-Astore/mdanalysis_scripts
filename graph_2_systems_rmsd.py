import numpy as np
import os 

import MDAnalysis as mda
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec

system_list=['_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_I37R_3_res_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_1_res_rmsd.dat']

def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

#cwd=os.getcwd()
#system_list=np.loadtxt('system_list',dtype=str)
#num_systems=len(system_list)

plt.figure(figsize=(10,5))
numblocks=1
sys_count=1
gridspec.GridSpec(2,1)
plt.subplot2grid((2,1),(0,0),colspan=1,rowspan=1)

for i in system_list:
    sys_name=i.split('_')[-4].upper()
    print (sys_name)
    data = np.loadtxt(str(i))
    x = data[:,1]*0.001
    y = data[:,3]
    plt.plot(x,y,label=sys_name)

plt.title('RMSD of the lasoo domain in the wild type and the I37R mutant')
plt.xlabel('time (ns)')
plt.ylabel('lasoo domain RMSD ($\AA$)')
plt.legend()

color_vec=['blue']*70
color_vec[36]='orange'
res_results=np.zeros((2,70))
row_num=0
for i in system_list:
    sys_name=i.split('_')[-4].upper()
    print (sys_name)
    data = np.loadtxt(str(i))
    resnum=len(data[0,3:-1])
    x = data[:,1]*0.001
    y = data[:,3]
    resnum=len(data[0,3:-1])

    data = np.loadtxt(str(i))
    resnum=len(data[0,3:-1])
    blocks=chunkIt(np.array(data[:,0]),numblocks)
    blocks=np.array(blocks).astype(int)

    width = 0.45
    inds = np.array(range(1,len(data[0,3:-1])+1))

    for j in range(numblocks):
        block=blocks[j]
        resrmsd=np.zeros(resnum)
        for k in range(resnum):
            resrmsd[k]=np.mean(data[block,4+k])
        #plt.bar(inds+j*width,resrmsd,width,color=color_vec)
        res_results[row_num,:]=resrmsd
        #plt.ylim([0,6])
        
#    plt.ylabel('Residue specific RMSDs $(\AA)$')
#    plt.xlabel('residue number')

    print (str(i))
    data = np.loadtxt(str(i))
    resnum=len(data[0,3:-1])
    x = data[:,1]*0.001
    y = data[:,3]
    sys_count=sys_count+1
    row_num=row_num+1

plt.subplot2grid((2,1),(1,0),colspan=1,rowspan=1)
separation=0.1
plt.bar(np.array([37])-width*0.49,1.05*np.amax(res_results),width,color='yellow')
plt.bar(np.array([37])+width*0.49,1.05*np.amax(res_results),width,color='yellow')
plt.bar(inds-width*0.5,res_results[0,:],width,)
plt.bar(inds+width*0.5,res_results[1,:],width,)
plt.xlim([0,70])
plt.ylim([0,1.05*np.amax(res_results)])
plt.title('Residue Specific Lasoo Domain RMSDs ')
plt.xlabel('residue number')
plt.ylabel('CA RMSD ($\AA$)')

plt.tight_layout()
plt.savefig('fig.pdf')
plt.show()
