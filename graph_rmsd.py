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

cwd=os.getcwd()
system_list=np.loadtxt('system_list',dtype=str)
num_systems=len(system_list)

plt.figure(figsize=(10,5))
numblocks=1
sys_count=1
gridspec.GridSpec(num_systems+1,1)
plt.subplot2grid((num_systems+1,1),(0,0),colspan=1,rowspan=1)

for i in system_list:
    os.chdir(i)
    sys_name=i.split('/')[-2].upper()
    print (sys_name)
    data = np.loadtxt('rmsd.dat')
    x = data[:,1]*0.001
    y = data[:,3]
    plt.plot(x,y,label=sys_name)
    plt.title('RMSD of the lasoo domain in the wild type and the I37R mutant')
    plt.xlabel('time (ns)')
    plt.ylabel('lasoo domain RMSD ($\AA$)')
    plt.legend()

color_vec=['blue']*70
color_vec[36]='orange'
for i in system_list:
    os.chdir(i)
    sys_name=i.split('/')[-2].upper()
    print (sys_name)
    data = np.loadtxt('rmsd.dat')
    resnum=len(data[0,3:-1])
    x = data[:,1]*0.001
    y = data[:,3]
    resnum=len(data[0,3:-1])

    plt.subplot2grid((num_systems+1,1),(sys_count,0),colspan=1,rowspan=1)
    data = np.loadtxt('rmsd.dat')
    resnum=len(data[0,3:-1])
    blocks=chunkIt(np.array(data[:,0]),numblocks)
    blocks=np.array(blocks).astype(int)

    width = 0.9
    inds = np.array(range(1,len(data[0,3:-1])+1))

    for j in range(numblocks):
        block=blocks[j]
        resrmsd=np.zeros(resnum)
        for i in range(resnum):
            resrmsd[i]=np.mean(data[block,4+i])
        plt.bar(inds+j*width,resrmsd,width,color=color_vec)
        #plt.ylim([0,6])
        
    plt.title(sys_name + ' residue specific lasoo domain RMSDs ')
    plt.ylabel('Residue specific RMSDs $(\AA)$')
    plt.xlabel('residue number')

    data = np.loadtxt('rmsd.dat')
    resnum=len(data[0,3:-1])
    x = data[:,1]*0.001
    y = data[:,3]
    sys_count=sys_count+1

os.chdir(cwd)
plt.tight_layout()
plt.savefig('fig.pdf')

