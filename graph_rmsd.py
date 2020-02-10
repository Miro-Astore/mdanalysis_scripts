import numpy as np

import MDAnalysis as mda
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec

plt.figure(figsize=(8,7))
numblocks=1
data=np.loadtxt('rmsd.dat')
num_steps_i37r=len(data[:,3])
timepoints=(num_steps_i37r)
def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

gridspec.GridSpec(3,1)
plt.subplot2grid((3,1),(0,0),colspan=1,rowspan=1)
data = np.loadtxt('../WT/open/rmsd.dat')
resnum=len(data[0,3:-1])
x = data[:,0]*0.4
y = data[:,3]
plt.plot(x,y,label='WT')
plt.title('RMSD of the lasoo domain in the wild type and the I37R mutant')
plt.xlabel('time (ns)')
plt.ylabel('lasoo domain RMSD ($\AA$)')
data = np.loadtxt('rmsd.dat')
resnum=len(data[0,3:-1])
x = data[:,0]*0.4
y = data[:,3]
plt.plot(x,y,label='I37R')
plt.legend()

plt.subplot2grid((3,1),(1,0),colspan=1,rowspan=1)
data = np.loadtxt('../WT/open/rmsd.dat')
resnum=len(data[0,3:-1])
blocks=chunkIt(range(timepoints),numblocks)

width = 0.9
inds = np.array(range(1,len(data[0,3:-1])+1))

for j in range(numblocks):
    block=blocks[j]
    resrmsd=np.zeros(resnum)
    for i in range(resnum):
        print (4+i)
        resrmsd[i]=np.mean(data[block,4+i])
    plt.bar(inds+j*width,resrmsd,width)
    plt.ylim([0,6])
plt.title('WT residue specific lasoo domain RMSDs (first 600ns)')
plt.ylabel('Residue specific RMSDs $(\AA)$')
plt.xlabel('residue number')

data = np.loadtxt('rmsd.dat')
resnum=len(data[0,3:-1])
x = data[:,0]*0.4
y = data[:,3]


plt.subplot2grid((3,1),(2,0),colspan=1,rowspan=1)
blocks=chunkIt(range(timepoints),numblocks)
for j in range(numblocks):
    block=blocks[j]
    resrmsd=np.zeros(resnum)
    for i in range(resnum):
        resrmsd[i]=np.mean(data[block,4+i])
    plt.bar(inds+j*width,resrmsd,width)
    plt.ylim([0,6])
plt.title('I37R residue specific lasoo domain RMSDs')
plt.ylabel('Residue specific RMSDs $(\AA)$')
plt.xlabel('residue number')


plt.tight_layout()
plt.savefig('fig.pdf')

