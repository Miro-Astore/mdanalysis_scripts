import MDAnalysis as mda
import IPython
import re 
import matplotlib.pyplot as plt
import numpy as np 


#data=np.loadtxt('../lastframe.sph',dtype=str)
#rad=np.array(data[:,-1],dtype=float)
#print (rad)
#z=np.array(data[:,8],dtype=float)[rad<12]
#rad=rad[rad<12]
#z=np.array(range(len(rad)))*0.2

u=mda.Universe('../R352Q_last.sph',topology_format='PDB')

z=u.atoms.positions[:,2]
rad=u.atoms.tempfactors
z, rad = (list(t) for t in zip(*sorted(zip(z, rad))))

plt.xlabel('distance into pore (z axis $\AA$)')
plt.ylabel('pore radius $\AA$')
plt.plot(z,rad,label='R352Q')


u=mda.Universe('../I37R_test_doubt.sph',topology_format='PDB')

z=u.atoms.positions[:,2]
z=list(z)
rad=u.atoms.tempfactors
rad=list(rad)
z, rad = (list(t) for t in zip(*sorted(zip(z, rad))))

plt.xlabel('distance into pore (z axis $\AA$)')
plt.ylabel('pore radius $\AA$')
plt.plot(z,rad,label='I37R')

u=mda.Universe('../wt_last.sph',topology_format='PDB')

z=u.atoms.positions[:,2]
rad=u.atoms.tempfactors
z, rad = (list(t) for t in zip(*sorted(zip(z, rad))))

plt.xlabel('distance into pore (z axis $\AA$)')
plt.ylabel('pore radius $\AA$')
plt.plot(z,rad,label='wt')
plt.legend()

plt.savefig('test.pdf')
