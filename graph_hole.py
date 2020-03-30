import MDAnalysis as mda
import re 
import matplotlib.pyplot as plt
import numpy as np 


data=np.loadtxt('../lastframe.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
print (rad)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2

plt.xlabel('distance into pore (z axis $\AA$)')
plt.ylabel('pore radius $\AA$')
plt.plot(z,rad)
plt.legend()

plt.savefig('test.pdf')
