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
file_list=['temp.sph','6msm_moved.sph']
legend=["cryo-em structure","open"]
y_max=0

for i in range(len(file_list)):

    u=mda.Universe(file_list[i],topology_format='PDB')

    z=u.atoms.positions[:,2]
    rad=u.atoms.tempfactors
    z, rad = (list(t) for t in zip(*sorted(zip(z, rad))))
    z=np.array(z)
    rad=np.array(rad)

    rad=rad[z>115]
    z=z[z>115]
    rad=rad[z<160]
    z=z[z<160]

    if y_max < np.max(rad):
        y_max=np.max(rad)
        plt.ylim([0,y_max*1.1])

    plt.xlabel('distance into pore (z axis $\AA$)')
    plt.ylabel('pore radius $\AA$')
    plt.plot(z,rad,label=legend[i])

    plt.legend()

#plt.savefig('test.pdf')
plt.title('Comparison Between Flooded Structure and Cryo EM Structure')
plt.show()
