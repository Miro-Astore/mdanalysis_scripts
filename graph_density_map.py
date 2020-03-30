import MDAnalysis.analysis.density as de 
import numpy as np 
import MDAnalysis as mda
import matplotlib.pyplot as plt

x=de.Density('wt_test.dx')
y=(x.centers())

#print(x.interpolated())

#FF = x.interpolated(p1[0],p1[1],p1[2])
u=mda.Universe('wt_last.sph',topology_format='PDB')
values=np.zeros(len(u.atoms))

for i in range(len(u.atoms)):
    pos=u.atoms.positions[i]
    values[i]=(x.interpolated(pos[0],pos[1],pos[2]))

plt.plot(range(len(u.atoms)),values,label='wt')

x=de.Density('_test.dx')
y=(x.centers())

#print(x.interpolated())

#FF = x.interpolated(p1[0],p1[1],p1[2])
u=mda.Universe('_test_doubt.sph',topology_format='PDB')
values=np.zeros(len(u.atoms))

for i in range(len(u.atoms)):
    pos=u.atoms.positions[i]
    values[i]=(x.interpolated(pos[0],pos[1],pos[2]))

plt.plot(range(len(u.atoms)),values,label='')

x=de.Density('_test.dx')
y=(x.centers())

u=mda.Universe('_last.sph',topology_format='PDB')
values=np.zeros(len(u.atoms))

for i in range(len(u.atoms)):
    pos=u.atoms.positions[i]
    values[i]=(x.interpolated(pos[0],pos[1],pos[2]))

plt.plot(range(len(u.atoms)),values,label='')

plt.legend()
plt.show()


