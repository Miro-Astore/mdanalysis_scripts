import MDAnalysis as mda
import sys
u=mda.Universe(sys.argv[1],sys.argv[2])
print(str(len(u.trajectory)*u.trajectory.dt*0.001) + 'ns in  ' + str(len(u.trajectory)) + ' frames')

