import MDAnalysis as mda
import sys
u=mda.Universe(sys.argv[1])
print(str(u.trajectory[-1].time * 0.001) + 'ns in  ' + str(u.trajectory[-1].frame) + ' frames')
