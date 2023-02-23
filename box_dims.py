import numpy as np 
import matplotlib.pyplot as plt 
import MDAnalysis as mda
import sys

print ('TOP ' + str(sys.argv[1]))
print ('traj ' + str(sys.argv[2]))

u = mda.Universe(sys.argv[1],sys.argv[2])

for ts in u.trajectory:
    print(u.dimensions[0:3])

