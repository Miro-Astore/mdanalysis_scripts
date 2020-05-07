import numpy as np
import matplotlib.pyplot as plt 

data=np.loadtxt('I37R_3_37R_823D.npy')
ymax=np.amax(data[:,1])
plt.ylim([0,ymax])
plt.plot(data[:,0],data[:,1])

plt.show()
