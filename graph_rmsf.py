import numpy as np
import matplotlib.pyplot as plt 

dat_files = ['/Users/mastore/temp/rmsf/WT/350K/1/elbow_rmsf.xvg','/Users/mastore/temp/rmsf/WT/350K/2/elbow_rmsf.xvg','/Users/mastore/temp/rmsf/WT/350K/3/elbow_rmsf.xvg']
start_res = 845
colors = ['tab:blue','green','orange']
labels = ['simulation 1', 'simulation 2', 'simulation 3']
y_max = 0

for i in range(len(dat_files)):
    data = np.loadtxt(dat_files[i])
    res = np.arange(start_res,start_res+len(data[:,0]),dtype=int)
    rmsf = data[:,1]*10
    print(np.max(rmsf))
    if y_max < np.max(rmsf):
        y_max = np.max(rmsf) 

    plt.plot(res,rmsf,color=colors[i],label=labels[i])

plt.ylim([0, 1.1*y_max])
plt.xlim([start_res-1, res[-1]+1])
plt.ylabel('WT R-domain Elbow RMSF $\AA$',size=16)
plt.xlabel('Amino Acid Number',size=16)
plt.legend()
plt.savefig('WT.pdf')
plt.show()
