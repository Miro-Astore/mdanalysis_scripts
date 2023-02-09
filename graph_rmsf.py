import numpy as np
import matplotlib.pyplot as plt 

dat_files = ['/mnt/ceph/users/mastore/apo_trpv1/rmsf_S5_A.xvg','/mnt/home/mastore/ceph/trpv1_soybean_nanodisc/rmsf_S5_A.xvg']
start_res = 562
colors = ['red','green']
labels = ['apo', 'toxin-bound']
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
plt.ylabel('S5 RMSF $\AA$',size=16)
plt.xlabel('Amino Acid Number',size=16)
plt.legend()
plt.savefig('S5_rmsf.pdf')
plt.show()
