import matplotlib.pyplot as plt
import os 
import numpy as np

os.chdir('..')

#label=['I37R_310K','I37R_350K','R352Q','WT']
label=['I37R','R352Q','WT']
label1=['I37R_310K_1','I37R_310K_2','I37R_310K_3']
label2=['R352Q_310K_1','R352Q_310K_2','R352Q_310K_3']
label3=['WT_310K_1','WT_310K_2','WT_310K_3']

filelist=['_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_I37R_3_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_2_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_1_rmsd.dat']
#filelist1=['_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_I37R_1_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_I37R_2_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_I37R_3_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_350K_I37R_2_rmsd.dat']

#filelist2=['_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_WT_1_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_WT_2_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_WT_3_rmsd.dat']

#filelist3=['_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_1_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_2_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_3_rmsd.dat']

max_rmsd=0

plt.figure()
for i in range(len(label)):

    data=np.loadtxt(filelist[i])
    rmsd=data[:,-1]
    t=data[:,1]*0.001
    plt.plot(t,rmsd,label=label[i])
    if max_rmsd <= np.amax(rmsd):
        max_rmsd = np.amax(rmsd)

plt.ylim([0,max_rmsd*1.1])
plt.legend()
plt.show()
#plt.savefig('ensemble_rmsds.pdf')
