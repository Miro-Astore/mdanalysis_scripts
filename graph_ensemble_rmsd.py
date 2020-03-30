import matplotlib.pyplot as plt
import os 
import numpy as np

os.chdir('..')
plt.subplots(nrows=2,ncols=2)

list1=['_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_I37R_1_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_I37R_2_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_I37R_3_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_350K_I37R_2_rmsd.dat']

list2=['_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_1_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_2_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_3_rmsd.dat']

list3=['_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_1_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_2_rmsd.dat','_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_3_rmsd.dat']

list4=['_scratch_r16_ma2374_gmx_cftr_2nd_round_350K_I37R_2_rmsd.dat']

list_all = list1+list2+list3+list4
max_rmsd=0

for rmsd in list_all:

    data=np.array(np.loadtxt(rmsd))
    print (data[0])
    data=np.array(data)
    temp_list=np.zeros(len(data))
    for i in range(len(data)):
        print (i)
        temp_list[i]=(data[i][-1])

    temp_max=np.amax(temp_list)
    if max_rmsd < temp_max:
        max_rmsd=temp_max
print(max_rmsd)




plt.savefig('ensemble_rmsds.pdf')
