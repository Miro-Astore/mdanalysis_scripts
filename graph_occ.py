import numpy as np
import matplotlib.pyplot as plt 



wt=np.load('_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_2_occ.npy')
R352Q=np.load('_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_R352Q_1_occ.npy')
inds=range(int(wt[0][-1]))
width=0.4

prop_wt=wt[1][0:-1]/(wt[1][-1])
#print(wt[1][0:-1][-1])
#print(len(prop_wt))
#print(len(inds))
plt.figure(figsize=(20,5))
prop_mut=R352Q[1][0:-1]/R352Q[1][-1]
inds=np.array([float(x) for x in inds ])
print((inds))
print((prop_wt))
plt.bar(inds-width*0.5,prop_wt,width,label='WT',color='tab:blue')
plt.xlim([-0.5,inds[-1]+0.5])
plt.bar(inds+width*0.5,prop_mut,width,label='R352Q',color='red')
res_list=wt[0][0:-1]
res_list=[int(x) for x in res_list]
plt.xticks(inds,res_list,fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Residue Number',fontsize=16)
plt.ylabel('Chloride Coordination Number',fontsize=16)
plt.ylim([0,1.2])
#plt.bar(inds+width*0.5,prop_wt,width,label='WT',color='tab:blue')
#plt.bar(inds+width*0.5,res_results[1,:],width,color='tab:blue',label='R352Q',color='red')
plt.tight_layout()
plt.show()
