import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

data_arr = [ 'flooded_wt/_rds_ma2374_flooding_CFTR_steer_open_candidate_1_pca_1_2_R334_E1126.dat', 'flooded_wt/_scratch_f91_ma2374_flooding_CFTR_steer_open_test_stability_R334_E1126_2.dat','flooded_wt/_scratch_f91_ma2374_flooding_CFTR_steer_open_test_stability_R334_E1126_5.dat']
legend_arr = [ '1', '2','3']
color_arr = ['tab:blue','green','red'] 

plt.figure(figsize=(10,5))
gridspec.GridSpec(len(data_arr),1)
temp_max=0

for i in range(len(data_arr)):
    print(i)
    data=np.loadtxt(data_arr[i])[0::,:]
    a_max=np.max(data[:,1])
    if temp_max < a_max:
        temp_max=a_max
    
    


for i in range(len(data_arr)):


    data=np.loadtxt(data_arr[i])[0::,:]


    data=np.array(data)
    data=data.reshape((len(data),2))

    plt.subplot2grid((len(data_arr),1),(i,0),colspan=1,rowspan=1)
    plt.plot(data[:,0],data[:,1],label=legend_arr[i],color=color_arr[i])
    plt.xlabel('time (ns)',fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(prop={'size':12})
    plt.ylim([0,temp_max])


plt.show()
