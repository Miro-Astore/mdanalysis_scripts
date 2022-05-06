import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

data_arr = [ 'flooded_wt/_rds_ma2374_flooding_CFTR_steer_open_candidate_1_pca_1_2_R334_E1126.dat', 'flooded_wt/_scratch_f91_ma2374_flooding_CFTR_steer_open_test_stability_R334_E1126_2.dat','flooded_wt/_scratch_f91_ma2374_flooding_CFTR_steer_open_test_stability_R334_E1126_5.dat']

cutoff = 500
touching_dist=4
touching_proportion = np.zeros(len(data_arr))

legend_arr = [ 'simulation 1', 'simulation 2','simulation 3']
color_arr = ['tab:blue','green','red'] 

plt.figure(figsize=(10,5))
temp_max=0

for i in range(len(data_arr)):
    print(i)
    data=np.loadtxt(data_arr[i])[0::2,:]
    a_max=np.max(data[:,1])
    if temp_max < a_max:
        temp_max=a_max
    
    


for i in range(len(data_arr)):


    data=np.loadtxt(data_arr[i])

    data=np.array(data)
    data=data.reshape((len(data),2))

    x=data[0:cutoff,0]
    y=data[0:cutoff,1]
    touching_proportion[i] = np.sum(y<touching_dist)/len(y)

    print(len(x))

    plt.plot(x,y,label=legend_arr[i],color=color_arr[i], linewidth=0.6)

touching_mean = round(np.mean(touching_proportion)*100,1)
touching_std = round(np.std(touching_proportion)*100,1)

plt.plot([0,cutoff],[touching_dist,touching_dist], 'k--', linewidth=1.0)

plt.ylabel('Distance ($\mathrm{\AA}$)',fontsize=16)
plt.xlabel('time (ns)',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(prop={'size':12})

plt.xlim([0,cutoff])
plt.ylim([0,temp_max])
text_string = "Salt bridge contact proportion " + str(touching_mean) + "$\mathrm{\pm}$ "  + str(touching_std) + "%"

plt.text (1.0,0.5,s=text_string,fontsize=16)
plt.tight_layout()
plt.savefig('wt_flooded_E1126_R334.pdf')
plt.show()
