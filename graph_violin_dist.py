import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

data_arr = [ 'ARD_distances/apo/ARD_distances.dat', 'ARD_distances/toxin/ARD_distances.dat']
legend_arr = ['apo','toxin bound']
color_arr = ['tab:blue','red']

plt.figure(figsize=(10,5))
temp_max=0
temp_min=np.inf

do_touching=False
cutoff = -1
for i in range(len(data_arr)):
    print(i)
    data=np.loadtxt(data_arr[i])[0:,:]
    a_max=np.max(data[:,1])
    a_min=np.min(data[:,1])
    if temp_max < a_max:
        temp_max=a_max
    if temp_min > a_min:
        temp_min = a_min
    
    


for i in range(len(data_arr)):


    data=np.loadtxt(data_arr[i])

    data=np.array(data)
    data=data.reshape((len(data),2))

    x=data[0:cutoff,0]
    y=data[0:cutoff,1]

    if do_touching:
        touching_proportion[i] = np.sum(y<touching_dist)/len(y)

    print(len(x))

    plt.violinplot(y, [i], showextrema=False)

if do_touching:
    touching_mean = round(np.mean(touching_proportion)*100,1)
    touching_std = round(np.std(touching_proportion)*100,1)

    plt.plot([0,cutoff],[touching_dist,touching_dist], 'k--', linewidth=1.0)

plt.ylabel('Distance ($\mathrm{\AA}$)',fontsize=16)
#plt.xlabel('System',fontsize=16)
plt.xticks([0,1], labels=['Apo','Toxin-Bound'],fontsize=14)
plt.yticks(fontsize=14)
plt.legend(prop={'size':12})

#plt.xlim([0,1000])
plt.ylim([temp_min*0.9,temp_max*1.1])
if do_touching:
    text_string = "Salt bridge contact proportion " + str(touching_mean) + "$\mathrm{\pm}$ "  + str(touching_std) + "%"

    plt.text (1.0,0.5,s=text_string,fontsize=16)
plt.tight_layout()
plt.savefig('violin_ARD_distances.pdf')
plt.show()
