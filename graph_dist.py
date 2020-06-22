import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

data_arr = [ '23_458.dat', '37_403.dat', '30_417.dat']
color_arr = ['orange'] * 3
legend_arr=['23_458', '37_403', '30_417']

plt.figure(figsize=(10,5))
gridspec.GridSpec(len(data_arr),1)
temp_max=0

for i in range(len(data_arr)):
    data=np.loadtxt(data_arr[i])[0::50,:]
    a_max=np.max(data[:,1])
    if temp_max < a_max:
        temp_max=a_max
    
    


for i in range(len(data_arr)):


    data=np.loadtxt(data_arr[i])[0::50,:]


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
