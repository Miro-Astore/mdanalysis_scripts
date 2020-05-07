import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

data=np.loadtxt('I37R_3_37R_823E.npy')
data2=np.loadtxt('I37R_3_37R_826E.npy')
data=data[0::50,:]
data2=data2[0::50,:]
data=[[float(row[0]),float(row[1])] for row in data if row[0] >= 800]
data2=[[float(row[0]),float(row[1])] for row in data2 if row[0] >= 800]
data=np.array(data)
data=data.reshape((len(data),2))
data2=np.array(data2)
data2=data2.reshape((len(data2),2))
plt.figure(figsize=(10,5))
gridspec.GridSpec(2,1)
plt.subplot2grid((2,1),(0,0),colspan=1,rowspan=1)
plt.plot(data[:,0],data[:,1],label='R37-E823',color='orange')
plt.xlabel('time (ns)',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0,20])
plt.ylabel("N-O distance ($\AA$)",fontsize=16)
plt.legend(prop={'size':12})
#plt.tight_layout()
plt.subplot2grid((2,1),(1,0),colspan=1,rowspan=1)
plt.plot(data2[:,0],data2[:,1],label='R37-E826',color='orange')
plt.xlabel('time (ns)',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0,20])
plt.legend(prop={'size':12})
#plt.tight_layout()
#plt.plot(data2[:][0],data2[:][1])

plt.ylabel("N-O distance ($\AA$)",fontsize=16)
plt.show()
