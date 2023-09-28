import numpy as np
from natsort import natsorted, ns
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import os, fnmatch
import re
print('imported')


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def find_hists_overlap(meta_file,out_file,fignum):
    file_list = np.loadtxt(meta_file,dtype=str)
    file_list = file_list[:,0]
    bins_dist = np.loadtxt(out_file)
    bins_dist = bins_dist[:,0]
    overlaps = np.zeros(len(file_list))
    last_hist = []
    mins = []

    for i in range(len(file_list)):
        dat = np.loadtxt(file_list[i])
        dat = np.array(dat)
        dat = dat[:,1]
        hists = np.histogram(dat,bins_dist,density=True)[0]
        
        if i != 0:
            mins = np.amin(np.array([hists,last_hist]),axis=0)
            overlaps[i] = np.sum(np.diff(bins_dist) * mins)
            #overlaps[i-1] = np.sum(np.diff(bins_dist) * mins) + overlaps[i-1]
            
        
        last_hist = hists
        norm_hists=(np.diff(bins_dist)*hists)
        x = (bins_dist[1:] + bins_dist[:-1]) / 2
        x = x[norm_hists!=0]
        norm_hists = norm_hists[norm_hists!=0]
        fignum.plot(x,norm_hists,'b-',linewidth=0.1)
        fignum.set_ylim([0,0.5])
    print (overlaps)

    

print ('starting')
wt_list = find ('out_wt_**.txt','.')
wt_list = natsorted (wt_list, key=lambda y:y.lower())
mut_list = find ('out_mut_**.txt','.')
wt_list=wt_list[-5:]
mut_list = natsorted (mut_list, key=lambda y:y.lower())
mut_list=mut_list[-5:]
print (wt_list)
num = np.zeros(len(wt_list))
for i in range(len(wt_list)):
    pat = re.compile("[0-9]+")
    m = pat.search(wt_list[i])
    num[i] = m.group()

    

num = np.sort(num)
#last_step = [wt_list[-1],mut_list[-1]]
#pat1 = re.compile("[0-9]+")
#pat2 = re.compile("[0-9]+")
#m1 = pat1.search(last_step[0])
#m2 = pat2.search(last_step[1])
#last_step1 = m1.group()
#last_step2 = m2.group()
#print (last_step1)
last_step=np.int(np.amax(num))
fig1= plt.figure(1,figsize=(8,5))

ax1 = fig1.add_subplot(2,2,1)
ax2 = fig1.add_subplot(2,2,2)
ax3 = fig1.add_subplot(2,2,3)
ax4 = fig1.add_subplot(2,2,4)
ax1.set_ylabel('Relative $\Delta$ G')
ax3.set_ylabel('Incidence (Normalised)')
ax3.set_xlabel('Distance from Y852 CA to S945L CA $\AA$')
ax4.set_xlabel('Distance from Y852 CA to S945L CA $\AA$')

for i in wt_list:
    print (i)
    time_step = i.replace('out_wt_','')
    time_step = time_step.replace('ns.txt','')
    time_step = time_step.replace('./','')
    print (time_step)
    data_wt = np.loadtxt ('out_wt_' + str(time_step) + 'ns.txt')
    data_mut = np.loadtxt ('out_mut_' + str(time_step) + 'ns.txt')
    data_wt = data_wt[~np.isnan(data_wt).any(axis=1)]
    data_mut = data_mut[~np.isnan(data_mut).any(axis=1)]
    data_wt = data_wt[~np.isinf(data_wt).any(axis=1)]
    data_mut = data_mut[~np.isinf(data_mut).any(axis=1)]
    mean_begin_wt = int(0.3*len(data_wt[:,0]))
    mean_begin_mut = int(0.3*len(data_mut[:,0]))
    mean_wt = np.mean(data_wt[0:mean_begin_wt,1])
    mean_mut = np.mean(data_mut[0:mean_begin_mut,1])
    #mut_y = data_mut[:,1]-np.amax(data_mut[:,1])
    #wt_y = data_wt[:,1]-np.amax(data_wt[:,1])
    #mut_x = data_mut[:,0][data_mut[:,0] > 14.5] 
    #mut_y = data_mut[:,1][data_mut[:,0] > 14.5] 
    mut_x = data_mut[:,0] 
    mut_y = data_mut[:,1] 
    print(mut_y)
    print(np.amin(mut_y))
    mut_y = mut_y-np.amin(mut_y)
    wt_y = data_wt[:,1]
    wt_y = wt_y - mean_wt
    mut_y = mut_y - mean_mut
    ax1.plot(data_wt[:,0],wt_y,label=str(i) + 'ns',linewidth=0.6)
    ax1.set_title('wt')
    ax2.plot(mut_x,mut_y,label=str(i) + 'ns',linewidth=0.6)
    ax2.set_title('mutant')
    #ax1.scatter(data_mut[:,0],data_mut[:,1])
    #ax1.errorbar(data_mut[:,0],data_mut[:,1],yerr=data_mut[:,2])


data_wt = np.loadtxt('out_wt.txt')
data_mut= np.loadtxt('out_mut.txt',)
data_wt = data_wt[~np.isnan(data_wt).any(axis=1)]
data_mut = data_mut[~np.isnan(data_mut).any(axis=1)]
data_wt = data_wt[~np.isinf(data_wt).any(axis=1)]
data_mut = data_mut[~np.isinf(data_mut).any(axis=1)]
wt_y = data_wt[:,1]-np.amax(data_wt[:,1])
mut_y = data_mut[:,1]-np.amax(data_mut[:,1])
#ax1.plot(data_wt[:,0],wt_y,label='tot')
#ax2.plot(data_mut[:,0],mut_y,label='tot')
ax1.legend(loc='best',fontsize=5)
ax2.legend(loc='best',fontsize=5)

#find_hists_overlap('meta_wt.txt','out_wt.txt',ax3)
#find_hists_overlap('meta_mut.txt','out_mut.txt',ax4)

fig1.savefig('test.pdf')

plt.figure(2,figsize=(8,5)) 


mean_begin_wt = int(0.9*len(data_wt[:,0]))
mean_begin_mut = int(0.9*len(data_mut[:,0]))
mean_wt = np.mean(data_wt[mean_begin_wt:,1])
mean_mut = np.mean(data_mut[mean_begin_mut:,1])

#wt_y = data_wt[:,1]-data_wt[-1,1]
#mut_y= data_mut[:,1]-data_mut[-1,1]
wt_y = data_wt[:,1]-mean_wt
mut_y= data_mut[:,1]-mean_mut
plt.plot(data_wt[:,0],wt_y,label='wt')
plt.plot(data_mut[:,0],mut_y,label='mut')
#plt.scatter(data_wt[:,0],wt_y)
#plt.scatter(data_mut[:,0],mut_y)
#plt.errorbar(data_wt[:,0],wt_y,yerr=data_wt[:,2])
#plt.errorbar(data_mut[:,0],mut_y,yerr=data_mut[:,2])
plt.legend(fontsize=6)
plt.xlabel('Distance from Y852 CA to S945L CA $\AA$')
plt.ylabel('Aligned $\Delta$ G kcal/mol') 
plt.savefig('second.pdf')
#plt.show()
