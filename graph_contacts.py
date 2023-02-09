import numpy as np
import itertools
from functools import reduce
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

root_dir = '/mnt/home/mastore/ceph/unbiased_trpv1_first_round_no_cter/trpv1_soybean_nanodisc/'

num_replicates = 10
files_arr = ['A_B_ARD_contacts.dat','B_C_ARD_contacts.dat','C_D_ARD_contacts.dat','D_A_ARD_contacts.dat']

legend_arr = [ 'contact 1', 'contact 2','contact 3','contact 4' ]
color_arr = ['tab:blue','pink','red', 'yellow'] 

plt.figure(figsize=(10,10))
temp_max=0

def factors(n):    
    return set(reduce(list.__add__, ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def choose_most_square (n) :
    target = np.sqrt(n)

facts = np.sort(list(factors(num_replicates)))


#determining rows and cols arrangement, more annoying than you'd think
num_plots_rows=0
num_plots_cols=0

#if a square number we'll have an odd number of factors so set number of rows and cols to the same thing 
if (len(facts) % 2) == 0:
    num_plots_rows=facts[int(len(facts)/2)]
    num_plots_cols=facts[int(len(facts)/2)-1]
else:
    num_plots_rows=facts[int((len(facts)-1)/2)]
    num_plots_cols=facts[int((len(facts)-1)/2)]

print(num_plots_rows)
print(num_plots_cols)

#make grid for viewing convergence
#plt.figure(figsize=(10,4))
#AX.update (wspace=0.4,hspace=0.8)

row_arr=range(num_plots_rows)
col_arr=range(num_plots_cols)

row_arr_t=list(itertools.chain.from_iterable(itertools.repeat(x, num_plots_cols) for x in row_arr))
row_arr=row_arr_t

col_arr=list(col_arr)*num_plots_rows

print(row_arr)
print(col_arr)

AX=gridspec.GridSpec(num_plots_rows,num_plots_cols)


for j in range(num_replicates):
    for i in range(len(files_arr)):
        curr_file = root_dir + str(j+1) + '/' + files_arr [i]

        data=np.loadtxt(curr_file)[0::,:]
        a_max=np.max(data[:,1])

        if temp_max < a_max:
            temp_max=a_max
    
    


for j in range(num_replicates):
    row_place=row_arr[j] 
    col_place=col_arr[j] 
    plt.subplot2grid((num_plots_rows,num_plots_cols),(row_place,col_place),colspan=1,rowspan=1)
    for i in range(len(files_arr)):

        curr_file = root_dir + str(j+1) + '/' + files_arr [i]

        data=np.loadtxt(curr_file)

        data=np.array(data)
        data=data.reshape((len(data),2))

        x=data[0:-1,0]
        y=data[0:-1,1]
        
        y = moving_average(y,10)

        x = x[-len(y):]

        plt.plot(x,y,label=legend_arr[i],color=color_arr[i], linewidth=2.0)
    plt.xlim([0,np.max(x)])
    plt.ylim([0,temp_max])

plt.ylabel('#ARD contacts',fontsize=16)
plt.xlabel('time (ns)',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(prop={'size':12})


plt.tight_layout()
plt.show()
