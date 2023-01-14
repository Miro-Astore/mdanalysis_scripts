import numpy as np 
import matplotlib.pyplot as plt 

z_pos = np.load('ion_z_position.npy')
ion_coord_num = np.load('chloride_coodination_number.npy')

z_max = np.amax(z_pos)
z_min = np.amin(z_pos)

bin_spacing = 0.2
n_bins = int(( z_max - z_min ) / 0.2 )

bins = np.linspace(z_min,z_max , n_bins )
bin_counts = np.zeros(len(bins))

for i in range(len(bins)-1): 

    zs = [ x for x in  range(len(z_pos)) if z_pos[x] > bins[i] and z_pos[x] < bins[i+1] ]
    bin_counts[i] = np.sum(ion_coord_num[zs])/len(zs)

print(bin_counts)

plt.plot(bins, bin_counts)
plt.xlim([10,40])

plt.show()
