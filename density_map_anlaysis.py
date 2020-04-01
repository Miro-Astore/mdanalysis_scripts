import MDAnalysis.analysis.density as de 
import numpy as np 
import MDAnalysis as mda
import matplotlib.pyplot as plt

def sample_cylinder_coords(cylinder_center,cylinder_radius, sample_points,rad_scale=0.9): 
    angles=np.array(2*np.pi/(sample_points-1) * np.arange(sample_points-1))
    points = np.zeros([sample_points,3])
    points[0,:] = cylinder_center
    for i in range(1,sample_points):
        x=points[0,0] + cylinder_radius*rad_scale*np.cos(angles[i-1])
        y=points[0,1] + cylinder_radius*rad_scale*np.sin(angles[i-1])
        z=points[0,2] #no change in elevation
        points[i,:] = [x,y,z] 
    return points

    

#method. Create a cylinder centered at the hole sampled site with the radius of the probe detected radius
#place 7 points in the cylinder, 1 in the center 6 at the edges (a little in radius=proberadius*0.90 so we don't clash with surrounding protein pore)
#mess with the above number of sample points to see if you have artefacts for sampling around the edges of the probe.
#average over readings in slice
#add up slices
#make sure you can graph density over z so you can look at it.

u=mda.Universe('../I37R_last.sph',topology_format='PDB')

interp_values=np.zeros(len(u.atoms))
cylinder_height=np.abs(u.atoms.positions[0][2]-u.atoms.positions[1][2])
spacing=0.2
stride=np.int(np.round(spacing/cylinder_height))
cylinder_height=stride*cylinder_height
stride_pos=u.atoms.positions[0:-1:stride]
stride_rad=u.atoms.tempfactors[0:-1:stride]
big_sample=np.zeros(len(stride_pos))

x=de.Density('../I37R_test_doubt.dx')
num_sample_points=10
y=(x.centers())

for i in range(len(u.atoms)):
    pos=u.atoms.positions[i]
    rad=u.atoms.tempfactors[i]
    interp_values[i]=(x.interpolated(pos[0],pos[1],pos[2]))



for i in range(len(stride_pos)):

    pos=stride_pos[i]
    rad=stride_rad[i]

    grid_sample_values=np.zeros(num_sample_points)
    sample_grid_points=sample_cylinder_coords(pos,rad,num_sample_points)


    for j in range(num_sample_points):
    #=x.interpolated(
        grid_sample_values[j] = x.interpolated(sample_grid_points[j][0],sample_grid_points[j][1],sample_grid_points[j][2])

    big_sample[i]=2*np.pi*rad**2*cylinder_height*np.mean(grid_sample_values)


z=u.atoms.positions[:,2]
sort_arr=np.argsort(z)
z=z[sort_arr]
interp_values=interp_values[sort_arr]

sort_arr=np.argsort(stride_pos[:,2])
stride_z=stride_pos[:,2][np.argsort(stride_pos[:,2])]
big_sample=big_sample[sort_arr]
paired_vol=np.array([stride_z,big_sample])
np.save('../I37R_data.npy',paired_vol)
#
#plt.plot(z,interp_values,label='linear')
#plt.legend()
#plt.figure()
#plt.plot(stride_z,big_sample,label='volume')
#plt.legend()
#plt.show()


