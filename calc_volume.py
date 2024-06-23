import MDAnalysis as mda
import tqdm
import argparse
import jax
import jax.numpy as  jnp
import pdb
import rlcompleter

pdb.Pdb.complete=rlcompleter.Completer(locals()).complete

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument ('--topology',type=str) 
parser.add_argument ('--trajectory',type=str) 
parser.add_argument ('--stride',type=int) 
parser.add_argument ('--spacing',default=0.3,type=float) 
parser.add_argument ('--n_batches',default = 10000,type=int) 
parser.add_argument ('--selection',default = 'all',type=str) 
parser.add_argument ('-o',type=str) 
args = parser.parse_args()

u = mda.Universe(args.topology, args.trajectory)

selected_atoms = u.select_atoms(args.selection)
volumes = np.zeros(len(range(u.trajectory.n_frames)[::args.stride]))
# Calculate grid dimensions


element_list = [first_letter[0] for first_letter in selected_atoms.names]  
size_dictionary = {'H':1.2, 'C':1.7, 'O':1.52, 'N':1.55 , 'S':1.80, 'P':1.8 }
#radius_array = jnp.array([3.26]*len(atoms_positions))
radius_array = jnp.array([size_dictionary[element] for element in element_list],dtype=jnp.float32)
max_atom_radius  = jnp.amax(radius_array)

@jax.jit
def calc_dist (grid_point, atom_positions,radius_array):
    tiled_array = jnp.array(jnp.tile(grid_point,(len(atom_positions),1)))
    distance_array = jnp.linalg.norm(atom_positions-tiled_array,axis=1)
    value = jnp.any(jnp.less(distance_array,radius_array))
    return value

#def calculate_volume(atoms_positions, radius):

atoms_positions = jnp.array(selected_atoms.positions,dtype=jnp.float32)

"""
Calculate the volume occupied by a set of atoms using a 3D grid approximation.

Parameters:
- atoms_positions (np.array): Array of shape (N, 3) containing positions of atoms.
- max_atom_radius (float): Radius within which to count grid cubes (Angstrom).

Returns:
- volume (float): Estimated volume occupied by the atoms (Angstrom^3).
"""

# Define the grid spacing (smaller than max_atom_radius for accuracy)
grid_spacing = args.spacing # angstrom
# Determine the bounding box for the atoms
x_min_coords = np.min(atoms_positions[:,0], axis=0)-max_atom_radius
x_max_coords = np.max(atoms_positions[:,0], axis=0)+max_atom_radius
y_min_coords = np.min(atoms_positions[:,1], axis=0)-max_atom_radius
y_max_coords = np.max(atoms_positions[:,1], axis=0)+max_atom_radius
z_min_coords = np.min(atoms_positions[:,2], axis=0)-max_atom_radius
z_max_coords = np.max(atoms_positions[:,2], axis=0)+max_atom_radius

n_x_grid_points = np.ceil(-(x_min_coords - x_max_coords)/grid_spacing ).astype(int)
n_y_grid_points = np.ceil(-(x_min_coords-x_max_coords)/grid_spacing ).astype(int)
n_z_grid_points = np.ceil(-(z_min_coords - z_max_coords)/grid_spacing ).astype(int)

grid_shape = np.array([n_x_grid_points, n_y_grid_points, n_z_grid_points])

x_coor = x_min_coords + (np.arange(n_x_grid_points)*grid_spacing)
y_coor = y_min_coords + (np.arange(n_y_grid_points)*grid_spacing)
z_coor = z_min_coords + (np.arange(n_z_grid_points)*grid_spacing)
X, Y, Z = np.meshgrid(x_coor, y_coor, z_coor, indexing='ij')

grid_points = jnp.float32(jnp.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=-1))


n_grid_points = np.prod(grid_shape)
print('There are ' + str(n_grid_points) + ' points in this grid.' )



# Create a grid to count contained grid cubes

n_batches = args.n_batches
batch_size = np.int64(n_grid_points/n_batches)

if (n_grid_points) % batch_size == 0:
    n_batches = ((n_grid_points)//batch_size)
else:
    n_batches = ((n_grid_points)//batch_size) + 1

batch = jax.vmap(calc_dist)
#bin_indexes=np.zeros()
for frame_index in range(u.trajectory.n_frames):
    running_volume = 0 
    u.trajectory[frame_index]

    for i in tqdm.tqdm(range(n_batches)):
        batch_grid_positions = grid_points[batch_size*i: (i+1) *batch_size]
        tiled_atoms = jnp.tile(atoms_positions, (len(batch_grid_positions),1,1))
        tiled_radius_array = jnp.tile(radius_array, (len(batch_grid_positions),1))

        result = batch(batch_grid_positions, tiled_atoms, tiled_radius_array)

        running_volume = running_volume + (jnp.sum(result) * grid_spacing ** 3 )

    volumes[frame_index] = running_volume

print(volumes)
np.save(args.o,volumes)
#for i in tqdm.tqdm(range(n_grid_points)):
#    current_grid_point = grid_points[i]
#    #print(current_grid_point)
#    #thing = [ i   for j in atoms_positions  if  np.linalg.norm(j-current_grid_point ) <  max_atom_radius ]
#    #found = any(np.linalg.norm(j - current_grid_point) < max_atom_radius for j in atoms_positions)
#    if calc_dist(current_grid_point,atoms_positions,radius_array) == True:
#        inside_grid[i]=True
#
#
#
#print(np.sum(inside_grid) * grid_spacing**3)
#    #inside_points[i] = [True if np.sum()
#            
#
##    # Mark grid cubes that are within max_atom_radius of any atom
##    for x in range(max(0, grid_coords[0] - r), min(grid_shape[0], grid_coords[0] + r + 1)):
##        for y in range(max(0, grid_coords[1] - r), min(grid_shape[1], grid_coords[1] + r + 1)):
##            for z in range(max(0, grid_coords[2] - r), min(grid_shape[2], grid_coords[2] + r + 1)):
##
##                x_dists = (atom_positions[:,0] - grid_coords[0]) * grid_spacing
##                print(x_dists)
##                y_dists = (y - grid_coords[1]) * grid_spacing
##                z_dists = (z - grid_coords[2]) * grid_spacing
##                norm_dist = calc_dist(x_dist,y_dist,z_dist, atom_number)
##                if norm_dist <= max_atom_radius:
##                    grid[x, y, z] = True
##    for atom_pos in atoms_positions:
##        print(grid_coords)
##        r = int(np.ceil(max_atom_radius / grid_spacing))
##
##
##    
##    # Calculate the volume based on the number of contained grid cubes
##    volume = np.sum(grid) * (grid_spacing ** 3)
##    
##    return volume 
