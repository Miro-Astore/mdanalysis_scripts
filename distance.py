from MDAnalysis.analysis.distances import dist
import MDAnalysis as mda
import MDAnalysis.analysis.rms

import numpy as np 
import sys

print ('psf ' + str(sys.argv[1]))
print ('traj ' + str(sys.argv[2]))
print ('sel1 ' + str(sys.argv[3]))
print ('sel2 ' + str(sys.argv[4]))
print ('out_name ' + str(sys.argv[5]))
u=mda.Universe(sys.argv[1],sys.argv[2])

file_out_name = str(sys.argv[5])

sel1txt=str(sys.argv[3])
sel2txt=str(sys.argv[4])

sel1=u.select_atoms(sel1txt)
sel2=u.select_atoms(sel2txt)
CAs1=sel1.select_atoms('name CA and protein')
CAs2=sel2.select_atoms('name CA and protein')
res_set={CAs1.resnames[0],CAs2.resnames[0]}
print(res_set)

dist_arr=np.zeros(u.trajectory.n_frames)

i=0

#if there is more than 1 residue in each selection do distances as normal.
if ((CAs1.n_atoms==CAs2.n_atoms==1))==False and (sel1.n_atoms>1) and (sel2.n_atoms>1):
    pass 
#else its single residue distances and we are going to want to do fancy things.
#want to check N-O distances for salt bridges. Should probably add one for lysine too .
else:
    #names based on charmm might need to change
    if (res_set=={'ARG','GLU'}):
        #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
        in_sel1=sel1.select_atoms('name OE1 or name NH1') 
        in_sel2=sel2.select_atoms('name OE1 or name NH1') 
        in_sel3=sel1.select_atoms('name OE2 or name NH2') 
        in_sel4=sel2.select_atoms('name OE2 or name NH2') 
        for ts in u.trajectory:
           dist1=dist(in_sel1,in_sel2) 
           dist2=dist(in_sel1,in_sel4) 
           dist3=dist(in_sel3,in_sel2) 
           dist4=dist(in_sel3,in_sel4) 
           dist_arr[i] =  np.amin ([dist1[-1],dist2[-1],dist3[-1],dist4[-1]])
           i=i+1
    elif (res_set=={'ARG','ASP'}):
        #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
        in_sel1=sel1.select_atoms('name OD1 or name NH1') 
        in_sel2=sel2.select_atoms('name OD1 or name NH1') 
        in_sel3=sel1.select_atoms('name OD2 or name NH2') 
        in_sel4=sel2.select_atoms('name OD2 or name NH2') 
        for ts in u.trajectory:
           dist1=dist(in_sel1,in_sel2) 
           dist2=dist(in_sel1,in_sel4) 
           dist3=dist(in_sel3,in_sel2) 
           dist4=dist(in_sel3,in_sel4) 
           dist_arr[i] =  np.amin ([dist1[-1],dist2[-1],dist3[-1],dist4[-1]])
           i=i+1
    elif (res_set=={'LYS','GLU'}):
        #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
        in_sel1=sel1.select_atoms('name OE1 or name NZ ') 
        in_sel2=sel2.select_atoms('name OE2 or name NZ ') 
        for ts in u.trajectory:
           dist1=dist(in_sel1,in_sel2) 
           dist2=dist(in_sel1,in_sel2) 
           dist_arr[i] =  np.amin ([dist1[-1],dist2[-1],dist3[-1],dist4[-1]])
           i=i+1
        pass
    elif (res_set=={'LYS','ASP'}):
        pass
print (dist_arr)
np.savetxt(file_out_name,dist_arr)
