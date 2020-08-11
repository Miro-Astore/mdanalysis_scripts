from MDAnalysis.analysis.distances import dist
import MDAnalysis as mda
import MDAnalysis.analysis.rms

import numpy as np 
import sys
#####
#USAGE: python mdanalysis_scripts/distance.py TOPFILE traj sel1 sel2 out_name
#
#TODO better detection of lysine bridges, trying to fix argenine detection
#####
print ('TOP ' + str(sys.argv[1]))
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

dist_arr=np.zeros([u.trajectory.n_frames,2])

i=0

#if there is more than 1 residue in each selection do distances things are weird and this code wont work properly.
if ((CAs1.n_atoms==CAs2.n_atoms==1))==False and (sel1.n_atoms>1) and (sel2.n_atoms>1):
    print("doing normal distance measurements")
    for ts in u.trajectory:
       dist1=dist(sel1,sel2) 
       dist_arr[i][1] =  dist1
       dist_arr[i][0] =  ts.time*0.001
       i=i+1
#else its single residue distances and we are going to want to do fancy things.
#want to check N-O distances for salt bridges. Should probably add one for lysine too .
else:
    #names based on charmm might need to change
    if (res_set=={'ARG','GLU'}):
        print("using ARG GLU salt bridge")
        #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
        #trying to do things efficienctly so we are going to combine the two selections, since we have already checked there is only one res in each selection and we can do funky things that would lose generality in other cases. 
        #
        #combine selections into separate residues and pick out their important atoms

        whole_sel_text=sel1txt + ' or ' + sel2txt
        whole_sel=u.select_atoms(whole_sel_text)

        in_sel1=whole_sel.select_atoms('name OE1') 
        in_sel2=whole_sel.select_atoms('name OE2') 
        in_sel3=whole_sel.select_atoms('name NE') 
        in_sel4=whole_sel.select_atoms('name NH1') 
        in_sel5=whole_sel.select_atoms('name NH2') 
        for ts in u.trajectory:
           dist1=dist(in_sel1,in_sel3) 
           dist2=dist(in_sel1,in_sel4) 
           dist3=dist(in_sel1,in_sel5) 
           dist4=dist(in_sel2,in_sel3) 
           dist5=dist(in_sel2,in_sel4) 
           dist6=dist(in_sel2,in_sel5) 
           dist_arr[i][1] =  np.amin ([dist1[-1],dist2[-1],dist3[-1],dist4[-1],dist5[-1],dist6[-1]])
           dist_arr[i][0] =  ts.time*0.001
           i=i+1
    elif (res_set=={'ARG','ASP'}):
        print("using ARG ASP salt bridge")
        #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
        whole_sel_text=sel1txt + ' or ' + sel2txt
        whole_sel=u.select_atoms(whole_sel_text)

        in_sel1=whole_sel.select_atoms('name OD1') 
        in_sel2=whole_sel.select_atoms('name OD2') 
        in_sel3=whole_sel.select_atoms('name NE') 
        in_sel4=whole_sel.select_atoms('name NH1') 
        in_sel5=whole_sel.select_atoms('name NH2') 
        for ts in u.trajectory:
           dist1=dist(in_sel1,in_sel3) 
           dist2=dist(in_sel1,in_sel4) 
           dist3=dist(in_sel1,in_sel5) 
           dist4=dist(in_sel2,in_sel3) 
           dist5=dist(in_sel2,in_sel4) 
           dist6=dist(in_sel2,in_sel5) 
           dist_arr[i][1] =  np.amin ([dist1[-1],dist2[-1],dist3[-1],dist4[-1],dist5[-1],dist6[-1]])
           dist_arr[i][0] =  ts.time*0.001
           i=i+1
    elif (res_set=={'LYS','GLU'}):
        print("using LYS GLU salt bridge")
        #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
        in_sel1=sel1.select_atoms('name OE1 or name NZ ') 
        in_sel2=sel2.select_atoms('name OE2 or name NZ ') 
        for ts in u.trajectory:
           dist1=dist(in_sel1,in_sel2) 
           dist2=dist(in_sel1,in_sel2) 
           dist_arr[i][1] =  np.amin ([dist1[-1],dist2[-1]])
           dist_arr[i][0] =  ts.time*0.001
           i=i+1
        pass
    elif (res_set=={'LYS','ASP'}):
        print("using LYS ASP salt bridge")
        whole_sel_text=sel1txt + ' or ' + sel2txt
        whole_sel=u.select_atoms(whole_sel_text)

        in_sel1=whole_sel.select_atoms('name OD1') 
        in_sel2=whole_sel.select_atoms('name OD2') 
        in_sel3=whole_sel.select_atoms('name NZ') 
        for ts in u.trajectory:
           dist1=dist(in_sel1,in_sel3) 
           dist2=dist(in_sel2,in_sel3) 
           dist_arr[i][1] =  np.amin ([dist1[-1],dist2[-1]])
           dist_arr[i][0] =  ts.time*0.001
           i=i+1
    elif (res_set=={'ASN','ARG'}):
        print("detected ASN ARG interaction")
        whole_sel_text=sel1txt + ' or ' + sel2txt
        whole_sel=u.select_atoms(whole_sel_text)

        in_sel1=whole_sel.select_atoms('name OD1') 
        in_sel2=whole_sel.select_atoms('name NE') 
        in_sel3=whole_sel.select_atoms('name NH1') 
        in_sel4=whole_sel.select_atoms('name NH2') 
        for ts in u.trajectory:
           dist1=dist(in_sel1,in_sel2) 
           dist2=dist(in_sel1,in_sel3) 
           dist3=dist(in_sel1,in_sel4) 
           dist_arr[i][1] =  np.amin ([dist1[-1],dist2[-1],dist3[-1]])
           dist_arr[i][0] =  ts.time*0.001
           i=i+1
    else:
        print("no special cases, doing normal residue to residue distances")
        for ts in u.trajectory:
           dist1=dist(sel1,sel2) 
           dist_arr[i][1] =  dist1[-1]
           dist_arr[i][0] =  ts.time*0.001
           i=i+1
print (dist_arr)
np.savetxt(file_out_name,dist_arr)
