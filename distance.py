from MDAnalysis.analysis.distances import dist
import MDAnalysis as mda

import numpy as np 
import sys
#####
#USAGE: python mdanalysis_scripts/distance.py TOPFILE traj out_file sel1 sel2 sel3 sel4 ...
#
#####
print ('TOP ' + str(sys.argv[1]))
print ('traj ' + str(sys.argv[2]))
print ('out_name ' + str(sys.argv[3]))
if len(sys.argv) % 2 != 0:
    raise ValueError ("wrong number of arguments, make sure you have topology, trajectory, output file and an even number of selection strings")

sel_string_pairs_list =  []

#make pairs of selection strings
for i in range(4,len(sys.argv),2):
    
    n_pairs = len(range(4,len(sys.argv),2))

    print('sel' + str (i) + ': ' + str(sys.argv[i]) +  ', sel' + str (i+1) + ': ' + str(sys.argv[i+1]) + ' ')
    sel_string_pairs_list.append([str(sys.argv[i]), str(sys.argv[i+1])])


#make universe object
u=mda.Universe(sys.argv[1],sys.argv[2])

file_out_name = str(sys.argv[3])

#make pairs of selection objects
sel_objects_pairs_list = [] 
for i in range(n_pairs):
    sel_pair_1 = u.select_atoms(sel_string_pairs_list[i][0])
    sel_pair_2 = u.select_atoms(sel_string_pairs_list[i][1])
    sel_objects_pairs_list.append([sel_pair_1,sel_pair_2])

#make a dictionary to save the type of selection we have for each pair
pair_type_dict = {}
for i in range(n_pairs):
    sel1 = sel_objects_pairs_list[i][0]
    sel2 = sel_objects_pairs_list[i][1]
    CAs1=sel1.select_atoms('name CA and protein')
    CAs2=sel2.select_atoms('name CA and protein')

    whole_sel_text= '(' + sel_string_pairs_list[i][0] + ')' + ' or ' + '(' +  sel_string_pairs_list[i][1] + ')'
    whole_sel=u.select_atoms(whole_sel_text)
    res_set={sel1.resnames[0],sel2.resnames[0]}
    if sel1.n_atoms == 1 and sel2.n_atoms == 1:
        pair_type_dict[i] = "simple_pair"
        continue
    elif (res_set=={'ARG','GLU'}) and (sel1.n_atoms > 1 and sel2.n_atoms > 1) and (CAs1.n_atoms == CAs2.n_atoms == 1):
        #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
        #trying to do things efficienctly so we are going to combine the two selections, since we have already checked there is only one res in each selection and we can do funky things that would lose generality in other cases. 
        #
        #combine selections into separate residues and pick out their important atoms

        in_sel1=whole_sel.select_atoms('name OE1') 
        in_sel2=whole_sel.select_atoms('name OE2') 
        in_sel3=whole_sel.select_atoms('name NE') 
        in_sel4=whole_sel.select_atoms('name NH1') 
        in_sel5=whole_sel.select_atoms('name NH2') 

        sel_objects_pairs_list [i] = [in_sel1, in_sel2, in_sel3, in_sel4, in_sel5] 

        pair_type_dict[i] = "RE_SB"
        continue
    elif (res_set=={'ARG','ASP'}) and (sel1.n_atoms > 1 and sel2.n_atoms > 1) and (CAs1.n_atoms == CAs2.n_atoms == 1):
        #similar situation to the above for the aspartate

        in_sel1=whole_sel.select_atoms('name OD1') 
        in_sel2=whole_sel.select_atoms('name OD2') 
        in_sel3=whole_sel.select_atoms('name NE') 
        in_sel4=whole_sel.select_atoms('name NH1') 
        in_sel5=whole_sel.select_atoms('name NH2') 

        sel_objects_pairs_list [i] = [in_sel1, in_sel2, in_sel3, in_sel4, in_sel5] 

        pair_type_dict[i] = "RD_SB"
        continue
    elif (res_set=={'LYS','GLU'}) and (sel1.n_atoms > 1 and sel2.n_atoms > 1) and (CAs1.n_atoms == CAs2.n_atoms == 1):

        in_sel1=whole_sel.select_atoms('name NZ') 
        in_sel2=whole_sel.select_atoms('name OE1') 
        in_sel3=whole_sel.select_atoms('name OE2') 
        sel_objects_pairs_list [i] = [in_sel1, in_sel2, in_sel3] 

        pair_type_dict[i] = "KE_SB"
        continue
    elif (res_set=={'LYS','ASP'}) and (sel1.n_atoms > 1 and sel2.n_atoms > 1) and (CAs1.n_atoms == CAs2.n_atoms == 1):

        in_sel1=whole_sel.select_atoms('name NZ') 
        in_sel2=whole_sel.select_atoms('name OD1') 
        in_sel3=whole_sel.select_atoms('name OD2') 
        sel_objects_pairs_list [i] = [in_sel1, in_sel2, in_sel3] 
        pair_type_dict[i] = "KD_SB"

    elif (res_set=={'ARG','ASN'}) and (sel1.n_atoms > 1 and sel2.n_atoms > 1) and (CAs1.n_atoms == CAs2.n_atoms == 1):

        in_sel1=whole_sel.select_atoms('name OD1') 
        in_sel2=whole_sel.select_atoms('name NE') 
        in_sel3=whole_sel.select_atoms('name NH1') 
        in_sel4=whole_sel.select_atoms('name NH2') 
        
        pair_type_dict[i] = "RN_SB"

        continue
    else:
        raise ValueError ('the selection pair: \"' + str(sel_string_pairs_list[i][0]) + '\", \"' + str(sel_string_pairs_list[i][1]) + '\" has a strange number of atoms. You should only have 1 atom in each selection string or a pair of residues that forms a salt bridge. Check the selection pairs.' )
        

print(pair_type_dict)
dist_arr=np.zeros([u.trajectory.n_frames,n_pairs+1])

frame=0

for ts in u.trajectory:
    dist_arr[frame][0] =  ts.time*0.001
    for i in range(n_pairs):

            
        #else its single residue distances and we are going to want to do fancy things.
        #want to check N-O distances for salt bridges. Should probably add one for lysine too .
            #names based on charmm might need to change

        if (pair_type_dict[i] == "simple_pair"):

            in_sel1=sel_objects_pairs_list[i][0]
            in_sel2=sel_objects_pairs_list[i][1]
            dist_arr[frame][i+1] = dist(in_sel1,in_sel2)[-1]

        elif (pair_type_dict[i] == "RE_SB" or pair_type_dict[i] == "RD_SB"):
            #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
            #trying to do things efficienctly so we are going to combine the two selections, since we have already checked there is only one res in each selection and we can do funky things that would lose generality in other cases. 

            #combine selections into separate residues and pick out their important atoms


            in_sel1=sel_objects_pairs_list[i][0]
            in_sel2=sel_objects_pairs_list[i][1]
            in_sel3=sel_objects_pairs_list[i][2]
            in_sel4=sel_objects_pairs_list[i][3]
            in_sel5=sel_objects_pairs_list[i][4]

            dist1=dist(in_sel1,in_sel3) 
            dist2=dist(in_sel1,in_sel4) 
            dist3=dist(in_sel1,in_sel5) 
            dist4=dist(in_sel2,in_sel3) 
            dist5=dist(in_sel2,in_sel4) 
            dist6=dist(in_sel2,in_sel5) 

            dist_arr[frame][i+1] =  np.amin ([dist1[-1],dist2[-1],dist3[-1],dist4[-1],dist5[-1],dist6[-1]])
            continue

        elif (pair_type_dict[i] == "KE_SB" or pair_type_dict[i] == "KD_SB"):
            #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
            in_sel1=sel_objects_pairs_list[i][0]
            in_sel2=sel_objects_pairs_list[i][1]
            in_sel3=sel_objects_pairs_list[i][2]

            dist1=dist(in_sel1,in_sel2) 
            dist2=dist(in_sel1,in_sel3) 

            dist_arr[frame][i+1] =  np.amin ([dist1[-1],dist2[-1]])
            continue


        elif (pair_type_dict[i] == "RN_SB" or pair_type_dict[i] == "RQ_SB"):

           dist1=dist(in_sel1,in_sel2) 
           dist2=dist(in_sel1,in_sel3) 
           dist3=dist(in_sel1,in_sel4) 
           dist_arr[frame][i+1] =  np.amin ([dist1[-1],dist2[-1],dist3[-1]])
           continue

    frame = frame + 1
#printing results and saving to file
header_str = 'time, '


print (dist_arr)
for i in range(4,len(sys.argv),2):
    header_str = header_str + '\"' + str(sys.argv[i]) +  '\"-\"' + str(sys.argv[i+1]) + '\", '
header_str = header_str [:-2]
print(header_str)
np.savetxt(file_out_name,dist_arr,header = header_str)
