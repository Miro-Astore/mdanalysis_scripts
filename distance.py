from MDAnalysis.analysis.distances import dist
import MDAnalysis as mda

import numpy as np 
import sys
import argparse
#####
#USAGE: python mdanalysis_scripts/distance.py TOPFILE traj out_file sel1 sel2 sel3 sel4 ...
#
#####

parser = argparse.ArgumentParser()

parser.add_argument('--topology', '-top', dest='top_file' , help='System topology file') 
parser.add_argument('--trajectory', '-traj', dest='traj_file', help='Simulation trajectory file') 
parser.add_argument('--stride', '-st', dest='stride' , default = 1, help='Stride through trajectory skipping this many frames.') 
parser.add_argument('--out', '-o', dest='out_file', default='distances.dat', help='Output file') 
parser.add_argument('--selections', '-s', dest='selection_strings', help='Selection strings with the syntax of MDAnalysis', nargs="+") 
parser.add_argument('--begin', '-b', dest='beginning_frame', help='Frame from which to begin analysis.',type=int) 
parser.add_argument('--end', '-e', dest='ending_frame', help='Frame at which to end analysis.',type=int) 
parser.add_argument('--first_n_frames', '-f', dest='first_n_frames', help='Analyse the first n frames of the trajectory.',type=int) 
parser.add_argument('--last_n_frames', '-l', dest='last_n_frames', help='Analyse the last n frames of the trajectory.',type=int) 

args = parser.parse_args()


#make pairs of selection strings
sel_string_pairs_list =  []

n_pairs = len(range(0,len(args.selection_strings),2))
for i in range(0,len(args.selection_strings),2):
    print ('pairs: ')
    print('\"' + args.selection_strings[i] + '\"-\"' + args.selection_strings[i+1] + '\"')
    sel_string_pairs_list.append([args.selection_strings[i], args.selection_strings[i+1]])

#make universe object
u=mda.Universe(args.top_file,args.traj_file)

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

#this was so much more work than it should have been. 
    
if args.first_n_frames is None and args.beginning_frame is None and args.ending_frame is None and args.last_n_frames is None : 

    print(1)
    first_frame = 0
    last_frame = -1

elif args.first_n_frames is not None and args.beginning_frame is None and args.ending_frame is None and args.last_n_frames is None : 

    print(2)
    first_frame = 0
    last_frame = args.first_n_frames

elif args.first_n_frames is None and args.beginning_frame is not None and args.ending_frame is None and args.last_n_frames is None : 

    print(3)
    first_frame = args.beginning_frame
    last_frame = -1
    
elif args.first_n_frames is None and args.beginning_frame is None and args.ending_frame is None and args.last_n_frames is not None : 

    print(4)
    first_frame = u.trajectory.n_frames - args.last_n_frames
    last_frame = -1

elif args.first_n_frames is None and args.beginning_frame is not None and args.ending_frame is not None and args.last_n_frames is None : 

    print(5)
    first_frame = args.beginning_frame
    last_frame = args.ending_frame

elif args.first_n_frames is None and args.beginning_frame is None and args.ending_frame is not None and args.last_n_frames is None : 

    print(6)
    first_frame = 0
    last_frame = args.ending_frame

    
elif args.first_n_frames is None and args.beginning_frame is None and args.ending_frame is not None and args.last_n_frames is not None : 

    print(7)
    first_frame = u.trajectory.n_frames - args.last_n_frames
    last_frame = args.ending_frame

elif args.first_n_frames is not None and args.beginning_frame is not None and args.ending_frame is None and args.last_n_frames is None : 
    
    print(8)
    first_frame = args.beginning_frame
    last_frame = args.first_n_frames

else:

    raise UserWarning('You need to make sure that the right set of beginning and ending flags are specified. Something is wrong. For example, have you tried to tell the script to analyse both the first m frames and the last m frames? This would create a disjoint interval. ')


dist_arr=np.zeros([len(u.trajectory[first_frame:last_frame:int(args.stride)]),n_pairs+1])
dist_array_row=0

for ts in u.trajectory[first_frame:last_frame:int(args.stride)]:
    dist_arr[dist_array_row][0] =  ts.time*0.001
    for i in range(n_pairs):

            
        #else its single residue distances and we are going to want to do fancy things.
        #want to check N-O distances for salt bridges. Should probably add one for lysine too .
            #names based on charmm might need to change

        if (pair_type_dict[i] == "simple_pair"):

            in_sel1=sel_objects_pairs_list[i][0]
            in_sel2=sel_objects_pairs_list[i][1]
            dist_arr[dist_array_row][i+1] = dist(in_sel1,in_sel2)[-1]

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

            dist_arr[dist_array_row][i+1] =  np.amin ([dist1[-1],dist2[-1],dist3[-1],dist4[-1],dist5[-1],dist6[-1]])
            continue

        elif (pair_type_dict[i] == "KE_SB" or pair_type_dict[i] == "KD_SB"):
            #since we dont know which selection was passed as the argenine we can just pair up atoms arbitrarily so long as the logic makes sense.
            in_sel1=sel_objects_pairs_list[i][0]
            in_sel2=sel_objects_pairs_list[i][1]
            in_sel3=sel_objects_pairs_list[i][2]

            dist1=dist(in_sel1,in_sel2) 
            dist2=dist(in_sel1,in_sel3) 

            dist_arr[dist_array_row][i+1] =  np.amin ([dist1[-1],dist2[-1]])
            continue


        elif (pair_type_dict[i] == "RN_SB" or pair_type_dict[i] == "RQ_SB"):

           dist1=dist(in_sel1,in_sel2) 
           dist2=dist(in_sel1,in_sel3) 
           dist3=dist(in_sel1,in_sel4) 
           dist_arr[dist_array_row][i+1] =  np.amin ([dist1[-1],dist2[-1],dist3[-1]])
           continue

    dist_array_row = dist_array_row + 1
#printing results and saving to file
header_str = 'time, '
print (dist_arr)
for i in range(0,len(args.selection_strings),2):
    header_str = header_str + '\"' + str(args.selection_strings[i]) +  '\"-\"' + str(args.selection_strings[i+1]) + '\", '
header_str = header_str [:-2]

np.savetxt(args.out_file,dist_arr,header = header_str)
