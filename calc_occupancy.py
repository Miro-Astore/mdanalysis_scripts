import MDAnalysis as mda
import pprint
import sys
import os 
from multiprocessing import Pool
import multiprocessing
import numpy as np
import glob, os

#currently this script works by parralelising over residues but this causes way too many disk acesses. we should chunk over the trajectory. 


#distribute frames equally over cores 
def distribute_frames_amongst_cores(num_frames, n_cores):
    #what the actual hell is this code
   #you've found the largest number that is less than num_frames and still divisible by n_cores
   #found the difference between n_frames and this number
   int_block=int(num_frames/n_cores)
   diff=num_frames-(int_block)*n_cores
   blocks=[int_block]*n_cores
   blocks=np.array(blocks)
   #now add 1 to some blocks to absorb the difference so we don't miss a bunch of frames at the end. 
   if diff!=0:
           blocks[0:diff]=blocks[0:diff]+1

   #these weird indexes are gonna do my head in why did i think this was going to be a good idea.
   inds=np.zeros([n_cores,2])
   inds[0,:]=[0,blocks[1]]
   for i in range(1,n_cores):
           inds[i,:]=[inds[i-1][1]+1,inds[i-1][1]+blocks[i]]
   inds[-1][-1]= inds[-1][-1]-1
   inds=[[int (x[0]),int(x[1])] for x in inds]

   return inds

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def track_occ(frame_arr):
    res_occ_temp = np.zeros(num_res) 
    u=mda.Universe(TOP,TRAJ)

    for us in u.trajectory[frame_arr[0]:frame_arr[1]]:
        print(u.trajectory.time)
        for i in range(num_res):
            seltext  = str("name CLA and around 5.0 (protein and resid " + str(res_list[i]) + ") " )
            cls_sel = u.select_atoms(seltext)
            res_occ_temp[i]=res_occ_temp[i]+len(cls_sel)
    #print(res_occ_temp[0])
    return res_occ_temp

        
if __name__ == '__main__':
    #res_list = [95, 134]
    res_list = np.array([95, 134, 153, 190, 248, 251, 254, 303, 334, 335, 347, 352, 370, 968, 975, 978, 1030, 1041, 1048, 1097, 1158, 1162, 1165])
    num_res=len(res_list)

    TOP=str(sys.argv[1])
    TRAJ=str(sys.argv[2])
    out_name=str(sys.argv[3])
    #TOP='../_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_1.psf'
    #TRAJ='../test.xtc'
    #out_name=str('outname')
    n_cores=multiprocessing.cpu_count()

    u=mda.Universe(TRAJ)
    num_frames=u.trajectory.n_frames
    blocks=distribute_frames_amongst_cores(num_frames,n_cores)

    ntasks=n_cores
    frame_inds=distribute_frames_amongst_cores(num_frames,ntasks)
    i=0
    res_occ=np.zeros([n_cores,num_res])
    with Pool(processes=n_cores) as p:

        x = (p.map(track_occ,frame_inds))
    res_occ=x
    res_occ=np.sum(res_occ,axis=0)
    res_occ=np.append(res_occ,num_frames)
    res_list=np.append(res_list,len(res_list))
    res_occ=np.array([res_list,res_occ])
    #res_occ=np.vstack(res_occ)

    print (res_occ)
    np.save(out_name,res_occ)

