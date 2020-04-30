import MDAnalysis as mda
import pprint
from multiprocessing import Pool
import numpy as np
import glob, os

#currently this script works by parralelising over residues but this causes way too many disk acesses. we should chunk over the trajectory. 

res_list = [95 134 153 190 248 251 254 303 334 335 347 352 968 975 978 1030 1041 1048 1097 1158 1162 1165]

#TOP=str(sys.argv[1])
#TRAJ=str(sys.argv[2])
TOP='../_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_1.psf'
TRAJ'../_scratch_r16_ma2374_gmx_cftr_2nd_round_310K_wt_1_sum.xtc'
n_cores=multiprocessing.cpu_count()

u=mda.Universe(TOP,TRAJ)

ntasks=n_cores

#distribute frames equally over cores 
def distribute_frames_amongst_cores(num_frames, n_cores):
   int_block=int(num_frames/n_cores)
   diff=len(arr)-(int_block)*n_cores
   blocks=[int_block]*n
   blocks=np.array(blocks)
   blocks[0:diff]=blocks[0:diff]+1
   return blocks

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def track_occ(frame_arr):
    cls = 0 
    for us in u.trajectory[frame_arr]:
        for res in res_list:
            seltext  = str("name CLA and around 5.0 (protein and resid " + str(resid) + ") " )
            cls_sel = u.select_atoms(seltext)
            cls=cls+len(cls_sel)
    return cls

dirs=np.loadtxt('dir_list',dtype=str)
for place in dirs:
    res_occ = []
        
    if __name__ == '__main__':
        for c in list(chunks(res_list,ntasks)):
            print (c)
            with Pool(ntasks,maxtasksperchild=1000) as p:
                res_occ = [res_occ, (p.map(track_occ,c))]
    print (res_occ)

