import mdtraj as md
import MDAnalysis as mda
import torch
import copy
import argparse
import sys
from torch_batch_svd import svd
import numpy as np
from tqdm import tqdm
from multiprocessing.pool import ThreadPool
import time
import pdb
import rlcompleter
pdb.Pdb.complete=rlcompleter.Completer(locals()).complete
import os

cuda = torch.device('cuda')

parser = argparse.ArgumentParser()

parser.add_argument('--topology', '-t', dest='topology' , help='topology file',type=str) 
parser.add_argument('--trajectory', '-j', dest='traj', help='trajectory file',type=str) 
parser.add_argument('--out', '-o', dest='out_file',default="rmsd_2d.npy", help='Output file',type=str) 
parser.add_argument('--stride', '-dt', dest='stride' , default=None, help='Stride through trajectory skipping this many frames.',type=int) 
parser.add_argument('--chunk', '-c', dest='chunk' , default=None, help='Chunk size for what to load into the gpu.',type=int) 
parser.add_argument('--selection', '-s', dest='selection_string', help='Selection string for fitting. Will be applied to both target and reference structures.',type=str, default="name CA") 

args = parser.parse_args()

chunk_size = args.chunk
traj_file = args.traj

#top_only = md.load_frame(traj_file,index = 0, top=args.topology)
top_only = mda.Universe(args.topology, n_atoms = 0)

#top_only = md.load('test_100_frames.xtc' , top = "prot.psf")
chunkTrajA = md.iterload(traj_file, top=args.topology, chunk = chunk_size, stride = args.stride)
print('loaded once')

#this is a clever way to copy the trajectory itterators so we dont have to load them from scratch all the time.
chunkTrajA_list = list(chunkTrajA)
chunkTrajA = (item for item in chunkTrajA_list)
chunkTrajB = (item for item in chunkTrajA_list)
chunkTraj_backup_list = list((item for item in chunkTrajA_list))
print('trajs copied')

print('copied')
print(chunkTrajB)


# use the atom selection from mdanalysis because it's a million times better than mdtraj's
#atom_indices = [a.index for a in top_only.topology.atoms if a.name == 'CA' and a.residue.resSeq > 105 ]
sel_object = top_only.select_atoms(args.selection_string)
atom_indices = sel_object.indices
#atom_indices = [a.index for a in top_only.topology.atoms if a.name == 'CA' and a.residue.resSeq > 105 ]
print(len(atom_indices),flush=True)
n_atoms = len(atom_indices)

print('selected',flush=True)
start_frame_index_A = 0
start_frame_index_B = 0
end_frame_index_A = 0
end_frame_index_B = 0

rmsd_results = np.empty([0,0],dtype=np.float32)
pass_index = 0
block_index = 0 

for tA in chunkTrajA:
    end_frame_index_A = end_frame_index_A  + tA.n_frames
    #rmsd_results = np.pad (rmsd_results, [(0,tA.n_frames),(0,0)] , mode='constant', constant_values = 0 )
    for tB in chunkTrajB:
        end_frame_index_B = end_frame_index_B  + tB.n_frames
        
        #$if pass_index > 0 :  
        #$    print('here')
        #$    end_frame_index_B = min(end_frame_index_B  + tB.n_frames, total_frames)
        #$else:
        #$    end_frame_index_B = end_frame_index_B  + tB.n_frames

        #$while block_index < pass_index:
        #$    print('frame skipped')
        #$    block_index = block_index + 1
        #$    start_frame_index_B = end_frame_index_B
        #$    end_frame_index_B = min(end_frame_index_B  + tB.n_frames, total_frames)
        #$    continue
        #$print(block_index)

        if pass_index  < 1 :
            rmsd_results = np.pad (rmsd_results, [(0,tB.n_frames),(0,tB.n_frames)] , mode='constant', constant_values = 0 )

        #print(np.shape(rmsd_results),flush=True)

        #tA = tB
        #curr_rmsd = np.zeros([chunk_size,chunk_size])
        if block_index >= pass_index : 
            print('compute block begins',flush=True)
            start = time.time()
            torch.cuda.empty_cache()

            posAtt = torch.from_numpy(tA.xyz[:,atom_indices,:]).cuda()
            posBtt = torch.from_numpy(tB.xyz[:,atom_indices,:]).cuda()

            posAtt -= torch.mean(posAtt, keepdim=True, dim = 1)
            posBtt -= torch.mean(posBtt, keepdim=True, dim = 1)
            Hs = torch.einsum('nji,mjk->nmik', posAtt, posBtt)
            u, s, v = svd(Hs.flatten(0,1))
            uT = u.transpose(1,2)
            vuT = torch.matmul(v, uT)
            d = torch.sign(torch.linalg.det(vuT)).cuda()
            Md = torch.eye(3)[None,:,:].repeat(Hs.shape[0] * Hs.shape[1], 1, 1).cuda()
            Md[:, 2, 2] = d
            R = v @ (Md @ uT)
            R = R.reshape(Hs.shape)
            posBttRotate = posBtt.matmul(R)
            rmsds = torch.sqrt(torch.sum((posAtt[:,None,:,:] - posBttRotate) ** 2, dim = (2,3)) / n_atoms).cpu()
            curr_rmsd = rmsds.numpy()* 10 
            #print('actual shape')
            #print(np.shape(curr_rmsd))
            ##print(np.shape(curr_rmsd))
            #print('expected shape')
            #print(np.shape(rmsd_results [start_frame_index_A:end_frame_index_A, start_frame_index_B:end_frame_index_B]))
            #print ('B indices')
            #print(start_frame_index_B)
            #print(end_frame_index_B)
            #print ('A indices')
            #print(start_frame_index_A)
            #print(end_frame_index_A)
            #np.savetxt('rmsd_parralel.txt',rmsd_results)
            #$mat_size = np.shape(rmsd_results)[0]
            rmsd_results [start_frame_index_A:end_frame_index_A, start_frame_index_B: end_frame_index_B] = curr_rmsd
            rmsd_results [start_frame_index_B:end_frame_index_B, start_frame_index_A:end_frame_index_A] = curr_rmsd.T
            end = time.time()
            print('compute block finished, took  ' + str( end - start),flush=True)
        else:
            #chunkTrajB = (item for item in chunkTraj_backup_list)
            
            print('here')
            pass

            
        start_frame_index_B = end_frame_index_B
        block_index = block_index + 1

    #since it's a generator we need to reload the trajb universe 
    #chunkTrajB = md.iterload(traj_file, top="prot.psf", chunk = chunk_size)
    print('here')
    chunkTrajB = (item for item in chunkTraj_backup_list)

    total_frames = np.shape(rmsd_results )[0]
    pass_index = pass_index + 1
    block_index = 0 

    start_frame_index_B = 0
    end_frame_index_B = 0

    start_frame_index_A = end_frame_index_A
print('saving')
np.save(args.out_file,rmsd_results)
