#usage: python mdanalysis/res_specific_rmsd.py TOPFILE TRAJFILE REF_STRUCTURE OUT_NAME.DATZZ
import numpy as np
import os 
import MDAnalysis as mda
import MDAnalysis.analysis.rms
import sys
import argparse
#system_list=np.loadtxt("system_list",dtype=str)
parser = argparse.ArgumentParser()
cwd=os.getcwd()
parser.add_argument ('--topology',dest='topology',help='topology',required=True)
parser.add_argument ('--trajectory',dest='trajectory',help='trajectory',required=True)
parser.add_argument ('--stride',dest='stride',help='stride of trajectory, skip this many frames each jump.',required=False,default=1,type=int)
parser.add_argument ('--sel_list',dest='selection_list',help='list of selections to collect the rmsd for.',nargs="+",required=True)
parser.add_argument ('--ref_sel',dest='reference_selection',help='list of selections to collect the rmsd for.',required=True)
parser.add_argument ('--ref_struct',dest='ref_struct',help='structure for reference',required=True)
parser.add_argument ('--ref_top',help='structure for reference')
parser.add_argument ('--o',dest='o',help='outfile')
args = parser.parse_args()

if args.ref_top == None:
    reference_topology = args.topology
else:
    reference_topology = args.ref_top


#f=open('system_list','r')

#f1=f.readlines()


#os.chdir ('..')
#print (i)
#files=i.split()
#psf=files[0]
#xtc=files[1]
u = mda.Universe(args.topology,args.trajectory)
r = mda.Universe(reference_topology,args.ref_struct)
R = MDAnalysis.analysis.rms.RMSD(u,r,select=args.reference_selection, groupselections=args.selection_list, filename=args.o)
R.run()

header_str = 'time, '
for i in range(0,len(args.selection_strings),2):
    header_str = header_str + '\"' + str(args.selection_strings[i]) +  '\"-\"' + str(args.selection_strings[i+1]) + '\", '
header_str = header_str [:-2]
np.savetxt(args.o, R.rmsd,header=header_str)
