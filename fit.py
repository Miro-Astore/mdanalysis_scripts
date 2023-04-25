import MDAnalysis as mda
import MDAnalysis.transformations.fit as fit
import sys 
import argparse
#####
#USAGE: python mdanalysis_scripts/fit.py -top TOPFILE -traj -o sel --stride n --begin --end 
#
#####

parser = argparse.ArgumentParser()

parser.add_argument('--ref-topology', '-rt', dest='ref_top' , help='reference topology file',type=str) 
parser.add_argument('--ref-structure', '-rs', dest='ref_struct', help='reference structure file',type=str) 
parser.add_argument('--target-topology', '-tt', dest='target_top' , help='target topology file',type=str) 
parser.add_argument('--target-simulation', '-ts', dest='target_sim', help='target trajectory file.',type=str) 
parser.add_argument('--out', '-o', dest='out_trajectory', help='Output trajectory file',type=str) 
parser.add_argument('--stride', '-st', dest='stride' , default = 1, help='Stride through trajectory skipping this many frames.',type=int) 
parser.add_argument('--selection', '-s', dest='selection_string', help='Selection string for fitting. Will be applied to both target and reference structures.',type=str) 
parser.add_argument('--begin', '-b', dest='beginning_frame', help='Frame from which to begin analysis.',type=int) 
parser.add_argument('--end', '-e', dest='ending_frame', help='Frame at which to end analysis.',type=int) 
parser.add_argument('--first_n_frames', '-f', dest='first_n_frames', help='Analyse the first n frames of the trajectory.',type=int) 
parser.add_argument('--last_n_frames', '-l', dest='last_n_frames', help='Analyse the last n frames of the trajectory.',type=int) 

args = parser.parse_args()


reference_universe = mda.Universe(args.ref_top,args.ref_struct)

target_universe = mda.Universe(args.target_top,args.target_traj)

selection_string = args.selection_string 

target_sel = target_universe.select_atoms(selection_string)
ref_sel = reference_universe.select_atoms(selection_string)

transform = fit.fit_rot_trans(target_sel, reference_sel)

target_universe.trajectory.add_transformations(transform)

with MDAnalysis.Writer(args.out_trajectory, target_universe.n_atoms) as W:
    for ts in u.trajectory:
        W.write(protein)
