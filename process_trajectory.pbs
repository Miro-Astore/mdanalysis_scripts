#!/bin/bash
#PBS -P r16
#PBS -q normal
#PBS -l ncpus=2
#PBS -l storage=scratch/r16
#PBS -m abe
#PBS -l walltime=48:00:00
#PBS -l mem=16GB
#PBS -l wd
#PBS -M yearlyboozefest@gmail.com

module load gromacs/2019.3
cd $PBS_O_WORKDIR

mkdir -p ../xtcs/
cd AWORKDIR
pwd
file_name=$(echo "AWORKDIR" | sed "s^\/^_^g")
cp ionized.psf /scratch/r16/ma2374/gmx_cftr_2nd_round/xtcs/$file_name.psf
#gmx trjcat -f memb_prod*.trr -o all.xtc -dt 1000
membs=$(ls -v  memb_prod*xtc | tail -n 1 )
echo $membs
echo -e "0\n" | gmx trjconv -f $membs -s memb0.tpr -o $PBS_O_WORKDIR/../xtcs/$file_name\_temp1.xtc -pbc whole
echo -e "1\n1\n0\n" | gmx trjconv -center yes  -f $PBS_O_WORKDIR/../xtcs/$file_name\_temp1.xtc -s memb0.tpr -o $PBS_O_WORKDIR/../xtcs/$file_name\_temp2.xtc -pbc cluster
echo -e "0\n" | gmx trjconv -f $PBS_O_WORKDIR/../xtcs/$file_name\_temp2.xtc  -s memb0.tpr -o $PBS_O_WORKDIR/../xtcs/$file_name.xtc -pbc res
echo -e "0\n" | gmx trjconv -f $PBS_O_WORKDIR/../xtcs/$file_name.xtc -dt 10000   -o $PBS_O_WORKDIR/../xtcs/$file_name\_sum.xtc 

source $HOME/.bashrc
python $PBS_O_WORKDIR/../xtcs/mdanalysis_scripts/res_specific_rmsd.py /scratch/r16/ma2374/gmx_cftr_2nd_round/xtcs/$file_name.psf $PBS_O_WORKDIR/../xtcs/$file_name.xtc $PBS_O_WORKDIR/../xtcs/6msm_prot.pdb
#python $PBS_O_WORKDIR/mdanalysis_scripts/grab_last_frame.py /scratch/r16/ma2374/gmx_cftr_2nd_round/xtcs/$file_name.psf  $PBS_O_WORKDIR/../xtcs/$file_name\_sum.xtc $file_name.pdb



rm $PBS_O_WORKDIR/../xtcs/$file_name\_temp1.xtc $PBS_O_WORKDIR/../xtcs/$file_name\_temp2.xtc
