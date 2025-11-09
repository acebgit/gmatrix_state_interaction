#!/bin/bash

############################################################
##             QCHEM: LAUNCHING CALCULATIONS
############################################################

#SBATCH --partition=xlong
#SBATCH --job-name=que_2
#SBATCH --mem=100gb
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=que_2.err
#SBATCH --time=8-00:00:00
#SBATCH --mail-user=antonio.cebreiro@dipc.org
#SBATCH --mail-type=END

export MODULEPATH=/scratch/abel/SOFTWARE/privatemodules:$MODULEPATH
module load qchem_group  # group or trunk
export QCHEM_CPUS=$SLURM_NTASKS

if [ ! $1 ]; then
 echo "An input file is required"
 exit
fi

export Project=${1:0:-4}
export runqchem=qchem
$runqchem -nt $SLURM_NTASKS $Project.inp $Project.out

if [ -d "$Project" ];then
   rm -r $Project
   mkdir $Project
 else
   mkdir $Project
fi
mv $Project.* scratch_data $Project/.
cp $Project/$Project.inp .

