#!/bin/bash
#SBATCH --mem=1G
#SBATCH -o logs/bash_%A_%a.out
#SBATCH -e logs/bash_%A_%a.err
#SBATCH -J bash
#SBATCH --array=1-3

input[1]='../1_data/AC_sorted'
input[2]='../1_data/AY30_sorted'
input[3]='../1_data/AC_downsampled_sorted'

ID=${SLURM_ARRAY_TASK_ID} 

egrep -v ".*_.*" ${input[$ID]}.ccf > ${input[$ID]}_clean.ccf









