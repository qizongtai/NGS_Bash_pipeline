#!/bin/bash
#SBATCH --mem=20G
#SBATCH -o logs/call_peaks_macs_%A_%a.out
#SBATCH -e logs/call_peaks_macs_%A_%a.err
#SBATCH -J call_peak_macs
#SBATCH --array=1-2
ID=${SLURM_ARRAY_TASK_ID} 

module load ccf_tools

declare -A win step
win[a]='500'
step[a]='250'
win[b]='1000'
step[b]='500'

for i in {a..b}; do
#for i in {a,b}; do #an alternative

exp[1]='../1_data/AC_sorted_clean.ccf'
out[1]="../3_output/AC_sorted_clean_macs_win${win[$i]}_step${step[$i]}_p0001.txt"
exp[2]='../1_data/AY30_sorted_clean.ccf'
out[2]="../3_output/AY30_sorted_clean_macs_win${win[$i]}_step${step[$i]}_p0001.txt"

#python $CCF_TOOLS/call_peaks_macs.py --help
#---without background insertions
#:<<'END'
python $CCF_TOOLS/call_peaks_macs.py    -e ${exp[$ID]} \
					-o ${out[$ID]} \
					-t /scratch/ref/rmlab/calling_card_ref/mouse/TTAA_mm10.ccf \
					-a /scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed \
					-pc 0.001 \
					--peak_finder_pvalue 0.01 \
					--window ${win[$i]} \
					--step ${step[$i]} \
					--pseudocounts 0.2
#END

#---with background insertions
:<<'END'
python $CCF_TOOLS/call_peaks_macs.py    -e ${exp[$ID]} \
					-b ${bg[$ID]} \
					-o ${out[$ID]} \
					-t /scratch/ref/rmlab/calling_card_ref/mouse/TTAA_mm10.ccf \
					-a /scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed \
					-pc 0.001 \
					--peak_finder_pvalue 0.01 \
					--window ${win[$i]} \
					--step ${step[$i]} \
					--pseudocounts 0.2
END

done
