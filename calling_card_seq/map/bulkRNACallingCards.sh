#!/bin/bash
#
#SBATCH --job-name=bulkRNA_CC
#SBATCH --output=logs/bulkRNA_CC_%A_%a.out
#SBATCH --error=logs/bulkRNA_CC_%A_%a.err
#SBATCH --array=1-120%20
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qizongtai@gmail.com

module load py-cutadapt
module load novoalign
module load samtools
module load bedtools2
module load drop-seq-tools

# Demultiplexes files. Needs a file of input filenames
# and a corresponding file of barcodes
filename=$( sed -n ${SLURM_ARRAY_TASK_ID}p manifest_R1.txt )
barcode=$( sed -n ${SLURM_ARRAY_TASK_ID}p 5prime_barcodes_R1.txt ) 
index1=$( sed -n ${SLURM_ARRAY_TASK_ID}p manifest_I7.txt )
index2=$( sed -n ${SLURM_ARRAY_TASK_ID}p manifest_I5.txt )
base=`basename $filename`
stem=$barcode"_""${base%%.*}"

out_file="../output_and_analysis/"$stem
out_sam=$out_file".sam"
out_bam=$out_file".bam"
out_map_sort_prefix=$out_file"_map_sort"
out_bam_map=$out_map_sort_prefix".bam"
out_bed=$out_map_sort_prefix".bed"
out_bedgraph=$out_map_sort_prefix".bedgraph"

SB_ITR="TAAGTGTATGTAAACTTCCGACTTCAACTGTA"
alt_SB_ITR="AAGTGTATGTAAACTTCCGACTTCAACTGTA"
PB_LTR="TTTACGCAGACTATCTTTCTAGGGTTAA"
Long_PB_LTR="GCGTCAAT"$PB_LTR

transposase=PB
genome=hg38

# Trim and demultiplex by LTR + GGTTAA insertion site (exact matching)
cutadapt \
    -g "^"$barcode$PB_LTR \
    -o $out_file".temp.fastq.gz" \
    --minimum-length 1 \
    --discard-untrimmed \
    -e 0 \
    --no-indels \
    "../raw/"$filename

# Trim any trailing Nextera adapters (allowing mismatches)
cutadapt \
    -a "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTNNNNNNNNNNTCTCGTATGCCGTCTTCTGCTTG" \
    -o $out_file".fastq.gz" \
    --minimum-length 1 \
    $out_file".temp.fastq.gz"

# Align the reads
novoalign \
    -d /scratch/ref/rmlab/novoalign_indexes/$genome/$genome.nvx \
    -f $out_file".fastq.gz" \
    -n 40 \
    -o SAM \
    -o SoftClip > $out_sam

# Convert to BAM
samtools view \
    -bS -h \
    $out_sam \
    -o $out_bam

# Filter only mapped reads, convert to BAM, and sort
samtools view \
    -bS -h -F 260 \
    $out_sam | \
samtools sort - -o $out_bam_map

## Tag reads with genes and features
#java -Xmx4g -jar $DROPSEQTOOLS_HOME/dropseq.jar \
    #TagReadWithGeneExon \
    #I=$out_bam_map \
    #O=$out_map_sort_prefix"_features.bam" \
    #ANNOTATIONS_FILE=/scratch/ref/rmlab/star_genomes/hg19/hg19.gtf \
    #TAG=GE

# Tag reads with barcodes
python3 TagBam.py \
    --tag XP:Z:$barcode \
    $out_bam_map \
    $out_map_sort_prefix"_tagged.bam"

# Tag reads wih i7 indexes
python3 TagBam.py \
    --tag XJ:Z:$index1 \
    $out_map_sort_prefix"_tagged.bam" \
    $out_map_sort_prefix"_tagged2.bam"

# Tag reads with i5 indexes
python3 TagBam.py \
    --tag XK:Z:$index2 \
    $out_map_sort_prefix"_tagged2.bam" \
    $out_map_sort_prefix"_tagged3.bam"

# Tag reads with transposon insertions
python3 TagBamWithTransposonInserts.py \
    --transposase $transposase \
    -f \
    $out_map_sort_prefix"_tagged3.bam" \
    /scratch/ref/rmlab/novoalign_indexes/$genome/$genome.2bit \
    $out_map_sort_prefix"_final.bam"

# Index the final BAM file
samtools index $out_map_sort_prefix"_final.bam"

# Clean up
rm $out_sam
rm $out_bam
rm $out_file".temp.fastq.gz"
rm $out_bam_map
rm $out_map_sort_prefix"_features.bam"
rm $out_map_sort_prefix"_tagged.bam"
rm $out_map_sort_prefix"_tagged2.bam"
rm $out_map_sort_prefix"_tagged3.bam"

# Get a list of unique insertions
python3 BamToCallingCard.py -b XP XJ -i $out_map_sort_prefix"_final.bam" -o $out_map_sort_prefix"_final.ccf" 
