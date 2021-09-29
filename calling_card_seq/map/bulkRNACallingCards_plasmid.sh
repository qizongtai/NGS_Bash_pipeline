#!/bin/bash
#
#SBATCH --job-name=bulkRNA_CC_plasmid
#SBATCH --output=logs/bulkRNA_CC_plasmid_%a.out
#SBATCH --error=logs/bulkRNA_CC_plasmid_%a.err
#SBATCH --array=1-2
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=16000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amoudgil@wustl.edu

# Demultiplexes files. Needs a file of input filenames
# and a corresponding file of barcodes
filename=$( sed -n ${SLURM_ARRAY_TASK_ID}p masterFile_R1.txt )
barcode=$( sed -n ${SLURM_ARRAY_TASK_ID}p 5prime_barcodes_R1.txt )
ref=$( sed -n ${SLURM_ARRAY_TASK_ID}p plasmids.txt )
base=`basename $filename`
stem="${base%%.*}"

module load cutadapt
module load novoalign
module load samtools
module load bedtools
module load java
module load drop-seq-tools

out_file="../output_and_analysis/plasmid_"$stem
out_sam=$out_file".sam"
out_bam=$out_file".bam"
out_map_sort_prefix=$out_file"_map_sort"
out_bam_map=$out_map_sort_prefix".bam"
out_bed=$out_map_sort_prefix".bed"
out_bedgraph=$out_map_sort_prefix".bedgraph"

# Trim and demultiplex by LTR + GGTTAA insertion site (exact matching)
cutadapt \
    -g "^"$barcode"GGTTAA" \
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
    -d $ref \
    -f $out_file".fastq.gz" \
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

# Index the BAM file
samtools index $out_bam_map

# Clean up
rm $out_sam
rm $out_file".temp.fastq.gz"

## Tag reads with genes and features
#java -Xmx4g -jar $DROPSEQTOOLS_HOME/dropseq.jar \
    #TagReadWithGeneExon \
    #I=$out_bam_map \
    #O=$out_map_sort_prefix"_tagged.bam" \
    #ANNOTATIONS_FILE=/scratch/rmlab/ref/star_genomes/hg19/hg19.gtf \
    #TAG=GE

## Tag reads with transposon insertions
#python TagBamWithTransposonInserts.py \
    #--transposase PB \
    #$out_map_sort_prefix"_tagged.bam" \
    #/scratch/rmlab/ref/star_genomes/hg19/hg19.2bit \
    #$out_map_sort_prefix"_final.bam"

## Index the final BAM file
#samtools index $out_map_sort_prefix"_final.bam"

## Get a list of unique insertions
#paste \
    #<(samtools view $out_map_sort_prefix"_final.bam" | egrep -o "XI:Z:\S*") \
    #<(samtools view $out_map_sort_prefix"_final.bam" | egrep -o "XZ:Z:\S*") \
    #<(samtools view $out_map_sort_prefix"_final.bam" | egrep -o "GS:Z:\S*") | \
    #sort -V | \
    #uniq | \
    #grep "XZ:Z:TTAA" > $out_map_sort_prefix"_final.txt"
