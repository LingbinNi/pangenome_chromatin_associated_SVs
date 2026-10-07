#!/bin/bash
# Written by Lingbin.
# Map nonhuman-primate Hi-C reads to Hap2, processing two replicates separately.
# Edit the input paths below before running.
# Uses SGE qsub and environment modules; adapt module names and SGE resources
# to your cluster. Running this script submits two mapping jobs per sample.
# Use a separate output root for each reference (CHM13, Hap1 and Hap2).
# Run: bash /path/to/map_nhp_hic_to_hap2.sh
#
# Read input: directories matching /path/to/clean_hic_reads/{A,J,P}*.
# Each sample directory must contain one paired FASTQ set per replicate:
# rep1: *_1_*_R1_*.clean.fastq.gz and *_1_*_R2_*.clean.fastq.gz;
# rep2: *_2_*_R1_*.clean.fastq.gz and *_2_*_R2_*.clean.fastq.gz.
# Each pattern must match one file. Use absolute paths without whitespace.
# Run from the Hap2 output root, after reference-index jobs have completed.
# Required per-sample reference: <sample>/database/reference/<sample>.hap2.fa
# with its BWA indexes, and <sample>/database/chr_size/<sample>_HAP2.chr.size.
# These reference files are prepared by prepare_nhp_hap2_references.sh.
#
# For each replicate: BWA-MEM alignment; pairtools parsing (MAPQ 40,
# walks-policy 5unique, max-inter-align-gap 30), sorting, deduplication and
# splitting; then coordinate-sort and index the BAM. Deduplication is performed
# independently within each replicate.
# Outputs under <sample>/rep1/ and <sample>/rep2/: <sample>.repN.mapped.pairs,
# <sample>.repN.mapped.PT.bam and its index, stats.txt, script.sh and job logs.

CURRENT_DIR=$(pwd)

for dir in $(ls -d /path/to/clean_hic_reads/{A,J,P}*)
do
sample=$(basename $dir)
echo $sample
cd ./$sample

#############################################
## Rep1
#############################################
mkdir rep1
cd rep1
echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -l disk_free=10G
#$ -pe serial 16
#$ -o output.log
#$ -e error.log

module load bedtools/2.31.1
module load deeptools/3.5.5
module load pairtools/1.0.3
module load bwa/0.7.17
module load samtools/1.19
module load preseq/3.2.0
module load miniconda/24.7.1

mkdir -p '$CURRENT_DIR'/'$sample'/rep1/temp

bwa mem -5SP -T0 -t16 '$CURRENT_DIR'/'$sample'/database/reference/'$sample'.hap2.fa '$dir'/*_1_*_R1_*.clean.fastq.gz '$dir'/*_1_*_R2_*.clean.fastq.gz| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path '$CURRENT_DIR'/'$sample'/database/chr_size/'$sample'_HAP2.chr.size | \
pairtools sort --tmpdir='$CURRENT_DIR'/'$sample'/rep1/temp/ --nproc 16|\
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt|\
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs '$sample'.rep1.mapped.pairs --output-sam -|\
samtools view -bS -@16 | samtools sort -@16 -T '$CURRENT_DIR'/'$sample'/rep1/temp/temp.bam -o '$sample'.rep1.mapped.PT.bam;
samtools index '$sample'.rep1.mapped.PT.bam' > script.sh
qsub script.sh
cd ..

#############################################
## Rep2
#############################################
mkdir rep2
cd rep2
echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -l disk_free=10G
#$ -pe serial 16
#$ -o output.log
#$ -e error.log

module load bedtools/2.31.1
module load deeptools/3.5.5
module load pairtools/1.0.3
module load bwa/0.7.17
module load samtools/1.19
module load preseq/3.2.0
module load miniconda/24.7.1

mkdir -p '$CURRENT_DIR'/'$sample'/rep2/temp

bwa mem -5SP -T0 -t16 '$CURRENT_DIR'/'$sample'/database/reference/'$sample'.hap2.fa '$dir'/*_2_*_R1_*.clean.fastq.gz '$dir'/*_2_*_R2_*.clean.fastq.gz| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path '$CURRENT_DIR'/'$sample'/database/chr_size/'$sample'_HAP2.chr.size | \
pairtools sort --tmpdir='$CURRENT_DIR'/'$sample'/rep2/temp/ --nproc 16|\
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt|\
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs '$sample'.rep2.mapped.pairs --output-sam -|\
samtools view -bS -@16 | samtools sort -@16 -T '$CURRENT_DIR'/'$sample'/rep2/temp/temp.bam -o '$sample'.rep2.mapped.PT.bam;
samtools index '$sample'.rep2.mapped.PT.bam' > script.sh
qsub script.sh
cd ..

cd ..
done
