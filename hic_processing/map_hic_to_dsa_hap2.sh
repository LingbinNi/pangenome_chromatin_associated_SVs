#!/bin/bash
#written by Lingbin

# Replace /path/to/hic_data with the directory containing sample folders.
# Adjust the sample glob and FASTQ filename pattern for your input files.
# Run from a Hap2-specific working directory containing these files:
#   <sample>/database/reference/ragtag.scaffold.fasta
#   <sample>/database/chr_size/<sample>_HAP2.chr.size
# Expected FASTQs: <hic_data>/<sample>/raw_data/hic/<prefix>_R[12]_001.fastq.gz
# Adjust SGE resource requests and module names for your cluster.

current_path=$(pwd)

for dir in $(ls -d /path/to/hic_data/{HG,NA}*)
do
sample=$(basename $dir)
cd $sample

for file in $(ls $dir/raw_data/hic/*_R1_*.fastq.gz)
do
filename=$(basename $file _R1_001.fastq.gz)
echo $filename
mkdir $filename
cd $filename

echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -pe serial 16
#$ -o output.log
#$ -e error.log

module load bedtools/2.31.1
module load deeptools/3.5.5
module load pairtools/1.0.3
module load bwa/0.7.17
module load samtools/1.19
module load preseq/3.2.0
module load miniconda/24.11.1

mkdir -p '$current_path'/'$sample'/'$filename'/temp

bwa mem -5SP -T0 -t16 '$current_path'/'$sample'/database/reference/ragtag.scaffold.fasta '$dir'/raw_data/hic/'$filename'_R1_001.fastq.gz '$dir'/raw_data/hic/'$filename'_R2_001.fastq.gz| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path '$current_path'/'$sample'/database/chr_size/'$sample'_HAP2.chr.size | \
pairtools sort --tmpdir='$current_path'/'$sample'/'$filename'/temp/ --nproc 16|\
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --keep-parent-id --output-dups dedup.dups.pairsam|\
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs '$filename'.mapped.pairs --output-sam -|\
samtools view -bS -@16 | samtools sort -@16 -T '$current_path'/'$sample'/'$filename'/temp/temp.bam -o '$filename'.mapped.PT.bam;

bgzip -@ 16 dedup.dups.pairsam
bgzip -@ 16 '$filenamebase'.mapped.pairs
' > script.sh
qsub script.sh

cd ..
done

cd ..
done
