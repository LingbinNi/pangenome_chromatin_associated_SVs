#!/bin/bash
#written by Lingbin

# Replace /path/to/hic_data and /path/to/chm13_database with your local paths.
# The reference database is expected to contain:
#   reference/T2T-CHM13v2.fasta (with BWA index files)
#   chr_size/T2T-CHM13v2.chr.size
# Expected FASTQs: <hic_data>/<sample>/<prefix>/<prefix>_R1*.fastq.gz and R2 mates.
# Adjust the sample glob, file names, SGE resources and modules for your system.
# Run from the directory in which sample output folders should be created.

current_path=$(pwd)
database_path="/path/to/chm13_database"

for dir in $(ls -d /path/to/hic_data/{HG,NA}*)
do
sample=$(basename $dir)
echo $sample
mkdir $sample
cd $sample

for file in $(ls $dir/*/*_R1*.fastq.gz)
do
filename1=$(basename $file)
filenamebase=${filename1%%_R1*}
filename2=${filename1/_R1/_R2}
echo $filenamebase
mkdir $filenamebase
cd $filenamebase

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

mkdir -p '$current_path'/'$sample'/'$filenamebase'/temp

bwa mem -5SP -T0 -t16 '$database_path'/reference/T2T-CHM13v2.fasta '$dir'/'$filenamebase'/'$filename1' '$dir'/'$filenamebase'/'$filename2'| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path '$database_path'/chr_size/T2T-CHM13v2.chr.size | \
pairtools sort --tmpdir='$current_path'/'$sample'/'$filenamebase'/temp/ --nproc 16|\
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --keep-parent-id --output-dups dedup.dups.pairsam|\
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs '$filenamebase'.mapped.pairs --output-sam -|\
samtools view -bS -@16 | samtools sort -@16 -T '$current_path'/'$sample'/'$filenamebase'/temp/temp.bam -o '$filenamebase'.mapped.PT.bam;

bgzip -@ 16 dedup.dups.pairsam
bgzip -@ 16 '$filenamebase'.mapped.pairs
' > script.sh
qsub script.sh

cd ..
done

cd ..
done
