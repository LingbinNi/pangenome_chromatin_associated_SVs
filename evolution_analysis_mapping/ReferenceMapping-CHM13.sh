#!/bin/bash
#written by Lingbin

DATABASE_DIR="/net/eichler/vol28/projects/hprc/nobackups/hic_hprc/01.ReferenceMapping-CHM13-HPRC-Y2-HPRC1/DataBase"
CURRENT_DIR=$(pwd)

for dir in $(ls -d /net/eichler/vol28/projects/hic_cohorts/nobackups/00.NHP/00.DataSeq/CleanData/Hi-C/{A,J,P}*)
do
sample=$(basename $dir)
echo $sample
mkdir ./$sample
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

bwa mem -5SP -T0 -t16 '$DATABASE_DIR'/reference/T2T-CHM13v2.fasta '$dir'/*_1_*_R1_*.clean.fastq.gz '$dir'/*_1_*_R2_*.clean.fastq.gz| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path '$DATABASE_DIR'/chr_size/T2T-CHM13v2.chr.size | \
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

bwa mem -5SP -T0 -t16 '$DATABASE_DIR'/reference/T2T-CHM13v2.fasta '$dir'/*_2_*_R1_*.clean.fastq.gz '$dir'/*_2_*_R2_*.clean.fastq.gz| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path '$DATABASE_DIR'/chr_size/T2T-CHM13v2.chr.size | \
pairtools sort --tmpdir='$CURRENT_DIR'/'$sample'/rep2/temp/ --nproc 16|\
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt|\
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs '$sample'.rep2.mapped.pairs --output-sam -|\
samtools view -bS -@16 | samtools sort -@16 -T '$CURRENT_DIR'/'$sample'/rep2/temp/temp.bam -o '$sample'.rep2.mapped.PT.bam;
samtools index '$sample'.rep2.mapped.PT.bam' > script.sh
qsub script.sh
cd ..

cd ..
done
