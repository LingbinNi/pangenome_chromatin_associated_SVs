#!/bin/bash
#written by Lingbin

current_path=$(pwd)
database_path="/net/eichler/vol28/projects/hprc/nobackups/hic_hprc/01.ReferenceMapping-CHM13-HPRC-Update/DataBase"

for dir in $(ls -d /net/eichler/vol28/projects/hic_cohorts/nobackups/00.DataSeq/01.RawData/HiC/HPRC-Y2-Files/HG00*)
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
rm '$filenamebase'.mapped.PT.bam

bgzip -@ 16 dedup.dups.pairsam
bgzip -@ 16 '$filenamebase'.mapped.pairs
' > script.sh
qsub script.sh

cd ..
done

cd ..
done
