#!/bin/bash
#written by Lingbin

for dir in $(ls -d /net/eichler/vol28/projects/hic_cohorts/nobackups/00.NHP/00.DataSeq/CleanData/Hi-C/{A,J,P}*)
do
sample=$(basename $dir)
mkdir -p ./$sample/database
cd ./$sample/database

echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=2G
#$ -pe serial 8
#$ -o output.log
#$ -e error.log

module load samtools/1.19
module load bwa/0.7.15

#reference
mkdir reference
cd reference
ln -s /net/eichler/vol28/projects/hic_cohorts/nobackups/00.NHP/00.DataBase/assemblies_old/'$sample'.hap2.fa .
samtools faidx '$sample'.hap2.fa
bwa index '$sample'.hap2.fa
cd ..

#chr_size
mkdir chr_size
cd chr_size
cut -f1,2 ../reference/*.fai > '$sample'_HAP2.chr.size
cd ..
' > script.sh
qsub script.sh

cd ../..
done
