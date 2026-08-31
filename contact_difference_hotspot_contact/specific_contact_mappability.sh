#!/bin/bash
#written by Lingbin

for file in $(ls /net/eichler/vol28/projects/hic_cohorts/nobackups/00.DataBase/scaffolds/{HG,NA}*.Hap*.scaffolds.fasta)
do
filename=$(basename $file)
sample=$(basename $file | cut -d'.' -f1)
hap=$(basename $file | cut -d'.' -f2)
echo $sample.$hap
mkdir $sample.$hap
cd $sample.$hap

echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -l disk_free=10G
#$ -pe serial 10
#$ -o output.log
#$ -e error.log

module load miniconda/24.7.1

################################################################
## Create an index for a genome 
################################################################
ln -s '$file' .
newmap index '$filename'

################################################################
## Find the minimum unique k-mer lengths for the genome using the index 
################################################################

newmap search --verbose --num-threads=10 --search-range=20:200 --output-directory=unique_lengths '$filename'

################################################################
## Convert the unique lengths to mappability tracks 
################################################################

newmap track --single-read=150.bed --multi-read=150.wig 150 unique_lengths/*.unique.uint8
' > script.sh
qsub script.sh

cd ..
done
