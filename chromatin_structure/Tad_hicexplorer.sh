#!/bin/bash
#written by Lingbin

matrix_dir="/net/eichler/vol28/projects/hgsvc/nobackups/analysis/genome_analysis/07.ContactMatrix-Hap1-HPRC-RagTag-Y2-HPRC1"

for dir in $(ls -d $matrix_dir/{HG,NA}*)
do
sample=$(basename $dir)
echo $sample
mkdir $sample
cd $sample

echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -pe serial 8
#$ -o '$sample'.o
#$ -e '$sample'.e

module load hicexplorer/3.7.3

################################################################
## tad
################################################################

hicFindTADs \
-m '$dir'/Merge/'$sample'.SelfMapping.Hap1.HPRC.25kb.corrected.h5 \
--outPrefix '$sample'_min75000_max250000_step25000_thres0.05_delta0.01_fdr \
--minDepth 75000 \
--maxDepth 250000 \
--step 25000 \
--thresholdComparisons 0.05 \
--delta 0.01 \
--correctForMultipleTesting fdr \
-p 8 ' > script.sh
qsub script.sh

cd ..
done
