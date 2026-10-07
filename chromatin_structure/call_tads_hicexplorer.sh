#!/bin/bash
#written by Lingbin

# Call TADs and boundaries with HiCExplorer 3.7.3 using 25 kb matrices.
# Replace matrix_dir with an absolute path on your system (without spaces).
# Input: <matrix_dir>/<sample>/Merge/<sample>.SelfMapping.Hap1.HPRC.25kb.corrected.h5
# Sample-level merging and matrix correction must be completed beforehand.
# For Hap2 or CHM13, edit both matrix_dir and the matrix filename below.
# Adjust the {HG,NA}* sample-directory pattern to match your sample names.
# Use a new TAD-output directory, separate from the loop-output directory.
# Adapt module loading and SGE resources to your cluster.
# Depth and step parameters are in bp; thresholdComparisons is a q-value cutoff with fdr.
# All original calling parameters and output prefixes are retained.

matrix_dir="/path/to/hap1_contact_matrices"

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
