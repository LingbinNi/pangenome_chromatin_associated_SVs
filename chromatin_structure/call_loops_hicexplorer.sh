#!/bin/bash
#written by Lingbin

# Call chromatin loops with HiCExplorer 3.7.3 using 10 kb matrices.
# Replace matrix_dir with an absolute path on your system (without spaces).
# Input: <matrix_dir>/<sample>/Merge/<sample>.SelfMapping.Hap1.HPRC.10kb.corrected.h5
# Sample-level merging and matrix correction must be completed beforehand.
# For Hap2 or CHM13, edit both matrix_dir and the matrix filename below.
# Adjust the {HG,NA}* sample-directory pattern to match your sample names.
# Run directly on a suitable compute node; this script does not submit SGE jobs.
# Use a new loop-output directory, separate from the TAD-output directory.
# Adapt module loading to your system while retaining the software version.
# Keep the original output prefix and .bedgraph extension for downstream use.

module load hicexplorer/3.7.3

matrix_dir="/path/to/hap1_contact_matrices"

for dir in $(ls -d $matrix_dir/{HG,NA}*)
do
sample=$(basename $dir)
echo $sample
mkdir $sample
cd $sample

################################################################
## loop
################################################################

hicDetectLoops \
-m $dir/Merge/$sample.SelfMapping.Hap1.HPRC.10kb.corrected.h5 \
-o $sample\_max2000000_window10_peak6_pvaluepre0.05_pvalue0.05.bedgraph \
--maxLoopDistance 2000000 \
--windowSize 10 \
--peakWidth 6 \
--pValuePreselection 0.05 \
--pValue 0.05 \
-t 12 \
> $sample.o 2> $sample.e

cd ..
done
