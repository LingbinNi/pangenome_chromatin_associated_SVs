#!/bin/bash
#written by Lingbin

module load hicexplorer/3.7.3

matrix_dir="/net/eichler/vol28/projects/hgsvc/nobackups/analysis/genome_analysis/07.ContactMatrix-Hap1-HPRC-RagTag-Y2-HPRC1-RawData"

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
