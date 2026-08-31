#!/bin/bash
#written by Lingbin

for file in $(ls /net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.ExpressionAnalysis/02.LongReadSequence-Ref/01.Align/*/*.transcripts.pbmm2.bam)
do
sample=$(basename $file .transcripts.pbmm2.bam)
echo $sample
mkdir $sample
cd $sample

echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -pe serial 1
#$ -o '$sample'.o
#$ -e '$sample'.e

module load isoquant/3.10.0

isoquant.py \
--reference /net/eichler/vol28/projects/hprc/nobackups/analysis/transcriptome_analysis/00.DataBase/references/CHM13/reference/T2T-CHM13v2.fasta \
--genedb /net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.ExpressionAnalysis/02.LongReadSequence-Ref/03.PigeonUpdate/DataBase/chm13v2.0_RefSeq_Liftoff_v5.2.sorted.gtf \
--data_type pacbio_ccs \
--bam '$file' \
-o ./
' > $sample.sh
qsub $sample.sh

cd ..
done
