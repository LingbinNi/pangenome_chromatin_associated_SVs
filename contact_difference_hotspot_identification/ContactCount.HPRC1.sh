#!/bin/bash
#written by Lingbin

cooler_dir="/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/cdh_identification_update/01.ContactMatrix"

for dir in $(ls -d /net/eichler/vol28/projects/hic_cohorts/nobackups/00.Analysis/mapping_difference_analysis_hprc1/HG*)
do
sample=$(basename $dir)
echo $sample
mkdir $sample
cd $sample

echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -l disk_free=10G
#$ -pe serial 1
#$ -o output.log
#$ -e error.log

module load pairix/0.3.8
module load cooltools/0.6.1
module load miniconda/24.7.1

cooler dump \
--join \
-o '$sample'.SelfMapping.Hap1.HPRC.25kb.Unique.bed \
'$cooler_dir'/'$sample'/'$sample'.SelfMapping.Hap1.HPRC.25kb.cool

cooler dump \
--join \
-o '$sample'.SelfMapping.Hap2.HPRC.25kb.Unique.bed \
'$cooler_dir'/'$sample'/'$sample'.SelfMapping.Hap2.HPRC.25kb.cool
' > script.sh
qsub script.sh

cd ..
done
