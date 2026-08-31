#!/bin/bash
#written by Lingbin

specific_dir="/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/cdh_identification_update/00.MappingDifference"
database_dir1="/net/eichler/vol28/projects/hprc/nobackups/hic_hprc/03.SelfMapping-Hap1-HPRC-RagTag-Y2-HPRC1"
database_dir2="/net/eichler/vol28/projects/hprc/nobackups/hic_hprc/04.SelfMapping-Hap2-HPRC-RagTag-Y2-HPRC1"

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
#$ -l disk_free=1G
#$ -pe serial 8
#$ -o output.log
#$ -e error.log

module load pairix/0.3.8
module load cooltools/0.6.1
module load miniconda/24.7.1
module load pairtools/1.0.3
export LC_ALL=C

################################################################
## Hap1
################################################################

ln -s '$specific_dir'/'$sample'/Merge/'$sample'.Hap1.unique.merge.pairs .
sort -k2,2 -k4,4 -k3,3n -k5,5n --parallel=8 -S 4G '$sample'.Hap1.unique.merge.pairs > '$sample'.Hap1.sorted.pairs
bgzip '$sample'.Hap1.sorted.pairs
pairix -p pairs '$sample'.Hap1.sorted.pairs.gz
cooler cload pairix -p 8 '$database_dir1'/'$sample'/database/chr_size/'$sample'_HAP1.chr.size:25000 '$sample'.Hap1.sorted.pairs.gz '$sample'.SelfMapping.Hap1.HPRC.25kb.cool
rm '$sample'.Hap1.sorted.pairs.gz*

################################################################
## Hap2
################################################################

ln -s '$specific_dir'/'$sample'/Merge/'$sample'.Hap2.unique.merge.pairs .
sort -k2,2 -k4,4 -k3,3n -k5,5n --parallel=8 -S 4G '$sample'.Hap2.unique.merge.pairs > '$sample'.Hap2.sorted.pairs
bgzip '$sample'.Hap2.sorted.pairs
pairix -p pairs '$sample'.Hap2.sorted.pairs.gz
cooler cload pairix -p 8 '$database_dir2'/'$sample'/database/chr_size/'$sample'_HAP2.chr.size:25000 '$sample'.Hap2.sorted.pairs.gz '$sample'.SelfMapping.Hap2.HPRC.25kb.cool
rm '$sample'.Hap2.sorted.pairs.gz*
' > script.sh
qsub script.sh

cd ..
done
