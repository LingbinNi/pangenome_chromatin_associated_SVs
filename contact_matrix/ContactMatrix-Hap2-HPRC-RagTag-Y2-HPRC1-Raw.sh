#!/bin/bash
#written by Lingbin

database_path="/net/eichler/vol28/projects/hprc/nobackups/hic_hprc/04.SelfMapping-Hap2-HPRC-RagTag-Y2-HPRC1"
mapping_path="/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/04.SelfMapping-Hap2-HPRC-RagTag-Y2-HPRC1-RawData"
hap="2"

for file in $(ls $mapping_path/*/*/*.mapped.pairs.gz)
do
sample=$(echo $file | awk -F'/' '{print $10}')
filename=$(echo $file | awk -F'/' '{print $11}')
mkdir -p ./$sample/$filename
cd ./$sample/$filename

echo $filename
echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -pe serial 8
#$ -o output.log
#$ -e error.log

module load pairix/0.3.8
module load cooltools/0.6.1
module load hicexplorer/3.7.3
module load miniconda/23.5.2

################################################################
## raw 
################################################################

ln -s '$file' '$filename'.mapped.pairs.gz

################################################################
## hic
################################################################

java -Xmx40000m -jar /net/eichler/vol28/home/lbni/software/nobackups/Juicer/Juicertools/juicer_tools_1.22.01.jar pre --threads 8 '$file' '$filename'.SelfMapping.Hap'$hap'.HPRC.hic '$database_path'/'$sample'/database/chr_size/'$sample'_HAP'$hap'.chr.size

################################################################
## cool
################################################################

pairix '$filename'.mapped.pairs.gz

reslist=(10 25 100)
for res in ${reslist[@]}
do
resbp=$(expr $res \* 1000)
cooler cload pairix -p 8 '$database_path'/'$sample'/database/chr_size/'$sample'_HAP'$hap'.chr.size:$resbp '$filename'.mapped.pairs.gz '$filename'.SelfMapping.Hap'$hap'.HPRC.$res\kb.cool
done

################################################################
## h5
################################################################

for file in $(ls '$filename'.SelfMapping.Hap'$hap'.HPRC.{10,25,100}kb.cool)
do
filename=$(basename $file .cool)
hicConvertFormat \
--matrices $file \
--outFileName $filename.h5 \
--inputFormat cool \
--outputFormat h5
done
' > script.sh
qsub script.sh

cd ../..
done
