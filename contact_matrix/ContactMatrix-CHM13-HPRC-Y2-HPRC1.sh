#!/bin/bash
#written by Lingbin

mapping_path="/net/eichler/vol28/projects/hprc/nobackups/hic_hprc/01.ReferenceMapping-CHM13-HPRC-Y2-HPRC1"

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

cp '$file' '$filename'.mapped.pairs.gz
gunzip '$filename'.mapped.pairs.gz
bgzip '$filename'.mapped.pairs

################################################################
## hic
################################################################

java -Xmx40000m -jar /net/eichler/vol28/home/lbni/software/nobackups/Juicer/Juicertools/juicer_tools_1.22.01.jar pre --threads 8 '$file' '$filename'.ReferenceMapping.CHM13.HPRC.hic '$mapping_path'/DataBase/chr_size/T2T-CHM13v2.chr.size

################################################################
## cool
################################################################

pairix '$filename'.mapped.pairs.gz

reslist=(10 25 100)
for res in ${reslist[@]}
do
resbp=$(expr $res \* 1000)
cooler cload pairix -p 8 '$mapping_path'/DataBase/chr_size/T2T-CHM13v2.chr.size:$resbp '$filename'.mapped.pairs.gz '$filename'.ReferenceMapping.CHM13.HPRC.$res\kb.cool
done

################################################################
## h5
################################################################

for file in $(ls '$filename'.ReferenceMapping.CHM13.HPRC.{10,25,100}kb.cool)
do
filename=$(basename $file .cool)
hicConvertFormat \
--matrices $file \
--outFileName $filename.h5 \
--inputFormat cool \
--outputFormat h5
done


################################################################
## corrected 
################################################################

################################################################
## hic
################################################################

ln -s '$filename'.ReferenceMapping.CHM13.HPRC.hic '$filename'.ReferenceMapping.CHM13.HPRC.corrected.hic

################################################################
## h5
################################################################

for file in $(ls *.h5)
do
filename=$(basename $file .h5)
### diagnostic_plot
hicCorrectMatrix diagnostic_plot \
--matrix $file \
--plotName $filename.png
### correct
hicCorrectMatrix correct \
--matrix $file \
--outFileName $filename.corrected.h5 \
--correctionMethod KR
done

################################################################
## cool
################################################################

for file in $(ls *.cool)
do
filename=$(basename $file .cool)
cooler balance --force -p 8 -c 10000 $file
ln -s $file $filename.corrected.cool
done' > script.sh
qsub script.sh

cd ../..
done
