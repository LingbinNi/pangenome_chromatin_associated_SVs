#!/bin/bash
#written by Lingbin

# Build contact matrices in T2T-CHM13 coordinates.
# Replace /path/to/... with absolute paths on your system (without spaces).
# Input layout: <mapping_path>/<sample>/<unit>/<unit>.mapped.pairs.gz
# Each input file is processed separately; this script does not merge lanes.
# Input pairs must be sorted; the local copy is recompressed with bgzip before indexing.
# For uncompressed upstream output, run: bgzip <unit>.mapped.pairs
# Run once from a new output directory, separate from mapping_path.
# Adjust module names and SGE resources; Java, bgzip and cooler must be available.
# Chromosome sizes: <mapping_path>/DataBase/chr_size/T2T-CHM13v2.chr.size
# .cool/.h5 resolutions: 10, 25 and 100 kb; .hic uses Juicer defaults.
# The original KR correction and cooler balancing commands are retained.

mapping_path="/path/to/chm13_mapping"

for file in $(ls $mapping_path/*/*/*.mapped.pairs.gz)
do
sample=$(basename "$(dirname "$(dirname "$file")")")
filename=$(basename "$(dirname "$file")")
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

java -Xmx40000m -jar /path/to/juicer_tools_1.22.01.jar pre --threads 8 '$file' '$filename'.ReferenceMapping.CHM13.HPRC.hic '$mapping_path'/DataBase/chr_size/T2T-CHM13v2.chr.size

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

# The corrected .hic filename is a symbolic link, not an additional normalization step.
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
# Balancing stores weights in the existing .cool; the corrected name is an alias.
cooler balance --force -p 8 -c 10000 $file
ln -s $file $filename.corrected.cool
done' > script.sh
qsub script.sh

cd ../..
done
