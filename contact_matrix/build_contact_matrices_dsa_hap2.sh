#!/bin/bash
#written by Lingbin

# Build contact matrices in donor-specific Hap2 coordinates.
# Replace /path/to/... with absolute paths on your system (without spaces).
# Input layout: <mapping_path>/<sample>/<unit>/<unit>.mapped.pairs.gz
# Each input file is processed separately; this script does not merge lanes.
# Input pairs must be sorted and BGZF-compressed for pairix indexing.
# For uncompressed upstream output, run: bgzip <unit>.mapped.pairs
# Run once from a new output directory, separate from mapping_path.
# Adjust module names and SGE resources; Java, bgzip and cooler must be available.
# Chromosome sizes: <database_path>/<sample>/database/chr_size/<sample>_HAP2.chr.size
# The chromosome-size file must match the assembly used for mapping.
# .cool/.h5 resolutions: 10, 25 and 100 kb; .hic uses Juicer defaults.
# This script creates raw .cool/.h5 matrices; it has no explicit balancing step.

database_path="/path/to/hap2_reference_data"
mapping_path="/path/to/hap2_mapping"
hap="2"

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

ln -s '$file' '$filename'.mapped.pairs.gz

################################################################
## hic
################################################################

java -Xmx40000m -jar /path/to/juicer_tools_1.22.01.jar pre --threads 8 '$file' '$filename'.SelfMapping.Hap'$hap'.HPRC.hic '$database_path'/'$sample'/database/chr_size/'$sample'_HAP'$hap'.chr.size

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
