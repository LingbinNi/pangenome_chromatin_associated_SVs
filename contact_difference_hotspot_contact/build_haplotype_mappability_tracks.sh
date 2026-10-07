#!/bin/bash
# written by Lingbin

# Build genome-wide mappability tracks independently for each haplotype assembly.
# Replace /path/to/haplotype_assemblies below with the absolute FASTA directory.
# Expected filenames: <sample>.Hap1.scaffolds.fasta or <sample>.Hap2.scaffolds.fasta.
# The current sample pattern is {HG,NA}*; edit it to match your sample names.
# Sample and haplotype names are taken from the first two dot-separated fields.
# Use paths and sample names without whitespace or shell-special characters.
# Each assembly is indexed, then searched for minimum unique lengths in 20:200 bp.
# The search uses 10 threads; tracks are generated for a read length of 150 bp.
# Output directory: <working_directory>/<sample>.<haplotype>/
# Outputs include the index, unique_lengths/*.unique.uint8, 150.bed and 150.wig.
# 150.bed contains single-read mappability scores; 150.wig contains multi-read scores.
# For downstream contact annotation, use the assemblies corresponding to the mapping references.
# Requires newmap on PATH in the job environment and SGE qsub.
# Adjust the module load and SGE resources for your cluster.
# Run in a new empty output directory: bash /path/to/build_haplotype_mappability_tracks.sh

for file in $(ls /path/to/haplotype_assemblies/{HG,NA}*.Hap*.scaffolds.fasta)
do
filename=$(basename $file)
sample=$(basename $file | cut -d'.' -f1)
hap=$(basename $file | cut -d'.' -f2)
echo $sample.$hap
mkdir $sample.$hap
cd $sample.$hap

echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -l disk_free=10G
#$ -pe serial 10
#$ -o output.log
#$ -e error.log

module load miniconda/24.7.1

################################################################
## Create an index for a genome 
################################################################
ln -s '$file' .
newmap index '$filename'

################################################################
## Find the minimum unique k-mer lengths for the genome using the index 
################################################################

newmap search --verbose --num-threads=10 --search-range=20:200 --output-directory=unique_lengths '$filename'

################################################################
## Convert the unique lengths to mappability tracks 
################################################################

newmap track --single-read=150.bed --multi-read=150.wig 150 unique_lengths/*.unique.uint8
' > script.sh
qsub script.sh

cd ..
done
