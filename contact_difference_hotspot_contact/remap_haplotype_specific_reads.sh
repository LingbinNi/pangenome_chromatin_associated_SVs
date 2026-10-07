#!/bin/bash
# written by Lingbin

# Remap reads from Hap1-specific contacts to the corresponding Hap2 assembly.
# Edit database_path, unique_reads_path and the raw FASTQ root in the loop below.
# Use absolute paths and sample names without whitespace or shell-special characters.
# Raw FASTQ layout: <raw_root>/<sample>/<unit>/<unit>_R1*.fastq.gz
# R2 is the matching filename with _R1 replaced by _R2.
# The unit directory name must equal the filename prefix before the first _R1.
# The current sample pattern is HG*; edit it to match your sample names.
# Read-ID input: <unique_reads_path>/<sample>/<unit>/Hap1.specific.dedup.readid
# Supply one readID per line, matching the read names in both paired FASTQ files.
# R1 and R2 inputs must contain corresponding read pairs in the same order.
# Hap2 FASTA: <database_path>/<sample>/database/reference/ragtag.scaffold.fasta
# Hap2 chromosome sizes: <database_path>/<sample>/database/chr_size/<sample>_HAP2.chr.size
# Prepare the FASTA BWA index and matching chromosome-size file before running.
# Each sample/unit is submitted as a separate SGE job in the current output directory.
# The job extracts R1/R2 with seqtk, then runs BWA and pairtools parse/sort/dedup/split.
# Outputs per unit: <unit>.mapped.pairs, stats.txt and dedup.dups.pairsam (uncompressed).
# The intermediate sorted BAM and extracted FASTQ files are removed after processing.
# Output layout: <working_directory>/<sample>/<unit>/
# Adjust module loads and SGE resources for your cluster; qsub must be available.
# Run in a new empty output directory: bash /path/to/remap_hap1_specific_reads_to_hap2.sh

current_path=$(pwd)
database_path="/path/to/hap2_reference_database"
unique_reads_path="/path/to/haplotype_specific_readids"

for dir in $(ls -d /path/to/raw_hic_fastq/HG*)
do
sample=$(basename $dir)
echo $sample
mkdir $sample
cd $sample

for file in $(ls $dir/*/*_R1*.fastq.gz)
do
filename1=$(basename $file)
filenamebase=${filename1%%_R1*}
filename2=${filename1/_R1/_R2}
echo $filenamebase
mkdir $filenamebase
cd $filenamebase

echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -pe serial 16
#$ -o output.log
#$ -e error.log

module load bedtools/2.31.1
module load deeptools/3.5.5
module load pairtools/1.0.3
module load bwa/0.7.17
module load samtools/1.19
module load preseq/3.2.0
module load miniconda/24.7.1
module load seqtk/1.4

################################################################
## Extract the selected paired reads
################################################################

seqtk subseq '$dir'/'$filenamebase'/'$filename1' '$unique_reads_path'/'$sample'/'$filenamebase'/Hap1.specific.dedup.readid > R1.extract.fastq
seqtk subseq '$dir'/'$filenamebase'/'$filename2' '$unique_reads_path'/'$sample'/'$filenamebase'/Hap1.specific.dedup.readid > R2.extract.fastq

################################################################
## Map to the Hap2 assembly and process contact pairs
################################################################

mkdir -p '$current_path'/'$sample'/'$filenamebase'/temp

bwa mem -5SP -T0 -t16 '$database_path'/'$sample'/database/reference/ragtag.scaffold.fasta R1.extract.fastq R2.extract.fastq| \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path '$database_path'/'$sample'/database/chr_size/'$sample'_HAP2.chr.size | \
pairtools sort --tmpdir='$current_path'/'$sample'/'$filenamebase'/temp/ --nproc 16|\
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --keep-parent-id --output-dups dedup.dups.pairsam|\
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs '$filenamebase'.mapped.pairs --output-sam -|\
samtools view -bS -@16 | samtools sort -@16 -T '$current_path'/'$sample'/'$filenamebase'/temp/temp.bam -o '$filenamebase'.mapped.PT.bam;
rm '$filenamebase'.mapped.PT.bam

################################################################
## Remove the extracted FASTQ files
################################################################

rm R*.extract.fastq
' > script.sh
qsub script.sh

cd ..
done

cd ..
done
