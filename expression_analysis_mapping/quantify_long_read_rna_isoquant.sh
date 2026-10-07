#!/bin/bash
# Written by Lingbin.
# Quantify genes and transcripts from aligned PacBio CCS RNA reads with IsoQuant.
# Requires Bash, SGE qsub and an environment providing IsoQuant 3.10.0.
# Edit the three absolute input paths below; adapt module loading and SGE resources
# to your computing environment.
# Run from a new output directory: bash /path/to/quantify_long_read_rna_isoquant.sh
# BAM input: /path/to/aligned_bams/<sample>/<sample>.transcripts.pbmm2.bam
# Supply coordinate-sorted, indexed BAMs aligned to the supplied reference genome.
# Reference FASTA, GTF annotation and BAMs must use matching assembly coordinates
# and chromosome names; the example filenames refer to T2T-CHM13 v2.0.
# Sample names are obtained by removing .transcripts.pbmm2.bam from each basename.
# Use unique sample names and paths without whitespace or shell-special characters.
# Each BAM creates <sample>/<sample>.sh and submits one SGE job with qsub.
# Jobs run in <sample>/ and write scheduler logs to <sample>.o and <sample>.e.
# IsoQuant outputs are written under <sample>/ using its default output naming.
# The command supplies reference annotation and retains default transcript-model
# construction, read assignment and quantification settings.

for file in $(ls /path/to/aligned_bams/*/*.transcripts.pbmm2.bam)
do
sample=$(basename $file .transcripts.pbmm2.bam)
echo $sample
mkdir $sample
cd $sample

echo '#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mfree=10G
#$ -pe serial 16
#$ -o '$sample'.o
#$ -e '$sample'.e

module load isoquant/3.10.0

isoquant.py \
--reference /path/to/T2T-CHM13v2.fasta \
--genedb /path/to/chm13v2.0_RefSeq_Liftoff_v5.2.sorted.gtf \
--data_type pacbio_ccs \
--bam '$file' \
-o ./
' > $sample.sh
qsub $sample.sh

cd ..
done
