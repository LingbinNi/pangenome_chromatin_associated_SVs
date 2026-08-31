#!/bin/bash
#written by Lingbin

module load bcftools

################################################################
## CSA-SVs 
################################################################

awk -F"\t" 'NR>1{print $5}' /net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/structure_associated_SV_identification/cdh_sv.d0.keep.tsv > csa-sv.names

################################################################
## CSA-SVs 
################################################################

{ printf "CHROM\tPOS\tID\tSVTYPE\tSVLEN\tPTR\tPPA\tGGO\tPPY\tPAB\n"
  bcftools query \
    -i 'ID=@csa-sv.names' \
    -s PTR,PPA,GGO,PPY,PAB \
    -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN[\t%GT]\n' \
    /net/eichler/vol28/projects/hgsvc/nobackups/analysis/genome_analysis/00.Analysis/contact_difference_enrichment_SVs/00.Jiadong-complete/complete_gt.vcf
} > csa-sv.gt

################################################################
## Names 
################################################################

awk -F"\t" 'NR>1{print $3}' Functional_INS_Ancestral.tab > Functional_INS_Ancestral.names
awk -F"\t" 'NR>1{print $3}' Functional_DEL_Ancestral.tab > Functional_DEL_Ancestral.names
awk -F"\t" 'NR>1{print $3}' Functional_INS_Derived.tab > Functional_INS_Derived.names
awk -F"\t" 'NR>1{print $3}' Functional_DEL_Derived.tab > Functional_DEL_Derived.names
