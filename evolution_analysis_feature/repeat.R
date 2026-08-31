#!/usr/bin/env Rscript
#written by Lingbin

library(GenomicRanges)
library(dplyr)
library(ggplot2)

repeat_dir="/net/eichler/vol28/projects/hic_cohorts/nobackups/00.DataBase/repeats/repeat_masker"
sv_dir="/net/eichler/vol28/projects/hgsvc/nobackups/analysis/genome_analysis/00.Analysis/contact_difference_enrichment_SVs/00.SVs-Bed-Filter"

################################################################
## data read 
################################################################

sv_ins=read.table("../../evolution_origin_analysis/Functional_INS_Ancestral.names",header=F)
sv_del=read.table("../../evolution_origin_analysis/Functional_DEL_Ancestral.names",header=F)
sv=rbind(sv_ins,sv_del)

################################################################
## sv match 
################################################################

samplelist=list.files(sv_dir,pattern="^(HG|NA)")
repeat_content_total=data.frame()

#for (i in 1:length(samplelist)){
#sample=samplelist[i]
sample="HG002"
#sv
hap1_ins=read.table(paste(sv_dir, sample, paste0(sample,".insdel.het.hap1.INS.pos.bed"),sep="/"), sep="\t", fill=T, header=F)
hap1_del=read.table(paste(sv_dir, sample, paste0(sample,".insdel.het.hap1.DEL.pos.bed"),sep="/"), sep="\t", fill=T, header=F)
hap2_ins=read.table(paste(sv_dir, sample, paste0(sample,".insdel.het.hap2.INS.pos.bed"),sep="/"), sep="\t", fill=T, header=F)
hap2_del=read.table(paste(sv_dir, sample, paste0(sample,".insdel.het.hap2.DEL.pos.bed"),sep="/"), sep="\t", fill=T, header=F)
hap1_sv=rbind(hap1_ins,hap2_del)
hap2_sv=rbind(hap2_ins,hap1_del)
#sv_filter
hap1_sv_filtered = hap1_sv[hap1_sv$V4 %in% sv$V1, ]
hap2_sv_filtered = hap2_sv[hap2_sv$V4 %in% sv$V1, ]
#sv_repeat
hap1_repeat=read.table(paste(repeat_dir, paste0(sample,".hap1.fasta_rm.bed"),sep="/"), header=F)
hap2_repeat=read.table(paste(repeat_dir, paste0(sample,".hap2.fasta_rm.bed"),sep="/"), header=F)

################################################################
## Hap1 
################################################################
#rename
colnames(hap1_sv_filtered) <- c("chr","start","end","svid")
colnames(hap1_repeat) <- c("chr","start","end","motif","len","strand","class","family","score","id")
#grange
sv_gr <- GRanges(
seqnames = hap1_sv_filtered$chr,
ranges = IRanges(
start = hap1_sv_filtered$start,
end = hap1_sv_filtered$end
),
svid = hap1_sv_filtered$svid
)
repeat_gr <- GRanges(
seqnames = hap1_repeat$chr,
ranges = IRanges(
start = hap1_repeat$start,
end = hap1_repeat$end
),
class = hap1_repeat$class
)
#hit
hits <- findOverlaps(sv_gr, repeat_gr)
#overlap_len
ov_width <- width(
pintersect(
sv_gr[queryHits(hits)],
repeat_gr[subjectHits(hits)]
)
)
#overlap_df
overlap_df <- data.frame(
svid = mcols(sv_gr)$svid[queryHits(hits)],
repeat_class = mcols(repeat_gr)$class[subjectHits(hits)],
overlap_bp = ov_width
)
#repeat_content
repeat_content <- overlap_df %>%
group_by(svid, repeat_class) %>%
summarise(total_bp = sum(overlap_bp), .groups="drop")
#sv_info
sv_len <- width(sv_gr)
sv_info <- data.frame(
svid = mcols(sv_gr)$svid,
sv_len = sv_len
)
#sv_content
repeat_content <- left_join(repeat_content, sv_info, by="svid")
repeat_content$percent <- repeat_content$total_bp / repeat_content$sv_len * 100
repeat_content$sample=paste(sample,"Hap1",sep=".")
repeat_content_total=rbind(repeat_content_total,repeat_content)

################################################################
## Hap2 
################################################################
#rename
colnames(hap2_sv_filtered) <- c("chr","start","end","svid")
colnames(hap2_repeat) <- c("chr","start","end","motif","len","strand","class","family","score","id")
#grange
sv_gr <- GRanges(
seqnames = hap2_sv_filtered$chr,
ranges = IRanges(
start = hap2_sv_filtered$start,
end = hap2_sv_filtered$end
),
svid = hap2_sv_filtered$svid
)
repeat_gr <- GRanges(
seqnames = hap2_repeat$chr,
ranges = IRanges(
start = hap2_repeat$start,
end = hap2_repeat$end
),
class = hap2_repeat$class
)
#hit
hits <- findOverlaps(sv_gr, repeat_gr)
#overlap_len
ov_width <- width(
pintersect(
sv_gr[queryHits(hits)],
repeat_gr[subjectHits(hits)]
)
)
#overlap_df
overlap_df <- data.frame(
svid = mcols(sv_gr)$svid[queryHits(hits)],
repeat_class = mcols(repeat_gr)$class[subjectHits(hits)],
overlap_bp = ov_width
)
#repeat_content
repeat_content <- overlap_df %>%
group_by(svid, repeat_class) %>%
summarise(total_bp = sum(overlap_bp), .groups="drop")
#sv_info
sv_len <- width(sv_gr)
sv_info <- data.frame(
svid = mcols(sv_gr)$svid,
sv_len = sv_len
)
#sv_content
repeat_content <- left_join(repeat_content, sv_info, by="svid")
repeat_content$percent <- repeat_content$total_bp / repeat_content$sv_len * 100
repeat_content$sample=paste(sample,"Hap2",sep=".")
repeat_content_total=rbind(repeat_content_total,repeat_content)

#}

################################################################
## Calculate 
################################################################

total_sv_bp <- repeat_content_total %>% distinct(sample, svid, sv_len) %>% summarise(total = sum(sv_len)) %>% pull(total)

repeat_summary_total <- repeat_content_total %>%
group_by(repeat_class) %>%
summarise(
total_repeat_bp = sum(total_bp),
.groups = "drop"
) %>%
mutate(
percent_of_all_sv =
total_repeat_bp / total_sv_bp * 100
) %>%
mutate(
percent_of_repeat =
total_repeat_bp /
sum(total_repeat_bp) * 100
)

write.table(repeat_summary_total,"repeat.tsv",row.names=F,quote=F,sep="\t")

################################################################
## Plot 
################################################################

pdf("repeat.svperc.pdf",height=10,width=10)
ggplot(repeat_summary_total,aes(x=reorder(repeat_class, percent_of_all_sv),y=percent_of_all_sv,fill=repeat_class)) +
geom_col() +
coord_flip() +
theme_bw() +
labs(x = "",y = "Percent of all SV (%)") +
theme(legend.position = "none")
dev.off()

pdf("repeat.repeatperc.pdf",height=10,width=10)
ggplot(repeat_summary_total,aes(x=reorder(repeat_class,percent_of_repeat),y=percent_of_repeat,fill=repeat_class)) +
geom_col() +
coord_flip() +
theme_bw() +
labs(x = "",y = "Percent of repeat (%)") +
theme(legend.position = "none")
dev.off()
