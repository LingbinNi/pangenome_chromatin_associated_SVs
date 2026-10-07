#!/usr/bin/env Rscript
# Written by Lingbin.
# Plot correlations between insulation at each SV start and nearby-gene TPM.
# Requires ggplot2 and dplyr. Set dir below to an absolute input directory path.
# Run from an existing empty output directory:
# Rscript /path/to/plot_insulation_expression_correlations.R
#
# Expression input: <dir>/<sv_id>/<sample>/AdjacentGenes.tab.
# SV directory names start with chr and use chromosome-start-TYPE-length IDs.
# Sample directory names start with HG or NA. Tables have a header containing
# gene_id and TPM; supply one row per gene/sample with numeric, nonmissing TPM.
# Score input: <dir>/PlotInsulation/Files/<sample>_min75000_max250000_step25000_thres0.05_delta0.01_fdr_score.bedgraph.
# Score files are headerless with four columns: chromosome, start, end, score.
# Use matching reference coordinates and chromosome names for SV IDs and scores.
# Each 0-based SV start must match one interval with start <= position < end,
# supplying one numeric, nonmissing score for the corresponding sample.
#
# For each SV, append the sample's score to its expression rows and combine
# samples. For correlation plots, exclude gene IDs beginning with novel_gene.
# Test genes with at least three rows and nonzero variation in both variables.
# Calculate two-sided Pearson and Spearman correlations using raw TPM values;
# display their correlation coefficients and individual test P values.
# Plot a linear regression fit with its confidence band.
# Outputs under the current directory: <sv_id>/<sv_id>.tab (combined input
# rows with insulation scores) and <sv_id>/scatter_<gene_id>.pdf per tested gene.

library(ggplot2)

dir="/path/to/adjacent_gene_results"
varlist=list.files(dir,pattern="^chr")

for (i in 1:length(varlist)){
  var=varlist[i]
  filelist=list.files(paste(dir,var,sep="/"),pattern="^(HG|NA)")

  dir.create(var)
  setwd(var)

  dat=data.frame()
  for (j in 1:length(filelist)){
    sample=filelist[j]
    file=paste(dir,var,sample,"AdjacentGenes.tab",sep="/")
    if (file.exists(file)) {
      data=read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      insulation=read.table(paste(dir,"PlotInsulation","Files",paste0(sample,"_min75000_max250000_step25000_thres0.05_delta0.01_fdr_score.bedgraph"),sep="/"), header = F, sep = "\t", stringsAsFactors = FALSE)
      colnames(insulation)=c("chr","start","end","insulation")
      tmp <- strsplit(var, "-")[[1]]
      chr <- tmp[1]
      pos <- as.numeric(tmp[2])
      idx <- which(insulation$chr == chr & insulation$start <= pos & insulation$end > pos)
      sample_score <- insulation$insulation[idx]
      data$insulation=sample_score
      dat=rbind(dat,data)
    }
  }

  write.table(dat,file=paste(var,"tab",sep="."),quote=F,row.names=F,sep="\t")

  library(ggplot2)
  library(dplyr)

  df=dat

  genes <- unique(df$gene_id)
  genes_filtered <- genes[!grepl("^novel_gene", genes)]

  for(g in genes_filtered) {
    df_gene <- df %>% filter(gene_id == g)

    # Require at least three rows and nonzero variation in both variables.
    if(nrow(df_gene) < 3 || sd(df_gene$insulation, na.rm = TRUE) == 0 || sd(df_gene$TPM, na.rm = TRUE) == 0) {
      next
    }

    # Calculate Pearson and Spearman correlations.
    cor_pearson <- cor.test(df_gene$insulation, df_gene$TPM, method = "pearson")
    cor_spearman <- cor.test(df_gene$insulation, df_gene$TPM, method = "spearman")

    r_pearson <- round(cor_pearson$estimate, 3)
    p_pearson <- signif(cor_pearson$p.value, 3)

    r_spearman <- round(cor_spearman$estimate, 3)
    p_spearman <- signif(cor_spearman$p.value, 3)

    # Assemble the correlation labels displayed in the plot.
    annot_text <- paste0(
      "Pearson r = ", r_pearson, ", p = ", p_pearson, "\n",
      "Spearman rho = ", r_spearman, ", p = ", p_spearman
    )

    p <- ggplot(df_gene, aes(x = insulation, y = TPM)) +
      geom_point(color = "steelblue", size = 2, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 0.8) +
      annotate("text",
               x = -Inf, y = Inf,
               label = annot_text,
               hjust = -0.1, vjust = 1.3,
               size = 3.5) +
      labs(title = paste0("Gene: ", g),
           subtitle = paste0("Pearson r = ", r_pearson, " (p = ", p_pearson, ");  ",
                             "Spearman rho = ", r_spearman, " (p = ", p_spearman, ")"),
           x = "Insulation Score", y = "TPM") +
      theme_bw()

    ggsave(filename = paste0("scatter_", g, ".pdf"), plot = p, width = 5.5, height = 5)
  }

  setwd("..")

}
