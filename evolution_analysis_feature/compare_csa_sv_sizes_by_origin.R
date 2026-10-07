#!/usr/bin/env Rscript
# Written by Lingbin.
# Compare CSA-SV sizes between Ancestral and Derived origin classes.
# Requires ggplot2. Supply the origin table produced by the origin-classification
# analysis; this script uses those existing classifications.
# Run from an existing empty output directory:
# Rscript /path/to/compare_csa_sv_sizes_by_origin.R /path/to/csa-sv.origin.tsv
# Without an argument, input defaults to evolution_origin_output/csa-sv.origin.tsv
# relative to the current working directory. All outputs use the current directory.
#
# Input: headered TSV with one row per SV and columns ID, SVTYPE, SVLEN and origin.
# SVTYPE values are INS or DEL; SVLEN contains numeric lengths in bp.
# Retain Ancestral and Derived rows, take abs(SVLEN), and keep nonmissing lengths > 0.
# Summarize counts and median lengths for all retained SVs and separately for INS/DEL.
# Compare the two origin groups using a two-sided, unpaired Wilcoxon rank-sum test
# on lengths in bp (exact = FALSE), requiring at least three SVs in each group.
# Report individual test P values directly; comparisons below this count have NA P.
# Plots use the natural logarithm of length and include all retained rows.
# Boxplots hide outlier points; violin plots show the distribution of log lengths.
#
# Outputs: Functional_<INS|DEL>_<Ancestral|Derived>.tab and sv_size_summary.tsv.
# For each SV type with data: Plot-<type>-Size-1.pdf (boxplot),
# Plot-<type>-Size-2.pdf (violin) and Plot-<type>-Size.pdf (both).
# Plot-Size-byType.pdf combines the two SV types in separate panels.

library(ggplot2)

################################################################
## Input
################################################################

args <- commandArgs(trailingOnly = TRUE)

# Usage:
#   Rscript compare_csa_sv_sizes_by_origin.R /path/to/csa-sv.origin.tsv
#
# If no argument is supplied, look in the local evolution output directory.
origin_file <- if (length(args) >= 1) {
  args[1]
} else {
  "evolution_origin_output/csa-sv.origin.tsv"
}

if (!file.exists(origin_file)) {
  stop("Cannot find origin file: ", origin_file)
}

################################################################
## Classification
################################################################

data <- read.table(
  origin_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("ID", "SVTYPE", "SVLEN", "origin") %in% colnames(data))) {
  stop("Input must contain ID, SVTYPE, SVLEN and origin columns")
}

# Use Ancestral and Derived SVs for the size comparison.
data <- data[data$origin %in% c("Ancestral", "Derived"), ]

data$variant_status <- factor(
  data$origin,
  levels = c("Ancestral", "Derived")
)

data$length <- abs(as.numeric(data$SVLEN))

data <- data[!is.na(data$length) & data$length > 0, ]

cat("Confidently assigned CSA-SVs:", nrow(data), "\n")
cat("Ancestral:", sum(data$variant_status == "Ancestral"), "\n")
cat("Derived:", sum(data$variant_status == "Derived"), "\n\n")

################################################################
## Output
################################################################

for (ty in c("INS", "DEL")) {
  for (st in c("Ancestral", "Derived")) {
    sub <- data[data$SVTYPE == ty & data$variant_status == st, ]
    tag <- paste0(ty, ".", st)
    write.table(
      sub,
      file.path(paste0("Functional_", ty, "_", st, ".tab")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    cat(sprintf("%-14s %6d\n", tag, nrow(sub)))
  }
}

################################################################
## Statistics
################################################################

med <- function(x) if (length(x)) median(x) else NA_real_

size_stats <- function(d, label) {
  la <- d$length[d$variant_status == "Ancestral"]
  ld <- d$length[d$variant_status == "Derived"]

  pw <- if (length(la) >= 3 && length(ld) >= 3) {
    wilcox.test(ld, la, exact = FALSE)$p.value
  } else {
    NA_real_
  }

  out <- data.frame(
    group = label,
    n_Ancestral = length(la),
    median_Ancestral_bp = med(la),
    n_Derived = length(ld),
    median_Derived_bp = med(ld),
    wilcoxon_P = pw
  )

  cat(sprintf(
    "%s | Ancestral n=%d median=%.0f bp | Derived n=%d median=%.0f bp | Wilcoxon P=%.3g\n",
    label, length(la), med(la), length(ld), med(ld), pw
  ))

  out
}

stat_all <- size_stats(data, "ALL")
stat_ins <- size_stats(data[data$SVTYPE == "INS", ], "INS")
stat_del <- size_stats(data[data$SVTYPE == "DEL", ], "DEL")

stat <- rbind(stat_all, stat_ins, stat_del)

write.table(
  stat,
  "sv_size_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

################################################################
## Plot
################################################################

big_theme <- theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 30),
    plot.subtitle = element_text(size = 18)
  )

make_plot <- function(d, kind, sub_lab) {
  p <- ggplot(d, aes(x = variant_status, y = log(length)))
  if (kind %in% c("violin", "both")) {
    p <- p + geom_violin(width = 0.5)
  }
  if (kind %in% c("box", "both")) {
    p <- p + geom_boxplot(
      aes(fill = variant_status),
      width = 0.1,
      outlier.shape = NA
    )
  }
  p + big_theme +
    labs(
      x = "",
      y = "Ln (SV size (bp))",
      subtitle = sub_lab
    )
}

for (ty in c("INS", "DEL")) {
  d <- data[data$SVTYPE == ty, ]
  if (nrow(d) == 0) {
    message("no data for ", ty)
    next
  }

  s <- stat[stat$group == ty, ]

  sub_lab <- sprintf(
    "%s: n = %d (Ancestral) vs %d (Derived); Wilcoxon P = %.3g",
    ty, s$n_Ancestral, s$n_Derived, s$wilcoxon_P
  )

  files <- c(
    box = sprintf("Plot-%s-Size-1.pdf", ty),
    violin = sprintf("Plot-%s-Size-2.pdf", ty),
    both = sprintf("Plot-%s-Size.pdf", ty)
  )

  for (kind in names(files)) {
    pdf(files[[kind]], width = 10, height = 10)
    print(make_plot(d, kind, sub_lab))
    dev.off()
  }
}

## Combined INS/DEL plot
pdf("Plot-Size-byType.pdf", width = 14, height = 8)
print(
  ggplot(data, aes(x = variant_status, y = log(length))) +
    geom_violin(width = 0.5) +
    geom_boxplot(
      aes(fill = variant_status),
      width = 0.1,
      outlier.shape = NA
    ) +
    facet_wrap(~SVTYPE) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 22)
    ) +
    labs(x = "", y = "Ln (SV size (bp))")
)
dev.off()

message("Done.")
