#!/usr/bin/env Rscript
#written by Lingbin

library(ggplot2)

################################################################
## Classification 
################################################################

data=read.table("../evolution_origin_analysis/csa-sv.gt",header=T)
last5=data[, (ncol(data)-4):ncol(data)]
has_variant=apply(last5,1,function(x) any(grepl("(^|[|/])1([|/]|$)", x)))
data$variant_status=ifelse(has_variant,"with_variant","without_variant")
data$variant_status <- factor(data$variant_status,levels = c("with_variant", "without_variant"),labels = c("Ancestral", "Derived"))
data$length=as.numeric(sub(".*-(\\d+)$", "\\1", data$ID))

################################################################
## Output 
################################################################

for (ty in c("INS", "DEL")) {
  for (st in c("Ancestral", "Derived")) {
    sub <- data[data$SVTYPE == ty & data$variant_status == st, ]
    tag <- paste0(ty, ".", st)
    write.table(sub, file.path(paste0("Functional_",ty,"_",st,".tab")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("%-14s %6d\n", tag, nrow(sub)))
  }
}

################################################################
## Plot
################################################################

big_theme <- theme_bw() +
theme(legend.position = "none",
axis.text = element_text(size = 30),
axis.title = element_text(size = 30),
plot.subtitle = element_text(size = 18))
 
med <- function(x) if (length(x)) median(x) else NA_real_
 
make_plot <- function(d, kind, sub_lab) {
  p <- ggplot(d, aes(x = variant_status, y = log(length)))
  if (kind %in% c("violin", "both")) p <- p + geom_violin(width = 0.5)
  if (kind %in% c("box", "both"))
    p <- p + geom_boxplot(aes(fill = variant_status), width = 0.1, outlier.shape = NA)
  p + big_theme + labs(x = "", y = "Ln (SV size (bp))", subtitle = sub_lab)
}
 
for (ty in c("INS", "DEL")) {
  d <- data[data$SVTYPE == ty, ]
  if (nrow(d) == 0) { message("no data for ", ty); next }
 
  la <- d$length[d$variant_status == "Ancestral"]
  ld <- d$length[d$variant_status == "Derived"]
  pw <- if (length(la) >= 3 && length(ld) >= 3)
          wilcox.test(la, ld, exact = FALSE)$p.value else NA_real_
 
  cat(sprintf("%s | Ancestral n=%d median=%.0f bp | Derived n=%d median=%.0f bp | Wilcoxon P=%.3g\n",
              ty, length(la), med(la), length(ld), med(ld), pw))
  sub_lab <- sprintf("%s:  n = %d (Ancestral) vs %d (Derived);  Wilcoxon P = %.3g",
                     ty, length(la), length(ld), pw)
 
  files <- c(box = sprintf("Plot-%s-Size-1.pdf", ty),
             violin = sprintf("Plot-%s-Size-2.pdf", ty),
             both = sprintf("Plot-%s-Size.pdf", ty))
  for (kind in names(files)) {
    pdf(files[[kind]], width = 10, height = 10)
    print(make_plot(d, kind, sub_lab))
    dev.off()
  }
}
 
## 合并版: 一张图两个 panel, 更适合放进正文/Extended Data
pdf("Plot-Size-byType.pdf", width = 14, height = 8)
print(
  ggplot(data, aes(x = variant_status, y = log(length))) +
    geom_violin(width = 0.5) +
    geom_boxplot(aes(fill = variant_status), width = 0.1, outlier.shape = NA) +
    facet_wrap(~SVTYPE) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text  = element_text(size = 20),
          axis.title = element_text(size = 20),
          strip.text = element_text(size = 22)) +
    labs(x = "", y = "Ln (SV size (bp))")
)
dev.off()
