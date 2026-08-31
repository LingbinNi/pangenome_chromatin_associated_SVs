#!/usr/bin/env Rscript
# =============================================================================
# CDH x SV 富集
#   配对: INS 同侧 (Hap1 <-> hap1 INS) / DEL 交叉 (Hap1 <-> hap2 DEL)
#   零模型: 保持宽度与所在染色体,在整条染色体内随机重投,nperm 次
# =============================================================================

suppressMessages(library(GenomicRanges))

nperm <- 1000
set.seed(1)
outdir <- "./"

cdhroot <- "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/contact_difference_hotspot_identification/04.ContactDifferenceHotspots"
svroot  <- "/net/eichler/vol28/projects/hgsvc/nobackups/analysis/genome_analysis/00.Analysis/contact_difference_enrichment_SVs/00.SVs-Bed-Filter"

VALID_CHR <- "^chr([1-9]|1[0-9]|2[0-2])_RagTag$"
LEVELS    <- c("region", "window")

############################################################
## Readers
############################################################

keep_valid <- function(gr) if (length(gr) == 0) gr else
  gr[grepl(VALID_CHR, as.character(seqnames(gr)))]

read_bed <- function(f) {
  if (!file.exists(f) || file.info(f)$size == 0) return(GRanges())
  d <- tryCatch(read.table(f, header = FALSE, sep = "\t", quote = "",
                           comment.char = "", fill = TRUE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || ncol(d) < 3) return(GRanges())
  d <- d[, 1:3]; names(d) <- c("chr", "start", "end")
  d$start <- suppressWarnings(as.integer(d$start))
  d$end   <- suppressWarnings(as.integer(d$end))
  d <- d[!is.na(d$chr) & !is.na(d$start) & !is.na(d$end), ]
  z <- which(d$end == d$start)            # INS 常为零宽度,补 1bp 而不是丢掉
  if (length(z)) d$end[z] <- d$start[z] + 1L
  d <- d[d$end > d$start, ]
  if (nrow(d) == 0) return(GRanges())
  keep_valid(GRanges(d$chr, IRanges(d$start + 1L, d$end)))
}

read_cdh <- function(f) {
  if (!file.exists(f) || file.info(f)$size == 0) return(GRanges())
  d <- tryCatch(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !all(c("Chr","Start","End") %in% names(d)))
    return(GRanges())
  keep_valid(GRanges(d$Chr, IRanges(as.integer(d$Start) + 1L, as.integer(d$End))))
}

############################################################
## Permutation
############################################################

chr_lengths <- function(...) {
  gl <- list(...)
  ch <- unlist(lapply(gl, function(g) as.character(seqnames(g))), use.names = FALSE)
  en <- unlist(lapply(gl, function(g) end(g)), use.names = FALSE)
  if (!length(ch)) return(numeric(0))
  tapply(en, ch, max, na.rm = TRUE)
}

randomize <- function(gr, chrlen) {
  ch <- as.character(seqnames(gr)); w <- width(gr)
  L  <- as.numeric(chrlen[ch])
  ok <- !is.na(L) & w <= L
  if (!any(ok)) return(GRanges())
  ch <- ch[ok]; w <- w[ok]; L <- L[ok]
  st <- floor(runif(length(w)) * (L - w + 1)) + 1
  GRanges(ch, IRanges(as.integer(st), width = w))
}

test_one <- function(cdh, feat, chrlen, fname, S, HAP, LV) {
  if (length(cdh) == 0 || length(feat) == 0)
    return(data.frame(sample = S, hap = HAP, level = LV, feature = fname,
                      n_CDH = length(cdh), n_SV = length(feat),
                      obs_hit = NA_integer_, obs_pct = NA_real_,
                      exp_hit = NA_real_, exp_sd = NA_real_, exp_pct = NA_real_,
                      enrichment = NA_real_, z_score = NA_real_,
                      p_enrich = NA_real_, stringsAsFactors = FALSE))
  feat <- reduce(feat, ignore.strand = TRUE)
  obs  <- sum(overlapsAny(cdh, feat, ignore.strand = TRUE))
  rand <- replicate(nperm, sum(overlapsAny(randomize(cdh, chrlen), feat,
                                           ignore.strand = TRUE)))
  rm_ <- mean(rand); rsd <- sd(rand); n <- length(cdh)
  data.frame(sample = S, hap = HAP, level = LV, feature = fname,
             n_CDH = n, n_SV = length(feat),
             obs_hit = obs, obs_pct = round(100 * obs / n, 2),
             exp_hit = round(rm_, 2), exp_sd = round(rsd, 2),
             exp_pct = round(100 * rm_ / n, 2),
             enrichment = round((obs + 0.5) / (rm_ + 0.5), 3),
             z_score = if (rsd > 0) round((obs - rm_) / rsd, 2) else NA_real_,
             p_enrich = (sum(rand >= obs) + 1) / (nperm + 1),
             stringsAsFactors = FALSE)
}

############################################################
## Main
############################################################

samples <- sort(basename(list.dirs(cdhroot, full.names = FALSE, recursive = FALSE)))
samples <- samples[grepl("^(HG|NA)", samples)]
message("Total samples found: ", length(samples))

res <- list()
for (S in samples) {
  for (HAP in c("Hap1", "Hap2")) {
    self  <- if (HAP == "Hap1") "hap1" else "hap2"
    other <- if (HAP == "Hap1") "hap2" else "hap1"

    ins <- read_bed(file.path(svroot, S, sprintf("%s.insdel.het.%s.INS.pos.bed", S, self)))
    del <- read_bed(file.path(svroot, S, sprintf("%s.insdel.het.%s.DEL.pos.bed", S, other)))
    if (length(ins) == 0 && length(del) == 0) next
    sv <- reduce(c(ins, del), ignore.strand = TRUE)

    for (LV in LEVELS) {
      cdh <- read_cdh(file.path(cdhroot, S, sprintf("CDH.%s.%s.%s.tsv", LV, S, HAP)))
      if (length(cdh) == 0) next
      chrlen <- chr_lengths(cdh, ins, del)
      message("  ", S, " ", HAP, " ", LV, "  nCDH=", length(cdh))
      for (fn in c("INS", "DEL", "SV_all")) {
        ft <- switch(fn, INS = ins, DEL = del, SV_all = sv)
        res[[length(res) + 1]] <- test_one(cdh, ft, chrlen, fn, S, HAP, LV)
      }
    }
  }
}

df <- do.call(rbind, res)
write.table(df, file.path(outdir, "CDH_SV_enrichment_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## Summary
############################################################

d2 <- df[!is.na(df$enrichment), ]
agg <- do.call(rbind, lapply(split(d2, list(d2$level, d2$feature, d2$hap), drop = TRUE),
  function(x) data.frame(
    level = x$level[1], feature = x$feature[1], hap = x$hap[1],
    n_sample_hap = nrow(x),
    median_n_CDH = median(x$n_CDH),
    median_obs_pct = round(median(x$obs_pct), 2),
    median_exp_pct = round(median(x$exp_pct), 2),
    median_enrichment = round(median(x$enrichment), 3),
    mean_enrichment = round(mean(x$enrichment), 3),
    median_z = round(median(x$z_score, na.rm = TRUE), 2),
    frac_p_lt_0.05 = round(mean(x$p_enrich < 0.05), 3),
    stringsAsFactors = FALSE)))
agg <- agg[order(agg$level, agg$feature, agg$hap), ]
write.table(agg, file.path(outdir, "CDH_SV_enrichment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n================ CDH x SV enrichment ================\n")
print(agg, row.names = FALSE)
cat("\n样本数:", length(unique(df$sample)),
    " 单倍型条目:", nrow(df) / (length(LEVELS) * 3), "\n")
message("Done.")
