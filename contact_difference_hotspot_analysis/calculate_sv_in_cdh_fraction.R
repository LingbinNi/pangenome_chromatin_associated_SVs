#!/usr/bin/env Rscript
# Calculate the number and fraction of heterozygous SV records overlapping CDHs.
# Requires GenomicRanges. Edit cdhroot, svroot and outdir before running.
# Run: Rscript /path/to/calculate_sv_in_cdh_fraction.R
# CDH input: <cdhroot>/<sample>/CDH.<region|window>.<sample>.<Hap1|Hap2>.tsv
# CDH tables have a header containing Chr, Start and End, as written by the CDH caller.
# SV input: <svroot>/<sample>/<sample>.insdel.het.<hap1|hap2>.<INS|DEL>.pos.bed
# Supply prefiltered heterozygous SV intervals as headerless BED with at least 3 columns.
# SV file pairing:
#   Hap1 CDHs: hap1 INS and hap2 DEL.
#     Both represent sequence intervals present in Hap1 but absent from Hap2,
#     expressed in Hap1 assembly coordinates.
#   Hap2 CDHs: hap2 INS and hap1 DEL.
#     Both represent sequence intervals present in Hap2 but absent from Hap1,
#     expressed in Hap2 assembly coordinates.
# All input intervals use 0-based starts and half-open ends; zero-width SV sites
# are represented by a 1-bp interval beginning at the supplied start.
# Samples match HG* or NA*; VALID_CHR selects chr1_RagTag through chr22_RagTag.
# Both region- and window-level CDH files are analyzed separately.
# Keep SV records separate; count each SV overlapping at least one CDH once.
# sv_hit is the number of overlapping SVs and sv_pct is their percentage.
# n_SV is the number of supplied SV records retained after chromosome/coordinate filtering.
# The denominator is SV records, unlike the CDH-side test_cdh_sv_enrichment.R.
# Each permutation independently relocates CDHs on their original chromosomes,
# retaining their lengths. Bounds are the maximum observed ends of CDH, INS and DEL
# intervals for the sample/haplotype and CDH level being tested.
# Use 1,000 permutations and seed 1; enrichment = (observed + 0.5)/(random_mean + 0.5).
# Upper-tail empirical P = (1 + number of random counts >= observed)/(nperm + 1).
# Empirical P values are reported directly. Summaries are grouped by level, feature and haplotype.
# Pooled summaries give medians and IQRs of per-haplotype percentages across Hap1 and Hap2;
# they do not divide summed SV hit counts by summed SV counts.
# Outputs: SV_in_CDH_fraction_per_sample.tsv, SV_in_CDH_fraction_summary.tsv
# and SV_in_CDH_fraction_pooled.tsv.
# Use an existing empty output directory.

suppressMessages(library(GenomicRanges))

nperm <- 1000
set.seed(1)
outdir <- "./"

cdhroot <- "/path/to/contact_difference_hotspots"
svroot  <- "/path/to/heterozygous_sv_intervals"

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
  z <- which(d$end == d$start)            # Represent a zero-width site as a 1-bp interval.
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

# Set each chromosome bound to the maximum end in the supplied interval sets.
chr_lengths <- function(...) {
  gl <- list(...)
  ch <- unlist(lapply(gl, function(g) as.character(seqnames(g))), use.names = FALSE)
  en <- unlist(lapply(gl, function(g) end(g)), use.names = FALSE)
  if (!length(ch)) return(numeric(0))
  tapply(en, ch, max, na.rm = TRUE)
}

# Independently sample positions while preserving each CDH chromosome and width.
randomize <- function(gr, chrlen) {
  ch <- as.character(seqnames(gr)); w <- width(gr)
  L  <- as.numeric(chrlen[ch])
  ok <- !is.na(L) & w <= L
  if (!any(ok)) return(GRanges())
  ch <- ch[ok]; w <- w[ok]; L <- L[ok]
  st <- floor(runif(length(w)) * (L - w + 1)) + 1
  GRanges(ch, IRanges(as.integer(st), width = w))
}

# Count SV records overlapping at least one CDH; each SV is counted once.
test_sv_side <- function(cdh, sv, chrlen, fname, S, HAP, LV) {
  n_sv <- length(sv)
  if (length(cdh) == 0 || n_sv == 0)
    return(data.frame(sample = S, hap = HAP, level = LV, feature = fname,
                      n_CDH = length(cdh), n_SV = n_sv,
                      sv_hit = NA_integer_, sv_pct = NA_real_,
                      exp_hit = NA_real_, exp_sd = NA_real_, exp_pct = NA_real_,
                      enrichment = NA_real_,
                      z_score = NA_real_, p_enrich = NA_real_,
                      stringsAsFactors = FALSE))

  obs  <- sum(overlapsAny(sv, cdh, ignore.strand = TRUE))
  rand <- replicate(nperm, sum(overlapsAny(sv, randomize(cdh, chrlen),
                                           ignore.strand = TRUE)))
  rm_ <- mean(rand); rsd <- sd(rand)

  data.frame(sample = S, hap = HAP, level = LV, feature = fname,
             n_CDH = length(cdh), n_SV = n_sv,
             sv_hit = obs, sv_pct = round(100 * obs / n_sv, 3),
             exp_hit = round(rm_, 2), exp_sd = round(rsd, 2),
             exp_pct = round(100 * rm_ / n_sv, 3),
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
    # Combine INS and DEL records without merging intervals.
    sv <- c(ins, del)

    for (LV in LEVELS) {
      cdh <- read_cdh(file.path(cdhroot, S, sprintf("CDH.%s.%s.%s.tsv", LV, S, HAP)))
      if (length(cdh) == 0) next
      chrlen <- chr_lengths(cdh, ins, del)
      message("  ", S, " ", HAP, " ", LV, "  nCDH=", length(cdh),
              "  nINS=", length(ins), "  nDEL=", length(del))
      for (fn in c("INS", "DEL", "SV_all")) {
        ft <- switch(fn, INS = ins, DEL = del, SV_all = sv)
        res[[length(res) + 1]] <- test_sv_side(cdh, ft, chrlen, fn, S, HAP, LV)
      }
    }
  }
}

df <- do.call(rbind, res)
write.table(df, file.path(outdir, "SV_in_CDH_fraction_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## Summary
############################################################

d2 <- df[!is.na(df$sv_pct), ]
agg <- do.call(rbind, lapply(split(d2, list(d2$level, d2$feature, d2$hap), drop = TRUE),
  function(x) data.frame(
    level = x$level[1], feature = x$feature[1], hap = x$hap[1],
    n_sample_hap    = nrow(x),
    median_n_SV     = median(x$n_SV),
    median_n_CDH    = median(x$n_CDH),
    median_sv_hit   = median(x$sv_hit),
    median_sv_pct   = round(median(x$sv_pct), 3),
    q25_sv_pct      = round(quantile(x$sv_pct, 0.25), 3),
    q75_sv_pct      = round(quantile(x$sv_pct, 0.75), 3),
    median_exp_pct  = round(median(x$exp_pct), 3),
    median_enrich   = round(median(x$enrichment), 3),
    median_z        = round(median(x$z_score, na.rm = TRUE), 2),
    frac_p_lt_0.05  = round(mean(x$p_enrich < 0.05), 3),
    stringsAsFactors = FALSE)))
agg <- agg[order(agg$level, agg$feature, agg$hap), ]
write.table(agg, file.path(outdir, "SV_in_CDH_fraction_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Summarize per-haplotype percentages across Hap1 and Hap2 together.
pooled <- do.call(rbind, lapply(split(d2, list(d2$level, d2$feature), drop = TRUE),
  function(x) data.frame(
    level = x$level[1], feature = x$feature[1],
    n_hap = nrow(x),
    median_sv_pct = round(median(x$sv_pct), 2),
    IQR_low  = round(quantile(x$sv_pct, 0.25), 2),
    IQR_high = round(quantile(x$sv_pct, 0.75), 2),
    median_exp_pct = round(median(x$exp_pct), 2),
    stringsAsFactors = FALSE)))
pooled <- pooled[order(pooled$level, pooled$feature), ]
write.table(pooled, file.path(outdir, "SV_in_CDH_fraction_pooled.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=========== fraction of het SVs intersecting CDHs ===========\n")
print(agg, row.names = FALSE)
cat("\n---- pooled across haplotypes ----\n")
print(pooled, row.names = FALSE)
message("Done.")
