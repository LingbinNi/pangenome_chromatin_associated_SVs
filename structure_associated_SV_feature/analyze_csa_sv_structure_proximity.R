#!/usr/bin/env Rscript

# Analyze CSA-SV proximity to TAD boundaries and loop anchors in haplotype assemblies.
# Requires data.table, GenomicRanges, IRanges and ggplot2.
# Edit the input paths, outdir and expected input counts below before running.
# Run: Rscript /path/to/analyze_csa_sv_structure_proximity.R
# both require a header containing sv_id. f_csa must be a subset of f_all.
# CSA contains IDs in f_csa; nonCSA contains the remaining IDs in f_all.
# ContactCount input: <contact_dir>/<sample>.insdel.het.<hap1|hap2>.<INS|DEL>.sum.ContactCount
# Each file has a header containing chr, start, end and variant_id.
# Provide one record per SV per sample/haplotype; contact-count columns are not used.
# variant_id must exactly match sv_id in the association tables (CHM13 locus IDs).
# ContactCount chr/start/end use the assembly coordinates of the haplotype being tested,
# matching its TAD and loop files. These supplied coordinates link each ID to the assembly.
# SV file pairing:
#   Hap1 structures: hap1 INS and hap2 DEL, representing intervals present in Hap1
#     but absent from Hap2, expressed in Hap1 assembly coordinates.
#   Hap2 structures: hap2 INS and hap1 DEL, representing intervals present in Hap2
#     but absent from Hap1, expressed in Hap2 assembly coordinates.
# TAD input: <tad_roots[[hap]]>/<sample>/*_boundaries.bed, headerless, at least 3 columns.
# Loop input: <loop_roots[[hap]]>/<sample>/*.bedgraph, headerless, at least 6 BEDPE columns.
# Supply one matching TAD file and one matching loop file per sample/haplotype.
# BED intervals are 0-based and half-open; zero-width SV sites become 1-bp intervals.
# Cohort inputs: <cdh_bed_dir>/<sample>.<Hap1|Hap2>.CDH.CHM13.bed filenames and
# a headered sv_carrier_file containing sample. CDH contents are not read here.
# Retain samples with both haplotype filenames and an entry in the carrier table,
# then restrict to sample subdirectories shared by the two TAD roots.
# Merge overlapping or adjacent TAD intervals and loop anchors within each feature set.
# Any_structure is the union of TAD boundaries and loop anchors.
# Distances are interval gaps to the nearest feature on the same chromosome;
# overlapping or adjacent intervals have distance 0; unmatched chromosomes have Inf.
# Thresholds are 0, 5, 10, 25, 50 and 100 kb; MAIN_DIST defaults to 25 kb.
# With do_perm=TRUE, run 200 permutations using seed 1. Independently relocate SV
# intervals on their original chromosomes, retaining lengths. Per-chromosome bounds
# are the maximum observed ends across SV, TAD and loop intervals for that sample/haplotype.
# For each sample/haplotype/feature, pool CSA and nonCSA intervals to calculate the
# random expected fraction. Both groups use this same expectation at each threshold.
# Enrichment = (observed fraction + 1e-6)/(mean random fraction + 1e-6).
# Summarize fractions across sample/haplotype groups and proximity per unique SV.
# Per-SV summaries use available carrier-haplotype records and Any_structure at MAIN_DIST:
# near_any: >=1 nearby record; near_majority: >=50%; near_most: >=80%.
# Outputs: CSA_structure_proximity_per_sample.tsv, CSA_structure_proximity_per_sv_obs.tsv.gz,
# summary_by_haplotype.tsv, summary_by_sv.tsv and summary_by_sv_overall.tsv.
# Figures: proximity_vs_distance.pdf, enrichment_boxplot.pdf and distance_ecdf.pdf.
# All TSVs have headers. Use a new output directory.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
})

############################################################
## Parameters
############################################################

# Expected input counts for the supplied dataset; edit for a different dataset.
# expected_n_samples is checked before restricting to samples in the TAD roots.
expected_n_csa <- 2495
expected_n_samples <- 176

set.seed(1)
nperm   <- 200          # Permutations per sample/haplotype/feature.
do_perm <- TRUE

outdir <- "./CSA_structure_proximity"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## Distance thresholds in bp; MAIN_DIST selects the per-SV summary threshold.
dist_list <- c(0, 5000, 10000, 25000, 50000, 100000)
MAIN_DIST <- 25000

contact_dir <- "/path/to/contact_counts"

# f25 locus-level CSA-SV results
f_all <- "/path/to/cdh_sv.local.f25000.tsv"
f_csa <- "/path/to/cdh_sv.local.f25000.keep.tsv"

# Inputs used to select paired sample IDs from CDH filenames and the SV carrier table.
cdh_bed_dir <- "/path/to/CDH_BED_CHM13"
sv_carrier_file <- "/path/to/SV_CHM13_hap_carrier.tsv"

tad_roots <- list(
  Hap1 = "/path/to/Hap1/Tad",
  Hap2 = "/path/to/Hap2/Tad"
)
loop_roots <- list(
  Hap1 = "/path/to/Hap1/Loop",
  Hap2 = "/path/to/Hap2/Loop"
)

############################################################
## Input readers
############################################################

find_one_file <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0) return(NA_character_)
  f[1]
}

read_tad_boundary <- function(file) {
  dt <- fread(file, header = FALSE)
  dt <- dt[!is.na(V1) & !is.na(V2) & !is.na(V3)][V3 > V2]
  if (nrow(dt) == 0) return(GRanges())
  reduce(GRanges(dt$V1, IRanges(dt$V2 + 1, dt$V3)), ignore.strand = TRUE)
}

read_loop_anchors <- function(file) {
  dt <- fread(file, header = FALSE)
  dt <- dt[!is.na(V1) & !is.na(V2) & !is.na(V3) & !is.na(V4) & !is.na(V5) & !is.na(V6)]
  dt <- dt[V3 > V2 & V6 > V5]
  if (nrow(dt) == 0) return(GRanges())
  a <- unique(rbind(dt[, .(chr = V1, start = V2, end = V3)],
                    dt[, .(chr = V4, start = V5, end = V6)]))
  reduce(GRanges(a$chr, IRanges(a$start + 1, a$end)), ignore.strand = TRUE)
}

## SV intervals on one haplotype: assembly coordinates linked to a CHM13 variant_id.
read_sv_on_hap <- function(sample, hap) {
  spec <- if (hap == "Hap1")
    data.table(f = c(sprintf("%s.insdel.het.hap1.INS.sum.ContactCount", sample),
                     sprintf("%s.insdel.het.hap2.DEL.sum.ContactCount", sample)),
               svtype = c("INS", "DEL"))
  else
    data.table(f = c(sprintf("%s.insdel.het.hap2.INS.sum.ContactCount", sample),
                     sprintf("%s.insdel.het.hap1.DEL.sum.ContactCount", sample)),
               svtype = c("INS", "DEL"))
  out <- list()
  for (i in seq_len(nrow(spec))) {
    p <- file.path(contact_dir, spec$f[i])
    if (!file.exists(p) || file.info(p)$size == 0) next
    d <- tryCatch(fread(p, header = TRUE), error = function(e) NULL)
    if (is.null(d) || nrow(d) == 0) next
    if (!all(c("chr", "start", "end", "variant_id") %in% names(d))) next
    out[[length(out) + 1]] <- data.table(
      chr = d$chr, start = as.numeric(d$start), end = as.numeric(d$end),
      variant_id = d$variant_id, svtype = spec$svtype[i])
  }
  if (!length(out)) return(NULL)
  r <- rbindlist(out)
  r <- r[!is.na(start) & !is.na(end) & end >= start]
  if (nrow(r) == 0) return(NULL)
  r
}

############################################################
## Permutation helper
############################################################

chr_len_from <- function(...) {
  gl <- list(...)
  ch <- unlist(lapply(gl, function(g) as.character(seqnames(g))), use.names = FALSE)
  en <- unlist(lapply(gl, function(g) end(g)), use.names = FALSE)
  if (!length(ch)) return(numeric(0))
  tapply(en, ch, max, na.rm = TRUE)
}

randomize_by_chr <- function(gr, chrlen) {
  ch <- as.character(seqnames(gr)); w <- width(gr)
  L <- as.numeric(chrlen[ch])
  ok <- !is.na(L) & w <= L
  if (!any(ok)) return(GRanges())
  ch <- ch[ok]; w <- w[ok]; L <- L[ok]
  st <- floor(runif(length(w)) * (L - w + 1)) + 1
  GRanges(ch, IRanges(as.integer(st), width = w))
}

## Nearest interval-gap distance; return Inf when no feature is available on that chromosome.
nearest_dist <- function(q, s) {
  if (length(q) == 0 || length(s) == 0) return(rep(Inf, length(q)))
  h <- distanceToNearest(q, s, ignore.strand = TRUE)
  d <- rep(Inf, length(q))
  d[queryHits(h)] <- mcols(h)$distance
  d
}

############################################################
## 1. CSA and nonCSA sets
############################################################

stopifnot(file.exists(f_all), file.exists(f_csa))
tested_ids <- unique(fread(f_all, select = "sv_id")$sv_id)
csa_ids    <- unique(fread(f_csa, select = "sv_id")$sv_id)
message("f25 local-CDH SV: ", length(tested_ids), "   CSA-SV: ", length(csa_ids))
stopifnot(length(csa_ids) == expected_n_csa)

############################################################
## 2. Analyze each sample and haplotype
############################################################

# Select paired sample IDs shared by the CDH filenames and the carrier table.
cdh_files <- list.files(cdh_bed_dir, pattern = "\\.CDH\\.CHM13\\.bed$")
cdh_haps <- sub("^(.+)\\.(Hap[12])\\.CDH\\.CHM13\\.bed$", "\\1_\\2", cdh_files)
cdh_haps <- sort(unique(cdh_haps))

sv_dt <- fread(sv_carrier_file, select = "sample")
sv_haps <- sort(unique(c(
  paste0(sv_dt$sample, "_Hap1"),
  paste0(sv_dt$sample, "_Hap2")
)))

common_haps <- intersect(cdh_haps, sv_haps)
samp_of <- sub("_Hap[12]$", "", common_haps)
samples_all <- sort(names(which(table(samp_of) == 2L)))
stopifnot(length(samples_all) == expected_n_samples)

# Restrict to the cohort represented by both TAD roots; run separate cohorts separately.
cohort_samples <- intersect(
  list.dirs(tad_roots$Hap1, full.names = FALSE, recursive = FALSE),
  list.dirs(tad_roots$Hap2, full.names = FALSE, recursive = FALSE)
)
samples <- sort(intersect(samples_all, cohort_samples))
message("Samples in the selected cohort: ", length(samples))

res_rows <- list()   # One row per sample/haplotype/feature/threshold/group.
sv_rows  <- list()   # One row per available SV/haplotype record for Any_structure at MAIN_DIST.

for (S in samples) {
  for (HAP in c("Hap1", "Hap2")) {

    sv <- read_sv_on_hap(S, HAP)
    if (is.null(sv)) next
    sv <- sv[variant_id %chin% tested_ids]                 # Retain IDs from the f25000 local-CDH association table.
    if (nrow(sv) == 0) next
    sv[, group := fifelse(variant_id %chin% csa_ids, "CSA", "nonCSA")]

    tad_f  <- find_one_file(file.path(tad_roots[[HAP]],  S), "_boundaries\\.bed$")
    loop_f <- find_one_file(file.path(loop_roots[[HAP]], S), "\\.bedgraph$")
    if (is.na(tad_f) || is.na(loop_f)) { message("skip (no TAD/loop): ", S, " ", HAP); next }
    tad  <- read_tad_boundary(tad_f)
    loop <- read_loop_anchors(loop_f)
    if (length(tad) == 0 && length(loop) == 0) next
    feats <- list(TAD_boundary = tad, Loop_anchor = loop,
                  Any_structure = reduce(c(tad, loop), ignore.strand = TRUE))

    gr_sv  <- GRanges(sv$chr, IRanges(sv$start + 1, pmax(sv$end, sv$start + 1)))
    chrlen <- chr_len_from(gr_sv, tad, loop)

    message("  ", S, " ", HAP, "  nSV=", nrow(sv), "  CSA=", sum(sv$group == "CSA"))

    for (fn in names(feats)) {
      ft <- feats[[fn]]
      if (length(ft) == 0) next
      dobs <- nearest_dist(gr_sv, ft)

      # Pool all selected SV intervals to estimate one random fraction shared by both groups.
      exp_frac <- NULL
      if (do_perm) {
        pm <- matrix(NA_real_, nrow = nperm, ncol = length(dist_list))
        for (k in seq_len(nperm)) {
          dr <- nearest_dist(randomize_by_chr(gr_sv, chrlen), ft)
          pm[k, ] <- vapply(dist_list, function(dd) mean(dr <= dd), numeric(1))
        }
        exp_frac <- colMeans(pm)
        exp_sd   <- apply(pm, 2, sd)
      }

      for (g in c("CSA", "nonCSA")) {
        ix <- which(sv$group == g)
        if (!length(ix)) next
        for (j in seq_along(dist_list)) {
          dd <- dist_list[j]
          nn <- sum(dobs[ix] <= dd)
          res_rows[[length(res_rows) + 1]] <- data.table(
            sample = S, hap = HAP, feature = fn, dist = dd, group = g,
            n_SV = length(ix), n_near = nn, frac_near = nn / length(ix),
            exp_frac = if (do_perm) exp_frac[j] else NA_real_,
            enrichment = if (do_perm) (nn / length(ix) + 1e-6) / (exp_frac[j] + 1e-6) else NA_real_,
            median_dist = median(dobs[ix][is.finite(dobs[ix])]))
        }
      }

      if (fn == "Any_structure")
        sv_rows[[length(sv_rows) + 1]] <- data.table(
          variant_id = sv$variant_id, svtype = sv$svtype, group = sv$group,
          sample = S, hap = HAP, dist = dobs,
          near = dobs <= MAIN_DIST)
    }
  }
}

res <- rbindlist(res_rows)
svd <- rbindlist(sv_rows)
stopifnot(nrow(res) > 0)
fwrite(res, file.path(outdir, "CSA_structure_proximity_per_sample.tsv"), sep = "\t")
fwrite(svd, file.path(outdir, "CSA_structure_proximity_per_sv_obs.tsv.gz"), sep = "\t")

############################################################
## 3. Summary A: medians across sample/haplotype fractions
############################################################

sumA <- res[, .(
  n_sample_hap   = .N,
  median_frac    = median(frac_near),
  IQR_frac       = paste(round(quantile(frac_near, c(.25, .75)), 3), collapse = "-"),
  median_exp     = median(exp_frac, na.rm = TRUE),
  median_enrich  = median(enrichment, na.rm = TRUE),
  median_dist_bp = median(median_dist, na.rm = TRUE)
), by = .(feature, dist, group)][order(feature, dist, group)]
fwrite(sumA, file.path(outdir, "summary_by_haplotype.tsv"), sep = "\t")

############################################################
## 4. Summary B: proximity per unique SV
############################################################
## For each SV, calculate the fraction of available carrier-haplotype records
## within MAIN_DIST of Any_structure; apply the any, >=50% and >=80% rules.

cons <- svd[, .(n_carrier = .N, n_near = sum(near),
                frac_carrier_near = mean(near),
                median_dist = median(dist[is.finite(dist)])),
            by = .(variant_id, svtype, group)]
cons[, near_any      := n_near >= 1]
cons[, near_majority := frac_carrier_near >= 0.5]
cons[, near_most     := frac_carrier_near >= 0.8]
fwrite(cons, file.path(outdir, "summary_by_sv.tsv"), sep = "\t")

sumB <- cons[, .(
  n_SV = .N,
  pct_near_any      = round(100 * mean(near_any), 1),
  pct_near_majority = round(100 * mean(near_majority), 1),
  pct_near_most     = round(100 * mean(near_most), 1),
  median_carrier    = median(n_carrier),
  median_dist_bp    = median(median_dist, na.rm = TRUE)
), by = .(group, svtype)]
sumB_all <- cons[, .(
  svtype = "ALL", n_SV = .N,
  pct_near_any      = round(100 * mean(near_any), 1),
  pct_near_majority = round(100 * mean(near_majority), 1),
  pct_near_most     = round(100 * mean(near_most), 1),
  median_carrier    = median(n_carrier),
  median_dist_bp    = median(median_dist, na.rm = TRUE)
), by = group]
sumB <- rbind(sumB_all, sumB, fill = TRUE)[order(group, svtype)]
fwrite(sumB, file.path(outdir, "summary_by_sv_overall.tsv"), sep = "\t")

cat("\n===== A. Fractions by sample/haplotype (feature x dist x group) =====\n")
print(sumA[dist == MAIN_DIST])
cat("\n===== B. Per-SV summary, Any_structure, dist <=", MAIN_DIST, "bp =====\n")
print(sumB)
cat("\npct_near_majority = percentage of SVs near a structure in at least half of their available carrier-haplotype records.\n")

############################################################
## 5. Figures
############################################################

## Figure 1: median proximity fractions across thresholds, with random expectations.
pd <- sumA[, .(feature, dist, group, obs = median_frac, exp = median_exp)]
pl <- melt(pd, id.vars = c("feature", "dist", "group"),
           variable.name = "kind", value.name = "frac")
pl[, lab := fifelse(kind == "exp", paste0(group, " (permuted)"), group)]
p1 <- ggplot(pl, aes(dist / 1000, frac, colour = lab, linetype = kind)) +
  geom_line() + geom_point(size = 1.5) +
  geom_vline(xintercept = MAIN_DIST / 1000, linetype = "dotted", colour = "grey40") +
  facet_wrap(~feature) +
  scale_linetype_manual(values = c(obs = "solid", exp = "dashed"), guide = "none") +
  theme_classic(base_size = 12) +
  labs(x = "distance threshold (kb)", y = "median fraction of SVs within threshold",
       colour = "", title = "SV proximity to chromatin architectural features")
ggsave(file.path(outdir, "proximity_vs_distance.pdf"), p1, width = 10, height = 4)

## Figure 2: enrichment per sample/haplotype at MAIN_DIST.
p2 <- ggplot(res[dist == MAIN_DIST & is.finite(enrichment)],
             aes(feature, enrichment, fill = group)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.3) +
  scale_y_log10() +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = "", y = "enrichment over permuted positions", fill = "",
       title = paste0("Enrichment at ", MAIN_DIST / 1000, " kb"))
ggsave(file.path(outdir, "enrichment_boxplot.pdf"), p2, width = 7, height = 4)

## Figure 3: ECDF of finite nearest-structure distances across SV/haplotype records.
## Distances of 0 are plotted at 1 bp on the logarithmic axis.
p3 <- ggplot(svd[is.finite(dist)], aes(pmax(dist, 1), colour = group)) +
  stat_ecdf(linewidth = 0.8) +
  geom_vline(xintercept = MAIN_DIST, linetype = "dotted", colour = "grey40") +
  scale_x_log10() +
  theme_classic(base_size = 12) +
  labs(x = "distance to nearest architectural feature (bp)",
       y = "cumulative fraction", colour = "",
       title = "Distance to nearest TAD boundary or loop anchor")
ggsave(file.path(outdir, "distance_ecdf.pdf"), p3, width = 6.5, height = 4)

message("Done. Output directory: ", normalizePath(outdir))
