#!/usr/bin/env Rscript
# written by Lingbin

# Test CDH overlap with TAD boundaries, loop anchors and their combined intervals.
# Requires data.table, GenomicRanges and IRanges; ggplot2 produces optional plots.
# Edit cdh_root, tad_roots, loop_roots and outdir before running.
# Run: Rscript /path/to/test_cdh_structure_enrichment.R
# CDH input: <cdh_root>/<sample>/CDH.region.<sample>.<Hap1|Hap2>.tsv
# CDH files contain a header with Chr, Start and End, as written by the CDH caller.
# TAD input: <tad_roots[[hap]]>/<sample>/*_boundaries.bed
#   Headerless BED with at least five columns: chr, start, end, name, score.
# Loop input: <loop_roots[[hap]]>/<sample>/*.bedgraph
#   Headerless paired intervals: chr1, start1, end1, chr2, start2, end2, ...
# Provide one matching TAD file and one matching loop file per sample/haplotype.
# Use matching sample/haplotype assembly coordinates and chromosome names in all inputs.
# All input intervals use 0-based starts and half-open ends.
# Sample directories currently match HG*; edit the sample filter for other names.
# The default analysis uses merged CDH regions; cdh_prefix can select CDH.window.
# Features are merged before testing; each CDH contributes at most one overlap hit.
# n_feature counts merged feature intervals; obs_overlap_events counts overlapping
# CDH-feature pairs, while obs_CDH_overlap counts CDHs with at least one overlap.
# Each permutation independently relocates CDHs on their original chromosomes,
# retaining their lengths. Per-chromosome bounds are the maximum observed ends
# across the CDH, TAD and loop-anchor inputs for that sample/haplotype.
# Use 1,000 permutations and seed 1; enrichment = (observed + 0.5)/(random_mean + 0.5).
# Upper-tail empirical P = (1 + number of random counts >= observed)/(nperm + 1).
# Empirical P values are reported directly; nearest-feature distances use interval gaps.
# Outputs: CDH_3D_structure_overlap_per_sample.tsv and CDH_3D_structure_overlap_summary.tsv.
# If ggplot2 is installed, two CDH_3D_structure_*_boxplot.pdf files are also written.
# Run in a new empty output directory.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

############################################################
## Parameters
############################################################

nperm <- 1000
set.seed(1)

outdir <- "./"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## CDH directory layout (output from call_contact_difference_hotspots.R):
##   cdh_root/HG002/CDH.region.HG002.Hap1.tsv
##   cdh_root/HG002/CDH.region.HG002.Hap2.tsv
## (CDH.window.*.tsv also present in the same directory; switch
##  cdh_prefix below if you want window-level instead of region-level.)
cdh_root <- "/path/to/contact_difference_hotspots"

## "CDH.region" (merged hotspots) or "CDH.window" (25-kb windows).
cdh_prefix <- "CDH.region"

## The input files are already filtered CDH calls. Keep NULL to use all rows.
## A numeric cutoff filters only an existing qval column (window-level input).
## The region-level min_qval column is left unchanged by this option.
qval_cutoff <- NULL

tad_roots <- list(
  Hap1 = "/path/to/hap1_tad_boundaries",
  Hap2 = "/path/to/hap2_tad_boundaries"
)

loop_roots <- list(
  Hap1 = "/path/to/hap1_loops",
  Hap2 = "/path/to/hap2_loops"
)

############################################################
## Reading functions
############################################################

read_cdh <- function(file, qval_cutoff = NULL) {
  dt <- fread(file)

  ## Region format:
  ##   Chr,Start,End,n_windows,sum_spec,sum_total,region_frac,max_OE,min_qval,peak_Start
  ## Window format:
  ##   Chr,Start,End,spec,shared,total,expected,theta_local,OE,pval,qval
  ## Other input tables use their first three columns as coordinates.
  if (all(c("Chr", "Start", "End") %in% names(dt))) {
    setnames(dt, old = c("Chr", "Start", "End"), new = c("chr", "start", "end"))
  } else {
    setnames(dt, old = names(dt)[1:3], new = c("chr", "start", "end"))
    if (ncol(dt) >= 4 && !("value" %in% names(dt))) {
      setnames(dt, old = names(dt)[4], new = "value")
    }
  }

  dt[, start := suppressWarnings(as.integer(start))]
  dt[, end := suppressWarnings(as.integer(end))]

  if (!is.null(qval_cutoff) && "qval" %in% names(dt)) {
    dt[, qval := suppressWarnings(as.numeric(qval))]
    dt <- dt[!is.na(qval) & qval <= qval_cutoff]
  }

  dt <- dt[!is.na(chr) & !is.na(start) & !is.na(end)]
  dt <- dt[end > start]

  if (nrow(dt) == 0) {
    return(GRanges())
  }

  ## Treat coordinates as BED-like: 0-based start, half-open end.
  gr <- GRanges(
    seqnames = dt$chr,
    ranges = IRanges(start = dt$start + 1, end = dt$end)
  )

  ## Store useful CDH metadata when present.
  if ("spec" %in% names(dt)) {
    mcols(gr)$value <- dt$spec
    mcols(gr)$spec <- dt$spec
  } else if ("value" %in% names(dt)) {
    mcols(gr)$value <- dt$value
  }
  if ("total" %in% names(dt)) mcols(gr)$total <- dt$total
  if ("frac" %in% names(dt)) mcols(gr)$frac <- dt$frac
  if ("qval" %in% names(dt)) mcols(gr)$qval <- dt$qval

  gr
}

read_tad_boundary <- function(file) {
  dt <- fread(file, header = FALSE)

  dt <- dt[!is.na(V1) & !is.na(V2) & !is.na(V3)]
  dt <- dt[V3 > V2]

  GRanges(
    seqnames = dt$V1,
    ranges = IRanges(start = dt$V2 + 1, end = dt$V3),
    name = dt$V4,
    score = suppressWarnings(as.numeric(dt$V5))
  )
}

read_loop_anchors <- function(file) {
  dt <- fread(file, header = FALSE)

  ## Expected format:
  ## chr1 start1 end1 chr2 start2 end2 [additional columns are not used]
  dt <- dt[!is.na(V1) & !is.na(V2) & !is.na(V3) &
             !is.na(V4) & !is.na(V5) & !is.na(V6)]
  dt <- dt[V3 > V2 & V6 > V5]

  anchor1 <- data.table(
    chr = dt$V1,
    start = dt$V2,
    end = dt$V3
  )

  anchor2 <- data.table(
    chr = dt$V4,
    start = dt$V5,
    end = dt$V6
  )

  anchors <- unique(rbind(anchor1, anchor2))
  anchors <- anchors[end > start]

  gr <- GRanges(
    seqnames = anchors$chr,
    ranges = IRanges(start = anchors$start + 1, end = anchors$end)
  )

  reduce(gr, ignore.strand = TRUE)
}

############################################################
## Helper functions
############################################################

find_one_file <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    return(NA_character_)
  }

  if (length(files) > 1) {
    warning("Multiple files found in ", dir, "; using first: ", basename(files[1]))
  }

  files[1]
}

# Set each chromosome bound to the maximum end across the supplied interval sets.
make_chr_lengths <- function(gr_list) {
  dt <- rbindlist(lapply(gr_list, function(gr) {
    if (length(gr) == 0) {
      return(data.table(chr = character(), end = integer()))
    }

    data.table(
      chr = as.character(seqnames(gr)),
      end = end(gr)
    )
  }))

  dt <- dt[!is.na(chr)]
  chr_dt <- dt[, .(chr_len = max(end, na.rm = TRUE)), by = chr]

  x <- chr_dt$chr_len
  names(x) <- chr_dt$chr

  x
}

# Independently sample a valid start for each CDH, preserving chromosome and width.
randomize_cdh_by_chr <- function(cdh, chr_len) {
  dt <- data.table(
    chr = as.character(seqnames(cdh)),
    width = width(cdh)
  )

  dt <- dt[chr %in% names(chr_len)]
  dt[, chr_len := chr_len[chr]]
  dt <- dt[width <= chr_len]

  if (nrow(dt) == 0) {
    return(GRanges())
  }

  dt[, max_start := as.integer(chr_len - width + 1)]
  dt[, rand_start := vapply(max_start, function(m) sample.int(m, 1), integer(1))]
  dt[, rand_end := rand_start + width - 1]

  GRanges(
    seqnames = dt$chr,
    ranges = IRanges(start = dt$rand_start, end = dt$rand_end)
  )
}

distance_summary <- function(cdh, feature) {
  if (length(cdh) == 0 || length(feature) == 0) {
    return(list(
      median_distance = NA_real_,
      mean_distance = NA_real_,
      p10_distance = NA_real_,
      p90_distance = NA_real_
    ))
  }

  hits <- distanceToNearest(cdh, feature, ignore.strand = TRUE)

  dist <- rep(NA_real_, length(cdh))
  dist[queryHits(hits)] <- mcols(hits)$distance

  list(
    median_distance = median(dist, na.rm = TRUE),
    mean_distance = mean(dist, na.rm = TRUE),
    p10_distance = as.numeric(quantile(dist, 0.10, na.rm = TRUE)),
    p90_distance = as.numeric(quantile(dist, 0.90, na.rm = TRUE))
  )
}

test_feature <- function(cdh, feature, feature_name, sample, hap, chr_len, nperm = 1000) {
  feature <- feature[as.character(seqnames(feature)) %in% as.character(unique(seqnames(cdh)))]
  feature <- reduce(feature, ignore.strand = TRUE)

  n_cdh <- length(cdh)
  n_feature <- length(feature)

  if (n_cdh == 0 || n_feature == 0) {
    return(data.table(
      sample = sample,
      hap = hap,
      feature = feature_name,
      n_CDH = n_cdh,
      n_feature = n_feature,
      obs_CDH_overlap = NA_integer_,
      obs_overlap_events = NA_integer_,
      obs_prop = NA_real_,
      random_mean = NA_real_,
      random_sd = NA_real_,
      random_prop_mean = NA_real_,
      enrichment = NA_real_,
      empirical_p = NA_real_,
      z_score = NA_real_,
      median_distance = NA_real_,
      mean_distance = NA_real_,
      p10_distance = NA_real_,
      p90_distance = NA_real_
    ))
  }

  obs_overlap_count <- sum(countOverlaps(cdh, feature, ignore.strand = TRUE) > 0)
  obs_overlap_events <- length(findOverlaps(cdh, feature, ignore.strand = TRUE))
  obs_prop <- obs_overlap_count / n_cdh

  rand_counts <- replicate(nperm, {
    rand_cdh <- randomize_cdh_by_chr(cdh, chr_len)
    sum(countOverlaps(rand_cdh, feature, ignore.strand = TRUE) > 0)
  })

  rand_mean <- mean(rand_counts)
  rand_sd <- sd(rand_counts)

  empirical_p <- (sum(rand_counts >= obs_overlap_count) + 1) / (nperm + 1)

  enrichment <- (obs_overlap_count + 0.5) / (rand_mean + 0.5)

  z_score <- ifelse(
    rand_sd > 0,
    (obs_overlap_count - rand_mean) / rand_sd,
    NA_real_
  )

  dist_sum <- distance_summary(cdh, feature)

  data.table(
    sample = sample,
    hap = hap,
    feature = feature_name,
    n_CDH = n_cdh,
    n_feature = n_feature,
    obs_CDH_overlap = obs_overlap_count,
    obs_overlap_events = obs_overlap_events,
    obs_prop = obs_prop,
    random_mean = rand_mean,
    random_sd = rand_sd,
    random_prop_mean = rand_mean / n_cdh,
    enrichment = enrichment,
    empirical_p = empirical_p,
    z_score = z_score,
    median_distance = dist_sum$median_distance,
    mean_distance = dist_sum$mean_distance,
    p10_distance = dist_sum$p10_distance,
    p90_distance = dist_sum$p90_distance
  )
}

############################################################
## Main function per sample/haplotype
############################################################

run_one_sample_hap <- function(sample, hap) {
  message("Processing ", sample, " ", hap)

  cdh_file <- file.path(
    cdh_root,
    sample,
    paste0(cdh_prefix, ".", sample, ".", hap, ".tsv")
  )

  tad_dir <- file.path(tad_roots[[hap]], sample)
  loop_dir <- file.path(loop_roots[[hap]], sample)

  tad_file <- find_one_file(tad_dir, "_boundaries\\.bed$")
  loop_file <- find_one_file(loop_dir, "\\.bedgraph$")

  if (!file.exists(cdh_file)) {
    warning("Missing CDH file: ", cdh_file)
    return(NULL)
  }

  if (is.na(tad_file) || !file.exists(tad_file)) {
    warning("Missing TAD boundary file for ", sample, " ", hap)
    return(NULL)
  }

  if (is.na(loop_file) || !file.exists(loop_file)) {
    warning("Missing loop file for ", sample, " ", hap)
    return(NULL)
  }

  cdh <- read_cdh(cdh_file, qval_cutoff = qval_cutoff)
  tad <- read_tad_boundary(tad_file)
  loop_anchor <- read_loop_anchors(loop_file)

  any_structure <- reduce(c(tad, loop_anchor), ignore.strand = TRUE)

  chr_len <- make_chr_lengths(list(cdh, tad, loop_anchor))

  rbindlist(list(
    test_feature(
      cdh = cdh,
      feature = tad,
      feature_name = "TAD_boundary",
      sample = sample,
      hap = hap,
      chr_len = chr_len,
      nperm = nperm
    ),
    test_feature(
      cdh = cdh,
      feature = loop_anchor,
      feature_name = "Loop_anchor",
      sample = sample,
      hap = hap,
      chr_len = chr_len,
      nperm = nperm
    ),
    test_feature(
      cdh = cdh,
      feature = any_structure,
      feature_name = "Any_3D_structure",
      sample = sample,
      hap = hap,
      chr_len = chr_len,
      nperm = nperm
    )
  ), fill = TRUE)
}

############################################################
## Traverse all HG* samples
############################################################

samples <- basename(list.dirs(cdh_root, full.names = TRUE, recursive = FALSE))
samples <- sort(samples[grepl("^HG", samples)])

## Keep samples with at least one Hap1/Hap2 CDH file.
samples <- samples[
  file.exists(file.path(cdh_root, samples, paste0(cdh_prefix, ".", samples, ".Hap1.tsv"))) |
    file.exists(file.path(cdh_root, samples, paste0(cdh_prefix, ".", samples, ".Hap2.tsv")))
]

message("Total samples found: ", length(samples))

all_results <- rbindlist(lapply(samples, function(s) {
  rbindlist(lapply(c("Hap1", "Hap2"), function(h) {
    tryCatch(
      run_one_sample_hap(s, h),
      error = function(e) {
        warning("Failed: ", s, " ", h, " -- ", conditionMessage(e))
        return(NULL)
      }
    )
  }), fill = TRUE)
}), fill = TRUE)

############################################################
## Write per-sample result
############################################################

per_sample_out <- file.path(outdir, "CDH_3D_structure_overlap_per_sample.tsv")
fwrite(all_results, per_sample_out, sep = "\t")

############################################################
## Summary by haplotype and feature
############################################################

summary_results <- all_results[
  !is.na(obs_prop),
  .(
    n_sample_hap = .N,
    mean_n_CDH = mean(n_CDH, na.rm = TRUE),
    median_n_CDH = median(n_CDH, na.rm = TRUE),

    mean_obs_prop = mean(obs_prop, na.rm = TRUE),
    median_obs_prop = median(obs_prop, na.rm = TRUE),

    mean_random_prop = mean(random_prop_mean, na.rm = TRUE),
    median_random_prop = median(random_prop_mean, na.rm = TRUE),

    mean_enrichment = mean(enrichment, na.rm = TRUE),
    median_enrichment = median(enrichment, na.rm = TRUE),

    mean_z_score = mean(z_score, na.rm = TRUE),
    median_z_score = median(z_score, na.rm = TRUE),

    n_empirical_p_lt_0.05 = sum(empirical_p < 0.05, na.rm = TRUE),
    frac_empirical_p_lt_0.05 = mean(empirical_p < 0.05, na.rm = TRUE),

    ## empirical_p is upper-tail: random >= observed.
    ## Count results whose upper-tail empirical P exceeds 0.95.
    n_empirical_p_gt_0.95 = sum(empirical_p > 0.95, na.rm = TRUE),
    frac_empirical_p_gt_0.95 = mean(empirical_p > 0.95, na.rm = TRUE),

    median_distance = median(median_distance, na.rm = TRUE),
    mean_distance = mean(mean_distance, na.rm = TRUE)
  ),
  by = .(hap, feature)
]

summary_out <- file.path(outdir, "CDH_3D_structure_overlap_summary.tsv")
fwrite(summary_results, summary_out, sep = "\t")

############################################################
## Optional simple plots
############################################################

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  p1 <- ggplot(all_results[!is.na(enrichment)],
               aes(x = feature, y = enrichment)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~hap) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "",
      y = "CDH enrichment over chromosome-matched random intervals",
      title = "CDH overlap enrichment at 3D genome structures"
    )

  ggsave(
    filename = file.path(outdir, "CDH_3D_structure_enrichment_boxplot.pdf"),
    plot = p1,
    width = 7,
    height = 4
  )

  p2 <- ggplot(all_results[!is.na(obs_prop)],
               aes(x = feature, y = obs_prop)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~hap) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "",
      y = "Fraction of CDHs overlapping structure",
      title = "Observed CDH overlap with 3D genome structures"
    )

  ggsave(
    filename = file.path(outdir, "CDH_3D_structure_overlap_fraction_boxplot.pdf"),
    plot = p2,
    width = 7,
    height = 4
  )
}

message("Done.")
message("Per-sample result: ", per_sample_out)
message("Summary result: ", summary_out)
