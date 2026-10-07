#!/usr/bin/env Rscript
# Identify recurrent CDH regions and test cohort-level recurrence by circular shifts.
# Requires R >= 3.6 and data.table: install.packages("data.table")
# Run: Rscript analyze_cdh_recurrence.R /path/to/CDH.CHM13.merged.tsv /path/to/output
# Alternatively, edit input_file and out_dir below.
# Input: projected CDH intervals in T2T-CHM13 v2.0 coordinates, with a tab-separated header:
#   chr, start, end, cdh_id, sample, haplotype, sample_haplotype
# Coordinates are 0-based, half-open; chromosome names must be chr1 through chr22.
# haplotype is Hap1 or Hap2; sample_haplotype is <sample>_Hap1 or <sample>_Hap2.
# Each donor must have both Hap1 and Hap2 records in the input.
# Prepare the projected intervals upstream; this script uses their supplied start/end bounds.
# Overlapping or adjacent intervals are merged within each chromosome and haplotype.
# Each haplotype contributes at most one unit of support at each genomic position.
# A recurrent core is a complete continuous interval supported by at least K=20 haplotypes.
# Retain cores with lengths from 12,500 (half of the CDH window) to 1,000,000 bp, inclusive, after forming complete cores.
# support_mean_bp is the base-pair-weighted mean haplotype support across a core.
# Chromosome lengths are embedded below and match the CHM13 reference used for projection.
# For each permutation, choose an independent circular offset for each donor/chromosome;
# that donor's Hap1 and Hap2 share the offset on that chromosome.
# Split intervals crossing the chromosome end into end/start segments, preserving covered length.
# Recurrent cores are evaluated in linear chromosome coordinates with the same rules in every run.
# Defaults: 999 permutations
# Report cohort-wide region counts and total recurrent bp, with upper-tail empirical P:
#   (1 + number of null values >= observed)/(B + 1).
# Outputs in out_dir: CDH_recurrent_core.tsv, CDH_recurrent_core.bed,
#   permutation_null.tsv, permutation_summary.tsv and run_parameters.txt.
# TSV files have headers; the BED file contains chr, start, end and region_id without a header.
# Use a new output directory. No plots are generated.

suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)

# ---------------- Input, output and analysis settings ----------------
input_file <- "/path/to/CDH.CHM13.merged.tsv"
out_dir <- "CDH_recurrence_results"
K <- 20L
min_length <- 12500 # half of the CDH window
max_length <- 1000000
B <- 999L
seed <- 2026100101L
# ------------------------------------------------

# Command-line arguments override the input path and output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) input_file <- args[1L]
if (length(args) >= 2L) out_dir <- args[2L]
stopifnot(K >= 1, min_length >= 1, max_length >= min_length, B >= 2)

# T2T-CHM13 v2.0 autosomal chromosome lengths for the projected input intervals.
chrom_sizes <- c(
  chr1=248387328, chr2=242696752, chr3=201105948, chr4=193574945,
  chr5=182045439, chr6=172126628, chr7=160567428, chr8=146259331,
  chr9=150617247, chr10=134758134, chr11=135127769, chr12=133324548,
  chr13=113566686, chr14=101161492, chr15=99753195, chr16=96330374,
  chr17=84276897, chr18=80542538, chr19=61707364, chr20=66210255,
  chr21=45090682, chr22=51324926
)

d <- fread(input_file)
required <- c("chr", "start", "end", "cdh_id", "sample", "haplotype", "sample_haplotype")
if (!all(required %in% names(d))) stop("Input must contain the seven required column names: ", paste(required, collapse="\t"))
d <- d[, ..required]
if (nrow(d) == 0L || anyNA(d)) stop("Input is empty or contains missing values.")
if (any(!d$chr %in% names(chrom_sizes))) stop("This analysis accepts only chr1 through chr22.")
if (!is.numeric(d$start) || !is.numeric(d$end)) stop("start/end must be numeric coordinates.")
d[, `:=`(start=as.numeric(start), end=as.numeric(end))]
if (any(!is.finite(d$start) | !is.finite(d$end) |
        d$start != floor(d$start) | d$end != floor(d$end) |
        d$start < 0 | d$end <= d$start | d$end > chrom_sizes[d$chr])) {
  stop("Coordinates must be BED integers: 0 <= start < end <= chromosome length.")
}
if (any(!d$haplotype %in% c("Hap1", "Hap2")) ||
    any(d$sample_haplotype != paste(d$sample, d$haplotype, sep="_"))) {
  stop("haplotype must be Hap1/Hap2; sample_haplotype must be sample_Hap1 or sample_Hap2.")
}
if (any(d[, uniqueN(haplotype), by=sample]$V1 != 2L)) {
  stop("Each donor must have both Hap1 and Hap2 records in the input.")
}
donors <- sort(unique(d$sample), method="radix")
chroms <- names(chrom_sizes)[names(chrom_sizes) %in% d$chr]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Merge overlapping or adjacent intervals within each haplotype and chromosome.
# This union limits each haplotype to one unit of support per position.
setorder(d, chr, sample_haplotype, start, end)
d[, block := cumsum(start > shift(cummax(end), fill=-1)),
  by=.(chr, sample, sample_haplotype)]
u <- d[, .(start=min(start), end=max(end)),
       by=.(chr, sample, sample_haplotype, block)]
u[, donor_id := match(sample, donors)]
by_chr <- setNames(lapply(chroms, function(ch) u[chr == ch]), chroms)

# Sweep interval endpoints, form complete continuous K-supported cores, then filter lengths.
# Keep BED coordinates; join adjacent qualifying segments and separate them at lower-support gaps.
get_cores <- function(starts, ends) {
  e <- data.table(pos=c(starts, ends),
                  delta=c(rep.int(1L, length(starts)), rep.int(-1L, length(ends))))
  e <- e[, .(delta=sum(delta)), by=pos]
  setorder(e, pos)
  e[, depth := cumsum(delta)]
  stopifnot(tail(e$depth, 1L) == 0L, all(e$depth >= 0L))
  e[, end := shift(pos, type="lead")]
  e <- e[!is.na(end)]
  if (!any(e$depth >= K)) {
    return(data.table(start=numeric(), end=numeric(), length_bp=numeric(),
                      support_min=integer(), support_max=integer(), support_mean_bp=numeric()))
  }
  e[, run := rleid(depth >= K)]
  cores <- e[depth >= K,
             .(start=min(pos), end=max(end), length_bp=sum(end-pos),
               support_min=min(depth), support_max=max(depth),
               support_mean_bp=sum(depth*(end-pos))/sum(end-pos)), by=run]
  cores <- cores[length_bp >= min_length & length_bp <= max_length]
  cores[, run := NULL]
  cores[]
}

# Generate the observed recurrent CDH catalog.
observed <- rbindlist(lapply(chroms, function(ch) {
  x <- by_chr[[ch]]
  ans <- get_cores(x$start, x$end)
  ans[, chr := rep(ch, .N)]
  ans
}), use.names=TRUE)
observed[, region_id := sprintf("CDH_CORE_%06d", seq_len(.N))]
setcolorder(observed, c("chr", "start", "end", "region_id", "length_bp",
                        "support_min", "support_max", "support_mean_bp"))
fwrite(observed, file.path(out_dir, "CDH_recurrent_core.tsv"), sep="\t")
fwrite(observed[, .(chr, start, end, region_id)],
       file.path(out_dir, "CDH_recurrent_core.bed"), sep="\t", col.names=FALSE)
obs_n <- nrow(observed)
obs_bp <- sum(observed$length_bp)
cat("Input: ", nrow(d), " intervals; ", length(donors), " donors; ",
    uniqueN(d$sample_haplotype), " haplotypes\n", sep="")
cat("Observed:", obs_n, "regions;", obs_bp, "bp\n")

# Draw independent donor/chromosome circular offsets; a donor shares the offset across both haplotypes.
# Uniform shifts on each complete chromosome preserve haplotype coverage and circular interval arrangement.
# Split wrapped intervals at the chromosome end and retain both segments.
# Evaluate final cores in linear coordinates, keeping chromosome-end and chromosome-start cores separate.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
null <- data.table(replicate=seq_len(B), n_regions=integer(B), total_bp=numeric(B))
for (b in seq_len(B)) {
  n <- 0L
  bp <- 0
  for (ch in chroms) {
    x <- by_chr[[ch]]
    L <- chrom_sizes[[ch]]
    offsets <- sample.int(L, length(donors), replace=TRUE) - 1L
    starts <- (x$start + offsets[x$donor_id]) %% L
    ends <- starts + (x$end - x$start)
    wrap <- ends > L
    shifted_starts <- c(starts, rep.int(0, sum(wrap)))
    shifted_ends <- c(pmin(ends, L), ends[wrap] - L)
    cores <- get_cores(shifted_starts, shifted_ends)
    n <- n + nrow(cores)
    bp <- bp + sum(cores$length_bp)
  }
  null[b, `:=`(n_regions=n, total_bp=bp)]
  if (b %% 100L == 0L || b == B) cat("Permutations:", b, "/", B, "\n")
}
fwrite(null, file.path(out_dir, "permutation_null.tsv"), sep="\t")

# Upper-tail P = (1 + number of null values >= observed)/(B + 1); SD is the sample standard deviation.
# Test cohort-wide region counts and total recurrent length.
summarize_null <- function(metric, observed_value) {
  v <- null[[metric]]
  data.table(metric=metric, observed=observed_value, B=B,
             null_mean=mean(v), null_sd=sd(v), null_min=min(v), null_max=max(v),
             fold_enrichment=if (mean(v) > 0) observed_value/mean(v) else NA_real_,
             n_null_ge_observed=sum(v >= observed_value),
             empirical_P=(1+sum(v >= observed_value))/(B+1))
}
result <- rbind(summarize_null("n_regions", obs_n), summarize_null("total_bp", obs_bp))
fwrite(result, file.path(out_dir, "permutation_summary.tsv"), sep="\t")
writeLines(c(paste("input:", normalizePath(input_file)),
             paste("K:", K), paste("min_length:", min_length), paste("max_length:", max_length),
             paste("B:", B), paste("seed:", seed), paste("RNG:", paste(RNGkind(), collapse=", ")),
             "Input length filter: none; length thresholds apply to final complete cores.",
             "Null: chromosome-wise, donor-paired circular shifts with uniform offsets on complete chromosomes.",
             capture.output(sessionInfo())), file.path(out_dir, "run_parameters.txt"))
print(result)
cat("Done. Output directory: ", normalizePath(out_dir), "\n", sep="")
