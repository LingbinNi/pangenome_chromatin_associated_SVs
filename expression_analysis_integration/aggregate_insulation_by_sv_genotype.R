#!/usr/bin/env Rscript
# Aggregate precomputed insulation scores by SV genotype.
# Requires data.table and ggplot2. Edit sv_id_file, genotype_file and tad_dirs.
# Run: Rscript /path/to/aggregate_insulation_by_sv_genotype.R
# out_base defaults to the current working directory; one subdirectory per SV.
#
# SV list: no header, one chromosome-start-type-length ID per line, with a
# 0-based start and a positive length in bp; interval end = start + length.
# Genotype TSV (.gz supported): sample, sv_id, hap1_carrier, hap2_carrier.
# Carrier columns contain 0, 1 or NA. Their sum defines dosage 0/1/2,
# corresponding to no_carrier (0/0), single_carrier (0/1), both_carrier (1/1).
# Exclude samples with missing dosage or no genotype record for the SV.
# Keep the first record for each sample/SV and the first score file per sample.
# Score input: recursively find *_fdr_score.bedgraph under each tad_dirs root;
# the immediate parent directory supplies the sample name for genotype matching.
# Supply headerless, four-column bedGraph files: chrom, start, end, score,
# with nonoverlapping 0-based, half-open intervals in matching CHM13 coordinates.
#
# Use a 25-kb grid from SV_start - 500000 through SV_start + length + 500000.
# Assign precomputed scores by interval containment; unmatched positions are NA.
# Express positions relative to the SV start. At each group/position, omit NA
# scores and calculate the mean, SD and SEM = SD/sqrt(n), where n is the number
# of nonmissing scores. Retain groups with N >= 2 matched, genotyped individuals;
# N is the total group size and n is the count at the individual grid position.
# Plot group means with SEM ribbons and mark the SV interval.
# Outputs per SV: <sv_id>_per_sample_grid.tsv,
# <sv_id>_insulation_profile_bygroup.tsv and <sv_id>_insulation_profile_bygroup.pdf.

suppressMessages({
  library(data.table)
  library(ggplot2)
})

## ------------------ Settings ------------------
sv_id_file <- "/path/to/sv_ids.txt"
genotype_file <- "/path/to/SV_CHM13_all_genotypes.tsv.gz"
out_base   <- "."
flank_bp   <- 500000
step_bp    <- 25000

tad_dirs <- c(
  "/path/to/cohort1/Tad",
  "/path/to/cohort2/Tad"
)

group_levels <- c("no_carrier", "single_carrier", "both_carrier")
group_colors <- c(no_carrier = "grey40", single_carrier = "darkorange", both_carrier = "firebrick")

## ------------------ 1. Parse SV IDs ------------------
parse_sv_id <- function(id_str) {
  parts <- strsplit(id_str, "-")[[1]]
  chrom  <- parts[1]
  start  <- as.numeric(parts[2])
  svtype <- parts[3]
  svlen  <- as.numeric(parts[4])
  end <- start + svlen
  list(chrom = chrom, start = start, end = end, svtype = svtype, svlen = svlen)
}

## ------------------ 2. Find sample bedGraph files ------------------
find_sample_files <- function(tad_dirs) {
  files <- unlist(lapply(tad_dirs, function(d) {
    list.files(d, pattern = "_fdr_score\\.bedgraph$",
                recursive = TRUE, full.names = TRUE)
  }))
  samples <- basename(dirname(files))
  dup <- duplicated(samples)
  if (any(dup)) {
    message("Multiple files found for these samples; keeping the first: ",
            paste(unique(samples[dup]), collapse = ", "))
  }
  setNames(files[!dup], samples[!dup])
}

sample_files <- find_sample_files(tad_dirs)
message("Found insulation score files for ", length(sample_files), " samples")
all_samples <- names(sample_files)

## ------------------ 3. Read all sample bedGraph files ------------------
message("Reading sample bedGraph files...")
sample_data <- lapply(sample_files, function(f) {
  dt <- fread(f, header = FALSE,
              col.names = c("chrom", "start", "end", "score"))
  setkey(dt, chrom, start, end)
  dt
})
names(sample_data) <- all_samples

## ------------------ 4. Read genotypes ------------------
message("Reading genotypes...")
geno <- fread(genotype_file, header = TRUE)

## Convert carrier columns to integers; missing values remain NA.
## A missing value in either carrier column gives a missing dosage.
geno[, hap1_carrier := suppressWarnings(as.integer(hap1_carrier))]
geno[, hap2_carrier := suppressWarnings(as.integer(hap2_carrier))]
geno[, dosage := hap1_carrier + hap2_carrier]   # Dosage 0 / 1 / 2; NA if either carrier value is missing

## ------------------ 5. Assign scores by interval containment ------------------
get_scores_at_grid <- function(bg_dt, chrom_, grid) {
  q <- data.table(chrom = chrom_, start = as.integer(grid), end = as.integer(grid) + 1L)
  setkey(q, chrom, start, end)
  res <- foverlaps(q, bg_dt, type = "within", mult = "first")
  res$score
}

## ------------------ 6. Process each SV ------------------
sv_ids <- readLines(sv_id_file)
sv_ids <- sv_ids[nzchar(sv_ids)]

for (cur_sv_id in sv_ids) {

  message("===> Processing ", cur_sv_id)
  sv <- parse_sv_id(cur_sv_id)

  sv_dir <- file.path(out_base, cur_sv_id)
  dir.create(sv_dir, showWarnings = FALSE, recursive = TRUE)

  ## --- Select genotype records for the current SV---
  geno_this <- geno[sv_id == cur_sv_id, .(sample, dosage)]
  geno_this <- unique(geno_this, by = "sample")   # Keep the first record for each sample

  if (nrow(geno_this) == 0) {
    message("  No genotype records for this SV; skipping.")
    next
  }

  ## --- Match genotypes to samples with insulation score files ---
  geno_sv <- merge(data.table(sample = all_samples), geno_this,
                   by = "sample", all.x = TRUE)

  ## --- Exclude samples with missing dosage or no genotype record---
  n_excluded <- geno_sv[is.na(dosage), .N]
  if (n_excluded > 0) {
    message("  Excluded ", n_excluded, " samples without genotype information for this SV: ",
            paste(geno_sv[is.na(dosage), sample], collapse = ", "))
  }
  geno_sv <- geno_sv[!is.na(dosage)]

  geno_sv[, group := factor(
    fifelse(dosage == 0, "no_carrier",
    fifelse(dosage == 1, "single_carrier", "both_carrier")),
    levels = group_levels
  )]

  group_counts <- geno_sv[, .N, by = group]
  message("Genotyped samples per group: ", paste(paste0(group_counts$group, "=", group_counts$N), collapse = ", "))

  region_start <- sv$start - flank_bp
  region_end   <- sv$end   + flank_bp
  grid <- seq(region_start, region_end, by = step_bp)

  ## --- Extract scores on the grid and attach genotype groups ---
  per_sample_list <- lapply(geno_sv$sample, function(s) {
    bg_dt <- sample_data[[s]]
    scores <- get_scores_at_grid(bg_dt, sv$chrom, grid)
    data.table(sample = s, rel_pos = grid - sv$start, score = scores)
  })
  dt_all <- rbindlist(per_sample_list)
  dt_all <- merge(dt_all, geno_sv[, .(sample, group)], by = "sample")

  fwrite(dt_all, file.path(sv_dir, paste0(cur_sv_id, "_per_sample_grid.tsv")), sep = "\t")

  ## --- Summarize nonmissing scores by group and grid position ---
  dt_summary <- dt_all[!is.na(score), .(
    mean_score = mean(score),
    sd_score   = sd(score),
    n          = .N
  ), by = .(group, rel_pos)]
  dt_summary[, sem := sd_score / sqrt(n)]
  dt_summary <- dt_summary[group_counts, on = "group"][N >= 2]
  setorder(dt_summary, group, rel_pos)

  fwrite(dt_summary, file.path(sv_dir, paste0(cur_sv_id, "_insulation_profile_bygroup.tsv")), sep = "\t")

  ## ------------------ Plot ------------------
  sv_start_kb <- 0
  sv_end_kb   <- (sv$end - sv$start) / 1000

  legend_labels <- setNames(
    paste0(group_levels, " (n=",
           sapply(group_levels, function(g) {
             v <- group_counts[group == g, N]
             if (length(v) == 0) 0 else v
           }), ")"),
    group_levels
  )

  p <- ggplot(dt_summary, aes(x = rel_pos / 1000, y = mean_score,
                               color = group, fill = group)) +
    annotate("rect", xmin = sv_start_kb, xmax = sv_end_kb,
             ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.08) +
    geom_ribbon(aes(ymin = mean_score - sem, ymax = mean_score + sem),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = c(sv_start_kb, sv_end_kb),
               linetype = "dashed", color = "grey30", linewidth = 0.4) +
    scale_color_manual(values = group_colors, labels = legend_labels, name = "Genotype") +
    scale_fill_manual(values = group_colors, labels = legend_labels, name = "Genotype") +
    labs(title = cur_sv_id,
         x = "Distance from SV start (kb)",
         y = "Mean insulation score \u00B1 SEM") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          legend.position = "bottom")

  ggsave(file.path(sv_dir, paste0(cur_sv_id, "_insulation_profile_bygroup.pdf")),
         p, width = 7, height = 4.5)
}

message("Done. Results are in the subdirectory for each SV.")