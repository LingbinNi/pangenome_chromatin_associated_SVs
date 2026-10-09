#!/usr/bin/env Rscript
# Plot local insulation versus expression on TPM and log2(TPM+1) scales.
# Requires R >= 4.0, data.table and ggplot2 >= 3.4.0.
# Rebuild mode also requires bash, gzip and awk.
# Edit EXPR_FILE, GENO_FILE and TAD_DIRS before running.
# Run: Rscript /path/to/plot_insulation_expression_correlations.R
# OUTDIR defaults to the current directory. Use an empty output directory.
# Expression TSV (.gz supported): sv_id, gene_id, sample, TPM; optional gene_name.
# The aliases var and sample_orig are accepted for sv_id and sample.
# Supply gene-level, nonnegative TPM. Zero TPM is retained.
# Genotype TSV (.gz supported): sv_id, sample, hap1_carrier, hap2_carrier.
# Carrier flags are 0, 1 or NA; their sum defines dosage 0/1/2.
# Genotypes come from GENO_FILE; an expression-table genotype column, if present,
# is used only for a consistency check. Missing genotypes are excluded.
# Outputs: Insulation_Expression.TPM.pdf/png and
# Insulation_Expression.log2TPM.pdf/png, plus individual gene panels.
# Tables: correlation_summary.tsv, sample_counts.tsv, sample_QC.tsv,
# plot_data.tsv, resolved_targets.tsv and per-gene insulation_expression.tsv.
# Rebuild mode also writes insulation_files_used.tsv.
# analysis_notes.txt and sessionInfo.txt record the analysis and R environment.

suppressPackageStartupMessages({
  for (pkg in c("data.table", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Please install package: ", pkg)
  }
  library(data.table)
  library(ggplot2)
})
if (utils::packageVersion("ggplot2") < "3.4.0") stop("ggplot2 >= 3.4.0 is required")

############################################################
## Settings
############################################################

EXPR_FILE <- "/path/to/AdjacentGenes.f25.tsv.gz"
GENO_FILE <- "/path/to/SV_CHM13_all_genotypes.tsv.gz"
OUTDIR <- "./"

TAD_DIRS <- c(
  "/path/to/cohort1/insulation_scores",
  "/path/to/cohort2/insulation_scores"
)

MODE <- "rebuild"                 # "rebuild" or "replot"
INTEGRATION_DIR <- "/path/to/insulation_expression_results"  # Replot input only
INSULATION_FILES <- NULL          # Optional TSV with sample and file columns
TARGET_FILE <- NULL               # Optional TSV with gene_name, sv_id, gene_id
OFFSET <- 0L                     # 0: use the supplied SV start; -1: convert 1-based POS

if (!MODE %in% c("rebuild", "replot")) stop("MODE must be rebuild or replot")
if (length(OFFSET) != 1L || is.na(OFFSET) || !OFFSET %in% c(0L, -1L)) {
  stop("OFFSET must be 0 or -1")
}
if (MODE == "replot" && OFFSET != 0L) stop("Coordinate offsets require rebuild mode")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Genotype colors used for the expression examples.
GT_LEVELS <- c("AA", "Aa", "aa")
GT_LABELS <- c(AA = "0/0", Aa = "0/1", aa = "1/1")
COLORS <- c(AA = "#706AB0", Aa = "#F37268", aa = "#56B4E9")
FIG_WIDTH <- 6.2             # Combined width for two gene panels, inches
FIG_HEIGHT <- 3.35
FONT_FAMILY <- "sans"
POINT_SIZE <- 1.35
POINT_ALPHA <- 0.72

############################################################
## Helpers
############################################################

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
out <- function(x, name) fwrite(x, file.path(OUTDIR, name), sep = "\t", na = "NA")
need <- function(d, cols, label) {
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(label, " lacks: ", paste(miss, collapse = ", "))
}
num <- function(x, label) {
  z <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & is.na(z))) stop("Non-numeric values in ", label)
  z
}
gt_dosage <- function(x) {
  fcase(x %in% c("0/0", "0|0", "AA"), 0L,
        x %in% c("0/1", "1/0", "0|1", "1|0", "Aa"), 1L,
        x %in% c("1/1", "1|1", "aa"), 2L, default = NA_integer_)
}
dedup_checked <- function(d, keys, values, label) {
  counts <- d[, lapply(.SD, uniqueN), by = keys, .SDcols = values]
  bad <- counts[rowSums(as.matrix(counts[, ..values]) > 1L) > 0L]
  if (nrow(bad)) {
    out(bad, paste0(safe_name(label), ".conflicting_duplicates.tsv"))
    stop("Conflicting duplicate observations in ", label)
  }
  unique(d, by = keys)
}

############################################################
## Targets
############################################################

targets <- data.table(gene_name = c("SESTD1", "SIGLEC5"),
  sv_id = c("chr2-179683834-DEL-16208", "chr19-54718506-DEL-16422"),
  gene_id = NA_character_)
target_file <- TARGET_FILE
if (!is.null(target_file)) {
  targets <- fread(target_file, colClasses = "character", na.strings = c("", "NA"))
  need(targets, c("gene_name", "sv_id"), "Targets")
  if (!"gene_id" %in% names(targets)) targets[, gene_id := NA_character_]
}
if (!nrow(targets) || anyNA(targets$sv_id) || anyNA(targets$gene_name)) stop("Invalid targets")
targets[, panel_index := .I]
loc <- rbindlist(lapply(unique(targets$sv_id), function(id) {
  a <- strsplit(id, "-", fixed = TRUE)[[1L]]
  if (length(a) != 4L || !a[3L] %in% c("DEL", "INS")) stop("Cannot parse SV ID: ", id)
  pos <- suppressWarnings(as.numeric(a[2L])) + OFFSET
  if (!is.finite(pos) || pos < 0 || pos != floor(pos)) stop("Invalid SV position: ", id)
  data.table(sv_id = id, chrom = a[1L], position = pos)
}))

############################################################
## Input reader
############################################################

# Read rows for the selected SVs; allow quoted headers and SV IDs.
read_sv_subset <- function(path, ids) {
  if (!file.exists(path)) stop("Input not found: ", path)
  for (cmd in c("bash", "gzip", "awk")) if (!nzchar(Sys.which(cmd))) stop("Missing program: ", cmd)
  key_file <- tempfile(fileext = ".txt")
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(key_file, tmp)), add = TRUE)
  if (any(grepl("[\r\n\t]", ids))) stop("Invalid SV identifier")
  writeLines(unique(ids), key_file)
  awk_code <- paste(
    'BEGIN { FS=OFS="\t" }',
    'function key(s) { sub(/\\r$/, "", s); sub(/^"/, "", s); sub(/"$/, "", s); return s }',
    'NR==FNR { ids[$1]=1; next }',
    'FNR==1 { for(i=1;i<=NF;i++) { if(key($i)=="sv_id") c=i; if(key($i)=="var") alt=i };',
    'if(!c) c=alt; if(!c) { print "Missing sv_id/var header" > "/dev/stderr"; exit 2 }; print; next }',
    '(key($c) in ids) { print }')
  source_cmd <- if (grepl("\\.gz$", path)) "gzip -cd" else "cat"
  pipeline <- paste(source_cmd, shQuote(path), "| awk", shQuote(awk_code), shQuote(key_file), "-")
  status <- system2(Sys.which("bash"), c("-o", "pipefail", "-c", shQuote(pipeline)), stdout = tmp)
  if (status != 0L) stop("Failed to extract input: ", path)
  d <- fread(tmp, na.strings = c("", "NA", "NaN", "."))
  if (!"sv_id" %in% names(d) && "var" %in% names(d)) setnames(d, "var", "sv_id")
  if (!"sample" %in% names(d) && "sample_orig" %in% names(d)) setnames(d, "sample_orig", "sample")
  d
}

############################################################
## Pair expression, genotypes and local insulation
############################################################

pair_list <- vector("list", nrow(targets))
if (MODE == "rebuild") {
  message("Reading expression and complete genotypes...")
  expr <- read_sv_subset(EXPR_FILE, targets$sv_id)
  geno <- read_sv_subset(GENO_FILE, targets$sv_id)
  need(expr, c("sv_id", "gene_id", "sample", "TPM"), "Expression table")
  need(geno, c("sv_id", "sample", "hap1_carrier", "hap2_carrier"), "Genotype table")
  if (!"gene_name" %in% names(expr)) expr[, gene_name := NA_character_]
  expr[, TPM := num(TPM, "TPM")]
  if (any(!is.na(expr$TPM) & (!is.finite(expr$TPM) | expr$TPM < 0))) stop("Invalid TPM")
  geno[, hap1_carrier := num(hap1_carrier, "hap1_carrier")]
  geno[, hap2_carrier := num(hap2_carrier, "hap2_carrier")]
  for (v in c("hap1_carrier", "hap2_carrier")) {
    if (any(!is.na(geno[[v]]) & !geno[[v]] %in% c(0, 1))) stop("Carrier flags must be 0, 1 or NA")
  }
  geno[, dosage := hap1_carrier + hap2_carrier]
  geno <- dedup_checked(geno[, .(sv_id, sample, dosage)], c("sv_id", "sample"), "dosage", "genotypes")

  for (i in seq_len(nrow(targets))) {
    t <- targets[i]
    if (!is.na(t$gene_id) && nzchar(t$gene_id)) {
      e <- expr[sv_id == t$sv_id & gene_id == t$gene_id]
    } else {
      e <- expr[sv_id == t$sv_id & (gene_id %in% t$gene_name | gene_name %in% t$gene_name)]
    }
    ids <- unique(e$gene_id)
    if (length(ids) != 1L || anyNA(ids)) {
      out(unique(expr[sv_id == t$sv_id, .(sv_id, gene_id, gene_name)]), paste0(t$gene_name, ".candidates.tsv"))
      stop("Need one exact SV-gene pair for ", t$gene_name, "; inspect candidates and set gene_id in TARGET_FILE if needed")
    }
    targets[i, gene_id := ids]
    # Use the complete genotype table; check expression-table genotypes if present.
    e[, expression_dosage := if ("genotype" %in% names(e)) gt_dosage(e$genotype) else NA_integer_]
    e <- e[, .(sv_id, gene_id, sample, TPM, expression_dosage)]
    e <- dedup_checked(e, c("sv_id", "gene_id", "sample"), c("TPM", "expression_dosage"), paste0(t$gene_name, "_expression"))
    e <- merge(e, geno, by = c("sv_id", "sample"), all.x = TRUE, sort = FALSE)
    if (any(!is.na(e$expression_dosage) & !is.na(e$dosage) & e$expression_dosage != e$dosage)) {
      stop("Expression/genotype tables disagree for ", t$gene_name)
    }
    e[, `:=`(plot_gene = t$gene_name, panel_index = i)]
    pair_list[[i]] <- e
  }
  pairs <- rbindlist(pair_list, fill = TRUE)
  if (anyNA(pairs$sample) || any(!nzchar(pairs$sample))) stop("Missing sample IDs")

  # Require one selected bedGraph per sample.
  manifest <- INSULATION_FILES
  if (!is.null(manifest)) {
    files <- fread(manifest, colClasses = "character")
    need(files, c("sample", "file"), "Insulation file manifest")
    files <- files[, .(sample, file)]
  } else {
    paths <- unlist(lapply(TAD_DIRS, function(d) {
      if (!dir.exists(d)) warning("TAD directory not found: ", d)
      list.files(d, pattern = "_fdr_score\\.bedgraph$", recursive = TRUE, full.names = TRUE)
    }), use.names = FALSE)
    files <- data.table(sample = basename(dirname(paths)), file = paths)
  }
  files <- unique(files[sample %in% pairs$sample])
  if (!nrow(files)) stop("No matching CHM13 insulation files; check TAD_DIRS or INSULATION_FILES")
  dup <- files[, .N, by = sample][N > 1L]
  if (nrow(dup)) {
    out(files[sample %in% dup$sample], "ambiguous_insulation_files.tsv")
    stop("Multiple bedGraphs for a sample. Select one per sample using INSULATION_FILES")
  }
  if (anyNA(files$file) || any(!file.exists(files$file))) stop("A selected bedGraph does not exist")
  out(files, "insulation_files_used.tsv")
  insulation_list <- vector("list", nrow(files))
  for (j in seq_len(nrow(files))) {
    if (j == 1L || j %% 20L == 0L || j == nrow(files)) message("Reading insulation: ", j, "/", nrow(files))
    f <- files[j]
    bg <- fread(f$file, header = FALSE, select = 1:4,
                col.names = c("chrom", "start", "end", "score"),
                na.strings = c("", "NA", "NaN", "nan", "."))
    bg <- bg[chrom %in% loc$chrom]
    bg[, `:=`(start = num(start, "bedGraph start"), end = num(end, "bedGraph end"),
               score = num(score, "bedGraph score"))]
    if (anyNA(bg$start) || anyNA(bg$end) ||
        any(!is.finite(bg$start) | !is.finite(bg$end) | bg$start < 0 | bg$end <= bg$start |
            bg$start != floor(bg$start) | bg$end != floor(bg$end))) stop("Invalid bedGraph coordinates: ", f$file)
    insulation_list[[j]] <- rbindlist(lapply(seq_len(nrow(loc)), function(k) {
      locus <- loc[k]
      # Select the half-open BED interval containing the query position.
      b <- unique(bg[chrom == locus$chrom & start <= locus$position & end > locus$position])
      if (nrow(b) > 1L) stop("Overlapping insulation bins at ", locus$sv_id, " in ", f$file)
      data.table(sample = f$sample, sv_id = locus$sv_id,
        insulation = if (nrow(b)) b$score else NA_real_,
        bin_start = if (nrow(b)) b$start else NA_real_,
        bin_end = if (nrow(b)) b$end else NA_real_, insulation_file = f$file)
    }))
  }
  ins <- rbindlist(insulation_list)
  pairs <- merge(pairs, ins, by = c("sample", "sv_id"), all.x = TRUE, sort = FALSE)
  pairs <- merge(pairs, loc, by = "sv_id", all.x = TRUE, sort = FALSE)
} else {
  message("Replotting the existing paired tables; insulation is not recomputed.")
  for (i in seq_len(nrow(targets))) {
    t <- targets[i]
    gid <- if (is.na(t$gene_id)) t$gene_name else t$gene_id
    path <- file.path(INTEGRATION_DIR, safe_name(t$sv_id), paste0(safe_name(gid), ".insulation_expression.tsv"))
    if (!file.exists(path)) stop("Paired table not found: ", path, "; use rebuild mode or explicit gene_id in TARGET_FILE")
    e <- fread(path, na.strings = c("", "NA", "NaN"))
    need(e, c("sample", "TPM", "insulation"), path)
    for (key in c("sv_id", "gene_id")) {
      expected <- if (key == "sv_id") t$sv_id else gid
      if (key %in% names(e) && any(is.na(e[[key]]) | e[[key]] != expected)) stop("Mismatched ", key, " in ", path)
    }
    if (!"dosage" %in% names(e)) {
      need(e, "gt_class", path)
      e[, dosage := gt_dosage(gt_class)]
    }
    e[, `:=`(TPM = num(TPM, "TPM"), insulation = num(insulation, "insulation"), dosage = num(dosage, "dosage"))]
    if ("gt_class" %in% names(e)) {
      gtd <- gt_dosage(e$gt_class)
      if (any(!is.na(gtd) & !is.na(e$dosage) & gtd != e$dosage)) stop("Genotype/dosage mismatch in ", path)
    }
    e <- dedup_checked(e, "sample", c("TPM", "insulation", "dosage"), paste0(t$gene_name, "_paired"))
    e <- e[, .(sample, TPM, insulation, dosage)]
    e[, `:=`(sv_id = t$sv_id, gene_id = gid, plot_gene = t$gene_name, panel_index = i,
               paired_source = path)]
    targets[i, gene_id := gid]
    pair_list[[i]] <- e
  }
  pairs <- rbindlist(pair_list, fill = TRUE)
}

############################################################
## Complete observations and output tables
############################################################

if (anyNA(pairs$sample) || any(!nzchar(pairs$sample))) stop("Missing sample IDs")
if (any(!is.na(pairs$TPM) & (!is.finite(pairs$TPM) | pairs$TPM < 0))) stop("Invalid TPM")
if (any(!is.na(pairs$dosage) & !pairs$dosage %in% 0:2)) stop("Invalid genotype dosage")
if (anyDuplicated(targets, by = c("sv_id", "gene_id"))) stop("Duplicate target associations")
out(targets, "resolved_targets.tsv")
pairs[, `:=`(missing_expression = !is.finite(TPM), missing_genotype = is.na(dosage),
               missing_insulation = !is.finite(insulation))]
pairs[, included := !missing_expression & !missing_genotype & !missing_insulation]
pairs[, gt_class := factor(GT_LEVELS[dosage + 1L], levels = GT_LEVELS)]
out(pairs, "sample_QC.tsv")
qc <- pairs[, .(n_expression_rows = .N, n_missing_expression = sum(missing_expression),
  n_missing_genotype = sum(missing_genotype), n_missing_insulation = sum(missing_insulation),
  n_complete = sum(included)), by = .(plot_gene, sv_id, gene_id)]
out(qc, "sample_counts.tsv")
print(qc)
d <- pairs[included == TRUE]
# Numerically stable log2(TPM + 1), retaining zero TPM.
d[, log2TPM := log1p(TPM) / log(2)]
setorder(d, panel_index, dosage, sample)
out(d, "plot_data.tsv")
for (i in seq_len(nrow(targets))) {
  t <- targets[i]
  z <- d[panel_index == i]
  if (nrow(z) < 3L) stop("Fewer than 3 complete samples for ", t$gene_name, "; inspect sample_QC.tsv")
  out(z, paste0(safe_name(t$gene_name), ".", safe_name(t$sv_id), ".insulation_expression.tsv"))
}

############################################################
## Correlation statistics
############################################################

cor_stats <- function(x, y) {
  if (length(x) < 3L || uniqueN(x) < 2L || uniqueN(y) < 2L) {
    return(list(spearman_rho = NA_real_, p_spearman = NA_real_, pearson_r = NA_real_, p_pearson = NA_real_))
  }
  sp <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE, alternative = "two.sided"))
  pe <- cor.test(x, y, method = "pearson", alternative = "two.sided")
  list(spearman_rho = unname(sp$estimate), p_spearman = sp$p.value,
       pearson_r = unname(pe$estimate), p_pearson = pe$p.value)
}
stats_list <- list()
for (i in seq_len(nrow(targets))) {
  t <- targets[i]
  z <- d[panel_index == i]
  for (mode in c("TPM", "log2TPM")) {
    st <- cor_stats(z$insulation, z[[mode]])
    stats_list[[length(stats_list) + 1L]] <- cbind(data.table(panel_index = i,
      gene_name = t$gene_name, gene_id = t$gene_id, sv_id = t$sv_id, scale = mode,
      n = nrow(z), n_AA = sum(z$dosage == 0), n_Aa = sum(z$dosage == 1), n_aa = sum(z$dosage == 2)), as.data.table(st))
  }
}
stats <- rbindlist(stats_list)
out(stats, "correlation_summary.tsv")
# Check that both scales preserve Spearman statistics on the same individuals.
for (i in seq_len(nrow(targets))) {
  a <- stats[panel_index == i & scale == "TPM"]
  b <- stats[panel_index == i & scale == "log2TPM"]
  if (!isTRUE(all.equal(a$spearman_rho, b$spearman_rho)) ||
      !isTRUE(all.equal(a$p_spearman, b$p_spearman))) stop("Spearman invariance check failed")
}

############################################################
## Plotting
############################################################

p_math <- function(p) {
  if (is.na(p)) return('"NA"')
  if (p == 0) return('"< 1e-300"')
  if (p >= 0.001) return(formatC(p, digits = 3, format = "fg"))
  exponent <- floor(log10(p))
  mantissa <- signif(p / 10^exponent, 2)
  if (mantissa >= 10) { mantissa <- mantissa / 10; exponent <- exponent + 1 }
  sprintf('%.2g %%*%% 10^{%d}', mantissa, exponent)
}
p_statement <- function(p, lhs = "italic(P)") {
  if (!is.na(p) && p == 0) paste(lhs, "< 10^{-300}") else paste(lhs, "==", p_math(p))
}
make_panel <- function(i, mode) {
  t <- targets[i]
  z <- copy(d[panel_index == i])
  z[, expression_value := get(mode)]
  st <- stats[panel_index == i & scale == mode]
  n_gt <- table(z$gt_class)
  legend_labels <- setNames(paste0(GT_LABELS[GT_LEVELS], " (n=", as.integer(n_gt), ")"), GT_LEVELS)
  rho_text <- if (is.na(st$spearman_rho)) '"NA"' else sprintf("%.2f", st$spearman_rho)
  annotation <- sprintf('atop("Spearman"~rho == %s, %s~~italic(n) == %d)',
                         rho_text, p_statement(st$p_spearman), st$n)
  top <- max(z$expression_value)
  top <- if (top > 0) top * 1.33 else 1
  p <- ggplot(z, aes(x = insulation, y = expression_value))
  if (uniqueN(z$insulation) >= 2L && uniqueN(z$expression_value) >= 2L) {
    p <- p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
                          color = "#333333", fill = "#BBBBBB", alpha = 0.30, linewidth = 0.55)
  }
  p + geom_point(aes(color = gt_class), size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = COLORS, limits = GT_LEVELS, labels = legend_labels, drop = FALSE,
                        name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.07, 0.07))) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    coord_cartesian(ylim = c(0, top)) +
    annotate("text", x = Inf, y = top * 0.97, label = annotation, parse = TRUE,
               hjust = 1.04, vjust = 1, size = 2.7, family = FONT_FAMILY) +
    labs(title = bquote(italic(.(t$gene_name))), subtitle = t$sv_id,
         x = "Local insulation score",
         y = if (mode == "TPM") "TPM" else expression(log[2](TPM + 1))) +
    theme_classic(base_size = 9, base_family = FONT_FAMILY) +
    theme(plot.title = element_text(hjust = 0.5, size = 11, margin = margin(b = 3)),
          plot.subtitle = element_text(hjust = 0.5, size = 6.8, color = "#444444", margin = margin(b = 5)),
          axis.text = element_text(size = 8, color = "black"),
          axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.3),
          axis.ticks.length = grid::unit(2, "pt"),
          axis.title.x = element_text(margin = margin(t = 5)),
          axis.title.y = element_text(margin = margin(r = 5)),
          legend.position = "bottom", legend.text = element_text(size = 7),
          legend.key.width = grid::unit(8, "pt"), legend.spacing.x = grid::unit(2, "pt"),
          plot.margin = margin(8, 9, 5, 7, "pt")) +
    guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 1.8)))
}
combine_panels <- function(plots) {
  grobs <- lapply(plots, ggplotGrob)
  if (length(unique(vapply(grobs, function(g) length(g$widths), integer(1)))) == 1L) {
    widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
    grobs <- lapply(grobs, function(g) { g$widths <- widths; g })
  }
  kids <- lapply(seq_along(grobs), function(i) {
    grid::grobTree(grobs[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
  })
  grid::gTree(children = do.call(grid::gList, kids),
              vp = grid::viewport(layout = grid::grid.layout(1, length(grobs))))
}
save_plot <- function(p, stem, width) {
  ggsave(file.path(OUTDIR, paste0(stem, ".pdf")), plot = p, device = grDevices::pdf,
         width = width, height = FIG_HEIGHT, useDingbats = FALSE, bg = "white")
  ggsave(file.path(OUTDIR, paste0(stem, ".png")), plot = p,
         width = width, height = FIG_HEIGHT, dpi = 600, bg = "white")
}
############################################################
## Save figures and analysis notes
############################################################

for (mode in c("TPM", "log2TPM")) {
  plots <- lapply(seq_len(nrow(targets)), make_panel, mode = mode)
  for (i in seq_along(plots)) {
    save_plot(plots[[i]], paste0(safe_name(targets$gene_name[i]), ".", safe_name(targets$sv_id[i]), ".", mode), FIG_WIDTH / 2)
  }
  save_plot(combine_panels(plots), paste0("Insulation_Expression.", mode), FIG_WIDTH * nrow(targets) / 2)
}
writeLines(c(
  "Each point represents one individual with nonmissing expression, SV genotype and local insulation.",
  "Local insulation is the CHM13 bedGraph score from the bin containing the specified SV start.",
  "Color indicates SV alternative-allele dosage: 0/0, 0/1 or 1/1.",
  "TPM and log2(TPM+1) plots use the same individuals; zero TPM is retained.",
  "Annotations report two-sided Spearman rho and nominal P values, without multiple-testing adjustment.",
  "Black lines show pooled ordinary least-squares fits; gray bands indicate 95% confidence intervals for the fitted mean.",
  "Pearson correlation is calculated separately on both scales in correlation_summary.tsv.",
  "No adjustment for genotype, population structure or other covariates is performed.",
  paste0("Mode: ", MODE), paste0("Expression: ", EXPR_FILE), paste0("Genotypes: ", GENO_FILE),
  paste0("Integration directory (replot only): ", INTEGRATION_DIR),
  if (MODE == "rebuild") paste0("SV ID coordinate offset used: ", OFFSET) else
    "Replot uses previously extracted insulation; the original extraction convention is unchanged."
), file.path(OUTDIR, "analysis_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
message("Done: ", normalizePath(OUTDIR))
print(stats)
