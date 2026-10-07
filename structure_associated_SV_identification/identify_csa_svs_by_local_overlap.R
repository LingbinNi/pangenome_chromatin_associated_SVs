#!/usr/bin/env Rscript
# Identify chromatin-structure-associated SVs (CSA-SVs) by local CDH overlap.
# Requires GenomicRanges, data.table and ggplot2.
# Edit cdh_bed_dir, sv_file and outdir below before running.
# Run: Rscript /path/to/identify_csa_svs_by_local_overlap.R
# CDH input: <cdh_bed_dir>/<sample>.<Hap1|Hap2>.CDH.CHM13.bed
# Supply headerless BED files with at least three columns: chr, start, end.
# CDH intervals use 0-based starts and half-open ends in CHM13 coordinates.
# SV input: a tab-separated table with the following header columns:
#   sample, chr, pos, sv_id, svtype, svlen, hap1_carrier, hap2_carrier
# Supply heterozygous INS/DEL records prepared upstream, with 0/1 carrier flags.
# pos is a 0-based SV position in the same CHM13 coordinates as the CDHs;
# chromosome names must match between the inputs. svlen is retained as metadata.
# Each sv_id must have consistent chr, pos, svtype and svlen across samples.
# The analyzed cohort contains samples present in the SV table with nonempty
# CDH records for both Hap1 and Hap2. All association counts refer to haplotypes.
# With flip_del=TRUE, INS uses the reported carrier haplotype and DEL uses
# the opposite haplotype of the same sample. These define effective carriers.
# For each SV, analyzed haplotypes absent from its effective carrier list are coded 0.
# Eligible SVs have at least one effective carrier in the analyzed cohort.
# A haplotype is CDH-positive if at least one of its own CDHs overlaps the SV query.
# Multiple overlaps for the same SV and haplotype contribute one CDH-positive count.
# Default flanks are 0, 5,000, 10,000, 25,000 and 50,000 bp, analyzed separately.
# flank=0 tests the 1-bp interval [pos, pos+1); each positive flank extends this
# interval on both sides, with its start clipped at chromosome coordinate 0.
# Test each SV using a two-sided Fisher test on the four haplotype counts.
# Report OR with a pseudocount of 1 added to each cell; formulas are given below.
# BH correction is applied separately for each flank, pooling INS and DEL loci.
# p_adj_overlap_only uses loci with at least one local CDH-positive haplotype.
# p_adj_all_eligible uses all eligible loci, assigning P=1 to loci with no local CDH.
# By default, p_adj is p_adj_all_eligible and defines the keep set.
# keep: OR > 2, p_adj < 0.05 and at least 3 effective carriers with a local CDH.
# keep2 additionally requires specificity > 0.5; its formula is given below.
# Per-flank outputs: cdh_sv.local.f<flank>.tsv, .keep.tsv and .keep2.tsv.
# The main TSV contains eligible SVs with at least one local CDH-positive haplotype.
# A nonempty keep set also produces enrichment.local.f<flank>.pdf.
# Combined summary: local_overlap_threshold_summary.tsv. All TSVs have headers.
# Flanks with no SV-local CDH overlaps are skipped.
# Use a new output directory.

suppressMessages({
    library(GenomicRanges)
    library(data.table)
    library(ggplot2)
})

################################################################
## Parameters / paths
################################################################

cdh_bed_dir <- "/path/to/CDH_BED_CHM13"

sv_file <- "/path/to/SV_CHM13_hap_carrier.tsv"

outdir <- "./local_overlap"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# FALSE reads the standard CDH BED files.
# TRUE also includes <sample>.<Hap1|Hap2>.CDH.CHM13.transloc.bed files.
include_transloc <- FALSE

# TRUE pairs INS with the reported carrier haplotype and DEL with its partner.
# FALSE uses the reported carrier haplotype for both SV types.
flip_del <- TRUE

# Analyze each listed flank separately; use c(0) to run only the SV position.
flank_list <- c(
    0,
    5000,
    10000,
    25000,
    50000
)

# keep requires OR > OR_cut, p_adj < padj_cut and joint-positive count >= min_a_cut.
# keep2 additionally requires specificity > spec_cut.
OR_cut    <- 2
padj_cut  <- 0.05
min_a_cut <- 3
spec_cut  <- 0.5

# Select the BH-adjusted P column used as p_adj for the keep/keep2 sets.
# "all_eligible": all eligible SV loci; "overlap_only": loci with local CDH occurrence.
primary_fdr_scope <- "all_eligible"

if (!primary_fdr_scope %in% c("all_eligible", "overlap_only")) {
    stop("primary_fdr_scope must be 'all_eligible' or 'overlap_only'")
}

################################################################
## 1. Read haplotype-specific CDH intervals
################################################################

pat <- if (include_transloc) {
    "\\.CDH\\.CHM13(\\.transloc)?\\.bed$"
} else {
    "\\.CDH\\.CHM13\\.bed$"
}

cdh_files <- list.files(
    cdh_bed_dir,
    pattern = pat,
    full.names = TRUE
)

stopifnot(length(cdh_files) > 0)

message(
    "CDH files: ", length(cdh_files),
    "  transloc files: ",
    sum(grepl("\\.transloc\\.bed$", cdh_files))
)

read_cdh <- function(f) {

    if (file.info(f)$size == 0) {
        return(NULL)
    }

    d <- tryCatch(
        fread(
            f,
            header = FALSE,
            select = 1:3,
            col.names = c("chr", "start", "end"),
            showProgress = FALSE
        ),
        error = function(e) {
            message("Failed to read: ", f, "  ", conditionMessage(e))
            NULL
        }
    )

    if (is.null(d) || nrow(d) == 0) {
        return(NULL)
    }

    hap <- sub(
        "^(.+)\\.(Hap[12])\\.CDH\\.CHM13(\\.transloc)?\\.bed$",
        "\\1_\\2",
        basename(f)
    )

    d[, hap := hap]

    d
}

cdh_dt <- rbindlist(
    lapply(cdh_files, read_cdh),
    use.names = TRUE,
    fill = TRUE
)

stopifnot(nrow(cdh_dt) > 0)

# BED is 0-based, half-open. Convert to GRanges 1-based closed coordinates.
gr_cdh_all <- GRanges(
    seqnames = cdh_dt$chr,
    ranges = IRanges(
        start = as.integer(cdh_dt$start) + 1L,
        end   = as.integer(cdh_dt$end)
    ),
    haplotype = cdh_dt$hap
)

cdh_haps <- sort(unique(cdh_dt$hap))

message(
    "CDH haplotypes: ", length(cdh_haps),
    "  total CDH intervals: ", length(gr_cdh_all)
)

################################################################
## 2. Read heterozygous SV carrier table
################################################################

sv_dt <- fread(
    sv_file,
    header = TRUE,
    showProgress = TRUE
)

required_sv_cols <- c(
    "sample",
    "chr",
    "pos",
    "sv_id",
    "svtype",
    "svlen",
    "hap1_carrier",
    "hap2_carrier"
)

missing_sv_cols <- setdiff(
    required_sv_cols,
    names(sv_dt)
)

if (length(missing_sv_cols) > 0) {
    stop(
        "SV table missing columns: ",
        paste(missing_sv_cols, collapse = ", ")
    )
}

sv_dt[, svtype := toupper(as.character(svtype))]
sv_dt[, pos := as.integer(pos)]
sv_dt[, svlen := abs(as.integer(svlen))]
sv_dt[, hap1_carrier := as.integer(hap1_carrier)]
sv_dt[, hap2_carrier := as.integer(hap2_carrier)]

bad_carrier <- sv_dt[
    is.na(hap1_carrier) |
    is.na(hap2_carrier) |
    !hap1_carrier %in% 0:1 |
    !hap2_carrier %in% 0:1
]

if (nrow(bad_carrier) > 0) {
    stop(
        "Carrier table contains invalid hap1/hap2 carrier values; n=",
        nrow(bad_carrier)
    )
}

# Same sv_id must represent one genomic event.
sv_def_check <- sv_dt[
    ,
    .(
        n_def = uniqueN(
            paste(chr, pos, svtype, svlen, sep = "|")
        )
    ),
    by = sv_id
]

if (any(sv_def_check$n_def > 1L)) {
    stop(
        "Some sv_id values have conflicting genomic definitions. Example: ",
        paste(
            head(
                sv_def_check[n_def > 1, sv_id],
                5
            ),
            collapse = ", "
        )
    )
}

sv_info <- sv_dt[
    ,
    .(
        chr    = chr[1L],
        pos    = pos[1L],
        svtype = svtype[1L],
        svlen  = svlen[1L]
    ),
    by = sv_id
]

sv_haps <- sort(
    unique(
        c(
            paste0(sv_dt$sample, "_Hap1"),
            paste0(sv_dt$sample, "_Hap2")
        )
    )
)

message(
    "SV haplotypes: ", length(sv_haps),
    "  unique SV loci: ", nrow(sv_info),
    "  DEL=", sum(sv_info$svtype == "DEL"),
    "  INS=", sum(sv_info$svtype == "INS")
)

################################################################
## 3. Common paired haplotypes
################################################################

common_haps <- intersect(
    cdh_haps,
    sv_haps
)

message(
    "Common haplotypes before paired-sample filter: ",
    length(common_haps)
)

stopifnot(length(common_haps) > 0)

samp_of <- sub(
    "_Hap[12]$",
    "",
    common_haps
)

both_ok <- names(
    which(
        table(samp_of) == 2L
    )
)

n_drop <- length(common_haps) -
    2L * length(both_ok)

if (n_drop > 0) {
    message(
        "Dropping ", n_drop,
        " haplotypes whose partner haplotype is missing."
    )
}

common_haps <- sort(
    common_haps[
        samp_of %in% both_ok
    ]
)

n_haps <- length(common_haps)

message(
    "Analyzed haplotypes: ", n_haps,
    "  samples: ", length(both_ok)
)

# Restrict CDHs to exactly the analyzed haplotypes.
gr_cdh_all <- gr_cdh_all[
    gr_cdh_all$haplotype %in% common_haps
]

################################################################
## 4. Effective carrier haplotypes
##
## With flip_del=TRUE: INS uses the reported carrier; DEL uses its partner.
################################################################

sub <- copy(sv_dt)

if (flip_del) {

    sub[
        ,
        `:=`(
            eff1 = fifelse(
                svtype == "DEL",
                hap2_carrier,
                hap1_carrier
            ),
            eff2 = fifelse(
                svtype == "DEL",
                hap1_carrier,
                hap2_carrier
            )
        )
    ]

} else {

    sub[
        ,
        `:=`(
            eff1 = hap1_carrier,
            eff2 = hap2_carrier
        )
    ]
}

car_long <- rbindlist(
    list(
        sub[
            eff1 == 1,
            .(
                sv_id,
                hap = paste0(sample, "_Hap1")
            )
        ],
        sub[
            eff2 == 1,
            .(
                sv_id,
                hap = paste0(sample, "_Hap2")
            )
        ]
    ),
    use.names = TRUE
)

car_long <- unique(
    car_long[
        hap %chin% common_haps
    ],
    by = c("sv_id", "hap")
)

sv_count_dt <- car_long[
    ,
    .(
        sv_count = uniqueN(hap)
    ),
    by = sv_id
]

# Eligible = the SV has at least one effective carrier haplotype
# among the analyzed common haplotypes.
eligible_info <- merge(
    sv_info,
    sv_count_dt,
    by = "sv_id",
    all = FALSE,
    sort = FALSE
)

message(
    "Eligible SV loci with >=1 effective carrier haplotype: ",
    nrow(eligible_info)
)

stopifnot(nrow(eligible_info) > 0)

################################################################
## 5. Fisher helper
################################################################

fisher_p <- function(a, b, c, d) {

    fisher.test(
        matrix(
            c(a, b, c, d),
            nrow = 2
        ),
        alternative = "two.sided"
    )$p.value
}

################################################################
## 6. Direct local-overlap analysis
################################################################

summary_list <- list()

for (flank_bp in flank_list) {

    tag <- paste0("f", flank_bp)

    message("")
    message("============================================================")
    message("SV-centered local flank = ", flank_bp, " bp")
    message("============================================================")

    ############################################################
    # 6A. Build one local query interval per eligible SV
    #
    # pos in carrier table is 0-based.
    # pos + 1 is the 1-based SV breakpoint/position.
    ############################################################

    point_1based <- eligible_info$pos + 1L

    query_start <- pmax(
        1L,
        point_1based - flank_bp
    )

    query_end <- point_1based + flank_bp

    gr_sv_local <- GRanges(
        seqnames = eligible_info$chr,
        ranges = IRanges(
            start = query_start,
            end = query_end
        ),
        sv_id = eligible_info$sv_id
    )

    ############################################################
    # 6B. Local CDH presence per SV x haplotype
    #
    # Mark each SV/haplotype once when that haplotype has at least one CDH
    # overlapping the SV query interval. Retain the haplotype identity of each CDH.
    ############################################################

    ov <- findOverlaps(
        gr_sv_local,
        gr_cdh_all,
        ignore.strand = TRUE
    )

    if (length(ov) == 0L) {

        warning(
            "No SV-local CDH overlaps at flank ",
            flank_bp,
            "; skipping."
        )

        next
    }

    cdh_long <- unique(
        data.table(
            sv_id = eligible_info$sv_id[
                queryHits(ov)
            ],
            hap = gr_cdh_all$haplotype[
                subjectHits(ov)
            ]
        ),
        by = c("sv_id", "hap")
    )

    cdh_count_dt <- cdh_long[
        ,
        .(
            cdh_count = uniqueN(hap)
        ),
        by = sv_id
    ]

    ############################################################
    # 6C. a = effective carrier AND local CDH
    ############################################################

    a_dt <- cdh_long[
        car_long,
        on = .(
            sv_id,
            hap
        ),
        nomatch = 0L
    ][
        ,
        .(
            w_cdh_w_sv = .N
        ),
        by = sv_id
    ]

    ############################################################
    # 6D. Assemble one 2x2 table per eligible SV
    # a = w_cdh_w_sv: local CDH present, effective carrier.
    # b = w_cdh_wo_sv: local CDH present, effective noncarrier.
    # c = wo_cdh_w_sv: local CDH absent, effective carrier.
    # d = wo_cdh_wo_sv: local CDH absent, effective noncarrier.
    # a + b + c + d equals the number of analyzed haplotypes.
    ############################################################

    result_all <- merge(
        eligible_info,
        cdh_count_dt,
        by = "sv_id",
        all.x = TRUE,
        sort = FALSE
    )

    result_all <- merge(
        result_all,
        a_dt,
        by = "sv_id",
        all.x = TRUE,
        sort = FALSE
    )

    result_all[
        is.na(cdh_count),
        cdh_count := 0L
    ]

    result_all[
        is.na(w_cdh_w_sv),
        w_cdh_w_sv := 0L
    ]

    result_all[
        ,
        `:=`(
            w_cdh_wo_sv =
                cdh_count - w_cdh_w_sv,

            wo_cdh_w_sv =
                sv_count - w_cdh_w_sv
        )
    ]

    result_all[
        ,
        wo_cdh_wo_sv :=
            n_haps -
            w_cdh_w_sv -
            w_cdh_wo_sv -
            wo_cdh_w_sv
    ]

    bad_table <- result_all[
        w_cdh_w_sv < 0 |
        w_cdh_wo_sv < 0 |
        wo_cdh_w_sv < 0 |
        wo_cdh_wo_sv < 0 |
        (
            w_cdh_w_sv +
            w_cdh_wo_sv +
            wo_cdh_w_sv +
            wo_cdh_wo_sv
        ) != n_haps
    ]

    if (nrow(bad_table) > 0) {
        stop(
            "Invalid 2x2 table(s) generated at flank ",
            flank_bp,
            ". Example SV: ",
            bad_table$sv_id[1]
        )
    }

    ############################################################
    # 6E. Fisher P values
    #
    # If cdh_count == 0, the table contains no CDH occurrence
    # and Fisher P is exactly 1, so we set it directly.
    ############################################################

    result_all[, p := 1.0]

    idx_test <- which(
        result_all$cdh_count > 0L
    )

    if (length(idx_test) > 0L) {

        result_all$p[idx_test] <- vapply(
            idx_test,
            function(i) {

                fisher_p(
                    result_all$w_cdh_w_sv[i],
                    result_all$w_cdh_wo_sv[i],
                    result_all$wo_cdh_w_sv[i],
                    result_all$wo_cdh_wo_sv[i]
                )
            },
            numeric(1)
        )
    }

    ############################################################
    # 6F. Effect size and specificity
    # OR = ((a+1)*(d+1))/((b+1)*(c+1)).
    # specificity = (a+1)/(a+c+2): the smoothed fraction of effective
    # carrier haplotypes with a local CDH.
    ############################################################

    result_all[
        ,
        OR :=
            (w_cdh_w_sv + 1) *
            (wo_cdh_wo_sv + 1) /
            (
                (w_cdh_wo_sv + 1) *
                (wo_cdh_w_sv + 1)
            )
    ]

    result_all[
        ,
        specificity :=
            (w_cdh_w_sv + 1) /
            (w_cdh_w_sv + wo_cdh_w_sv + 2)
    ]

    ############################################################
    # 6G. Multiple-testing correction
    #
    # overlap-only:
    #   SVs with >=1 local CDH-positive haplotype form the BH universe.
    #
    # all-eligible:
    #   all eligible SV loci are in the BH universe;
    #   no-CDH loci contribute P=1.
    ############################################################

    result_all[, p_adj_overlap_only := 1.0]

    if (length(idx_test) > 0L) {
        result_all$p_adj_overlap_only[idx_test] <-
            p.adjust(
                result_all$p[idx_test],
                method = "BH"
            )
    }

    result_all[
        ,
        p_adj_all_eligible :=
            p.adjust(
                p,
                method = "BH"
            )
    ]

    if (primary_fdr_scope == "all_eligible") {

        result_all[
            ,
            p_adj := p_adj_all_eligible
        ]

    } else {

        result_all[
            ,
            p_adj := p_adj_overlap_only
        ]
    }

    result_all[
        ,
        flank_bp := flank_bp
    ]

    ############################################################
    # 6H. Output SVs with at least one local CDH-positive haplotype
    #
    # SVs with no local CDH contribute P=1 to the all-eligible BH correction.
    # Output rows contain the subset with cdh_count > 0.
    ############################################################

    result <- result_all[
        cdh_count > 0L
    ]

    setcolorder(
        result,
        c(
            "chr",
            "pos",
            "sv_id",
            "svtype",
            "svlen",
            "flank_bp",
            "w_cdh_w_sv",
            "w_cdh_wo_sv",
            "wo_cdh_w_sv",
            "wo_cdh_wo_sv",
            "cdh_count",
            "sv_count",
            "OR",
            "p",
            "p_adj_overlap_only",
            "p_adj_all_eligible",
            "p_adj",
            "specificity"
        )
    )

    setorder(
        result,
        p_adj,
        -OR
    )

    ############################################################
    # 6I. Significant sets
    ############################################################

    keep <- result[
        OR > OR_cut &
        p_adj < padj_cut &
        w_cdh_w_sv >= min_a_cut
    ]

    keep2 <- keep[
        specificity > spec_cut
    ]

    ############################################################
    # 6J. Output
    ############################################################

    fwrite(
        result,
        file.path(
            outdir,
            paste0(
                "cdh_sv.local.",
                tag,
                ".tsv"
            )
        ),
        sep = "\t"
    )

    fwrite(
        keep,
        file.path(
            outdir,
            paste0(
                "cdh_sv.local.",
                tag,
                ".keep.tsv"
            )
        ),
        sep = "\t"
    )

    fwrite(
        keep2,
        file.path(
            outdir,
            paste0(
                "cdh_sv.local.",
                tag,
                ".keep2.tsv"
            )
        ),
        sep = "\t"
    )

    ############################################################
    # 6K. Diagnostic plot
    ############################################################

    if (nrow(keep) > 0L) {

        pp <- ggplot(
            keep,
            aes(
                x = cdh_count,
                y = sv_count
            )
        ) +
            geom_point(
                aes(
                    size = log2(OR),
                    color = -log10(p_adj)
                ),
                alpha = 0.75
            ) +
            scale_color_gradient(
                low = "lightblue",
                high = "red"
            ) +
            scale_size(
                range = c(2, 9)
            ) +
            facet_wrap(
                ~ svtype
            ) +
            theme_classic(
                base_size = 13
            ) +
            labs(
                x = "Haplotypes with local CDH",
                y = "Effective SV-carrier haplotypes",
                color = "-log10(FDR)",
                size = "log2(OR)",
                title = paste0(
                    "Direct local CSA-SV association; flank = ",
                    flank_bp,
                    " bp"
                )
            )

        ggsave(
            file.path(
                outdir,
                paste0(
                    "enrichment.local.",
                    tag,
                    ".pdf"
                )
            ),
            pp,
            width = 10,
            height = 6
        )
    }

    ############################################################
    # 6L. Summary
    ############################################################

    summary_list[[tag]] <- data.table(
        flank_bp = flank_bp,
        n_analyzed_haplotypes = n_haps,
        n_eligible_svs = nrow(result_all),
        n_svs_with_local_cdh = nrow(result),
        n_sig_keep = nrow(keep),
        n_sig_keep2 = nrow(keep2),
        keep_DEL = sum(keep$svtype == "DEL"),
        keep_INS = sum(keep$svtype == "INS"),
        keep2_DEL = sum(keep2$svtype == "DEL"),
        keep2_INS = sum(keep2$svtype == "INS"),
        primary_fdr_scope = primary_fdr_scope,
        include_transloc = include_transloc
    )

    message(
        "Eligible SVs = ", nrow(result_all),
        "; local-CDH SVs = ", nrow(result),
        "; keep = ", nrow(keep),
        "; keep2 = ", nrow(keep2)
    )
}

################################################################
## 7. Summary
################################################################

if (length(summary_list) == 0L) {
    stop("No flank analysis completed.")
}

summary_df <- rbindlist(
    summary_list,
    use.names = TRUE,
    fill = TRUE
)

fwrite(
    summary_df,
    file.path(
        outdir,
        "local_overlap_threshold_summary.tsv"
    ),
    sep = "\t"
)

print(summary_df)

message(
    "Done. Output directory: ",
    normalizePath(outdir)
)
