#!/usr/bin/env Rscript
# Test associations between CSA-SV genotype and gene-level expression.
# Requires R, data.table and the gzip command in PATH.
# Run from the directory containing AdjacentGenes.f25.tsv.gz:
#   Rscript /path/to/test_csa_sv_expression_associations.R
# Edit ROOT, INFILE or OUTDIR below to use different input/output locations.
# Input: a gzip-compressed, tab-separated table with the header columns:
#   sv_id, gene_id, gene_name, TPM, genotype, dosage, sample, svtype, svlen
# Supply the merged candidate SV-gene table prepared upstream using the 500-kb
# cis window and gene-level IsoQuant TPM values. This script tests the supplied pairs.
# Each record describes one individual's expression and genotype for an SV-gene pair.
# Genotype labels: AA = 0/0, Aa = 0/1, aa = 1/1; phased equivalents and 1/0 are accepted.
# Calculate dosage 0/1/2 from genotype and check agreement with nonmissing numeric
# entries in the supplied dosage column. Use the calculated dosage for trend tests.
# Exclude unrecognized/missing genotypes, missing TPM, empty/missing gene_id and
# gene_id values beginning with novel_gene. Zero TPM values are retained.
# Keep the first retained record for each sv_id/gene_id/sample combination.
# Test pairs with at least two observed genotype classes, each containing at least
# MIN_GROUP_N individuals; the default is 1. Use untransformed TPM for all tests.
# Primary test: Kruskal-Wallis across the observed genotype classes.
# Additional test: Spearman correlation between genotype dosage and TPM, using
# exact=FALSE when both variables have at least two distinct values.
# Apply Benjamini-Hochberg correction across the SV-gene result table separately
# for Kruskal-Wallis P values (q_kw) and Spearman P values (q_trend).
# Report sample counts and median TPM for each genotype class alongside test results.
# Outputs under association_f25/ by default; all tables have headers:
#   eQTL_scan_all.tsv: all pairs passing the genotype-class requirements.
#   eQTL_scan_sig_KW.tsv: primary significant associations with q_kw < 0.05.
#   eQTL_scan_sig_trend.tsv: associations with q_trend < 0.05.
#   eQTL_scan_sig_any.tsv: union of the KW- and trend-significant associations.
# Use a new output directory. Summary counts are printed to the console.

suppressMessages({
    library(data.table)
})


############################################################
# Paths
############################################################

ROOT <- normalizePath(getwd())

INFILE <- file.path(
    ROOT,
    "AdjacentGenes.f25.tsv.gz"
)

OUTDIR <- file.path(
    ROOT,
    "association_f25"
)

dir.create(
    OUTDIR,
    showWarnings = FALSE,
    recursive = TRUE
)


############################################################
# Parameters
############################################################

# Minimum number of individuals in every genotype class present for a pair.
# At least two genotype classes must be represented.
MIN_GROUP_N <- 1L


############################################################
# Load
############################################################

cat("Loading merged table...\n")

dat <- fread(
    cmd = paste(
        "gzip -cd",
        shQuote(INFILE)
    ),
    sep = "\t",
    header = TRUE,
    na.strings = c("NA", "")
)

setDT(dat)


required <- c(
    "sv_id",
    "gene_id",
    "gene_name",
    "TPM",
    "genotype",
    "dosage",
    "sample",
    "svtype",
    "svlen"
)

missing_cols <- setdiff(
    required,
    names(dat)
)

if (length(missing_cols) > 0) {
    stop(
        "Missing columns: ",
        paste(
            missing_cols,
            collapse = ","
        )
    )
}


cat(
    sprintf(
        "Raw rows: %d\n",
        nrow(dat)
    )
)

cat(
    sprintf(
        "Samples: %d\n",
        uniqueN(dat$sample)
    )
)

cat(
    sprintf(
        "SVs: %d\n",
        uniqueN(dat$sv_id)
    )
)


############################################################
# Genotype class
############################################################

dat[, gt_class := fcase(

    genotype %chin% c(
        "0|0",
        "0/0"
    ),
    "AA",

    genotype %chin% c(
        "0|1",
        "1|0",
        "0/1",
        "1/0"
    ),
    "Aa",

    genotype %chin% c(
        "1|1",
        "1/1"
    ),
    "aa",

    default = NA_character_
)]


dat[, dosage_calc := fcase(

    gt_class == "AA", 0,

    gt_class == "Aa", 1,

    gt_class == "aa", 2,

    default = NA_real_
)]


############################################################
# QC existing dosage column
############################################################

dat[, dosage_file :=
      suppressWarnings(
          as.numeric(dosage)
      )]


n_dosage_mismatch <- dat[
    !is.na(dosage_file) &
    !is.na(dosage_calc) &
    dosage_file != dosage_calc,
    .N
]


if (n_dosage_mismatch > 0) {

    stop(
        "Dosage mismatch detected: ",
        n_dosage_mismatch,
        " rows"
    )
}


cat("\nGenotype distribution:\n")

print(
    dat[
        ,
        .N,
        by = .(
            genotype,
            gt_class
        )
    ][
        order(-N)
    ]
)


############################################################
# TPM / gene QC
############################################################

dat[, TPM :=
      suppressWarnings(
          as.numeric(TPM)
      )]


dat[, is_novel :=
      grepl(
          "^novel_gene",
          gene_id
      )]


############################################################
# Analysis input
############################################################

d <- dat[
    !is.na(gt_class) &
    !is.na(TPM) &
    !is_novel &
    !is.na(gene_id) &
    gene_id != ""
]


############################################################
# Remove duplicate sample observations
############################################################

setorder(
    d,
    sv_id,
    gene_id,
    sample
)


dup_n <- d[
    duplicated(
        d,
        by = c(
            "sv_id",
            "gene_id",
            "sample"
        )
    ),
    .N
]


cat(
    sprintf(
        "\nDuplicate SV-gene-sample rows removed: %d\n",
        dup_n
    )
)


d <- unique(
    d,
    by = c(
        "sv_id",
        "gene_id",
        "sample"
    )
)


cat(
    sprintf(
        "Rows entering association test: %d\n",
        nrow(d)
    )
)

cat(
    sprintf(
        "SV x gene combinations: %d\n",
        uniqueN(
            d[
                ,
                .(
                    sv_id,
                    gene_id
                )
            ]
        )
    )
)


############################################################
# Helper
############################################################

first_nonempty <- function(x) {

    x <- unique(
        x[
            !is.na(x) &
            nzchar(x)
        ]
    )

    if (length(x) == 0) {
        return("")
    }

    x[1]
}


############################################################
# Statistical test for one SV-gene pair
############################################################

stat_one <- function(
    TPM,
    gt_class,
    dosage
) {

    gt_factor <- factor(
        gt_class,
        levels = c(
            "AA",
            "Aa",
            "aa"
        )
    )


    n <- table(
        gt_factor
    )


    present <- n[
        n > 0
    ]


    ########################################################
    # Need at least two observed genotype classes
    ########################################################

    if (length(present) < 2) {
        return(NULL)
    }


    if (any(
        present < MIN_GROUP_N
    )) {
        return(NULL)
    }


    ########################################################
    # Kruskal-Wallis
    ########################################################

    p_kw <- tryCatch(

        kruskal.test(
            TPM,
            gt_factor
        )$p.value,

        error = function(e)
            NA_real_
    )


    ########################################################
    # Spearman dosage trend
    ########################################################

    sp <- NULL


    if (
        length(unique(dosage)) >= 2 &&
        length(unique(TPM)) >= 2
    ) {

        sp <- tryCatch(

            suppressWarnings(
                cor.test(
                    dosage,
                    TPM,
                    method = "spearman",
                    exact = FALSE
                )
            ),

            error = function(e)
                NULL
        )
    }


    ########################################################
    # Median TPM
    ########################################################

    median_group <- function(cl) {

        x <- TPM[
            gt_class == cl
        ]

        if (length(x) == 0) {
            return(NA_real_)
        }

        median(
            x,
            na.rm = TRUE
        )
    }


    list(

        n_total =
            length(TPM),

        n_groups =
            length(present),

        n_AA =
            as.integer(
                n[["AA"]]
            ),

        n_Aa =
            as.integer(
                n[["Aa"]]
            ),

        n_aa =
            as.integer(
                n[["aa"]]
            ),

        med_AA =
            median_group("AA"),

        med_Aa =
            median_group("Aa"),

        med_aa =
            median_group("aa"),

        p_kw =
            p_kw,

        rho =
            if (is.null(sp))
                NA_real_
            else
                unname(
                    sp$estimate
                ),

        p_trend =
            if (is.null(sp))
                NA_real_
            else
                sp$p.value
    )
}


############################################################
# Run tests
############################################################

cat("\nRunning SV-gene tests...\n")


res <- d[
    ,
    {

        s <- stat_one(
            TPM,
            gt_class,
            dosage_calc
        )

        if (is.null(s)) {
            NULL
        } else {
            s
        }

    },
    by = .(
        sv_id,
        gene_id
    )
]


############################################################
# Metadata
############################################################

meta <- d[
    ,
    .(
        gene_name =
            first_nonempty(
                gene_name
            ),

        svtype =
            first(
                svtype
            ),

        svlen =
            first(
                svlen
            )
    ),
    by = .(
        sv_id,
        gene_id
    )
]


res <- merge(
    res,
    meta,
    by = c(
        "sv_id",
        "gene_id"
    ),
    all.x = TRUE,
    sort = FALSE
)


############################################################
# Global BH correction, applied separately to KW and Spearman P values
############################################################

res[
    ,
    q_kw :=
        p.adjust(
            p_kw,
            method = "BH"
        )
]


res[
    ,
    q_trend :=
        p.adjust(
            p_trend,
            method = "BH"
        )
]


############################################################
# Sort
############################################################

setorder(
    res,
    q_kw,
    q_trend
)


############################################################
# Output all
############################################################

fwrite(
    res,
    file.path(
        OUTDIR,
        "eQTL_scan_all.tsv"
    ),
    sep = "\t"
)


############################################################
# Primary significant associations:
# KW + global BH
############################################################

sig_kw <- res[
    q_kw < 0.05
]


fwrite(
    sig_kw,
    file.path(
        OUTDIR,
        "eQTL_scan_sig_KW.tsv"
    ),
    sep = "\t"
)


############################################################
# Spearman trend-significant associations
############################################################

sig_trend <- res[
    q_trend < 0.05
]


fwrite(
    sig_trend,
    file.path(
        OUTDIR,
        "eQTL_scan_sig_trend.tsv"
    ),
    sep = "\t"
)


############################################################
# Union of KW- and trend-significant associations
############################################################

sig_any <- res[
    q_kw < 0.05 |
    q_trend < 0.05
]


fwrite(
    sig_any,
    file.path(
        OUTDIR,
        "eQTL_scan_sig_any.tsv"
    ),
    sep = "\t"
)


############################################################
# Summary
############################################################

cat("\n========================================\n")

cat(
    sprintf(
        "Tested SV-gene pairs: %d\n",
        nrow(res)
    )
)


cat(
    sprintf(
        "KW FDR < 0.05: %d associations, %d genes, %d SVs\n",
        nrow(sig_kw),
        uniqueN(sig_kw$gene_id),
        uniqueN(sig_kw$sv_id)
    )
)


cat(
    sprintf(
        "Trend FDR < 0.05: %d\n",
        nrow(sig_trend)
    )
)


cat(
    sprintf(
        "Either FDR < 0.05: %d\n",
        nrow(sig_any)
    )
)


cat("\nKW significant by SV type:\n")

print(
    sig_kw[
        ,
        .N,
        by = svtype
    ]
)


cat("========================================\n")


print(
    head(
        sig_kw,
        20
    )
)