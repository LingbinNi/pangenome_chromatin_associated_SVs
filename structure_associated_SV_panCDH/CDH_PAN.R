#!/usr/bin/env Rscript
# ==================================================================
# CDH 区域合并  —  从 cdh_sv_fisher_merge.R 抽出的合并部分
#
# 逻辑与原脚本完全一致：
#   读 CDH bed → 限定共同单倍型(且样本成对) → 按各距离阈值 reduce
#   合并 → 剔除超大区域 → 输出
# 不做 SV 配对矩阵、不做 Fisher、不出富集图
#
# 额外输出：至少 min_hap_support 条单倍型支持的区域 (原全量输出保持不变)
# ==================================================================
suppressMessages({
  library(GenomicRanges)
  library(data.table)
})

################################################################
# 参数 / 路径
################################################################
cdh_bed_dir <- "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/structure_associated_SV_pan/01.CDH_CHM13"
sv_file     <- "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/structure_associated_SV_identification/00.SV_CHM13/SV_CHM13_hap_carrier.tsv"

outdir <- "./"
dir.create(outdir, showWarnings = FALSE)

include_transloc <- TRUE                              # TRUE 则 *.transloc.bed 也算 CDH
dist_list        <- c(0, 5000, 10000, 25000, 50000)   # 合并距离阈值 (bp)
max_region_bp    <- 1000000                           # 区域宽度上限, Inf = 不限制

# TRUE  = 与原脚本一致: 合并前先限定到「与 SV 表共有、且样本两条 hap 都在」的单倍型
# FALSE = 不读 SV 表, 直接合并全部 CDH 单倍型 (合并结果会与原脚本不同)
restrict_to_sv_haps <- TRUE

# 每个合并区域统计有多少条单倍型在其中有 CDH (即原脚本 build_cdh_matrix 的 rowSums)
count_haps <- TRUE

# 至少多少条单倍型支持才输出到 recurrent 文件 (需要 count_haps = TRUE)
min_hap_support <- 20

################################################################
# 1. 读 CDH bed  ->  单倍型 GRanges
################################################################
pat <- if (include_transloc) "\\.CDH\\.CHM13(\\.transloc)?\\.bed$" else "\\.CDH\\.CHM13\\.bed$"
cdh_files <- list.files(cdh_bed_dir, pattern = pat, full.names = TRUE)
stopifnot(length(cdh_files) > 0)
message("CDH bed 文件数: ", length(cdh_files),
        "  (transloc: ", sum(grepl("\\.transloc\\.bed$", cdh_files)), ")")

read_cdh <- function(f) {
  if (file.info(f)$size == 0) return(NULL)            # transloc bed 可能为空
  d <- tryCatch(fread(f, header = FALSE, select = 1:3,
                      col.names = c("chr", "start", "end")),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  hap <- sub("^(.+)\\.(Hap[12])\\.CDH\\.CHM13(\\.transloc)?\\.bed$", "\\1_\\2", basename(f))
  d[, hap := hap]
  d
}

cdh_dt <- rbindlist(lapply(cdh_files, read_cdh))
stopifnot(nrow(cdh_dt) > 0)

# bed 是 0-based 半开区间, 转 1-based 供 GRanges 使用
gr_cdh_all <- GRanges(
  seqnames  = cdh_dt$chr,
  ranges    = IRanges(start = cdh_dt$start + 1L, end = cdh_dt$end),
  haplotype = cdh_dt$hap
)
cdh_haps <- sort(unique(cdh_dt$hap))
message("CDH 单倍型数: ", length(cdh_haps), "   CDH 区间总数: ", length(gr_cdh_all))

################################################################
# 2. 限定共同单倍型 (与原脚本第 2-3 步相同)
################################################################
if (restrict_to_sv_haps) {
  sv_dt <- fread(sv_file, header = TRUE)
  sv_haps <- sort(unique(c(paste0(sv_dt$sample, "_Hap1"),
                           paste0(sv_dt$sample, "_Hap2"))))
  message("SV 单倍型数: ", length(sv_haps))

  common_haps <- intersect(cdh_haps, sv_haps)
  message("共同单倍型数: ", length(common_haps))
  stopifnot(length(common_haps) > 0)

  # 要求样本两条 hap 都在 (与原脚本一致, 保证区域定义域相同)
  samp_of <- sub("_Hap[12]$", "", common_haps)
  both_ok <- names(which(table(samp_of) == 2L))
  n_drop  <- length(common_haps) - 2L * length(both_ok)
  if (n_drop > 0)
    message("注意: ", n_drop, " 条单倍型因同一样本另一条缺失而被剔除")
  common_haps <- sort(common_haps[samp_of %in% both_ok])
  message("成对后的共同单倍型数: ", length(common_haps),
          "  (样本数 ", length(both_ok), ")")

  gr_cdh_all <- gr_cdh_all[gr_cdh_all$haplotype %in% common_haps]
} else {
  common_haps <- cdh_haps
  message("跳过 SV 单倍型限定, 直接使用全部 ", length(common_haps), " 条单倍型")
}
message("参与合并的 CDH 区间数: ", length(gr_cdh_all))

################################################################
# 3. 每个区域的单倍型计数 (原 build_cdh_matrix 的 rowSums)
################################################################
count_hap_per_region <- function(merged, gr_all, haps) {
  cnt <- integer(length(merged))
  hv  <- gr_all$haplotype
  for (h in haps) {
    gr_h <- gr_all[hv == h]
    if (length(gr_h) == 0) next
    idx <- unique(queryHits(findOverlaps(merged, gr_h)))
    cnt[idx] <- cnt[idx] + 1L
  }
  cnt
}

################################################################
# 4. 扫描距离阈值, 合并并输出
################################################################
summary_list <- list()

for (dist in dist_list) {
  tag <- paste0("d", dist)
  message("== merge distance = ", dist, " bp ==")

  merged <- reduce(gr_cdh_all, min.gapwidth = dist + 1L)   # dist+1 等价 bedtools merge -d dist
  n_before <- length(merged)

  nbig <- 0L
  if (is.finite(max_region_bp)) {
    nbig <- sum(width(merged) > max_region_bp)
    if (nbig > 0) message("  剔除超大区域: ", nbig, " 个 (> ", max_region_bp, " bp)")
    merged <- merged[width(merged) <= max_region_bp]
  }
  if (length(merged) == 0) { message("  no region left, skip"); next }

  message("  区域数: ", length(merged),
          "   中位宽度: ", median(width(merged)),
          "   最大宽度: ", max(width(merged)))

  out <- data.frame(
    chr       = as.character(seqnames(merged)),
    start     = start(merged) - 1L,            # 输出回 0-based, 与输入 bed 一致
    end       = end(merged),
    region_id = paste0(tag, "_", seq_along(merged)),
    region_bp = width(merged),
    stringsAsFactors = FALSE
  )

  n_hap_med <- NA_real_
  if (count_haps) {
    out$n_hap <- count_hap_per_region(merged, gr_cdh_all, common_haps)
    n_hap_med <- median(out$n_hap)
    message("  区域内单倍型数中位: ", n_hap_med)
  }

  write.table(out, file.path(outdir, paste0("cdh_merged.", tag, ".bed")),
              row.names = FALSE, quote = FALSE, sep = "\t")

  # ---- 至少 min_hap_support 条单倍型支持的区域 (额外输出, 不影响上面的全量文件) ----
  n_keep <- NA_integer_; keep_bp <- NA_real_; keep_med_bp <- NA_real_
  if (count_haps) {
    keep   <- out[out$n_hap >= min_hap_support, ]
    n_keep <- nrow(keep)
    message("  >=", min_hap_support, " 单倍型支持的区域: ", n_keep,
            " / ", nrow(out), " (", sprintf("%.1f%%", 100 * n_keep / nrow(out)), ")")
    if (n_keep > 0) {
      keep_bp     <- sum(as.numeric(keep$region_bp))
      keep_med_bp <- median(keep$region_bp)
      message("    中位宽度: ", keep_med_bp,
              "   单倍型数中位: ", median(keep$n_hap),
              "   最大: ", max(keep$n_hap))
    }
    write.table(keep,
                file.path(outdir, paste0("cdh_merged.", tag, ".minhap", min_hap_support, ".bed")),
                row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    message("  count_haps = FALSE, 跳过单倍型支持筛选")
  }

  summary_list[[tag]] <- data.frame(
    dist                    = dist,
    n_regions_raw           = n_before,
    n_dropped_oversize      = nbig,
    n_regions               = length(merged),
    median_region_bp        = as.numeric(median(width(merged))),
    max_region_bp           = as.numeric(max(width(merged))),
    total_bp                = sum(as.numeric(width(merged))),
    median_n_hap            = n_hap_med,
    min_hap_support         = min_hap_support,
    n_regions_minhap        = n_keep,
    median_region_bp_minhap = keep_med_bp,
    total_bp_minhap         = keep_bp
  )
}

################################################################
# 5. 阈值比较汇总
################################################################
summary_df <- do.call(rbind, summary_list)
write.table(summary_df, file.path(outdir, "merge_summary.tsv"),
            row.names = FALSE, quote = FALSE, sep = "\t")
print(summary_df, row.names = FALSE)

message("Done. 输出目录: ", normalizePath(outdir))