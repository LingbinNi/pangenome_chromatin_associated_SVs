#!/usr/bin/env Rscript
# ==================================================================
# CDH (CHM13) vs SV (CHM13) 关联分析  —  Fisher 检验
#
# 检验单位 = 单倍型 (sample_Hap1 / sample_Hap2)
#   - CDH  : 每个单倍型一个 bed, 内容是 CHM13 坐标上的 CDH 区间
#            (含 .transloc.bed —— 易位也算作 CDH)
#   - SV   : 单倍型 carrier 表 (hap1_carrier / hap2_carrier), 只含杂合 SV
#
# 配对规则:
#   INS  同侧 —— 插入在哪条 hap, 就和哪条 hap 的 CDH 配对
#   DEL  交叉 —— 缺失侧少了片段, 保留序列的另一条 hap 才该有 CDH
#
# 流程:
#   1. 读所有单倍型的 CDH bed, 在 CHM13 上按距离阈值合并成 CDH 区域
#   2. 读 SV carrier 表, 按 svtype 决定有效携带单倍型, 建 SV × 单倍型 0/1 矩阵
#   3. CDH 区域 与 SV 位点取 overlap, 对每个 (CDH区域, SV) 组合,
#      在共同单倍型上做 2x2 Fisher 检验
#   4. 扫描不同合并距离阈值, 比较结果
# ==================================================================

suppressMessages({
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
})

################################################################
## 参数 / 路径  (按需修改)
################################################################

cdh_bed_dir <- "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/structure_associated_SV_pan/01.CDH_CHM13"
sv_file     <- "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/structure_associated_SV_identification/00.SV_CHM13/SV_CHM13_hap_carrier.tsv"

outdir <- "./"
dir.create(outdir, showWarnings = FALSE)

include_transloc <- TRUE                       # TRUE 则把 *.transloc.bed 也算作 CDH
flip_del         <- TRUE                       # TRUE 则 DEL 交叉配对 (INS 始终同侧)
dist_list        <- c(0, 5000, 10000, 25000, 50000)   # 合并距离阈值 (bp)

# 合并后区域的宽度上限 (bp)。Inf = 不限制。
# 加入 transloc 后极易出现跨越大半条染色体的超大合并区域,
# 一个区域吞掉成千上万个 SV, 会主导整张结果表 —— 需要时设成 1e6。
max_region_bp <- 100000

# 显著性筛选阈值
OR_cut    <- 2
padj_cut  <- 0.05
min_a_cut <- 3       # w_cdh_w_sv (同时有 CDH 有 SV 的单倍型数) 下限
spec_cut  <- 0.5     # keep2 额外要求 specificity

################################################################
## 1. 读 CDH bed  ->  单倍型 GRanges
################################################################
# 文件名: {sample}.Hap{1,2}.CDH.CHM13[.transloc].bed  (4 列: chr start end name)

pat <- if (include_transloc) "\\.CDH\\.CHM13(\\.transloc)?\\.bed$" else "\\.CDH\\.CHM13\\.bed$"
cdh_files <- list.files(cdh_bed_dir, pattern = pat, full.names = TRUE)
stopifnot(length(cdh_files) > 0)
message("CDH bed 文件数: ", length(cdh_files),
        "  (transloc: ", sum(grepl("\\.transloc\\.bed$", cdh_files)), ")")

read_cdh <- function(f) {
  if (file.info(f)$size == 0) return(NULL)          # transloc bed 可能为空
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

# bed 是 0-based 半开区间, 转 1-based (start + 1) 供 GRanges 使用
gr_cdh_all <- GRanges(
  seqnames  = cdh_dt$chr,
  ranges    = IRanges(start = cdh_dt$start + 1L, end = cdh_dt$end),
  haplotype = cdh_dt$hap
)
cdh_haps <- sort(unique(cdh_dt$hap))
message("CDH 单倍型数: ", length(cdh_haps), "   CDH 区间总数: ", length(gr_cdh_all))

################################################################
## 2. 读 SV carrier 表  ->  SV 信息 + 单倍型集合
################################################################
# 表头: sample chr pos sv_id svtype svlen hap1_carrier hap2_carrier

sv_dt <- fread(sv_file, header = TRUE)
sv_dt[, svtype := toupper(svtype)]

# 每个 SV 取一行坐标 (0-based pos -> 1-based)
sv_info <- sv_dt[, .(chr = chr[1L], pos = pos[1L], svtype = svtype[1L]), by = sv_id]
gr_sv <- GRanges(
  seqnames = sv_info$chr,
  ranges   = IRanges(start = sv_info$pos + 1L, end = sv_info$pos + 1L)
)

sv_haps <- sort(unique(c(paste0(sv_dt$sample, "_Hap1"),
                         paste0(sv_dt$sample, "_Hap2"))))
message("SV 单倍型数: ", length(sv_haps),
        "   SV 位点数: ", nrow(sv_info),
        "  (DEL=", sum(sv_info$svtype == "DEL"),
        "  INS=", sum(sv_info$svtype == "INS"), ")")

################################################################
## 3. 共同单倍型 (两边都有的才纳入检验)
################################################################
# !! 假设: CDH 的 Hap1/Hap2 与 SV 的 hap1/hap2 相位一致 (同一套 phased assembly)

common_haps <- intersect(cdh_haps, sv_haps)
message("共同单倍型数: ", length(common_haps))
stopifnot(length(common_haps) > 0)

# DEL 交叉配对要求样本的两条 hap 都在: 否则翻转会把携带者移到被丢弃的一侧
samp_of <- sub("_Hap[12]$", "", common_haps)
both_ok <- names(which(table(samp_of) == 2L))
n_drop  <- length(common_haps) - 2L * length(both_ok)
if (n_drop > 0)
  message("注意: ", n_drop, " 条单倍型因同一样本另一条缺失而被剔除 (DEL 交叉需要成对)")
common_haps <- sort(common_haps[samp_of %in% both_ok])
message("成对后的共同单倍型数: ", length(common_haps),
        "  (样本数 ", length(both_ok), ")")

# CDH 区域定义与检验都限定在共同单倍型上
gr_cdh_all <- gr_cdh_all[gr_cdh_all$haplotype %in% common_haps]

################################################################
## 4. 预建 SV 存在矩阵 (只针对可能与 CDH overlap 的 SV, 用最大阈值取超集)
################################################################

merged_max <- reduce(gr_cdh_all, min.gapwidth = max(dist_list) + 1L)
if (is.finite(max_region_bp)) merged_max <- merged_max[width(merged_max) <= max_region_bp]
sv_super   <- unique(sv_info$sv_id[subjectHits(findOverlaps(merged_max, gr_sv))])
message("与 CDH overlap 的 SV 数 (最大阈值): ", length(sv_super))
stopifnot(length(sv_super) > 0)

sub <- sv_dt[sv_id %chin% sv_super]

# ---- 有效携带单倍型: INS 同侧, DEL 交叉 ----
if (flip_del) {
  sub[, `:=`(eff1 = fifelse(svtype == "DEL", hap2_carrier, hap1_carrier),
             eff2 = fifelse(svtype == "DEL", hap1_carrier, hap2_carrier))]
} else {
  sub[, `:=`(eff1 = hap1_carrier, eff2 = hap2_carrier)]
}

car_long <- rbindlist(list(
  sub[eff1 == 1, .(sv_id, hap = paste0(sample, "_Hap1"))],
  sub[eff2 == 1, .(sv_id, hap = paste0(sample, "_Hap2"))]
))
car_long <- car_long[hap %chin% common_haps]

SV_M <- matrix(0L, nrow = length(sv_super), ncol = length(common_haps),
               dimnames = list(sv_super, common_haps))
SV_M[cbind(match(car_long$sv_id, sv_super),
           match(car_long$hap,   common_haps))] <- 1L
sv_rowsum <- rowSums(SV_M)

################################################################
## 5. 函数: 合并区域的 单倍型存在矩阵
################################################################

build_cdh_matrix <- function(merged, gr_all, haps) {
  M <- matrix(0L, nrow = length(merged), ncol = length(haps),
              dimnames = list(NULL, haps))
  hv <- gr_all$haplotype
  for (j in seq_along(haps)) {
    gr_h <- gr_all[hv == haps[j]]
    if (length(gr_h) == 0) next
    M[unique(queryHits(findOverlaps(merged, gr_h))), j] <- 1L
  }
  M
}

fisher_p <- function(a, b, c, d) fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value

################################################################
## 6. 扫描距离阈值, Fisher 检验
################################################################

n <- length(common_haps)
summary_list <- list()

for (dist in dist_list) {
  tag <- paste0("d", dist)
  message("== merge distance = ", dist, " bp ==")

  merged <- reduce(gr_cdh_all, min.gapwidth = dist + 1L)   # dist+1 等价 bedtools merge -d dist
  if (is.finite(max_region_bp)) {
    nbig <- sum(width(merged) > max_region_bp)
    if (nbig > 0) message("  剔除超大区域: ", nbig, " 个 (> ", max_region_bp, " bp)")
    merged <- merged[width(merged) <= max_region_bp]
  }
  if (length(merged) == 0) { message("  no region left, skip"); next }
  message("  区域数: ", length(merged),
          "   中位宽度: ", median(width(merged)),
          "   最大宽度: ", max(width(merged)))

  CDH_M  <- build_cdh_matrix(merged, gr_cdh_all, common_haps)
  cdh_rowsum <- rowSums(CDH_M)

  ov <- findOverlaps(merged, gr_sv)
  if (length(ov) == 0) { message("  no CDH-SV overlap, skip"); next }

  reg_idx <- queryHits(ov)
  sv_idx  <- subjectHits(ov)
  sv_ids  <- sv_info$sv_id[sv_idx]

  keep_pair <- sv_ids %chin% rownames(SV_M)
  reg_idx <- reg_idx[keep_pair]; sv_idx <- sv_idx[keep_pair]; sv_ids <- sv_ids[keep_pair]
  npair <- length(sv_ids)
  if (npair == 0) { message("  no valid pair, skip"); next }
  message("  CDH-SV 组合数: ", npair)

  # 2x2 计数: a=有CDH有SV, b=有CDH无SV, c=无CDH有SV, d=无CDH无SV
  # (分块向量化, 比逐对循环快一到两个数量级)
  a <- integer(npair)
  chunk <- 20000L
  for (i0 in seq(1L, npair, by = chunk)) {
    ix <- i0:min(i0 + chunk - 1L, npair)
    a[ix] <- rowSums(CDH_M[reg_idx[ix], , drop = FALSE] &
                     SV_M[sv_ids[ix],  , drop = FALSE])
  }
  b <- cdh_rowsum[reg_idx] - a
  c <- sv_rowsum[sv_ids]  - a
  d <- n - a - b - c

  result <- data.frame(
    chr_cdh   = as.character(seqnames(merged))[reg_idx],
    start_cdh = start(merged)[reg_idx],
    end_cdh   = end(merged)[reg_idx],
    region_bp = width(merged)[reg_idx],
    sv_id     = sv_ids,
    chr_sv    = sv_info$chr[sv_idx],
    pos_sv    = sv_info$pos[sv_idx],
    svtype    = sv_info$svtype[sv_idx],
    w_cdh_w_sv = a, w_cdh_wo_sv = b, wo_cdh_w_sv = c, wo_cdh_wo_sv = d,
    stringsAsFactors = FALSE
  )
  result$cdh_count   <- a + b
  result$sv_count    <- a + c
  result$OR          <- (a + 1) * (d + 1) / ((b + 1) * (c + 1))   # +1 平滑
  result$p           <- mapply(fisher_p, a, b, c, d)
  result$p_adj       <- p.adjust(result$p, method = "BH")
  result$specificity <- (a + 1) / (a + c + 2)   # SV 携带者里同时有 CDH 的比例

  write.table(result, file.path(outdir, paste0("cdh_sv.", tag, ".tsv")),
              row.names = FALSE, quote = FALSE, sep = "\t")

  keep  <- subset(result, OR > OR_cut & p_adj < padj_cut & w_cdh_w_sv >= min_a_cut)
  keep2 <- subset(keep, specificity > spec_cut)
  write.table(keep,  file.path(outdir, paste0("cdh_sv.", tag, ".keep.tsv")),
              row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(keep2, file.path(outdir, paste0("cdh_sv.", tag, ".keep2.tsv")),
              row.names = FALSE, quote = FALSE, sep = "\t")

  if (nrow(keep) > 0) {
    p <- ggplot(keep, aes(cdh_count, sv_count)) +
      geom_point(aes(size = log2(OR), color = -log10(p_adj)), alpha = 0.8) +
      scale_color_gradient(low = "lightblue", high = "red") +
      scale_size(range = c(2, 10)) +
      facet_wrap(~svtype) +
      theme_classic(base_size = 14) +
      labs(x = "haplotypes with CDH", y = "haplotypes with SV",
           color = "-log10(FDR)", size = "log2(OR)",
           title = paste0("CDH-SV enrichment (merge d = ", dist, ")"))
    ggsave(file.path(outdir, paste0("enrichment.", tag, ".pdf")), p, width = 11, height = 6)
  }

  summary_list[[tag]] <- data.frame(
    dist             = dist,
    n_cdh_regions    = length(merged),
    median_region_bp = as.numeric(median(width(merged))),
    max_region_bp    = as.numeric(max(width(merged))),
    n_tested_pairs   = npair,
    n_sig_keep       = nrow(keep),
    n_sig_keep2      = nrow(keep2),
    keep_DEL         = sum(keep$svtype  == "DEL"),
    keep_INS         = sum(keep$svtype  == "INS"),
    keep2_DEL        = sum(keep2$svtype == "DEL"),
    keep2_INS        = sum(keep2$svtype == "INS")
  )
}

################################################################
## 7. 阈值比较汇总
################################################################

summary_df <- do.call(rbind, summary_list)
write.table(summary_df, file.path(outdir, "threshold_summary.tsv"),
            row.names = FALSE, quote = FALSE, sep = "\t")
print(summary_df)

long <- do.call(rbind, lapply(
  c("n_cdh_regions", "n_tested_pairs", "n_sig_keep", "n_sig_keep2"),
  function(v) data.frame(dist = summary_df$dist, metric = v, value = summary_df[[v]])
))
pdf(file.path(outdir, "threshold_summary.pdf"), width = 8, height = 6)
print(
  ggplot(long, aes(dist, value, color = metric)) +
    geom_line() + geom_point(size = 2) +
    theme_classic(base_size = 14) +
    labs(x = "merge distance (bp)", y = "count", title = "Effect of CDH merge distance")
)
dev.off()

message("Done. 输出目录: ", normalizePath(outdir))
