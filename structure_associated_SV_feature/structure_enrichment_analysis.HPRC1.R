#!/usr/bin/env Rscript
# ==================================================================
# CSA-SV 与染色质结构特征 (TAD boundary / loop anchor) 的邻近性
#
# 坐标桥梁:
#   CSA-SV 定义在 CHM13 (variant_id = chr-pos-TYPE-len)
#   TAD/loop 定义在每条单倍型自己的 RagTag 组装
#   -> ContactCount 文件同时含 RagTag 坐标与 CHM13 variant_id, 用它做映射
#
# 文件 -> 单倍型 (接触/坐标始终在"保留序列"一侧):
#   Hap1 的结构  <-  hap1.INS  +  hap2.DEL
#   Hap2 的结构  <-  hap2.INS  +  hap1.DEL
#
# 两个对照 (只报 30% 而不给背景是没有意义的):
#   (1) nonCSA : 落在 CDH 内但未入选 CSA 的 SV —— 直接可比
#   (2) 置换   : SV 位置在同一染色体内随机重投, 给出富集倍数
# ==================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
})

############################################################
## Parameters
############################################################

set.seed(1)
nperm   <- 200          # 置换次数 (每 sample x hap x feature 一次 distanceToNearest)
do_perm <- TRUE

outdir <- "CSA_structure_proximity"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## 距离阈值 (bp)。正文用的是 25000,其余为敏感性分析
dist_list <- c(0, 5000, 10000, 25000, 50000, 100000)
MAIN_DIST <- 25000

contact_dir <- "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/structure_associated_SV_feature/contact_enrichment_analysis/02.ContactSum"

fisher_dir <- "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/structure_associated_SV_identification"
TAG        <- "d0"
csa_level  <- "keep"          # "keep2" 或 "keep"

tad_roots <- list(
  Hap1 = "/net/eichler/vol28/projects/hic_cohorts/nobackups/11.StructureIdentification-Hap1-HPRC-RagTag-Y2-HPRC1/Tad",
  Hap2 = "/net/eichler/vol28/projects/hic_cohorts/nobackups/12.StructureIdentification-Hap2-HPRC-RagTag-Y2-HPRC1/Tad"
)
loop_roots <- list(
  Hap1 = "/net/eichler/vol28/projects/hic_cohorts/nobackups/11.StructureIdentification-Hap1-HPRC-RagTag-Y2-HPRC1/Loop",
  Hap2 = "/net/eichler/vol28/projects/hic_cohorts/nobackups/12.StructureIdentification-Hap2-HPRC-RagTag-Y2-HPRC1/Loop"
)

############################################################
## Readers  (TAD / loop 沿用参考脚本)
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

## 某条单倍型上的 SV: RagTag 坐标 + CHM13 variant_id
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

## 到最近特征的距离 (无匹配染色体的记 Inf)
nearest_dist <- function(q, s) {
  if (length(q) == 0 || length(s) == 0) return(rep(Inf, length(q)))
  h <- distanceToNearest(q, s, ignore.strand = TRUE)
  d <- rep(Inf, length(q))
  d[queryHits(h)] <- mcols(h)$distance
  d
}

############################################################
## 1. CSA / nonCSA 集合
############################################################

f_all <- file.path(fisher_dir, sprintf("cdh_sv.%s.tsv", TAG))
f_csa <- file.path(fisher_dir, sprintf("cdh_sv.%s.%s.tsv", TAG, csa_level))
stopifnot(file.exists(f_all), file.exists(f_csa))
tested_ids <- unique(fread(f_all, select = "sv_id")$sv_id)
csa_ids    <- unique(fread(f_csa, select = "sv_id")$sv_id)
message("CDH 内被检验的 SV: ", length(tested_ids), "   CSA-SV: ", length(csa_ids))

############################################################
## 2. 逐 sample x hap
############################################################

samples <- sort(unique(sub("\\.insdel\\.het\\.hap[12]\\.(INS|DEL)\\.sum\\.ContactCount$", "",
                           list.files(contact_dir, pattern = "\\.sum\\.ContactCount$"))))
message("样本数: ", length(samples))

res_rows <- list()   # 逐 sample x hap x feature x dist x group
sv_rows  <- list()   # 逐 SV (主阈值, Any_structure)

for (S in samples) {
  for (HAP in c("Hap1", "Hap2")) {

    sv <- read_sv_on_hap(S, HAP)
    if (is.null(sv)) next
    sv <- sv[variant_id %chin% tested_ids]                 # 只看落在 CDH 内的
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

      # 置换期望 (两组共用同一批随机位置的距离分布, 按组抽取相应数量)
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
## 3. 汇总 A: 逐单倍型的比例, 再取样本间中位数
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
## 4. 汇总 B: 以 SV 为单位 (正文那个 "30% of CSA-SVs" 的分母)
############################################################
## 每个 CSA-SV 在其携带单倍型中有多少比例是"邻近结构的",
## 再按共识阈值判定该 SV 是否算 near。

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

cat("\n===== A. 逐单倍型比例 (feature x dist x group) =====\n")
print(sumA[dist == MAIN_DIST])
cat("\n===== B. 以 SV 为单位, Any_structure, dist <=", MAIN_DIST, "bp =====\n")
print(sumB)
cat("\npct_near_majority = 在过半数携带单倍型中邻近结构的 SV 占比 (建议用这个报告)\n")

############################################################
## 5. 图
############################################################

## 图1: 比例 vs 距离阈值, CSA vs nonCSA vs 置换期望
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
  labs(x = "distance threshold (kb)", y = "fraction of SVs within threshold",
       colour = "", title = "SV proximity to chromatin architectural features")
ggsave(file.path(outdir, "proximity_vs_distance.pdf"), p1, width = 10, height = 4)

## 图2: 主阈值下的富集倍数 (逐单倍型)
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

## 图3: 到最近结构的距离分布 (ECDF)
p3 <- ggplot(svd[is.finite(dist)], aes(pmax(dist, 1), colour = group)) +
  stat_ecdf(linewidth = 0.8) +
  geom_vline(xintercept = MAIN_DIST, linetype = "dotted", colour = "grey40") +
  scale_x_log10() +
  theme_classic(base_size = 12) +
  labs(x = "distance to nearest architectural feature (bp)",
       y = "cumulative fraction", colour = "",
       title = "Distance to nearest TAD boundary or loop anchor")
ggsave(file.path(outdir, "distance_ecdf.pdf"), p3, width = 6.5, height = 4)

message("Done. 输出目录: ", normalizePath(outdir))
