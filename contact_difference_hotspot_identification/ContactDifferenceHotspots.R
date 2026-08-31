#!/usr/bin/env Rscript

# =============================================================================
# CDH PRODUCTION CALLER  (final, locked method)
#   Model = M2 : NB with COVERAGE-DEPENDENT DISPERSION
#           spec ~ ns(log1p(shared), df) + Chr           (mean model)
#           dispersion sigma = 1/theta varies with coverage (NBI),
#           estimated per coverage stratum (robust) + interpolated per window.
#   Output = per sample x haplotype CDH lists at BOTH window and region level.
#
# 为什么是这个模型(方法定稿的完整依据):
#   - 密度法(绝对 spec 计数)而非比例法:大 SV 使 spec/total 饱和,比例法失明;
#     覆盖度用回归(shared 协变量)处理,不用除法。比例法 Beta-Binomial cov_ratio=0.20
#     (热点全挤低覆盖),M2 cov_ratio~1.06(无覆盖度混杂)。
#   - 覆盖度依赖离散修正了单一 theta 的唯一真短板(theta_chr 0.3-5.7):
#     lambda 0.74->0.95,并按覆盖度重新校准(高覆盖去假阳性、低覆盖找回真信号)。
#   - 三方对决裁决:M2 用约一半数量拿到 >= M1 的富集(SV 窗口 5.1->6.0、区域 8.9->9.0;
#     3D TAD 1.24->1.54、loop 1.26->1.70),证明砍掉的是假阳性,M2 更干净。
#
# 判定:CDH = FDR(BH) q < Q_CUTOFF  AND  OE(obs/exp) >= EFFECT  AND  spec >= MIN_SPEC。
# 区域:相邻(间隔<=REGION_GAP)CDH 窗口合并为一个热点区域(生物学上诚实的计数)。
#
# 参数固定(依据前面的稳健性/敏感性分析):
#   数据 25kb + filter25kb;MIN_TOTAL=50;MIN_SPEC=10;Q=0.01;EFFECT(OE)=2;
#   NS_DF=4(df 2-8 对结果无实质影响);离散分层 N_STRATA=20。
#
# 注意(诚实声明):本 caller 保证热点在群体层面显著富集真生物学信号(SV/3D),
#   但不保证每个热点都是真信号,也未针对比对/组装假象设计。若需排除比对噪音,
#   应另做"比对易错区(segdup/低mappability/gap)富集检验"。
# =============================================================================

suppressPackageStartupMessages({ library(MASS); library(splines) })

## ===========================================================================
## Paths
## ===========================================================================
spec_dir  <- "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/contact_difference_hotspot_identification/03.ContactPosition-Specific"
total_dir <- "/net/eichler/vol28/projects/hgsvc/nobackups/hic_hgsvc/00.Analysis/contact_difference_hotspot_identification/03.ContactPosition-Total"
anno_base <- "/net/eichler/vol28/projects/hic_cohorts/nobackups/00.DataBase/centromere"
OUTDIR    <- "./"

## ===========================================================================
## Locked parameters
## ===========================================================================
THR        <- "25kb"       # 距离过滤档
Q_CUTOFF   <- 0.01
EFFECT     <- 2            # OE = obs/exp 门槛
MIN_TOTAL  <- 50
MIN_SPEC   <- 10
NS_DF      <- 4
N_REFIT    <- 2
TRIM_UPPER <- 0.995
FIT_MAX    <- 40000
N_STRATA   <- 20          # 覆盖度依赖离散的分层数
REGION_GAP <- 0           # 区域合并间隔(0 = 仅相接窗口)
WIN        <- 25000
VALID_CHR_PATTERN <- "^chr([1-9]|1[0-9]|2[0-2])_RagTag$"   # 常染色体(排除性染色体)

N_SUBSET   <- Inf         # 全样本;测试可设 2
set.seed(1)

## ===========================================================================
## I/O helpers
## ===========================================================================
valid_chr <- function(x) grepl(VALID_CHR_PATTERN, x)
get_acro <- function(sample, hap) {
  cand <- if (hap=="Hap1") c(file.path(anno_base,paste0(sample,".pat.acro.bed")), file.path(anno_base,paste0(sample,".hap1.acro.bed")))
          else            c(file.path(anno_base,paste0(sample,".mat.acro.bed")), file.path(anno_base,paste0(sample,".hap2.acro.bed")))
  for (f in cand) if (file.exists(f)) return(f); NULL
}
read_pos <- function(path, nm) {
  d <- read.table(path, header=FALSE, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)[,1:4]
  colnames(d) <- c("Chr","Start","End",nm); d$Start<-as.numeric(d$Start); d$End<-as.numeric(d$End); d[[nm]]<-as.numeric(d[[nm]]); d
}
remove_acro <- function(d, af) {
  if (is.null(af)) { message("    (no acro BED for masking)"); return(d) }
  a <- read.table(af, header=FALSE, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)[,1:3]; colnames(a) <- c("Chr","Start","End")
  drop <- rep(FALSE, nrow(d)); for (i in seq_len(nrow(a))) drop <- drop | (d$Chr==a$Chr[i] & d$Start<a$End[i] & d$End>a$Start[i]); d[!drop,]
}
build_table <- function(sample, hap) {
  sf <- sprintf("%s/%s/%s.SelfMapping.%s.HPRC.25kb.Unique.position.filter%s.position", spec_dir, sample, sample, hap, THR)
  tf <- sprintf("%s/%s/%s.SelfMapping.%s.HPRC.25kb.Total.position.filter%s.position",  total_dir, sample, sample, hap, THR)
  if (!file.exists(sf)||!file.exists(tf)) stop("missing spec/total file")
  spec <- read_pos(sf,"spec"); total <- read_pos(tf,"total")
  d <- merge(total, spec, by=c("Chr","Start","End"), all.x=TRUE, sort=FALSE); d$spec[is.na(d$spec)] <- 0
  bad <- which(d$spec>d$total); if (length(bad)) { message(sprintf("    note: %d windows spec>total dropped", length(bad))); d <- d[-bad,] }
  d$shared <- d$total - d$spec
  d <- d[valid_chr(d$Chr) & (d$End-d$Start==WIN),]; d <- remove_acro(d, get_acro(sample,hap))
  if (nrow(d) < 1000) stop(sprintf("too few autosomal windows: %d", nrow(d))); d
}
lambda_gc <- function(pv){ p <- pv[is.finite(pv)&pv>0&pv<=1]; if (length(p)<100) return(NA_real_)
  qchisq(median(p),1,lower.tail=FALSE)/qchisq(0.5,1,lower.tail=FALSE) }
strat_sample <- function(chr,n){ ib<-split(seq_along(chr),chr); nt<-length(chr)
  unlist(lapply(ib,function(ix){k<-min(length(ix),max(1L,round(n*length(ix)/nt))); if(length(ix)<=k) ix else sample(ix,k)}),use.names=FALSE) }

## ===========================================================================
## Final caller (M2): mean model + coverage-dependent dispersion
## ===========================================================================
call_cdh <- function(d) {
  dm <- d[d$total>=MIN_TOTAL,]; dm$Chr <- droplevels(factor(dm$Chr))
  form <- if (nlevels(dm$Chr)>1) spec ~ ns(log1p(shared),df=NS_DF)+Chr else spec ~ ns(log1p(shared),df=NS_DF)

  # ---- mean model: robust NB fit with iterative upper-tail trimming ----
  fi <- if (nrow(dm)>FIT_MAX) strat_sample(dm$Chr,FIT_MAX) else seq_len(nrow(dm)); fd <- dm[fi,]; keep <- rep(TRUE,nrow(fd)); fit <- NULL
  for (it in seq_len(N_REFIT)) {
    fit <- tryCatch(MASS::glm.nb(form,data=fd[keep,],link=log,control=glm.control(maxit=100,epsilon=1e-8)),error=function(e) NULL)
    if (is.null(fit)) { form <- if (nlevels(dm$Chr)>1) spec~log1p(shared)+Chr else spec~log1p(shared); fit <- MASS::glm.nb(form,data=fd[keep,],link=log) }
    muf<-predict(fit,newdata=fd,type="response"); pe<-(fd$spec-muf)/sqrt(pmax(muf+muf^2/fit$theta,1e-12))
    ct<-as.numeric(quantile(pe,TRIM_UPPER,na.rm=TRUE,type=8)); nk<-is.finite(pe)&pe<=ct; if (identical(nk,keep)) break; keep<-nk
  }
  mu <- predict(fit,newdata=dm,type="response"); theta_global <- fit$theta

  # ---- coverage-dependent dispersion: sigma(coverage) via robust strata + interpolation ----
  r  <- ((dm$spec - mu)^2 - mu) / mu^2                       # per-window dispersion contribution; E[r]=sigma=1/theta (NBI)
  br <- unique(quantile(mu, seq(0,1,length.out=N_STRATA+1), na.rm=TRUE)); br[1]<- -Inf; br[length(br)]<-Inf
  b  <- cut(mu, br, include.lowest=TRUE)
  centers <- tapply(mu, b, median)
  sig_s   <- tapply(seq_along(r), b, function(ix){ v <- r[ix]; v <- v[is.finite(v)]
    v <- v[v <= quantile(v, TRIM_UPPER, na.rm=TRUE)]; max(mean(v, na.rm=TRUE), 1e-6) })   # trim hotspots -> robust
  ok <- is.finite(centers) & is.finite(sig_s)
  sig_i <- pmax(approx(centers[ok], sig_s[ok], xout=mu, rule=2)$y, 1e-6)                   # smooth per-window sigma
  theta_i <- 1/sig_i

  # ---- one-sided NB tail p (per-window theta) + BH-FDR ----
  pv <- pnbinom(dm$spec-1, size=theta_i, mu=mu, lower.tail=FALSE); pv[!is.finite(pv)] <- 1; pv <- pmax(pmin(pv,1),0)
  dm$expected <- mu; dm$theta_local <- round(theta_i,3)
  dm$OE <- ifelse(mu>0, dm$spec/mu, Inf); dm$pval <- pv; dm$qval <- p.adjust(pv,"BH")
  dm$is_cdh <- dm$qval<Q_CUTOFF & dm$OE>=EFFECT & dm$spec>=MIN_SPEC

  list(d=dm, theta_global=theta_global, lambda=lambda_gc(pv), n_tested=nrow(dm),
       med_total=median(dm$total), mu_baseline=sum(dm$spec)/sum(dm$total))
}

# merge adjacent CDH windows (gap<=REGION_GAP) into regions; aggregate stats
merge_regions <- function(cd) {
  if (nrow(cd)==0) return(data.frame())
  x <- cd[order(cd$Chr, cd$Start),]; rows <- list()
  Chr<-x$Chr[1]; S<-x$Start[1]; E<-x$End[1]; idx<-1L
  flush <- function(sub) data.frame(Chr=sub$Chr[1], Start=min(sub$Start), End=max(sub$End),
      n_windows=nrow(sub), sum_spec=sum(sub$spec), sum_total=sum(sub$total),
      max_OE=round(max(sub$OE),2), min_qval=signif(min(sub$qval),3),
      peak_Start=sub$Start[which.max(sub$OE)])
  cur <- x[1,,drop=FALSE]
  if (nrow(x)>=2) for (i in 2:nrow(x)) {
    if (x$Chr[i]==Chr && x$Start[i] <= E+REGION_GAP) { E <- max(E, x$End[i]); cur <- rbind(cur, x[i,]) }
    else { rows[[length(rows)+1]] <- flush(cur); Chr<-x$Chr[i]; S<-x$Start[i]; E<-x$End[i]; cur <- x[i,,drop=FALSE] }
  }
  rows[[length(rows)+1]] <- flush(cur)
  reg <- do.call(rbind, rows); reg$region_frac <- round(reg$sum_spec/reg$sum_total,4); reg
}

## ===========================================================================
## Main
## ===========================================================================
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
samplelist <- list.files(spec_dir, pattern="^(HG|NA)")
samplelist <- samplelist[file.info(file.path(spec_dir, samplelist))$isdir %in% TRUE]
if (is.finite(N_SUBSET)) samplelist <- head(samplelist, N_SUBSET)
if (length(samplelist)==0) stop("no samples found")
cat(sprintf("Calling CDH on %d samples (25kb, filter%s, M2 coverage-dependent dispersion)\n\n", length(samplelist), THR))

summ <- data.frame()
for (sample in samplelist) {
  for (hap in c("Hap1","Hap2")) {
    message("Processing ", sample, " ", hap)
    res <- tryCatch(call_cdh(build_table(sample, hap)),
                    error=function(e){ message("  SKIP: ",conditionMessage(e)); NULL })
    if (is.null(res)) next
    dm <- res$d; win <- dm[dm$is_cdh, ]
    reg <- merge_regions(win[, c("Chr","Start","End","spec","total","OE","qval")])

    so <- file.path(OUTDIR, sample); dir.create(so, showWarnings=FALSE)
    # window-level CDH
    write.table(win[, c("Chr","Start","End","spec","shared","total","expected","theta_local","OE","pval","qval")],
                file.path(so, sprintf("CDH.window.%s.%s.tsv", sample, hap)), sep="\t", quote=FALSE, row.names=FALSE)
    # region-level CDH (merged)
    if (nrow(reg)>0)
      write.table(reg[, c("Chr","Start","End","n_windows","sum_spec","sum_total","region_frac","max_OE","min_qval","peak_Start")],
                  file.path(so, sprintf("CDH.region.%s.%s.tsv", sample, hap)), sep="\t", quote=FALSE, row.names=FALSE)

    summ <- rbind(summ, data.frame(sample, hap,
      n_windows_tested=res$n_tested, med_total=round(res$med_total),
      theta_global=round(res$theta_global,2), lambda_GC=round(res$lambda,3),
      baseline_frac=round(res$mu_baseline,5),
      n_CDH_window=nrow(win), n_CDH_region=ifelse(nrow(reg)>0,nrow(reg),0),
      median_region_windows=ifelse(nrow(reg)>0, median(reg$n_windows), NA),
      median_OE_CDH=round(ifelse(nrow(win)>0, median(win$OE), NA),2),
      cov_ratio=round(ifelse(nrow(win)>0, mean(win$total)/mean(dm$total), NA),2)))
  }
}

if (nrow(summ)==0) stop("nothing produced")
write.table(summ, file.path(OUTDIR,"ContactDifferenceHotspots.summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ===========================================================================
## Report
## ===========================================================================
med <- function(x) median(x, na.rm=TRUE)
cat("\n============================================================\n")
cat("CDH FINAL (M2) — per sample x hap\n")
cat("============================================================\n")
print(summ, row.names=FALSE)

cat("\n---- cohort summary (median across sample x hap) ----\n")
cat(sprintf("  lambda_GC        : %.3f   (calibration; ~1 good)\n", med(summ$lambda_GC)))
cat(sprintf("  cov_ratio        : %.2f   (~1 = no coverage confound)\n", med(summ$cov_ratio)))
cat(sprintf("  n_CDH windows    : median=%.0f  range=%d-%d\n", med(summ$n_CDH_window), min(summ$n_CDH_window), max(summ$n_CDH_window)))
cat(sprintf("  n_CDH regions    : median=%.0f  range=%d-%d\n", med(summ$n_CDH_region), min(summ$n_CDH_region), max(summ$n_CDH_region)))
cat(sprintf("  median OE of CDH : %.2f\n", med(summ$median_OE_CDH)))
for (h in c("Hap1","Hap2")) { s <- summ[summ$hap==h,]
  cat(sprintf("  %s regions: mean=%.0f median=%.0f range=%d-%d\n", h, mean(s$n_CDH_region), median(s$n_CDH_region), min(s$n_CDH_region), max(s$n_CDH_region))) }
by <- split(summ, summ$sample); rr <- c()
for (s in by) if (nrow(s)==2 && min(s$n_CDH_region)>0) rr <- c(rr, max(s$n_CDH_region)/min(s$n_CDH_region))
if (length(rr)) cat(sprintf("  Hap1/Hap2 region-count comparability: median ratio=%.2f\n", median(rr)))

cat("\nOutput per sample dir:\n")
cat("  CDH.window.<sample>.<hap>.tsv  — window-level hotspots (Chr,Start,End,spec,shared,total,expected,theta_local,OE,pval,qval)\n")
cat("  CDH.region.<sample>.<hap>.tsv  — merged regions (Chr,Start,End,n_windows,sum_spec,sum_total,region_frac,max_OE,min_qval,peak_Start)\n")
cat(sprintf("Plus %s/ContactDifferenceHotspots.summary.tsv\n", OUTDIR))
