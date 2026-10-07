#!/usr/bin/env Rscript
# written by Lingbin

# Call contact difference hotspots (CDHs) separately for each sample and haplotype.
# Requires R with MASS and splines. Edit spec_dir, total_dir, anno_base and OUTDIR.
# Run: Rscript /path/to/call_contact_difference_hotspots.R
# Input files (four tab-separated columns, no header: Chr, Start, End, count):
#   <spec_dir>/<sample>/<sample>.SelfMapping.<hap>.HPRC.25kb.Unique.position.filter25kb.position
#   <total_dir>/<sample>/<sample>.SelfMapping.<hap>.HPRC.25kb.Total.position.filter25kb.position
# Here <hap> is Hap1 or Hap2; sample directories start with HG or NA.
# Prepare both inputs with count_long_range_cis_contacts.pl and matching bins/filters.
# Total counts include both shared and haplotype-specific contacts for that haplotype.
# Use unique window coordinates and nonnegative integer counts from the same reference.
# Start/End are 0-based, half-open coordinates. Missing specific windows are set to zero.
# The original filename convention, including .position.filter25kb.position, is retained.
# Acrocentric masks are headerless BED3 files in the corresponding assembly coordinates.
# Hap1 mask lookup: <sample>.pat.acro.bed, then <sample>.hap1.acro.bed.
# Hap2 mask lookup: <sample>.mat.acro.bed, then <sample>.hap2.acro.bed.
# The first existing file is used; when neither exists, no acrocentric mask is applied.
# Retain full 25-kb autosomal windows matching VALID_CHR_PATTERN; shared = total - spec.
# Fit a log-link negative-binomial mean model using a spline of log1p(shared) and Chr.
# Fitting uses up to two iterations with upper-tail residual trimming and seed 1.
# Estimate dispersion in predicted-count strata, then linearly interpolate per window.
# Test P(X >= observed specific count); apply BH within each sample and haplotype.
# Tested windows require total >= 50; CDHs require q < 0.01, OE >= 2 and spec >= 10.
# Merge adjacent CDH windows on the same chromosome with REGION_GAP = 0.
# Per-sample outputs: CDH.window.<sample>.<hap>.tsv and CDH.region.<sample>.<hap>.tsv.
# Region summaries include the minimum window q value and maximum window OE.
# peak_Start is the start of the window with the largest OE within each region.
# Region files are written when CDH regions are present; output tables contain headers.
# The cohort summary is written to ContactDifferenceHotspots.summary.tsv.
# Run in a new empty output directory; adjust paths and chromosome naming for your data.

suppressPackageStartupMessages({ library(MASS); library(splines) })

## ===========================================================================
## Paths
## ===========================================================================
spec_dir  <- "/path/to/specific_contact_counts"
total_dir <- "/path/to/total_contact_counts"
anno_base <- "/path/to/acrocentric_masks"
OUTDIR    <- "./"

## ===========================================================================
## Analysis parameters
## ===========================================================================
THR        <- "25kb"       # Distance-filter label.
Q_CUTOFF   <- 0.01
EFFECT     <- 2            # Observed/expected count threshold.
MIN_TOTAL  <- 50
MIN_SPEC   <- 10
NS_DF      <- 4
N_REFIT    <- 2
TRIM_UPPER <- 0.995
FIT_MAX    <- 40000
N_STRATA   <- 20          # Number of predicted-count strata for dispersion estimation.
REGION_GAP <- 0           # Maximum gap for merging CDH windows (0 = touching or overlapping).
WIN        <- 25000
VALID_CHR_PATTERN <- "^chr([1-9]|1[0-9]|2[0-2])_RagTag$"   # Autosomal chromosome names used by the mapping references.

N_SUBSET   <- Inf         # Process all samples; set a finite value to process a subset.
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
# Summary statistic: upper-tail chi-square transform of the median p value / chi-square null median.
lambda_gc <- function(pv){ p <- pv[is.finite(pv)&pv>0&pv<=1]; if (length(p)<100) return(NA_real_)
  qchisq(median(p),1,lower.tail=FALSE)/qchisq(0.5,1,lower.tail=FALSE) }
strat_sample <- function(chr,n){ ib<-split(seq_along(chr),chr); nt<-length(chr)
  unlist(lapply(ib,function(ix){k<-min(length(ix),max(1L,round(n*length(ix)/nt))); if(length(ix)<=k) ix else sample(ix,k)}),use.names=FALSE) }

## ===========================================================================
## CDH caller: mean model and coverage-dependent dispersion
## ===========================================================================
call_cdh <- function(d) {
  dm <- d[d$total>=MIN_TOTAL,]; dm$Chr <- droplevels(factor(dm$Chr))
  form <- if (nlevels(dm$Chr)>1) spec ~ ns(log1p(shared),df=NS_DF)+Chr else spec ~ ns(log1p(shared),df=NS_DF)

  # ---- Mean model: NB fitting with iterative upper-tail residual trimming ----
  fi <- if (nrow(dm)>FIT_MAX) strat_sample(dm$Chr,FIT_MAX) else seq_len(nrow(dm)); fd <- dm[fi,]; keep <- rep(TRUE,nrow(fd)); fit <- NULL
  for (it in seq_len(N_REFIT)) {
    fit <- tryCatch(MASS::glm.nb(form,data=fd[keep,],link=log,control=glm.control(maxit=100,epsilon=1e-8)),error=function(e) NULL)
    if (is.null(fit)) { form <- if (nlevels(dm$Chr)>1) spec~log1p(shared)+Chr else spec~log1p(shared); fit <- MASS::glm.nb(form,data=fd[keep,],link=log) }
    muf<-predict(fit,newdata=fd,type="response"); pe<-(fd$spec-muf)/sqrt(pmax(muf+muf^2/fit$theta,1e-12))
    ct<-as.numeric(quantile(pe,TRIM_UPPER,na.rm=TRUE,type=8)); nk<-is.finite(pe)&pe<=ct; if (identical(nk,keep)) break; keep<-nk
  }
  mu <- predict(fit,newdata=dm,type="response"); theta_global <- fit$theta

  # ---- Dispersion: estimate sigma in predicted-count strata and interpolate ----
  r  <- ((dm$spec - mu)^2 - mu) / mu^2                       # per-window dispersion contribution; E[r]=sigma=1/theta (NBI)
  br <- unique(quantile(mu, seq(0,1,length.out=N_STRATA+1), na.rm=TRUE)); br[1]<- -Inf; br[length(br)]<-Inf
  b  <- cut(mu, br, include.lowest=TRUE)
  centers <- tapply(mu, b, median)
  sig_s   <- tapply(seq_along(r), b, function(ix){ v <- r[ix]; v <- v[is.finite(v)]
    v <- v[v <= quantile(v, TRIM_UPPER, na.rm=TRUE)]; max(mean(v, na.rm=TRUE), 1e-6) })   # Exclude the upper tail within each stratum.
  ok <- is.finite(centers) & is.finite(sig_s)
  sig_i <- pmax(approx(centers[ok], sig_s[ok], xout=mu, rule=2)$y, 1e-6)                   # Linearly interpolate per-window sigma.
  theta_i <- 1/sig_i

  # ---- One-sided NB tail probability (per-window theta) and BH adjustment ----
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
cat(sprintf("Calling CDH on %d samples (25kb, filter%s, coverage-dependent dispersion)\n\n", length(samplelist), THR))

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

    # cov_ratio = mean total count in CDH windows / mean total count in tested windows.
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
cat("CDH results — per sample x haplotype\n")
cat("============================================================\n")
print(summ, row.names=FALSE)

cat("\n---- cohort summary (median across sample x hap) ----\n")
cat(sprintf("  lambda_GC        : %.3f\n", med(summ$lambda_GC)))
cat(sprintf("  cov_ratio        : %.2f\n", med(summ$cov_ratio)))
cat(sprintf("  n_CDH windows    : median=%.0f  range=%d-%d\n", med(summ$n_CDH_window), min(summ$n_CDH_window), max(summ$n_CDH_window)))
cat(sprintf("  n_CDH regions    : median=%.0f  range=%d-%d\n", med(summ$n_CDH_region), min(summ$n_CDH_region), max(summ$n_CDH_region)))
cat(sprintf("  median OE of CDH : %.2f\n", med(summ$median_OE_CDH)))
for (h in c("Hap1","Hap2")) { s <- summ[summ$hap==h,]
  cat(sprintf("  %s regions: mean=%.0f median=%.0f range=%d-%d\n", h, mean(s$n_CDH_region), median(s$n_CDH_region), min(s$n_CDH_region), max(s$n_CDH_region))) }
by <- split(summ, summ$sample); rr <- c()
for (s in by) if (nrow(s)==2 && min(s$n_CDH_region)>0) rr <- c(rr, max(s$n_CDH_region)/min(s$n_CDH_region))
if (length(rr)) cat(sprintf("  Larger/smaller haplotype region-count ratio: median=%.2f\n", median(rr)))

cat("\nOutput per sample dir:\n")
cat("  CDH.window.<sample>.<hap>.tsv  — window-level hotspots (Chr,Start,End,spec,shared,total,expected,theta_local,OE,pval,qval)\n")
cat("  CDH.region.<sample>.<hap>.tsv  — merged regions (Chr,Start,End,n_windows,sum_spec,sum_total,region_frac,max_OE,min_qval,peak_Start)\n")
cat(sprintf("Plus %s/ContactDifferenceHotspots.summary.tsv\n", OUTDIR))
