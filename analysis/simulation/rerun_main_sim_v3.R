# rerun_main_sim_v3.R -- main simulation under engine v0.3.0 (model-implied centering, parametric reference
# with direct-dependence confirmation, heterogeneity test power for group-specific pairs, Q3* and UVA comparators).
# 162 conditions x 100 reps, B = 2,000. Chunk across terminals: --chunk=1/3, 2/3, 3/3 ; restart with --skip-done.
#   nohup caffeinate -i Rscript rerun_main_sim_v3.R --chunk=1/3 --skip-done > log_sim_v3_1.txt 2>&1 &
suppressPackageStartupMessages({library(biasawareLID); library(lavaan); library(Matrix)})
a <- commandArgs(TRUE); getarg <- function(k,d){h<-grep(paste0("^--",k,"="),a,value=TRUE); if(length(h)) sub(paste0("^--",k,"="),"",h[1]) else d}
REPS <- as.integer(getarg("reps","100")); B <- as.integer(getarg("B","2000")); chunk <- getarg("chunk","1/1"); SKIP <- any(a=="--skip-done")
OUT <- "sim_v3"; dir.create(OUT, showWarnings = FALSE)
grid <- expand.grid(N = c(250,500,1000), p = c(10,20,30), rho = c(.1,.3,.5), k = c(2,3), loadings = c("equal","spread"), stringsAsFactors = FALSE)
grid$cid <- seq_len(nrow(grid)); ck <- as.integer(strsplit(chunk,"/")[[1]])
if (ck[2] > 1) grid <- grid[split(seq_len(nrow(grid)), cut(seq_len(nrow(grid)), ck[2], labels = FALSE))[[ck[1]]], ]
ap <- function(path, df){ h<-!file.exists(path); suppressWarnings(write.table(df,path,sep=",",row.names=FALSE,col.names=h,append=!h)) }
done <- character(0); f <- file.path(OUT,"replication_level.csv"); if (SKIP && file.exists(f)) { d<-read.csv(f); done <- paste(d$cid,d$rep) }
has_qgraph <- requireNamespace("qgraph", quietly = TRUE)
uva_sel <- function(X, cut) { if (!has_qgraph) return(NULL)
  R <- as.matrix(Matrix::nearPD(cor(X), corr = TRUE)$mat); g <- qgraph::EBICglasso(R, n = nrow(X)); A <- abs(g); diag(A) <- 0
  k <- rowSums(A); L <- A %*% A; W <- (L + A)/(outer(k,k,pmin) + 1 - A); diag(W) <- 0; W[upper.tri(W)] > cut }
for (r in seq_len(nrow(grid))) { cd <- grid[r,]
  for (rep in seq_len(REPS)) { if (paste(cd$cid,rep) %in% done) next
    lam <- if (cd$loadings == "equal") .7 else seq(.4, .85, length.out = cd$p)
    s <- simulate_lid_data(N = cd$N, p = cd$p, rho = cd$rho, n_groups = cd$k, loading = lam, seed = 100000L*cd$cid + rep)
    t0 <- proc.time()[3]
    res <- tryCatch(bias_aware_lid(s$data, s$group, B = B, seed = 5000000L + 100000L*cd$cid + rep, centering = "model",
                    null = "parametric", direct_test = TRUE, heterogeneity_test = TRUE), error = function(e) NULL)
    if (is.null(res)) { cat("FAIL", cd$cid, rep, "\n", file = file.path(OUT,"failures.log"), append = TRUE); next }
    e <- res$edges; key <- paste(e$item_i, e$item_j); tr <- s$truth; tkey <- paste(tr$item_i, tr$item_j)
    planted <- key %in% tkey; gs <- key %in% tkey[tr$class == "group_specific"]; cs <- key %in% tkey[tr$class == "consistent"]
    tp <- sum(e$final & planted); fp <- sum(e$final & !planted); tpd <- sum(e$direct_BH & planted, na.rm = TRUE); fpd <- sum(e$direct_BH & !planted, na.rm = TRUE)
    het <- res$heterogeneity; hkey <- if (!is.null(het)) paste(het$item_i, het$item_j) else character(0)
    het_gs <- if (length(hkey)) sum(het$het_BH & hkey %in% tkey[tr$class=="group_specific"]) else 0
    het_cs_false <- if (length(hkey)) sum(het$het_BH & hkey %in% tkey[tr$class=="consistent"]) else 0
    q3 <- abs(e$Q3_star_signed); u25 <- uva_sel(s$data, .25); u20 <- uva_sel(s$data, .20)
    ap(f, data.frame(cid = cd$cid, rep = rep, N = cd$N, p = cd$p, rho = cd$rho, k = cd$k, loadings = cd$loadings,
      n_bh = res$n_bh, n_E0 = res$n_E0, n_final = res$n_final, beta = res$beta_star,
      recall = tp/sum(planted), precision = ifelse(tp+fp > 0, tp/(tp+fp), NA), fdr_uncond = ifelse(tp+fp > 0, fp/(tp+fp), 0),
      recall_direct = tpd/sum(planted), fdr_direct = ifelse(tpd+fpd > 0, fpd/(tpd+fpd), 0), n_direct = tpd+fpd,
      recall_gs = ifelse(sum(gs) > 0, sum(e$final & gs)/sum(gs), NA), recall_cs = ifelse(sum(cs) > 0, sum(e$final & cs)/sum(cs), NA),
      gs_in_E0 = sum(e$E0 & gs), het_power_gs = ifelse(sum(e$E0 & gs) > 0, het_gs/sum(e$E0 & gs), NA), het_false_cs = het_cs_false, n_tested = length(hkey),
      q3_20_recall = sum(q3 > .20 & planted)/sum(planted), q3_20_fdr = ifelse(sum(q3 > .20) > 0, sum(q3 > .20 & !planted)/sum(q3 > .20), 0),
      q3_15_recall = sum(q3 > .15 & planted)/sum(planted), q3_15_fdr = ifelse(sum(q3 > .15) > 0, sum(q3 > .15 & !planted)/sum(q3 > .15), 0),
      uva25_recall = if (is.null(u25)) NA else sum(u25 & planted)/sum(planted), uva20_recall = if (is.null(u20)) NA else sum(u20 & planted)/sum(planted),
      uva20_fdr = if (is.null(u20)) NA else ifelse(sum(u20) > 0, sum(u20 & !planted)/sum(u20), 0),
      runtime = round(proc.time()[3]-t0,1)))
  }
  cat(sprintf("[%s] cond %d/%d (id %d) done\n", format(Sys.time(),"%H:%M:%S"), r, nrow(grid), cd$cid))
}
cat("chunk finished\n")
