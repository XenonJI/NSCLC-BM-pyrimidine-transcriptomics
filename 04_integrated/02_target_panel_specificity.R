# Target-panel specificity calibration
# Benchmarks the fixed five-gene module against non-overlapping proliferation/replication
# comparators and 10,000 proliferation-matched five-gene sets constructed without anatomical-site labels.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260818)

ROOT <- "D:/LUAD_LM_PYRIMIDINE"
OUT <- file.path(
  ROOT, "results", "major_revision",
  "E9A_target_panel_specificity"
)
CACHE <- file.path(
  ROOT, "results", "major_revision",
  "E2B_target_gene_cache"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

required <- c("data.table", "clubSandwich")
missing_required <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_required) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_required, collapse = ", "),
    "\nInstall once with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_required), collapse = ", "),
    "), dependencies = NA)"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(clubSandwich)
})

cat("\n============================================================\n")
cat("STEP E9A: TARGET-PANEL SPECIFICITY\n")
cat("============================================================\n")

# ----------------------------------------------------------------------
# 1. Fixed gene definitions
# ----------------------------------------------------------------------
module_genes <- c(
  "DHFR", "DHODH", "SHMT1", "TYMS", "UMPS"
)

s_genes_original <- c(
  "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2",
  "MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2",
  "RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7",
  "POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1",
  "TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1",
  "CHAF1B","BRIP1","E2F8"
)

g2m_genes_original <- c(
  "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80",
  "CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A",
  "SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E",
  "TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK",
  "CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8",
  "ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5",
  "CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"
)

# Exactly the same non-overlap rule used in E2B.
s_genes <- setdiff(s_genes_original, module_genes)
g2m_genes <- setdiff(g2m_genes_original, module_genes)
cycle_pool <- unique(c(s_genes, g2m_genes))

if (!identical(intersect(s_genes_original, module_genes), "TYMS")) {
  stop("Unexpected S-phase/module overlap.")
}
if (length(intersect(g2m_genes_original, module_genes)) != 0L) {
  stop("Unexpected G2M/module overlap.")
}

# Fixed replication-machinery comparator chosen BEFORE outcome inspection.
# All genes come from the pre-existing E2B S-phase cache.
dna_replication_genes <- c(
  "MCM2","MCM4","MCM5","MCM6",
  "PCNA","FEN1","GINS2","PRIM1",
  "RFC2","RPA2","POLD3",
  "CDC45","CDC6","POLA1",
  "TIPIN","DSCC1","CLSPN","CHAF1B"
)

nested_panels <- list(
  five_gene_fixed = module_genes,
  TYMS_UMPS_only = c("TYMS", "UMPS"),
  minus_TYMS = setdiff(module_genes, "TYMS"),
  minus_UMPS = setdiff(module_genes, "UMPS"),
  minus_TYMS_and_UMPS = setdiff(module_genes, c("TYMS", "UMPS"))
)

comparator_panels <- list(
  fixed_five_gene = module_genes,
  nonoverlap_S_phase = s_genes,
  nonoverlap_G2M = g2m_genes,
  DNA_replication_machinery = dna_replication_genes
)

# ----------------------------------------------------------------------
# 2. Helpers
# ----------------------------------------------------------------------
find_exact_file <- function(basename_target) {
  roots <- c(
    file.path(ROOT, "results", "major_revision"),
    file.path(ROOT, "results", "reviewer_defense"),
    file.path(ROOT, "results"),
    ROOT
  )

  hits <- character(0)

  for (rr in roots) {
    if (!dir.exists(rr)) next
    hh <- list.files(
      rr,
      recursive = TRUE,
      full.names = TRUE,
      pattern = paste0("^", gsub("\\.", "\\\\.", basename_target), "$")
    )
    hits <- unique(c(hits, hh))
  }

  if (length(hits) == 0L) {
    stop("Required file not found: ", basename_target)
  }

  ord <- order(
    !grepl("[/\\\\]major_revision[/\\\\]", hits),
    !grepl("[/\\\\]reviewer_defense[/\\\\]", hits),
    nchar(hits)
  )
  hits[ord][1]
}

cr2_score <- function(dat, label, formula_rhs = "site_BM") {
  fml <- as.formula(
    paste("score ~", formula_rhs)
  )
  fit <- lm(fml, data = dat)

  tt <- as.data.frame(
    clubSandwich::coef_test(
      fit,
      vcov = "CR2",
      cluster = dat$corrected_patient_id,
      test = "Satterthwaite",
      alternative = "two-sided"
    )
  )

  coef_col <- grep(
    "^Coef\\.?$|^coefficient$|^term$",
    names(tt),
    ignore.case = TRUE,
    value = TRUE
  )
  if (length(coef_col) >= 1L) {
    tt$term_internal <- as.character(tt[[coef_col[1]]])
  } else {
    tt$term_internal <- rownames(tt)
  }

  rr <- tt[tt$term_internal == "site_BM", , drop = FALSE]
  if (nrow(rr) != 1L) {
    stop("Could not extract site_BM from: ", label)
  }

  pick_num <- function(patterns) {
    nms <- names(rr)
    for (pat in patterns) {
      hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
      if (length(hit) >= 1L) return(as.numeric(rr[[hit[1]]]))
    }
    NA_real_
  }

  beta <- unname(coef(fit)["site_BM"])
  se <- pick_num(c("^SE$", "std"))
  df <- pick_num(c("df_Satt", "d.f"))
  p <- pick_num(c("p_Satt", "p-val", "p_value", "^p$"))

  crit <- if (is.finite(df)) qt(0.975, df = df) else qnorm(0.975)

  data.table(
    analysis = label,
    model = paste(deparse(fml), collapse = ""),
    n_obs = nrow(dat),
    n_clusters = uniqueN(dat$corrected_patient_id),
    n_BM = sum(dat$site_BM == 1L),
    n_PT = sum(dat$site_BM == 0L),
    estimate_site_BM = beta,
    CR2_SE = se,
    Satterthwaite_df = df,
    two_sided_p_CR2 = p,
    CI95_low_CR2 = beta - crit * se,
    CI95_high_CR2 = beta + crit * se
  )
}

score_from_genes <- function(zmat, genes) {
  gg <- intersect(genes, rownames(zmat))
  if (length(gg) == 0L) stop("No genes available for requested panel.")
  colMeans(zmat[gg, , drop = FALSE])
}

safe_spearman <- function(x, y) {
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    return(NA_real_)
  }
  suppressWarnings(
    cor(x, y, method = "spearman", use = "complete.obs")
  )
}

# ----------------------------------------------------------------------
# 3. Load the already-generated compact E2B caches
# ----------------------------------------------------------------------
if (!dir.exists(CACHE)) {
  stop(
    "E2B cache directory is missing:\n", CACHE,
    "\nRun STEP_E2B V2 first. This E9A script will NOT rescan giant raw matrices."
  )
}

cache_files <- list.files(
  CACHE,
  pattern = "_E2B_targets\\.rds$",
  full.names = TRUE
)

if (length(cache_files) != 40L) {
  stop(
    "Expected 40 E2B target-gene cache files; found ",
    length(cache_files), "."
  )
}

cycle_file <- find_exact_file(
  "STEP_E2B_04_specimen_cellcycle_composition.tsv"
)
cycle_dt <- fread(cycle_file, showProgress = FALSE)

needed_cycle <- c(
  "sample", "corrected_patient_id", "cohort",
  "cycling_like_fraction"
)
if (!all(needed_cycle %in% names(cycle_dt))) {
  stop("Unexpected E2B cycling table schema.")
}
cycle_dt <- cycle_dt[, ..needed_cycle]

# One specimen-level pseudobulk vector per cache.
pb_list <- vector("list", length(cache_files))
meta_list <- vector("list", length(cache_files))

for (i in seq_along(cache_files)) {
  obj <- readRDS(cache_files[i])

  if (!all(c(
    "sample", "cohort", "corrected_patient_id",
    "libsize", "counts"
  ) %in% names(obj))) {
    stop("Unexpected E2B cache object: ", basename(cache_files[i]))
  }

  cnt <- obj$counts
  lib <- as.numeric(obj$libsize)

  if (
    is.null(rownames(cnt)) ||
    any(!is.finite(lib)) ||
    any(lib <= 0)
  ) {
    stop("Invalid cache content: ", basename(cache_files[i]))
  }

  total_lib <- sum(lib)
  gene_count <- rowSums(cnt)

  log2cpm <- log2(
    gene_count / total_lib * 1e6 + 1
  )

  pb_list[[i]] <- data.table(
    gene = rownames(cnt),
    sample = as.character(obj$sample),
    log2CPM = as.numeric(log2cpm),
    pseudobulk_count = as.numeric(gene_count),
    total_library_size = total_lib
  )

  meta_list[[i]] <- data.table(
    sample = as.character(obj$sample),
    corrected_patient_id = as.character(obj$corrected_patient_id),
    cohort_cache = as.character(obj$cohort)
  )

  rm(obj, cnt)
  if (i %% 10L == 0L) gc()
}

pb_long <- rbindlist(pb_list)
meta <- unique(rbindlist(meta_list))

meta <- merge(
  meta,
  cycle_dt,
  by = c("sample", "corrected_patient_id"),
  all.x = TRUE,
  sort = FALSE
)

if (
  nrow(meta) != 40L ||
  any(!is.finite(meta$cycling_like_fraction))
) {
  stop("Failed to merge cycling fraction onto 40 cached specimens.")
}

meta[, site_BM := as.integer(cohort == "Brain_metastasis")]

if (sum(meta$site_BM) != 30L || sum(meta$site_BM == 0L) != 10L) {
  stop("Unexpected E9A cohort counts.")
}

fwrite(
  pb_long,
  file.path(OUT, "STEP_E9A_01_target_panel_specimen_pseudobulk_long.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 4. Construct gene x specimen expression matrix
# ----------------------------------------------------------------------
expr_wide <- dcast(
  pb_long,
  gene ~ sample,
  value.var = "log2CPM"
)

gene_names <- expr_wide$gene
expr_mat <- as.matrix(expr_wide[, -1])
rownames(expr_mat) <- gene_names

# Match columns to metadata.
meta <- meta[match(colnames(expr_mat), sample)]
if (anyNA(meta$sample)) stop("Expression/meta sample order failed.")

# Remove invariant genes only from secondary comparator calculations.
gene_mean <- rowMeans(expr_mat)
gene_sd <- apply(expr_mat, 1, sd)
usable <- is.finite(gene_sd) & gene_sd > 0

zmat <- sweep(
  sweep(expr_mat[usable, , drop = FALSE], 1, gene_mean[usable], "-"),
  1,
  gene_sd[usable],
  "/"
)

# ----------------------------------------------------------------------
# 5. Gene-level cycle coupling and site-effect audit
# ----------------------------------------------------------------------
gene_metrics <- rbindlist(lapply(
  rownames(zmat),
  function(g) {
    x <- as.numeric(expr_mat[g, ])
    z <- as.numeric(zmat[g, ])

    data.table(
      gene = g,
      gene_class = fifelse(
        g %chin% module_genes, "fixed_five_gene",
        fifelse(
          g %chin% s_genes, "nonoverlap_S_phase",
          fifelse(
            g %chin% g2m_genes, "nonoverlap_G2M",
            "other"
          )
        )
      ),
      mean_log2CPM = mean(x),
      sd_log2CPM = sd(x),
      spearman_rho_with_cycling_fraction =
        safe_spearman(x, meta$cycling_like_fraction),
      site_beta_standardized =
        mean(z[meta$site_BM == 1L]) -
        mean(z[meta$site_BM == 0L])
    )
  }
))

gene_metrics[
  ,
  abs_rho_cycling := abs(
    spearman_rho_with_cycling_fraction
  )
]

fwrite(
  gene_metrics,
  file.path(OUT, "STEP_E9A_02_gene_cycle_and_site_metrics.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 6. Fixed comparator panels
# ----------------------------------------------------------------------
panel_results <- list()

for (nm in names(comparator_panels)) {
  gg <- intersect(
    comparator_panels[[nm]],
    rownames(zmat)
  )

  if (length(gg) < 2L) {
    stop("Too few usable genes in comparator panel: ", nm)
  }

  sc <- score_from_genes(zmat, gg)

  dd <- copy(meta)
  dd[, score := sc]

  rr0 <- cr2_score(
    dd,
    label = paste0(nm, ": site only"),
    formula_rhs = "site_BM"
  )
  rr0[, `:=`(
    panel = nm,
    n_genes = length(gg),
    genes = paste(gg, collapse = ";"),
    spearman_score_vs_cycling =
      safe_spearman(sc, dd$cycling_like_fraction)
  )]

  rr1 <- cr2_score(
    dd,
    label = paste0(nm, ": site + cycling"),
    formula_rhs = "site_BM + cycling_like_fraction"
  )
  rr1[, `:=`(
    panel = nm,
    n_genes = length(gg),
    genes = paste(gg, collapse = ";"),
    spearman_score_vs_cycling =
      safe_spearman(sc, dd$cycling_like_fraction)
  )]

  panel_results[[nm]] <- rbind(rr0, rr1, fill = TRUE)
}

panel_results <- rbindlist(panel_results, fill = TRUE)

fwrite(
  panel_results,
  file.path(OUT, "STEP_E9A_03_fixed_comparator_panels_CR2.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 7. Nested five-gene decomposition
# ----------------------------------------------------------------------
nested_results <- list()

for (nm in names(nested_panels)) {
  gg <- intersect(
    nested_panels[[nm]],
    rownames(zmat)
  )

  sc <- score_from_genes(zmat, gg)

  dd <- copy(meta)
  dd[, score := sc]

  rr0 <- cr2_score(
    dd,
    label = paste0(nm, ": site only"),
    formula_rhs = "site_BM"
  )
  rr0[, `:=`(
    panel = nm,
    n_genes = length(gg),
    genes = paste(gg, collapse = ";"),
    spearman_score_vs_cycling =
      safe_spearman(sc, dd$cycling_like_fraction)
  )]

  rr1 <- cr2_score(
    dd,
    label = paste0(nm, ": site + cycling"),
    formula_rhs = "site_BM + cycling_like_fraction"
  )
  rr1[, `:=`(
    panel = nm,
    n_genes = length(gg),
    genes = paste(gg, collapse = ";"),
    spearman_score_vs_cycling =
      safe_spearman(sc, dd$cycling_like_fraction)
  )]

  nested_results[[nm]] <- rbind(
    rr0, rr1, fill = TRUE
  )
}

nested_results <- rbindlist(
  nested_results,
  fill = TRUE
)

fwrite(
  nested_results,
  file.path(OUT, "STEP_E9A_04_nested_pyrimidine_panels_CR2.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 8. Proliferation-matched random 5-gene null
#
# Matching features:
#   1) absolute correlation with cycling fraction;
#   2) mean log2CPM;
#   3) SD of log2CPM.
#
# Site labels are NOT used to define candidate neighborhoods.
# ----------------------------------------------------------------------
candidate_genes <- intersect(
  cycle_pool,
  gene_metrics[
    is.finite(abs_rho_cycling) &
      is.finite(mean_log2CPM) &
      is.finite(sd_log2CPM),
    gene
  ]
)

if (length(candidate_genes) < 50L) {
  stop(
    "Too few usable non-overlapping cycle genes for matched null: ",
    length(candidate_genes)
  )
}

module_metrics <- gene_metrics[
  gene %chin% module_genes
][match(module_genes, gene)]

if (
  nrow(module_metrics) != 5L ||
  any(!is.finite(module_metrics$abs_rho_cycling))
) {
  stop("Five module genes are not all usable for matching.")
}

cand <- gene_metrics[
  gene %chin% candidate_genes
]

feature_names <- c(
  "abs_rho_cycling",
  "mean_log2CPM",
  "sd_log2CPM"
)

# Standardize matching features using candidate-pool location/scale.
mu_f <- vapply(
  feature_names,
  function(v) mean(cand[[v]]),
  numeric(1)
)
sd_f <- vapply(
  feature_names,
  function(v) sd(cand[[v]]),
  numeric(1)
)

if (any(!is.finite(sd_f)) || any(sd_f == 0)) {
  stop("Invalid matching-feature scale.")
}

scaled_candidate <- sapply(
  feature_names,
  function(v) (cand[[v]] - mu_f[v]) / sd_f[v]
)

scaled_module <- sapply(
  feature_names,
  function(v) (module_metrics[[v]] - mu_f[v]) / sd_f[v]
)

K_NEIGHBORS <- min(20L, nrow(cand))

nearest_list <- vector(
  "list",
  length(module_genes)
)
names(nearest_list) <- module_genes

neighbor_audit <- list()

for (j in seq_along(module_genes)) {
  d2 <- rowSums(
    sweep(
      scaled_candidate,
      2,
      scaled_module[j, ],
      "-"
    )^2
  )

  oo <- order(d2)
  take <- oo[seq_len(K_NEIGHBORS)]

  nearest_list[[j]] <- cand$gene[take]

  neighbor_audit[[j]] <- data.table(
    module_gene = module_genes[j],
    candidate_gene = cand$gene[take],
    rank = seq_along(take),
    squared_feature_distance = d2[take],
    module_abs_rho_cycling =
      module_metrics$abs_rho_cycling[j],
    candidate_abs_rho_cycling =
      cand$abs_rho_cycling[take],
    module_mean_log2CPM =
      module_metrics$mean_log2CPM[j],
    candidate_mean_log2CPM =
      cand$mean_log2CPM[take],
    module_sd_log2CPM =
      module_metrics$sd_log2CPM[j],
    candidate_sd_log2CPM =
      cand$sd_log2CPM[take]
  )
}

neighbor_audit <- rbindlist(neighbor_audit)

fwrite(
  neighbor_audit,
  file.path(
    OUT,
    "STEP_E9A_05_proliferation_matched_candidate_neighborhoods.tsv"
  ),
  sep = "\t"
)

# Observed target-panel five-gene score.
obs_score <- score_from_genes(
  zmat,
  module_genes
)
obs_beta <- (
  mean(obs_score[meta$site_BM == 1L]) -
  mean(obs_score[meta$site_BM == 0L])
)
obs_rho_cycle <- safe_spearman(
  obs_score,
  meta$cycling_like_fraction
)

B_NULL <- 10000L
null_out <- vector("list", B_NULL)

for (b in seq_len(B_NULL)) {
  # Randomize allocation order so one module gene does not always get
  # first choice in overlapping candidate neighborhoods.
  ord <- sample(seq_along(module_genes))

  chosen <- rep(NA_character_, length(module_genes))
  used <- character(0)

  for (jj in ord) {
    pool <- setdiff(
      nearest_list[[jj]],
      used
    )

    if (length(pool) == 0L) {
      # Transparent fallback remains within the fixed cycle-gene pool
      # and is independent of anatomical site.
      pool <- setdiff(candidate_genes, used)
    }

    pick <- sample(pool, 1L)
    chosen[jj] <- pick
    used <- c(used, pick)
  }

  if (
    anyNA(chosen) ||
    length(unique(chosen)) != 5L
  ) {
    next
  }

  sc <- score_from_genes(
    zmat,
    chosen
  )

  beta <- (
    mean(sc[meta$site_BM == 1L]) -
    mean(sc[meta$site_BM == 0L])
  )

  null_out[[b]] <- data.table(
    iteration = b,
    gene1 = chosen[1],
    gene2 = chosen[2],
    gene3 = chosen[3],
    gene4 = chosen[4],
    gene5 = chosen[5],
    site_beta_standardized_score = beta,
    spearman_score_vs_cycling =
      safe_spearman(
        sc,
        meta$cycling_like_fraction
      )
  )
}

null_dt <- rbindlist(
  null_out,
  fill = TRUE
)

if (nrow(null_dt) < 9500L) {
  stop(
    "Too few valid matched random sets: ",
    nrow(null_dt)
  )
}

empirical_two_sided_p <- (
  1 +
  sum(
    abs(null_dt$site_beta_standardized_score) >=
      abs(obs_beta)
  )
) / (1 + nrow(null_dt))

empirical_positive_tail_p <- (
  1 +
  sum(
    null_dt$site_beta_standardized_score >=
      obs_beta
  )
) / (1 + nrow(null_dt))

obs_percentile <- mean(
  null_dt$site_beta_standardized_score <
    obs_beta
)

null_summary <- data.table(
  analysis =
    "5-gene score vs proliferation-matched 5-gene cycle-marker null",
  n_specimens = nrow(meta),
  n_patient_clusters =
    uniqueN(meta$corrected_patient_id),
  n_valid_random_sets = nrow(null_dt),
  observed_site_beta = obs_beta,
  observed_spearman_score_vs_cycling =
    obs_rho_cycle,
  null_beta_median =
    median(null_dt$site_beta_standardized_score),
  null_beta_q025 =
    quantile(
      null_dt$site_beta_standardized_score,
      0.025,
      names = FALSE
    ),
  null_beta_q975 =
    quantile(
      null_dt$site_beta_standardized_score,
      0.975,
      names = FALSE
    ),
  observed_beta_percentile =
    obs_percentile,
  empirical_two_sided_p =
    empirical_two_sided_p,
  empirical_positive_tail_p_descriptive =
    empirical_positive_tail_p,
  null_score_cycle_rho_median =
    median(
      null_dt$spearman_score_vs_cycling,
      na.rm = TRUE
    ),
  null_score_cycle_rho_q025 =
    quantile(
      null_dt$spearman_score_vs_cycling,
      0.025,
      na.rm = TRUE,
      names = FALSE
    ),
  null_score_cycle_rho_q975 =
    quantile(
      null_dt$spearman_score_vs_cycling,
      0.975,
      na.rm = TRUE,
      names = FALSE
    )
)

fwrite(
  null_dt,
  file.path(
    OUT,
    "STEP_E9A_06_proliferation_matched_random_set_null_distribution.tsv"
  ),
  sep = "\t"
)

fwrite(
  null_summary,
  file.path(
    OUT,
    "STEP_E9A_07_proliferation_matched_random_set_null_summary.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 9. Matching-balance audit
# ----------------------------------------------------------------------
balance <- rbindlist(list(
  data.table(
    quantity = "abs rho(gene, cycling fraction)",
    observed_module_median =
      median(module_metrics$abs_rho_cycling),
    candidate_pool_median =
      median(cand$abs_rho_cycling)
  ),
  data.table(
    quantity = "mean log2CPM",
    observed_module_median =
      median(module_metrics$mean_log2CPM),
    candidate_pool_median =
      median(cand$mean_log2CPM)
  ),
  data.table(
    quantity = "SD log2CPM",
    observed_module_median =
      median(module_metrics$sd_log2CPM),
    candidate_pool_median =
      median(cand$sd_log2CPM)
  )
))

fwrite(
  balance,
  file.path(
    OUT,
    "STEP_E9A_08_matching_feature_balance.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 10. Guardrails / status
# ----------------------------------------------------------------------
five_site <- nested_results[
  panel == "five_gene_fixed" &
    grepl("site only", analysis, fixed = TRUE)
]

two_gene_site <- nested_results[
  panel == "TYMS_UMPS_only" &
    grepl("site only", analysis, fixed = TRUE)
]

three_gene_site <- nested_results[
  panel == "minus_TYMS_and_UMPS" &
    grepl("site only", analysis, fixed = TRUE)
]

specificity_flag <- if (
  empirical_two_sided_p < 0.05
) {
  "TARGET_PANEL_SPECIFICITY_SUPPORTED"
} else {
  "TARGET_PANEL_SPECIFICITY_NOT_ESTABLISHED"
}

guardrails <- c(
  "status=PASS_E9A_TARGET_PANEL_SPECIFICITY_COMPLETE",
  "primary_endpoint_replaced=FALSE",
  "gene_set_changed=FALSE",
  "site_labels_used_for_matching=FALSE",
  "target_panel_not_genome_wide=TRUE",
  paste0(
    "observed_5gene_site_beta=",
    signif(obs_beta, 6)
  ),
  paste0(
    "observed_5gene_cycle_rho=",
    signif(obs_rho_cycle, 6)
  ),
  paste0(
    "matched_null_two_sided_p=",
    signif(empirical_two_sided_p, 6)
  ),
  paste0(
    "observed_beta_percentile=",
    signif(obs_percentile, 6)
  ),
  paste0(
    "two_gene_TYMS_UMPS_site_beta_CR2=",
    signif(two_gene_site$estimate_site_BM, 6),
    "; p=",
    signif(two_gene_site$two_sided_p_CR2, 6)
  ),
  paste0(
    "three_gene_background_site_beta_CR2=",
    signif(three_gene_site$estimate_site_BM, 6),
    "; p=",
    signif(three_gene_site$two_sided_p_CR2, 6)
  ),
  paste0("specificity_decision=", specificity_flag),
  "wording_guardrail=Even if the matched-null test is significant, call the result pyrimidine-related rather than proof of pyrimidine-specific metabolic reprogramming.",
  "next_step_guardrail=If target-panel specificity is supported and a stronger pathway-specific claim is desired, perform a separate genome-wide pseudobulk specificity analysis. If it is not supported, do not escalate to a stronger pyrimidine-specific title."
)

writeLines(
  guardrails,
  file.path(OUT, "STEP_E9A_COMPLETE.txt")
)

# ----------------------------------------------------------------------
# 11. Console summary
# ----------------------------------------------------------------------
cat("\n============================================================\n")
cat("STEP E9A COMPLETE\n")
cat("============================================================\n\n")

cat("Fixed comparator panels:\n")
print(
  panel_results[
    grepl("site only", analysis, fixed = TRUE),
    .(
      panel, n_genes,
      estimate_site_BM,
      CI95_low_CR2,
      CI95_high_CR2,
      two_sided_p_CR2,
      spearman_score_vs_cycling
    )
  ]
)

cat("\nNested pyrimidine panels:\n")
print(
  nested_results[
    grepl("site only", analysis, fixed = TRUE),
    .(
      panel, n_genes,
      estimate_site_BM,
      CI95_low_CR2,
      CI95_high_CR2,
      two_sided_p_CR2,
      spearman_score_vs_cycling
    )
  ]
)

cat("\nProliferation-matched random-set null:\n")
print(null_summary)

cat("\nSpecificity decision: ", specificity_flag, "\n", sep = "")
cat("\nStatus: PASS_E9A_TARGET_PANEL_SPECIFICITY_COMPLETE\n")

cat("\nKey output files:\n")
cat(file.path(OUT, "STEP_E9A_02_gene_cycle_and_site_metrics.tsv"), "\n")
cat(file.path(OUT, "STEP_E9A_03_fixed_comparator_panels_CR2.tsv"), "\n")
cat(file.path(OUT, "STEP_E9A_04_nested_pyrimidine_panels_CR2.tsv"), "\n")
cat(file.path(OUT, "STEP_E9A_05_proliferation_matched_candidate_neighborhoods.tsv"), "\n")
cat(file.path(OUT, "STEP_E9A_07_proliferation_matched_random_set_null_summary.tsv"), "\n")
cat(file.path(OUT, "STEP_E9A_08_matching_feature_balance.tsv"), "\n")
cat(file.path(OUT, "STEP_E9A_COMPLETE.txt"), "\n")
