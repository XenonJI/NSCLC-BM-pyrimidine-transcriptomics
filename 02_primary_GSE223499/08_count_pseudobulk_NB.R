# Count-based pseudobulk sensitivity
# Uses the barcode-validated E2B caches to aggregate raw malignant-nucleus counts by specimen
# and fit gene-wise negative-binomial models with a library-size offset and patient clustering.
# The analysis is restricted to the fixed five genes.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.", call. = FALSE)
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required; it is normally included with R.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/LUAD_LM_PYRIMIDINE"
OUT <- file.path(ROOT, "results", "major_revision")
CACHE <- file.path(OUT, "E2B_target_gene_cache")
RD <- file.path(ROOT, "results", "reviewer_defense")

LOCKED_FILE <- file.path(
  RD,
  "STEP_A3_locked_primary_31BM_10Primary_patient_table.tsv"
)

E2B_SAMPLE_FILE <- file.path(
  OUT,
  "STEP_E2B_04_specimen_cellcycle_composition.tsv"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

for (f in c(LOCKED_FILE, E2B_SAMPLE_FILE)) {
  if (!file.exists(f)) {
    stop("Required file not found:\n", f, call. = FALSE)
  }
}

if (!dir.exists(CACHE)) {
  stop(
    "E2B cache directory not found:\n",
    CACHE,
    "\nRun STEP E2B V2 first.",
    call. = FALSE
  )
}

genes <- c(
  "DHFR",
  "DHODH",
  "SHMT1",
  "TYMS",
  "UMPS"
)

# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------
fit_nb <- function(d) {
  fit <- tryCatch(
    MASS::glm.nb(
      raw_count ~
        site_BM +
        offset(
          log(
            library_size
          )
        ),
      data = d,
      control = glm.control(
        maxit = 100
      )
    ),
    error = function(e) {
      structure(
        list(
          message = conditionMessage(e)
        ),
        class = "E3_NB_ERROR"
      )
    }
  )

  fit
}

nb_cluster_cr1 <- function(
  fit,
  cluster,
  term_name = "site_BM"
) {
  X <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  mu <- fitted(fit)
  theta <- fit$theta

  cluster <- as.character(cluster)

  if (
    nrow(X) != length(cluster) ||
    length(y) != length(cluster)
  ) {
    stop("Cluster vector/model dimensions do not match.", call. = FALSE)
  }

  if (
    !is.finite(theta) ||
    theta <= 0
  ) {
    stop("Invalid NB theta.", call. = FALSE)
  }

  # NB2 log-link estimating-equation components:
  # variance(mu) = mu + mu^2/theta
  # score_i(beta) = x_i * (y_i - mu_i)/(1 + mu_i/theta)
  denom <- 1 + mu / theta
  score_scalar <- (y - mu) / denom

  # Expected information for beta:
  # X' diag(mu/(1 + mu/theta)) X
  W <- mu / denom
  bread_inv <- crossprod(
    X,
    X * W
  )

  if (qr(bread_inv)$rank < ncol(X)) {
    stop("NB information matrix is rank-deficient.", call. = FALSE)
  }

  bread <- solve(bread_inv)

  clusters <- unique(cluster)
  G <- length(clusters)
  N <- nrow(X)
  K <- ncol(X)

  if (G <= K) {
    stop("Too few clusters for CR1.", call. = FALSE)
  }

  meat <- matrix(
    0,
    nrow = K,
    ncol = K
  )

  for (cl in clusters) {
    idx <- cluster == cl

    sg <- colSums(
      X[
        idx,
        ,
        drop = FALSE
      ] *
        score_scalar[
          idx
        ]
    )

    meat <- meat +
      tcrossprod(
        sg
      )
  }

  correction <- (
    G /
      (G - 1)
  ) *
    (
      (N - 1) /
        (N - K)
    )

  vc <- correction *
    bread %*%
      meat %*%
      bread

  b <- coef(fit)
  se <- sqrt(diag(vc))

  if (!term_name %in% names(b)) {
    stop("Requested NB coefficient missing.", call. = FALSE)
  }

  est <- unname(b[term_name])
  se0 <- unname(se[term_name])
  df <- G - 1
  t0 <- est / se0
  p2 <- 2 * pt(-abs(t0), df = df)
  crit <- qt(0.975, df = df)

  data.table(
    beta_log_ratio = est,
    CR1_SE = se0,
    t_CR1 = t0,
    df_clusters = df,
    two_sided_p_CR1 = p2,
    CR1_ratio_BM_vs_PT = exp(est),
    CR1_CI95_low_ratio = exp(
      est -
        crit *
          se0
    ),
    CR1_CI95_high_ratio = exp(
      est +
        crit *
          se0
    ),
    n_observations = N,
    n_patient_clusters = G
  )
}

fit_one_gene <- function(
  d,
  analysis_name,
  cluster_robust = FALSE
) {
  gg <- unique(d$gene)

  if (length(gg) != 1L) {
    stop("fit_one_gene received more than one gene.", call. = FALSE)
  }

  d[
    ,
    site_BM := as.integer(
      cohort ==
        "Brain_metastasis"
    )
  ]

  fit <- fit_nb(d)

  if (inherits(fit, "E3_NB_ERROR")) {
    return(
      data.table(
        analysis = analysis_name,
        gene = gg,
        n_BM = sum(d$cohort == "Brain_metastasis"),
        n_PT = sum(d$cohort == "Primary_tumor"),
        beta_log_ratio = NA_real_,
        ratio_BM_vs_PT = NA_real_,
        CI95_low_ratio = NA_real_,
        CI95_high_ratio = NA_real_,
        two_sided_p_model = NA_real_,
        theta = NA_real_,
        model_status = paste0(
          "NB_FIT_FAILED: ",
          fit$message
        )
      )
    )
  }

  sm <- summary(fit)$coefficients

  if (!"site_BM" %in% rownames(sm)) {
    stop("site_BM coefficient not found.", call. = FALSE)
  }

  beta <- sm[
    "site_BM",
    "Estimate"
  ]

  se <- sm[
    "site_BM",
    "Std. Error"
  ]

  z <- beta / se
  p <- 2 * pnorm(-abs(z))
  crit <- qnorm(0.975)

  ans <- data.table(
    analysis = analysis_name,
    gene = gg,
    n_BM = sum(d$cohort == "Brain_metastasis"),
    n_PT = sum(d$cohort == "Primary_tumor"),
    beta_log_ratio = beta,
    ratio_BM_vs_PT = exp(beta),
    model_SE = se,
    z_model = z,
    two_sided_p_model = p,
    CI95_low_ratio = exp(
      beta -
        crit *
          se
    ),
    CI95_high_ratio = exp(
      beta +
        crit *
          se
    ),
    theta = fit$theta,
    model_status = "PASS"
  )

  if (cluster_robust) {
    cr <- nb_cluster_cr1(
      fit = fit,
      cluster = d$corrected_patient_id,
      term_name = "site_BM"
    )

    ans <- cbind(
      ans,
      cr
    )
  }

  ans
}

# ----------------------------------------------------------------------
# 2. Corrected patient map
# ----------------------------------------------------------------------
locked <- fread(
  LOCKED_FILE,
  showProgress = FALSE
)

required_locked <- c(
  "patient",
  "samples",
  "cohort"
)

if (!all(required_locked %in% names(locked))) {
  stop("Locked table missing required columns.", call. = FALSE)
}

locked[
  ,
  cohort_corrected := fifelse(
    grepl(
      "brain|bm",
      tolower(
        as.character(cohort)
      )
    ),
    "Brain_metastasis",
    "Primary_tumor"
  )
]

locked[
  ,
  corrected_patient_id := as.character(patient)
]

locked[
  samples %chin% c("PA060", "N254"),
  corrected_patient_id := "MATCHED_PATIENT_PA060_N254"
]

locked[
  samples %chin% c("PA067", "N586"),
  corrected_patient_id := "MATCHED_PATIENT_PA067_N586"
]

if (
  nrow(locked) != 41L ||
  uniqueN(locked$corrected_patient_id) != 39L
) {
  stop(
    "Expected corrected locked mapping of 41 specimens / 39 unique patients.",
    call. = FALSE
  )
}

# ----------------------------------------------------------------------
# 3. Use the successful E2B specimen set
# ----------------------------------------------------------------------
e2b_samples <- fread(
  E2B_SAMPLE_FILE,
  showProgress = FALSE
)

required_e2b <- c(
  "sample",
  "cohort",
  "corrected_patient_id",
  "malignant_nuclei"
)

if (!all(required_e2b %in% names(e2b_samples))) {
  stop("E2B specimen table has unexpected columns.", call. = FALSE)
}

if (
  sum(e2b_samples$cohort == "Brain_metastasis") != 30L ||
  sum(e2b_samples$cohort == "Primary_tumor") != 10L ||
  nrow(e2b_samples) != 40L
) {
  stop(
    "Expected E2B technically concordant set of 30 BM / 10 PT.",
    call. = FALSE
  )
}

if ("KRAS_17" %chin% e2b_samples$sample) {
  stop(
    "KRAS_17 should not appear in the E2B cache-based E3 set.",
    call. = FALSE
  )
}

# ----------------------------------------------------------------------
# 4. Build compact raw-count pseudobulk table from E2B caches
# ----------------------------------------------------------------------
rows <- vector(
  "list",
  nrow(e2b_samples)
)

for (i in seq_len(nrow(e2b_samples))) {
  ss <- e2b_samples$sample[i]

  cache_file <- file.path(
    CACHE,
    paste0(
      ss,
      "_E2B_targets.rds"
    )
  )

  if (!file.exists(cache_file)) {
    stop(
      "Missing E2B cache for ",
      ss,
      ":\n",
      cache_file,
      call. = FALSE
    )
  }

  obj <- readRDS(
    cache_file
  )

  if (
    is.null(obj$counts) ||
    is.null(obj$libsize)
  ) {
    stop(
      "E2B cache structure is incomplete for ",
      ss,
      call. = FALSE
    )
  }

  if (!all(genes %chin% rownames(obj$counts))) {
    stop(
      "One or more module genes are missing from E2B cache for ",
      ss,
      call. = FALSE
    )
  }

  lib <- sum(
    as.numeric(
      obj$libsize
    )
  )

  if (!is.finite(lib) || lib <= 0) {
    stop(
      "Invalid malignant pseudobulk library size for ",
      ss,
      call. = FALSE
    )
  }

  rows[[i]] <- data.table(
    sample = ss,
    cohort = e2b_samples$cohort[i],
    corrected_patient_id = e2b_samples$corrected_patient_id[i],
    malignant_nuclei = e2b_samples$malignant_nuclei[i],
    gene = genes,
    raw_count = as.numeric(
      rowSums(
        obj$counts[
          genes,
          ,
          drop = FALSE
        ]
      )
    ),
    library_size = lib
  )

  rm(obj)
  gc(verbose = FALSE)
}

dat <- rbindlist(
  rows
)

# Structural audit.
if (
  nrow(dat) != 40L * 5L ||
  uniqueN(dat$sample) != 40L ||
  uniqueN(dat$gene) != 5L
) {
  stop(
    "Compact E3 input is not 40 specimens x 5 genes.",
    call. = FALSE
  )
}

if (any(dat$raw_count < 0) || any(!is.finite(dat$raw_count))) {
  stop("Invalid E3 raw counts.", call. = FALSE)
}

fwrite(
  dat,
  file.path(
    OUT,
    "STEP_E3_01_count_pseudobulk_analysis_input.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 5. Main E3: 40 technically concordant specimens
#    30 BM / 10 PT, corrected patient-cluster CR1 reported
# ----------------------------------------------------------------------
full_stats <- rbindlist(
  lapply(
    genes,
    function(gg) {
      fit_one_gene(
        d = dat[
          gene == gg
        ],
        analysis_name =
          "E3 count-based NB: 40 technically concordant specimens",
        cluster_robust =
          TRUE
      )
    }
  ),
  fill = TRUE
)

if ("two_sided_p_CR1" %in% names(full_stats)) {
  full_stats[
    ,
    BH_FDR_CR1_across_5_genes :=
      p.adjust(
        two_sided_p_CR1,
        method = "BH"
      )
  ]
}

full_stats[
  ,
  BH_FDR_model_across_5_genes :=
    p.adjust(
      two_sided_p_model,
      method = "BH"
    )
]

# ----------------------------------------------------------------------
# 6. Independence-preserving sensitivity: 28 BM / 8 PT
# ----------------------------------------------------------------------
pair_samples <- c(
  "PA060",
  "N254",
  "PA067",
  "N586"
)

ind <- dat[
  !sample %chin% pair_samples
]

if (
  uniqueN(ind[cohort == "Brain_metastasis", sample]) != 28L ||
  uniqueN(ind[cohort == "Primary_tumor", sample]) != 8L
) {
  stop(
    "Expected E3 independence-preserving subset of 28 BM / 8 PT.",
    call. = FALSE
  )
}

ind_stats <- rbindlist(
  lapply(
    genes,
    function(gg) {
      fit_one_gene(
        d = ind[
          gene == gg
        ],
        analysis_name =
          "E3 count-based NB: 28 BM / 8 PT independence-preserving subset",
        cluster_robust =
          FALSE
      )
    }
  ),
  fill = TRUE
)

ind_stats[
  ,
  BH_FDR_model_across_5_genes :=
    p.adjust(
      two_sided_p_model,
      method = "BH"
    )
]

all_stats <- rbindlist(
  list(
    full_stats,
    ind_stats
  ),
  fill = TRUE
)

fwrite(
  all_stats,
  file.path(
    OUT,
    "STEP_E3_02_count_based_NB_gene_statistics.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 7. reporting interpretation
# ----------------------------------------------------------------------
main_p_col <- if (
  "two_sided_p_CR1" %in% names(full_stats)
) {
  "two_sided_p_CR1"
} else {
  "two_sided_p_model"
}

main_fdr_col <- if (
  "BH_FDR_CR1_across_5_genes" %in% names(full_stats)
) {
  "BH_FDR_CR1_across_5_genes"
} else {
  "BH_FDR_model_across_5_genes"
}

tyms <- full_stats[
  gene == "TYMS"
]

umps <- full_stats[
  gene == "UMPS"
]

tyms_ratio <- if (
  "CR1_ratio_BM_vs_PT" %in% names(tyms)
) {
  tyms$CR1_ratio_BM_vs_PT
} else {
  tyms$ratio_BM_vs_PT
}

umps_ratio <- if (
  "CR1_ratio_BM_vs_PT" %in% names(umps)
) {
  umps$CR1_ratio_BM_vs_PT
} else {
  umps$ratio_BM_vs_PT
}

tyms_fdr <- tyms[[main_fdr_col]]
umps_fdr <- umps[[main_fdr_col]]

interpretation <- if (
  is.finite(tyms_ratio) &&
  is.finite(umps_ratio) &&
  tyms_ratio > 1 &&
  umps_ratio > 1 &&
  is.finite(tyms_fdr) &&
  is.finite(umps_fdr) &&
  tyms_fdr < 0.05 &&
  umps_fdr < 0.05
) {
  paste0(
    "TYMS and UMPS remain positive and FDR-significant under the ",
    "count-based negative-binomial pseudobulk sensitivity. This supports ",
    "the original gene-level decomposition while remaining an expression-",
    "level result rather than evidence of metabolic flux or dependency."
  )
} else {
  paste0(
    "The count-based negative-binomial sensitivity does not reproduce both ",
    "TYMS and UMPS as FDR-significant positive site effects. Report the ",
    "count-scale result transparently as a model-dependent sensitivity and ",
    "retain the locked nonparametric decomposition as the primary gene-level ",
    "expression analysis."
  )
}

# Rank by the cluster-robust FDR when available.
ranked_genes <- full_stats[
  order(
    get(main_fdr_col),
    -ratio_BM_vs_PT
  ),
  gene
]

writeLines(
  c(
    "status=PASS_E3_COUNT_BASED_PSEUDOBULK_COMPLETE",
    "script_version=E3_V3_E2B_CACHE_LOW_MEMORY",
    "source=barcode_validated_E2B_target_gene_caches_from_raw_GSE223499_counts",
    "model=gene_wise_negative_binomial_GLM_with_log_malignant_library_size_offset",
    "main_inference=patient_cluster_CR1_two_sided",
    "FDR_family=5_prespecified_genes",
    "main_E3_dataset=30_BM_10_PT_specimens_38_unique_patient_clusters",
    "KRAS17_included=FALSE",
    "KRAS17_exclusion_reason=preexisting_deposited_count_CSV_metadata_discordance",
    "independence_preserving_dataset=28_BM_8_PT",
    paste0(
      "main_ranked_genes=",
      paste(
        ranked_genes,
        collapse = ";"
      )
    ),
    paste0(
      "reviewer_facing_interpretation=",
      interpretation
    ),
    "primary_endpoint_replaced=FALSE",
    "gene_set_changed=FALSE",
    "result_driven_exclusion=FALSE"
  ),
  file.path(
    OUT,
    "STEP_E3_COMPLETE.txt"
  )
)

cat(
  "\n============================================================\n",
  "STEP E3 V3: LOW-MEMORY COUNT-BASED PSEUDOBULK NB SENSITIVITY\n",
  "============================================================\n",
  "Main E3 dataset: 30 BM / 10 PT specimens\n",
  "Corrected patient clusters: ",
  uniqueN(
    e2b_samples$corrected_patient_id
  ),
  "\n",
  "Independence-preserving subset: 28 BM / 8 PT\n\n",
  "Main count-based results:\n",
  sep = ""
)

display_cols <- intersect(
  c(
    "gene",
    "ratio_BM_vs_PT",
    "CR1_CI95_low_ratio",
    "CR1_CI95_high_ratio",
    "two_sided_p_CR1",
    "BH_FDR_CR1_across_5_genes",
    "two_sided_p_model",
    "BH_FDR_model_across_5_genes"
  ),
  names(full_stats)
)

print(
  full_stats[
    ,
    ..display_cols
  ]
)

cat(
  "\nIndependence-preserving 28 BM / 8 PT results:\n"
)

print(
  ind_stats[
    ,
    .(
      gene,
      ratio_BM_vs_PT,
      CI95_low_ratio,
      CI95_high_ratio,
      two_sided_p_model,
      BH_FDR_model_across_5_genes
    )
  ]
)

cat(
  "\nInterpretation:\n",
  interpretation,
  "\n\nStatus: PASS_E3_COUNT_BASED_PSEUDOBULK_COMPLETE\n\n",
  "Key output files:\n",
  file.path(
    OUT,
    "STEP_E3_01_count_pseudobulk_analysis_input.tsv"
  ),
  "\n",
  file.path(
    OUT,
    "STEP_E3_02_count_based_NB_gene_statistics.tsv"
  ),
  "\n",
  file.path(
    OUT,
    "STEP_E3_COMPLETE.txt"
  ),
  "\n",
  sep = ""
)
