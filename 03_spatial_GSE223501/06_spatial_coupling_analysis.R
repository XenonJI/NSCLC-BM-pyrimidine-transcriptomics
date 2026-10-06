# Spatial pyrimidine-proliferation coupling
# Aggregates the fixed targets on the prespecified grid, computes raw and adjusted within-specimen
# associations, and treats each spatial specimen as the inferential unit. Invariant genes retain
# a zero z contribution so that fixed module denominators are preserved.

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/LUAD_LM_PYRIMIDINE"

INVENTORY_FILE <- file.path(
  ROOT,
  "results/reviewer_defense",
  "STEP_C1_01_GSE223501_spatial_file_inventory.tsv"
)

TARGET_DIR <- file.path(
  ROOT,
  "data/processed/GSE223501_spatial_targets"
)

DEPTH_DIR <- file.path(
  ROOT,
  "data/processed/GSE223501_spatial_library_depth"
)

OUT_DIR <- file.path(
  ROOT,
  "results/reviewer_defense"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

GRID_WIDTH <- 400
MIN_BEADS_PER_GRID <- 250L
MIN_ELIGIBLE_GRIDS <- 20L

PYRIMIDINE_GENES <- c(
  "TYMS",
  "UMPS"
)

PROLIFERATION_GENES <- c(
  "MKI67",
  "PCNA",
  "TOP2A",
  "UBE2C",
  "CENPF",
  "BIRC5",
  "CCNB1"
)

EPITHELIAL_GENES <- c(
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "MUC1"
)

TARGET_GENES <- c(
  PYRIMIDINE_GENES,
  PROLIFERATION_GENES,
  EPITHELIAL_GENES
)

banner <- function(x) {
  cat(
    "\n========================================\n",
    x,
    "\n========================================\n",
    sep = ""
  )
}

safe_stub <- function(x) {
  gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x
  )
}

zscore_keep_fixed_gene <- function(x) {
  s <- sd(x)

  if (
    is.na(s) ||
      s <= 0
  ) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }

  (x - mean(x)) /
    s
}

partial_spearman_rank_residual <- function(
  y,
  x,
  epithelial,
  log_library
) {
  dat <- data.frame(
    y_rank = rank(
      y,
      ties.method = "average"
    ),
    x_rank = rank(
      x,
      ties.method = "average"
    ),
    epithelial_rank = rank(
      epithelial,
      ties.method = "average"
    ),
    library_rank = rank(
      log_library,
      ties.method = "average"
    )
  )

  # If epithelial score is constant, including it would create
  # a redundant covariate. Keep the locked adjustment concept by
  # including it only when it has rank variation; library depth is
  # always included if variable.
  rhs <- character()

  if (
    sd(
      dat$epithelial_rank
    ) >
      0
  ) {
    rhs <- c(
      rhs,
      "epithelial_rank"
    )
  }

  if (
    sd(
      dat$library_rank
    ) >
      0
  ) {
    rhs <- c(
      rhs,
      "library_rank"
    )
  }

  if (length(rhs) == 0L) {
    return(
      cor(
        dat$y_rank,
        dat$x_rank,
        method = "pearson"
      )
    )
  }

  form_y <- as.formula(
    paste(
      "y_rank ~",
      paste(
        rhs,
        collapse = " + "
      )
    )
  )

  form_x <- as.formula(
    paste(
      "x_rank ~",
      paste(
        rhs,
        collapse = " + "
      )
    )
  )

  fit_y <- lm(
    form_y,
    data = dat
  )

  fit_x <- lm(
    form_x,
    data = dat
  )

  ry <- residuals(
    fit_y
  )

  rx <- residuals(
    fit_x
  )

  if (
    sd(ry) <= 0 ||
      sd(rx) <= 0
  ) {
    return(
      NA_real_
    )
  }

  cor(
    ry,
    rx,
    method = "pearson"
  )
}

exact_signed_rank_signflip <- function(x) {
  x <- as.numeric(x)
  x <- x[
    is.finite(x)
  ]

  x <- x[
    abs(x) >
      .Machine$double.eps^0.5
  ]

  n <- length(x)

  if (n == 0L) {
    return(
      list(
        n_nonzero = 0L,
        statistic_signed_rank = NA_real_,
        p_two_sided = NA_real_
      )
    )
  }

  if (n > 20L) {
    stop(
      "exact sign-flip仅设计用于n<=20。"
    )
  }

  ranks <- rank(
    abs(x),
    ties.method = "average"
  )

  obs <- sum(
    ranks *
      sign(x)
  )

  sign_grid <- as.matrix(
    expand.grid(
      rep(
        list(
          c(
            -1,
            1
          )
        ),
        n
      )
    )
  )

  perm_stats <- as.numeric(
    sign_grid %*%
      ranks
  )

  p_two <- mean(
    abs(
      perm_stats
    ) >=
      abs(obs) -
        1e-12
  )

  list(
    n_nonzero = n,
    statistic_signed_rank = obs,
    p_two_sided = p_two
  )
}

banner("STEP C5B v2：固定基因集 + invariant-gene neutral contribution")

if (!file.exists(INVENTORY_FILE)) {
  stop(
    "缺少STEP C1 inventory：\n",
    INVENTORY_FILE
  )
}

inv <- fread(
  INVENTORY_FILE,
  data.table = TRUE,
  showProgress = FALSE
)

if (nrow(inv) != 14L) {
  stop(
    "STEP C1 inventory不是14个样本。"
  )
}

# ------------------------------------------------------------
# 1. Write explicit v2 analysis lock BEFORE rerun
# ------------------------------------------------------------

writeLines(
  c(
    paste0(
      "created_at=",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    "version=C5B_v2",
    "reason_for_v2=PA025 UBE2C had zero within-specimen variance after locked aggregation",
    "fixed_gene_set_changed=FALSE",
    "invariant_gene_rule=z contribution set to 0; gene retained in fixed module denominator",
    "apply_rule_to_all_14_specimens=TRUE",
    "partial_v1_results_for_inference=DISCARD",
    "dataset=GSE223501",
    "modality=Slide-seq",
    "inferential_unit=spatial specimen",
    "n_spatial_specimens=14",
    paste0(
      "grid_width=",
      GRID_WIDTH
    ),
    paste0(
      "minimum_beads_per_eligible_grid=",
      MIN_BEADS_PER_GRID
    ),
    "grid_origin=per specimen min(x), min(y)",
    "normalization=log2(CPM+1) using summed all-gene grid library size",
    "standardization=within-specimen gene-wise z-score; invariant gene contributes 0",
    paste0(
      "pyrimidine_genes=",
      paste(
        PYRIMIDINE_GENES,
        collapse = ","
      )
    ),
    paste0(
      "proliferation_nonoverlap_genes=",
      paste(
        PROLIFERATION_GENES,
        collapse = ","
      )
    ),
    paste0(
      "epithelial_content_genes=",
      paste(
        EPITHELIAL_GENES,
        collapse = ","
      )
    ),
    "module_denominators_fixed=2,7,5",
    "unadjusted_association=within-specimen Spearman rho",
    "adjusted_association=rank-residual partial Spearman controlling epithelial score and log10(grid library size+1)",
    "across_specimen_test=exact two-sided signed-rank sign-flip test versus rho=0",
    "BM_vs_primary_expression_endpoint=NOT_TESTED",
    "grid_level_p_values=NOT_COMPUTED",
    "interpretation=orthogonal-modality support for proliferation coupling only"
  ),
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_00_ANALYSIS_LOCK.txt"
  )
)

grid_list <- vector(
  "list",
  nrow(inv)
)

sample_corr_list <- vector(
  "list",
  nrow(inv)
)

gene_qc_list <- vector(
  "list",
  nrow(inv)
)

for (i in seq_len(nrow(inv))) {

  sample_name <- inv$official_paper_id[i]
  cohort <- inv$cohort[i]

  target_file <- file.path(
    TARGET_DIR,
    paste0(
      safe_stub(
        sample_name
      ),
      "_target_counts_coords.rds"
    )
  )

  depth_file <- file.path(
    DEPTH_DIR,
    paste0(
      safe_stub(
        sample_name
      ),
      "_library_depth.rds"
    )
  )

  cat(
    "\n[",
    i,
    "/14] ",
    sample_name,
    " (",
    cohort,
    ")\n",
    sep = ""
  )

  if (!file.exists(target_file)) {
    stop(
      "缺少target RDS：\n",
      target_file
    )
  }

  if (!file.exists(depth_file)) {
    stop(
      "缺少depth RDS：\n",
      depth_file
    )
  }

  d <- readRDS(
    target_file
  )

  depth <- readRDS(
    depth_file
  )

  required_target <- c(
    "bead_key",
    "x",
    "y",
    TARGET_GENES
  )

  missing_target <- setdiff(
    required_target,
    names(d)
  )

  if (length(missing_target) > 0L) {
    stop(
      sample_name,
      " target RDS缺少字段：",
      paste(
        missing_target,
        collapse = ", "
      )
    )
  }

  if (
    nrow(d) !=
      nrow(depth)
  ) {
    stop(
      sample_name,
      " target/depth行数不一致。"
    )
  }

  if (!identical(
    as.character(
      d$bead_key
    ),
    as.character(
      depth$bead_key
    )
  )) {
    stop(
      sample_name,
      " target/depth bead_key不一致。"
    )
  }

  x_min <- min(
    d$x
  )

  y_min <- min(
    d$y
  )

  bead_dt <- data.table(
    gx = floor(
      (d$x - x_min) /
        GRID_WIDTH
    ),
    gy = floor(
      (d$y - y_min) /
        GRID_WIDTH
    ),
    library_size =
      depth$library_size
  )

  for (g in TARGET_GENES) {
    bead_dt[
      ,
      (g) :=
        d[[g]]
    ]
  }

  grid_dt <- bead_dt[
    ,
    c(
      list(
        n_beads = .N,
        grid_library_size =
          sum(
            library_size
          )
      ),
      lapply(
        .SD,
        sum
      )
    ),
    by = .(
      gx,
      gy
    ),
    .SDcols =
      TARGET_GENES
  ]

  all_grids <- nrow(
    grid_dt
  )

  eligible_dt <- grid_dt[
    n_beads >=
      MIN_BEADS_PER_GRID
  ]

  n_eligible <- nrow(
    eligible_dt
  )

  if (
    n_eligible <
      MIN_ELIGIBLE_GRIDS
  ) {
    stop(
      sample_name,
      " eligible grids不足：",
      n_eligible
    )
  }

  # Normalize target pseudobulk counts by ALL-GENE pseudobulk library size.
  for (g in TARGET_GENES) {
    eligible_dt[
      ,
      (paste0(
        "logCPM_",
        g
      )) :=
        log2(
          (
            get(g) /
              grid_library_size *
              1e6
          ) +
            1
        )
    ]
  }

  gene_qc <- rbindlist(
    lapply(
      TARGET_GENES,
      function(g) {
        x <- eligible_dt[[
          paste0(
            "logCPM_",
            g
          )
        ]]

        s <- sd(x)

        data.table(
          sample = sample_name,
          cohort = cohort,
          gene = g,
          n_eligible_grids =
            length(x),
          mean_logCPM =
            mean(x),
          sd_logCPM =
            s,
          invariant =
            is.na(s) ||
              s <= 0,
          fraction_zero_raw_grid_count =
            mean(
              eligible_dt[[g]] ==
                0
            ),
          total_raw_count_in_eligible_grids =
            sum(
              eligible_dt[[g]]
            )
        )
      }
    )
  )

  gene_qc_list[[i]] <- gene_qc

  for (g in TARGET_GENES) {
    eligible_dt[
      ,
      (paste0(
        "z_",
        g
      )) :=
        zscore_keep_fixed_gene(
          get(
            paste0(
              "logCPM_",
              g
            )
          )
        )
    ]
  }

  # Fixed denominators are preserved because invariant genes remain columns
  # whose contribution is exactly 0.
  eligible_dt[
    ,
    pyrimidine_score :=
      rowSums(
        .SD
      ) /
        length(
          PYRIMIDINE_GENES
        ),
    .SDcols =
      paste0(
        "z_",
        PYRIMIDINE_GENES
      )
  ]

  eligible_dt[
    ,
    proliferation_score :=
      rowSums(
        .SD
      ) /
        length(
          PROLIFERATION_GENES
        ),
    .SDcols =
      paste0(
        "z_",
        PROLIFERATION_GENES
      )
  ]

  eligible_dt[
    ,
    epithelial_score :=
      rowSums(
        .SD
      ) /
        length(
          EPITHELIAL_GENES
        ),
    .SDcols =
      paste0(
        "z_",
        EPITHELIAL_GENES
      )
  ]

  eligible_dt[
    ,
    log10_grid_library :=
      log10(
        grid_library_size +
          1
      )
  ]

  # Module-level variation must exist for correlation to be defined.
  if (
    sd(
      eligible_dt$pyrimidine_score
    ) <=
      0
  ) {
    stop(
      sample_name,
      " pyrimidine score无空间变异。"
    )
  }

  if (
    sd(
      eligible_dt$proliferation_score
    ) <=
      0
  ) {
    stop(
      sample_name,
      " proliferation score无空间变异。"
    )
  }

  rho_raw <- suppressWarnings(
    cor(
      eligible_dt$pyrimidine_score,
      eligible_dt$proliferation_score,
      method = "spearman",
      use = "complete.obs"
    )
  )

  rho_partial <- partial_spearman_rank_residual(
    y =
      eligible_dt$pyrimidine_score,
    x =
      eligible_dt$proliferation_score,
    epithelial =
      eligible_dt$epithelial_score,
    log_library =
      eligible_dt$log10_grid_library
  )

  invariant_pyr <- gene_qc[
    gene %in%
      PYRIMIDINE_GENES &
      invariant,
    gene
  ]

  invariant_prolif <- gene_qc[
    gene %in%
      PROLIFERATION_GENES &
      invariant,
    gene
  ]

  invariant_epi <- gene_qc[
    gene %in%
      EPITHELIAL_GENES &
      invariant,
    gene
  ]

  sample_corr_list[[i]] <- data.table(
    sample = sample_name,
    cohort = cohort,
    grid_width =
      GRID_WIDTH,
    min_beads_per_grid =
      MIN_BEADS_PER_GRID,
    all_nonempty_grids =
      all_grids,
    eligible_grids =
      n_eligible,
    eligible_beads =
      sum(
        eligible_dt$n_beads
      ),
    eligible_bead_fraction =
      sum(
        eligible_dt$n_beads
      ) /
        sum(
          grid_dt$n_beads
        ),
    invariant_pyrimidine_n =
      length(
        invariant_pyr
      ),
    invariant_pyrimidine_genes =
      if (
        length(
          invariant_pyr
        ) == 0L
      ) {
        ""
      } else {
        paste(
          invariant_pyr,
          collapse = ","
        )
      },
    invariant_proliferation_n =
      length(
        invariant_prolif
      ),
    invariant_proliferation_genes =
      if (
        length(
          invariant_prolif
        ) == 0L
      ) {
        ""
      } else {
        paste(
          invariant_prolif,
          collapse = ","
        )
      },
    invariant_epithelial_n =
      length(
        invariant_epi
      ),
    invariant_epithelial_genes =
      if (
        length(
          invariant_epi
        ) == 0L
      ) {
        ""
      } else {
        paste(
          invariant_epi,
          collapse = ","
        )
      },
    variable_pyrimidine_n =
      length(
        PYRIMIDINE_GENES
      ) -
        length(
          invariant_pyr
        ),
    variable_proliferation_n =
      length(
        PROLIFERATION_GENES
      ) -
        length(
          invariant_prolif
        ),
    variable_epithelial_n =
      length(
        EPITHELIAL_GENES
      ) -
        length(
          invariant_epi
        ),
    rho_raw_spearman =
      rho_raw,
    rho_partial_spearman =
      rho_partial
  )

  eligible_dt[
    ,
    `:=`(
      sample =
        sample_name,
      cohort =
        cohort,
      grid_width =
        GRID_WIDTH,
      grid_id =
        paste(
          gx,
          gy,
          sep = "_"
        )
    )
  ]

  keep_cols <- c(
    "sample",
    "cohort",
    "grid_width",
    "grid_id",
    "gx",
    "gy",
    "n_beads",
    "grid_library_size",
    "log10_grid_library",
    TARGET_GENES,
    paste0(
      "logCPM_",
      TARGET_GENES
    ),
    paste0(
      "z_",
      TARGET_GENES
    ),
    "pyrimidine_score",
    "proliferation_score",
    "epithelial_score"
  )

  grid_list[[i]] <- eligible_dt[
    ,
    ..keep_cols
  ]

  cat(
    "  eligible grids=",
    n_eligible,
    "/",
    all_grids,
    " | bead retention=",
    round(
      100 *
        sample_corr_list[[i]]$
          eligible_bead_fraction,
      1
    ),
    "%\n",
    "  invariant pyr/prolif/epi=",
    length(
      invariant_pyr
    ),
    "/",
    length(
      invariant_prolif
    ),
    "/",
    length(
      invariant_epi
    ),
    "\n",
    "  raw Spearman rho=",
    round(
      rho_raw,
      4
    ),
    "\n",
    "  adjusted partial rho=",
    round(
      rho_partial,
      4
    ),
    "\n",
    sep = ""
  )

  rm(
    d,
    depth,
    bead_dt,
    grid_dt,
    eligible_dt,
    gene_qc
  )

  gc(
    verbose = FALSE
  )
}

grid_scores <- rbindlist(
  grid_list,
  fill = TRUE
)

sample_corr <- rbindlist(
  sample_corr_list,
  fill = TRUE
)

gene_qc <- rbindlist(
  gene_qc_list,
  fill = TRUE
)

setorder(
  sample_corr,
  cohort,
  sample
)

setorder(
  gene_qc,
  cohort,
  sample,
  gene
)

# ------------------------------------------------------------
# 2. Specimen-level inference
# ------------------------------------------------------------

raw_test <- exact_signed_rank_signflip(
  sample_corr$rho_raw_spearman
)

partial_test <- exact_signed_rank_signflip(
  sample_corr$rho_partial_spearman
)

inference <- rbindlist(
  list(
    data.table(
      association = "raw_spearman",
      n_specimens =
        sum(
          is.finite(
            sample_corr$rho_raw_spearman
          )
        ),
      n_positive =
        sum(
          sample_corr$rho_raw_spearman >
            0,
          na.rm = TRUE
        ),
      fraction_positive =
        mean(
          sample_corr$rho_raw_spearman >
            0,
          na.rm = TRUE
        ),
      median_rho =
        median(
          sample_corr$rho_raw_spearman,
          na.rm = TRUE
        ),
      mean_rho =
        mean(
          sample_corr$rho_raw_spearman,
          na.rm = TRUE
        ),
      min_rho =
        min(
          sample_corr$rho_raw_spearman,
          na.rm = TRUE
        ),
      max_rho =
        max(
          sample_corr$rho_raw_spearman,
          na.rm = TRUE
        ),
      signed_rank_statistic =
        raw_test$
          statistic_signed_rank,
      exact_two_sided_p =
        raw_test$
          p_two_sided
    ),
    data.table(
      association = "partial_spearman_epithelial_and_library_adjusted",
      n_specimens =
        sum(
          is.finite(
            sample_corr$rho_partial_spearman
          )
        ),
      n_positive =
        sum(
          sample_corr$rho_partial_spearman >
            0,
          na.rm = TRUE
        ),
      fraction_positive =
        mean(
          sample_corr$rho_partial_spearman >
            0,
          na.rm = TRUE
        ),
      median_rho =
        median(
          sample_corr$rho_partial_spearman,
          na.rm = TRUE
        ),
      mean_rho =
        mean(
          sample_corr$rho_partial_spearman,
          na.rm = TRUE
        ),
      min_rho =
        min(
          sample_corr$rho_partial_spearman,
          na.rm = TRUE
        ),
      max_rho =
        max(
          sample_corr$rho_partial_spearman,
          na.rm = TRUE
        ),
      signed_rank_statistic =
        partial_test$
          statistic_signed_rank,
      exact_two_sided_p =
        partial_test$
          p_two_sided
    )
  )
)

cohort_descriptive <- sample_corr[
  ,
  .(
    n_specimens = .N,
    median_raw_rho =
      median(
        rho_raw_spearman,
        na.rm = TRUE
      ),
    positive_raw =
      sum(
        rho_raw_spearman >
          0,
        na.rm = TRUE
      ),
    median_partial_rho =
      median(
        rho_partial_spearman,
        na.rm = TRUE
      ),
    positive_partial =
      sum(
        rho_partial_spearman >
          0,
        na.rm = TRUE
      ),
    specimens_with_any_invariant_pyrimidine =
      sum(
        invariant_pyrimidine_n >
          0
      ),
    specimens_with_any_invariant_proliferation =
      sum(
        invariant_proliferation_n >
          0
      ),
    specimens_with_any_invariant_epithelial =
      sum(
        invariant_epithelial_n >
          0
      )
  ),
  by = cohort
]

# ------------------------------------------------------------
# 3. Outputs
# ------------------------------------------------------------

fwrite(
  grid_scores,
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_01_grid_level_scores.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_qc,
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_02_gene_variability_QC.tsv"
  ),
  sep = "\t"
)

fwrite(
  sample_corr,
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_03_sample_level_spatial_correlations.tsv"
  ),
  sep = "\t"
)

fwrite(
  inference,
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_04_sample_level_inference.tsv"
  ),
  sep = "\t"
)

fwrite(
  cohort_descriptive,
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_05_cohort_descriptive_correlations.tsv"
  ),
  sep = "\t"
)

status <- if (
  nrow(sample_corr) ==
    14L &&
    all(
      is.finite(
        sample_corr$
          rho_raw_spearman
      )
    ) &&
    all(
      is.finite(
        sample_corr$
          rho_partial_spearman
      )
    )
) {
  "PASS_SPATIAL_ORTHOGONAL_VALIDATION_ANALYSIS_V2"
} else {
  "REVIEW_REQUIRED"
}

writeLines(
  c(
    paste0(
      "completed_at=",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "status=",
      status
    ),
    "fixed_gene_set_changed=FALSE",
    "invariant_gene_rule=zero neutral z contribution with fixed denominator",
    "partial_v1_results_used_for_inference=FALSE",
    paste0(
      "specimens=",
      nrow(
        sample_corr
      )
    ),
    paste0(
      "specimens_with_invariant_pyrimidine=",
      sum(
        sample_corr$
          invariant_pyrimidine_n >
          0
      )
    ),
    paste0(
      "specimens_with_invariant_proliferation=",
      sum(
        sample_corr$
          invariant_proliferation_n >
          0
      )
    ),
    paste0(
      "specimens_with_invariant_epithelial=",
      sum(
        sample_corr$
          invariant_epithelial_n >
          0
      )
    ),
    paste0(
      "raw_median_rho=",
      signif(
        inference[
          association ==
            "raw_spearman",
          median_rho
        ],
        8
      )
    ),
    paste0(
      "raw_exact_two_sided_p=",
      signif(
        inference[
          association ==
            "raw_spearman",
          exact_two_sided_p
        ],
        8
      )
    ),
    paste0(
      "adjusted_median_partial_rho=",
      signif(
        inference[
          association ==
            "partial_spearman_epithelial_and_library_adjusted",
          median_rho
        ],
        8
      )
    ),
    paste0(
      "adjusted_exact_two_sided_p=",
      signif(
        inference[
          association ==
            "partial_spearman_epithelial_and_library_adjusted",
          exact_two_sided_p
        ],
        8
      )
    ),
    "inferential_unit=spatial specimen",
    "grid_level_p_values=FALSE",
    "BM_vs_primary_inferential_comparison=FALSE",
    "interpretation_boundary=orthogonal support for proliferation coupling only"
  ),
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_COMPLETE.txt"
  )
)

banner("STEP C5B v2结果")

cat(
  "\nSample-level spatial correlations:\n"
)

print(
  sample_corr
)

cat(
  "\nAcross-specimen inference:\n"
)

print(
  inference
)

cat(
  "\nCohort descriptive summary:\n"
)

print(
  cohort_descriptive
)

cat(
  "\nStatus: ",
  status,
  "\n\n",
  "Key output files:\n",
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_02_gene_variability_QC.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_03_sample_level_spatial_correlations.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_04_sample_level_inference.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C5Bv2_05_cohort_descriptive_correlations.tsv"
  ),
  "\n",
  sep = ""
)
