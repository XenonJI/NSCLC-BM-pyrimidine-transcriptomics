# Spatial aggregation-scale sensitivity
# Repeats the prespecified spatial association calculations at adjacent coordinate-only grid
# widths while keeping the target genes, normalization, and specimen-level inference unchanged.

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

GRID_WIDTHS <- c(
  300,
  400,
  500
)

PRIMARY_GRID_WIDTH <- 400
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

  rhs <- character()

  if (
    sd(
      dat$epithelial_rank
    ) > 0
  ) {
    rhs <- c(
      rhs,
      "epithelial_rank"
    )
  }

  if (
    sd(
      dat$library_rank
    ) > 0
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

  fit_y <- lm(
    as.formula(
      paste(
        "y_rank ~",
        paste(
          rhs,
          collapse = " + "
        )
      )
    ),
    data = dat
  )

  fit_x <- lm(
    as.formula(
      paste(
        "x_rank ~",
        paste(
          rhs,
          collapse = " + "
        )
      )
    ),
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
        statistic = NA_real_,
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

  list(
    n_nonzero = n,
    statistic = obs,
    p_two_sided = mean(
      abs(
        perm_stats
      ) >=
        abs(obs) -
          1e-12
    )
  )
}

banner("STEP C5C：空间聚合尺度敏感性分析")

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

writeLines(
  c(
    paste0(
      "created_at=",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    "analysis_type=reviewer-defense spatial aggregation scale sensitivity",
    paste0(
      "primary_locked_grid_width=",
      PRIMARY_GRID_WIDTH
    ),
    paste0(
      "sensitivity_grid_widths=",
      paste(
        setdiff(
          GRID_WIDTHS,
          PRIMARY_GRID_WIDTH
        ),
        collapse = ","
      )
    ),
    paste0(
      "all_evaluated_grid_widths=",
      paste(
        GRID_WIDTHS,
        collapse = ","
      )
    ),
    paste0(
      "minimum_beads_per_grid=",
      MIN_BEADS_PER_GRID
    ),
    "candidate_widths_identified_before_expression_correlations=TRUE",
    "fixed_gene_set_changed=FALSE",
    "invariant_gene_rule=zero neutral z contribution with fixed denominator",
    "inferential_unit=spatial specimen",
    "BM_vs_primary_inferential_comparison=FALSE"
  ),
  file.path(
    OUT_DIR,
    "STEP_C5C_00_SCALE_SENSITIVITY_LOCK.txt"
  )
)

results_list <- list()
idx <- 0L

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

  cat(
    "\n",
    sample_name,
    " (",
    cohort,
    ")\n",
    sep = ""
  )

  for (w in GRID_WIDTHS) {

    bead_dt <- data.table(
      gx = floor(
        (d$x - x_min) /
          w
      ),
      gy = floor(
        (d$y - y_min) /
          w
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

    eligible_dt <- grid_dt[
      n_beads >=
        MIN_BEADS_PER_GRID
    ]

    if (
      nrow(
        eligible_dt
      ) <
        MIN_ELIGIBLE_GRIDS
    ) {
      stop(
        sample_name,
        " width=",
        w,
        " eligible grids不足。"
      )
    }

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

    invariant_pyr <- sum(
      vapply(
        PYRIMIDINE_GENES,
        function(g) {
          sd(
            eligible_dt[[
              paste0(
                "logCPM_",
                g
              )
            ]]
          ) <=
            0
        },
        logical(1)
      )
    )

    invariant_prolif <- sum(
      vapply(
        PROLIFERATION_GENES,
        function(g) {
          sd(
            eligible_dt[[
              paste0(
                "logCPM_",
                g
              )
            ]]
          ) <=
            0
        },
        logical(1)
      )
    )

    invariant_epi <- sum(
      vapply(
        EPITHELIAL_GENES,
        function(g) {
          sd(
            eligible_dt[[
              paste0(
                "logCPM_",
                g
              )
            ]]
          ) <=
            0
        },
        logical(1)
      )
    )

    idx <- idx + 1L

    results_list[[idx]] <- data.table(
      sample = sample_name,
      cohort = cohort,
      grid_width = w,
      primary_locked_width =
        w ==
          PRIMARY_GRID_WIDTH,
      nonempty_grids =
        nrow(
          grid_dt
        ),
      eligible_grids =
        nrow(
          eligible_dt
        ),
      eligible_bead_fraction =
        sum(
          eligible_dt$n_beads
        ) /
          sum(
            grid_dt$n_beads
          ),
      invariant_pyrimidine_n =
        invariant_pyr,
      invariant_proliferation_n =
        invariant_prolif,
      invariant_epithelial_n =
        invariant_epi,
      rho_raw_spearman =
        rho_raw,
      rho_partial_spearman =
        rho_partial
    )

    cat(
      "  width=",
      w,
      " | grids=",
      nrow(
        eligible_dt
      ),
      " | raw rho=",
      round(
        rho_raw,
        4
      ),
      " | adjusted rho=",
      round(
        rho_partial,
        4
      ),
      "\n",
      sep = ""
    )

    rm(
      bead_dt,
      grid_dt,
      eligible_dt
    )

    gc(
      verbose = FALSE
    )
  }

  rm(
    d,
    depth
  )

  gc(
    verbose = FALSE
  )
}

sample_results <- rbindlist(
  results_list,
  fill = TRUE
)

setorder(
  sample_results,
  grid_width,
  cohort,
  sample
)

# ------------------------------------------------------------
# Width-level specimen inference
# ------------------------------------------------------------

inference_list <- list()
j <- 0L

for (w in GRID_WIDTHS) {

  x <- sample_results[
    grid_width ==
      w
  ]

  raw_test <- exact_signed_rank_signflip(
    x$rho_raw_spearman
  )

  partial_test <- exact_signed_rank_signflip(
    x$rho_partial_spearman
  )

  j <- j + 1L
  inference_list[[j]] <- data.table(
    grid_width = w,
    primary_locked_width =
      w ==
        PRIMARY_GRID_WIDTH,
    association = "raw_spearman",
    n_specimens = nrow(x),
    n_positive =
      sum(
        x$rho_raw_spearman >
          0
      ),
    fraction_positive =
      mean(
        x$rho_raw_spearman >
          0
      ),
    median_rho =
      median(
        x$rho_raw_spearman
      ),
    mean_rho =
      mean(
        x$rho_raw_spearman
      ),
    min_rho =
      min(
        x$rho_raw_spearman
      ),
    max_rho =
      max(
        x$rho_raw_spearman
      ),
    exact_two_sided_p =
      raw_test$
        p_two_sided
  )

  j <- j + 1L
  inference_list[[j]] <- data.table(
    grid_width = w,
    primary_locked_width =
      w ==
        PRIMARY_GRID_WIDTH,
    association = "partial_spearman_epithelial_and_library_adjusted",
    n_specimens = nrow(x),
    n_positive =
      sum(
        x$rho_partial_spearman >
          0
      ),
    fraction_positive =
      mean(
        x$rho_partial_spearman >
          0
      ),
    median_rho =
      median(
        x$rho_partial_spearman
      ),
    mean_rho =
      mean(
        x$rho_partial_spearman
      ),
    min_rho =
      min(
        x$rho_partial_spearman
      ),
    max_rho =
      max(
        x$rho_partial_spearman
      ),
    exact_two_sided_p =
      partial_test$
        p_two_sided
  )
}

inference <- rbindlist(
  inference_list,
  fill = TRUE
)

setorder(
  inference,
  grid_width,
  association
)

# ------------------------------------------------------------
# Per-specimen sign robustness across widths
# ------------------------------------------------------------

sign_robustness <- sample_results[
  ,
  .(
    cohort = cohort[1],
    raw_positive_all_3 =
      all(
        rho_raw_spearman >
          0
      ),
    adjusted_positive_all_3 =
      all(
        rho_partial_spearman >
          0
      ),
    raw_min_rho =
      min(
        rho_raw_spearman
      ),
    raw_max_rho =
      max(
        rho_raw_spearman
      ),
    adjusted_min_rho =
      min(
        rho_partial_spearman
      ),
    adjusted_max_rho =
      max(
        rho_partial_spearman
      )
  ),
  by = sample
]

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------

fwrite(
  sample_results,
  file.path(
    OUT_DIR,
    "STEP_C5C_01_scale_sensitivity_sample_correlations.tsv"
  ),
  sep = "\t"
)

fwrite(
  inference,
  file.path(
    OUT_DIR,
    "STEP_C5C_02_scale_sensitivity_inference.tsv"
  ),
  sep = "\t"
)

fwrite(
  sign_robustness,
  file.path(
    OUT_DIR,
    "STEP_C5C_03_per_specimen_sign_robustness.tsv"
  ),
  sep = "\t"
)

status <- if (
  nrow(sample_results) ==
    14L *
      length(
        GRID_WIDTHS
      ) &&
    all(
      is.finite(
        sample_results$
          rho_raw_spearman
      )
    ) &&
    all(
      is.finite(
        sample_results$
          rho_partial_spearman
      )
    )
) {
  "PASS_SPATIAL_SCALE_SENSITIVITY"
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
    paste0(
      "primary_locked_grid_width=",
      PRIMARY_GRID_WIDTH
    ),
    paste0(
      "evaluated_grid_widths=",
      paste(
        GRID_WIDTHS,
        collapse = ","
      )
    ),
    paste0(
      "minimum_beads_per_grid=",
      MIN_BEADS_PER_GRID
    ),
    paste0(
      "specimens_raw_positive_all_3_widths=",
      sum(
        sign_robustness$
          raw_positive_all_3
      ),
      "/14"
    ),
    paste0(
      "specimens_adjusted_positive_all_3_widths=",
      sum(
        sign_robustness$
          adjusted_positive_all_3
      ),
      "/14"
    ),
    "primary_inference_remains_STEP_C5Bv2_width400=TRUE",
    "sensitivity_results_do_not_replace_primary=TRUE"
  ),
  file.path(
    OUT_DIR,
    "STEP_C5C_COMPLETE.txt"
  )
)

banner("STEP C5C结果")

cat(
  "\nWidth-level inference:\n"
)

print(
  inference
)

cat(
  "\nPer-specimen sign robustness:\n"
)

print(
  sign_robustness
)

cat(
  "\nStatus: ",
  status,
  "\n\n",
  "Key output files:\n",
  file.path(
    OUT_DIR,
    "STEP_C5C_02_scale_sensitivity_inference.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C5C_03_per_specimen_sign_robustness.tsv"
  ),
  "\n",
  sep = ""
)
