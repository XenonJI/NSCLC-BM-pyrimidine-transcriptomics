# Spatial robustness analyses
# Evaluates duplicate-payload handling, an epithelial-enriched proxy subset, toroidal-shift
# autocorrelation-preserving nulls, and manuscript spatial maps. GSE223501 is used as same-study
# spatial context rather than an independent validation cohort.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/LUAD_LM_PYRIMIDINE"
RD <- file.path(ROOT, "results", "reviewer_defense")
OUT <- file.path(ROOT, "results", "major_revision")
FIGDIR <- file.path(OUT, "spatial_maps_E5")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

GRID_FILE <- file.path(
  RD,
  "STEP_C5Bv2_01_grid_level_scores.tsv"
)

CORR_FILE <- file.path(
  RD,
  "STEP_C5Bv2_03_sample_level_spatial_correlations.tsv"
)

for (f in c(GRID_FILE, CORR_FILE)) {
  if (!file.exists(f)) {
    stop("Required prior spatial file not found:\n", f, call. = FALSE)
  }
}

DUP_A <- "PA056"
DUP_B <- "KRAS_11"
PRIMARY_WIDTH <- 400L
MIN_OVERLAP_FRACTION <- 0.50
MIN_OVERLAP_GRIDS <- 20L
EPITHELIAL_QUANTILE <- 0.75

# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------
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
    stop("Exact sign-flip helper is limited to n<=20.", call. = FALSE)
  }

  r <- rank(
    abs(x),
    ties.method = "average"
  )

  obs <- sum(
    r *
      sign(x)
  )

  signs <- as.matrix(
    expand.grid(
      rep(
        list(
          c(-1, 1)
        ),
        n
      )
    )
  )

  null <- as.numeric(
    signs %*%
      r
  )

  list(
    n_nonzero = n,
    statistic = obs,
    p_two_sided = mean(
      abs(null) >=
        abs(obs) -
        1e-12
    )
  )
}

rank_residual_pair <- function(
  y,
  x,
  epithelial,
  log_library
) {
  z <- data.frame(
    y_rank = rank(y, ties.method = "average"),
    x_rank = rank(x, ties.method = "average"),
    epithelial_rank = rank(epithelial, ties.method = "average"),
    library_rank = rank(log_library, ties.method = "average")
  )

  rhs <- character(0)

  if (
    is.finite(sd(z$epithelial_rank)) &&
    sd(z$epithelial_rank) > 0
  ) {
    rhs <- c(rhs, "epithelial_rank")
  }

  if (
    is.finite(sd(z$library_rank)) &&
    sd(z$library_rank) > 0
  ) {
    rhs <- c(rhs, "library_rank")
  }

  if (length(rhs) == 0L) {
    return(
      list(
        y_resid = z$y_rank - mean(z$y_rank),
        x_resid = z$x_rank - mean(z$x_rank)
      )
    )
  }

  fy <- as.formula(
    paste(
      "y_rank ~",
      paste(rhs, collapse = " + ")
    )
  )

  fx <- as.formula(
    paste(
      "x_rank ~",
      paste(rhs, collapse = " + ")
    )
  )

  list(
    y_resid = residuals(lm(fy, data = z)),
    x_resid = residuals(lm(fx, data = z))
  )
}

partial_spearman <- function(
  y,
  x,
  epithelial,
  log_library
) {
  rr <- rank_residual_pair(
    y = y,
    x = x,
    epithelial = epithelial,
    log_library = log_library
  )

  if (
    sd(rr$y_resid) <= 0 ||
    sd(rr$x_resid) <= 0
  ) {
    return(NA_real_)
  }

  cor(
    rr$y_resid,
    rr$x_resid,
    method = "pearson"
  )
}

summarize_rhos <- function(
  x,
  scenario,
  association
) {
  x <- x[
    is.finite(x)
  ]

  tst <- exact_signed_rank_signflip(x)

  data.table(
    scenario = scenario,
    association = association,
    n_specimens = length(x),
    n_positive = sum(x > 0),
    fraction_positive = mean(x > 0),
    median_rho = median(x),
    mean_rho = mean(x),
    min_rho = min(x),
    max_rho = max(x),
    signed_rank_statistic = tst$statistic,
    exact_two_sided_p = tst$p_two_sided
  )
}

toroidal_shift_null <- function(
  d,
  value_y,
  value_x,
  min_overlap_fraction = 0.50,
  min_overlap_grids = 20L
) {
  d <- copy(d)

  xmin <- min(d$gx)
  xmax <- max(d$gx)
  ymin <- min(d$gy)
  ymax <- max(d$gy)

  nx <- xmax - xmin + 1L
  ny <- ymax - ymin + 1L

  if (nx < 2L || ny < 2L) {
    stop("Spatial extent too small for toroidal shifts.", call. = FALSE)
  }

  key_y <- d[
    ,
    .(
      gx,
      gy,
      y_value = get(value_y)
    )
  ]

  base_x <- d[
    ,
    .(
      gx,
      gy,
      x_value = get(value_x)
    )
  ]

  observed <- cor(
    d[[value_y]],
    d[[value_x]],
    method = "pearson",
    use = "complete.obs"
  )

  out <- list()
  k <- 0L

  for (dx in 0:(nx - 1L)) {
    for (dy in 0:(ny - 1L)) {
      if (dx == 0L && dy == 0L) next

      shifted <- copy(base_x)

      shifted[
        ,
        gx := (
          (gx - xmin + dx) %% nx
        ) + xmin
      ]

      shifted[
        ,
        gy := (
          (gy - ymin + dy) %% ny
        ) + ymin
      ]

      m <- merge(
        key_y,
        shifted,
        by = c("gx", "gy"),
        all = FALSE,
        sort = FALSE
      )

      min_n <- max(
        min_overlap_grids,
        ceiling(
          min_overlap_fraction *
            nrow(d)
        )
      )

      if (nrow(m) < min_n) next

      rho <- suppressWarnings(
        cor(
          m$y_value,
          m$x_value,
          method = "pearson",
          use = "complete.obs"
        )
      )

      if (!is.finite(rho)) next

      k <- k + 1L

      out[[k]] <- data.table(
        dx = dx,
        dy = dy,
        overlap_grids = nrow(m),
        overlap_fraction = nrow(m) / nrow(d),
        shifted_rho = rho
      )
    }
  }

  null <- if (length(out) > 0L) {
    rbindlist(out)
  } else {
    data.table()
  }

  if (nrow(null) == 0L) {
    return(
      list(
        observed = observed,
        null = null,
        p_two_sided = NA_real_,
        null_median = NA_real_,
        null_q025 = NA_real_,
        null_q975 = NA_real_
      )
    )
  }

  p2 <- (
    1 +
      sum(
        abs(null$shifted_rho) >=
          abs(observed) -
          1e-12
      )
  ) /
    (
      nrow(null) +
        1
    )

  list(
    observed = observed,
    null = null,
    p_two_sided = p2,
    null_median = median(null$shifted_rho),
    null_q025 = quantile(null$shifted_rho, 0.025),
    null_q975 = quantile(null$shifted_rho, 0.975)
  )
}

map_colors <- function(z, lim = 2.5) {
  pal <- hcl.colors(
    101,
    palette = "Blue-Red 2"
  )

  zz <- pmax(
    -lim,
    pmin(
      lim,
      z
    )
  )

  ii <- floor(
    (zz + lim) /
      (2 * lim) *
      100
  ) + 1L

  pal[
    pmax(
      1L,
      pmin(
        101L,
        ii
      )
    )
  ]
}

plot_spatial_panel <- function(
  d,
  score_col,
  main,
  lim = 2.5
) {
  z <- d[[score_col]]

  plot(
    d$gx,
    d$gy,
    pch = 15,
    cex = 1.9,
    col = map_colors(z, lim = lim),
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    family = "sans"
  )

  box(
    bty = "l"
  )
}

# ----------------------------------------------------------------------
# 2. Load locked grid-level spatial scores
# ----------------------------------------------------------------------
grid <- fread(
  GRID_FILE,
  showProgress = FALSE
)

corr <- fread(
  CORR_FILE,
  showProgress = FALSE
)

required_grid <- c(
  "sample",
  "cohort",
  "grid_width",
  "gx",
  "gy",
  "pyrimidine_score",
  "proliferation_score",
  "epithelial_score",
  "log10_grid_library"
)

if (!all(required_grid %in% names(grid))) {
  stop("Grid-level spatial file lacks required columns.", call. = FALSE)
}

if (
  unique(grid$grid_width) != PRIMARY_WIDTH
) {
  stop(
    "E5 expects the locked width=400 grid-level score file.",
    call. = FALSE
  )
}

# Conservative corrected set.
retain12 <- setdiff(
  unique(grid$sample),
  c(DUP_A, DUP_B)
)

grid12 <- grid[
  sample %chin% retain12
]

corr12 <- corr[
  sample %chin% retain12
]

if (
  uniqueN(corr12$sample) != 12L
) {
  stop("Expected conservative n=12 spatial specimens.", call. = FALSE)
}

# ----------------------------------------------------------------------
# 3. Duplicate-payload retention sensitivity
# ----------------------------------------------------------------------
scenarios <- list(
  exclude_both_primary_n12 =
    corr[
      !sample %chin% c(DUP_A, DUP_B)
    ],
  retain_PA056_only_n13 =
    corr[
      sample != DUP_B
    ],
  retain_KRAS11_only_n13 =
    corr[
      sample != DUP_A
    ],
  retain_both_descriptive_n14 =
    corr
)

dup_stats <- rbindlist(
  lapply(
    names(scenarios),
    function(nm) {
      dd <- scenarios[[nm]]

      rbindlist(
        list(
          summarize_rhos(
            dd$rho_raw_spearman,
            scenario = nm,
            association = "raw_spearman"
          ),
          summarize_rhos(
            dd$rho_partial_spearman,
            scenario = nm,
            association = "partial_spearman_epithelial_and_library_adjusted"
          )
        )
      )
    }
  ),
  fill = TRUE
)

fwrite(
  dup_stats,
  file.path(
    OUT,
    "STEP_E5_01_duplicate_payload_sensitivity.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 4. Epithelial-enriched proxy subset
# ----------------------------------------------------------------------
epi_rows <- list()

for (ss in sort(retain12)) {
  d <- grid12[
    sample == ss
  ]

  q <- quantile(
    d$epithelial_score,
    EPITHELIAL_QUANTILE,
    na.rm = TRUE
  )

  de <- d[
    epithelial_score >= q
  ]

  if (nrow(de) < MIN_OVERLAP_GRIDS) {
    # Include ties / next-highest bins deterministically until at least 20.
    setorder(
      d,
      -epithelial_score,
      gx,
      gy
    )

    de <- d[
      seq_len(
        min(
          max(
            MIN_OVERLAP_GRIDS,
            ceiling(
              (1 - EPITHELIAL_QUANTILE) *
                nrow(d)
            )
          ),
          nrow(d)
        )
      )
    ]
  }

  raw_rho <- suppressWarnings(
    cor(
      de$pyrimidine_score,
      de$proliferation_score,
      method = "spearman",
      use = "complete.obs"
    )
  )

  adj_rho <- partial_spearman(
    y = de$pyrimidine_score,
    x = de$proliferation_score,
    epithelial = de$epithelial_score,
    log_library = de$log10_grid_library
  )

  epi_rows[[length(epi_rows) + 1L]] <- data.table(
    sample = ss,
    cohort = d$cohort[1],
    all_eligible_grids = nrow(d),
    epithelial_enriched_grids = nrow(de),
    epithelial_quantile_rule = EPITHELIAL_QUANTILE,
    epithelial_threshold = q,
    rho_raw_spearman = raw_rho,
    rho_partial_spearman = adj_rho
  )
}

epi_corr <- rbindlist(
  epi_rows
)

fwrite(
  epi_corr,
  file.path(
    OUT,
    "STEP_E5_02_epithelial_enriched_proxy_correlations.tsv"
  ),
  sep = "\t"
)

epi_inference <- rbindlist(
  list(
    summarize_rhos(
      epi_corr$rho_raw_spearman,
      scenario = "epithelial_enriched_top_quartile_proxy_n12",
      association = "raw_spearman"
    ),
    summarize_rhos(
      epi_corr$rho_partial_spearman,
      scenario = "epithelial_enriched_top_quartile_proxy_n12",
      association = "partial_spearman_epithelial_and_library_adjusted"
    )
  )
)

fwrite(
  epi_inference,
  file.path(
    OUT,
    "STEP_E5_03_epithelial_enriched_proxy_inference.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 5. Toroidal-shift spatial null
# ----------------------------------------------------------------------
shift_summary <- list()
shift_all <- list()

for (ii in seq_along(sort(retain12))) {
  ss <- sort(retain12)[ii]

  cat(
    "[E5 spatial shift ",
    ii,
    "/12] ",
    ss,
    "\n",
    sep = ""
  )

  d <- copy(
    grid12[
      sample == ss
    ]
  )

  # Raw Spearman = Pearson correlation of within-specimen ranks.
  d[
    ,
    raw_y_rank :=
      rank(
        pyrimidine_score,
        ties.method = "average"
      )
  ]

  d[
    ,
    raw_x_rank :=
      rank(
        proliferation_score,
        ties.method = "average"
      )
  ]

  rr <- rank_residual_pair(
    y = d$pyrimidine_score,
    x = d$proliferation_score,
    epithelial = d$epithelial_score,
    log_library = d$log10_grid_library
  )

  d[
    ,
    adjusted_y_resid :=
      rr$y_resid
  ]

  d[
    ,
    adjusted_x_resid :=
      rr$x_resid
  ]

  raw_null <- toroidal_shift_null(
    d,
    value_y = "raw_y_rank",
    value_x = "raw_x_rank",
    min_overlap_fraction = MIN_OVERLAP_FRACTION,
    min_overlap_grids = MIN_OVERLAP_GRIDS
  )

  adj_null <- toroidal_shift_null(
    d,
    value_y = "adjusted_y_resid",
    value_x = "adjusted_x_resid",
    min_overlap_fraction = MIN_OVERLAP_FRACTION,
    min_overlap_grids = MIN_OVERLAP_GRIDS
  )

  shift_summary[[length(shift_summary) + 1L]] <- rbindlist(
    list(
      data.table(
        sample = ss,
        cohort = d$cohort[1],
        association = "raw_spearman_rank_surface",
        eligible_grids = nrow(d),
        admissible_nonzero_shifts = nrow(raw_null$null),
        observed_rho = raw_null$observed,
        toroidal_null_median = raw_null$null_median,
        toroidal_null_q025 = raw_null$null_q025,
        toroidal_null_q975 = raw_null$null_q975,
        toroidal_two_sided_p = raw_null$p_two_sided
      ),
      data.table(
        sample = ss,
        cohort = d$cohort[1],
        association = "partial_rank_residual_surface",
        eligible_grids = nrow(d),
        admissible_nonzero_shifts = nrow(adj_null$null),
        observed_rho = adj_null$observed,
        toroidal_null_median = adj_null$null_median,
        toroidal_null_q025 = adj_null$null_q025,
        toroidal_null_q975 = adj_null$null_q975,
        toroidal_two_sided_p = adj_null$p_two_sided
      )
    )
  )

  if (nrow(raw_null$null) > 0L) {
    x <- copy(raw_null$null)
    x[
      ,
      `:=`(
        sample = ss,
        cohort = d$cohort[1],
        association = "raw_spearman_rank_surface"
      )
    ]
    shift_all[[length(shift_all) + 1L]] <- x
  }

  if (nrow(adj_null$null) > 0L) {
    x <- copy(adj_null$null)
    x[
      ,
      `:=`(
        sample = ss,
        cohort = d$cohort[1],
        association = "partial_rank_residual_surface"
      )
    ]
    shift_all[[length(shift_all) + 1L]] <- x
  }
}

shift_summary <- rbindlist(
  shift_summary,
  fill = TRUE
)

shift_all <- rbindlist(
  shift_all,
  fill = TRUE
)

fwrite(
  shift_summary,
  file.path(
    OUT,
    "STEP_E5_04_toroidal_shift_null_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  shift_all,
  file.path(
    OUT,
    "STEP_E5_05_toroidal_shift_null_all_shifts.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 6. Coordinate-scale transparency
# ----------------------------------------------------------------------
scale_audit <- grid12[
  ,
  .(
    grid_width_deposited_coordinate_units = unique(grid_width)[1],
    n_eligible_grids = .N,
    gx_min = min(gx),
    gx_max = max(gx),
    gy_min = min(gy),
    gy_max = max(gy),
    approximate_x_span_coordinate_units =
      (max(gx) - min(gx) + 1L) *
        PRIMARY_WIDTH,
    approximate_y_span_coordinate_units =
      (max(gy) - min(gy) + 1L) *
        PRIMARY_WIDTH
  ),
  by = .(
    sample,
    cohort
  )
]

scale_audit[
  ,
  physical_unit_statement :=
    "400 refers to deposited coordinate units used for grid aggregation; no micron conversion is asserted without an explicit calibration source"
]

fwrite(
  scale_audit,
  file.path(
    OUT,
    "STEP_E5_06_coordinate_scale_transparency.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 7. Spatial map atlas: all conservative n=12 specimens
# ----------------------------------------------------------------------
atlas_pdf <- file.path(
  FIGDIR,
  "STEP_E5_spatial_map_atlas_n12.pdf"
)

pdf(
  atlas_pdf,
  width = 9.2,
  height = 3.4,
  family = "sans",
  useDingbats = FALSE
)

op <- par(
  mfrow = c(1, 3),
  mar = c(1.0, 1.0, 2.2, 0.5),
  oma = c(0, 0, 1.2, 0)
)

for (ss in sort(retain12)) {
  d <- grid12[
    sample == ss
  ]

  plot_spatial_panel(
    d,
    "pyrimidine_score",
    "TYMS–UMPS"
  )

  plot_spatial_panel(
    d,
    "proliferation_score",
    "Proliferation"
  )

  plot_spatial_panel(
    d,
    "epithelial_score",
    "Epithelial content"
  )

  mtext(
    paste0(
      ss,
      " | ",
      gsub("_", " ", d$cohort[1]),
      " | width=400 coordinate units"
    ),
    outer = TRUE,
    side = 3,
    line = 0.1,
    cex = 0.85,
    font = 2
  )
}

par(op)
dev.off()

# Representative pair is chosen by eligible-grid count ONLY.
rep_candidates <- corr12[
  ,
  .SD[
    which.max(
      eligible_grids
    )
  ],
  by = cohort
]

rep_samples <- rep_candidates$sample

rep_rule <- data.table(
  cohort = rep_candidates$cohort,
  sample = rep_candidates$sample,
  eligible_grids = rep_candidates$eligible_grids,
  selection_rule =
    "largest eligible-grid count within cohort; spatial correlation magnitude not used"
)

fwrite(
  rep_rule,
  file.path(
    OUT,
    "STEP_E5_07_representative_map_selection.tsv"
  ),
  sep = "\t"
)

rep_pdf <- file.path(
  FIGDIR,
  "STEP_E5_representative_spatial_maps.pdf"
)

pdf(
  rep_pdf,
  width = 9.0,
  height = 5.8,
  family = "sans",
  useDingbats = FALSE
)

op <- par(
  mfrow = c(2, 3),
  mar = c(1.0, 1.0, 2.2, 0.5),
  oma = c(0, 0, 0.5, 0)
)

for (ss in rep_samples) {
  d <- grid12[
    sample == ss
  ]

  plot_spatial_panel(
    d,
    "pyrimidine_score",
    paste0(
      ss,
      " — TYMS–UMPS"
    )
  )

  plot_spatial_panel(
    d,
    "proliferation_score",
    paste0(
      gsub("_", " ", d$cohort[1]),
      " — Proliferation"
    )
  )

  plot_spatial_panel(
    d,
    "epithelial_score",
    "Epithelial content"
  )
}

par(op)
dev.off()

# ----------------------------------------------------------------------
# 8. Overall reporting summary
# ----------------------------------------------------------------------
primary_dup <- dup_stats[
  scenario == "exclude_both_primary_n12" &
    association == "partial_spearman_epithelial_and_library_adjusted"
]

epi_adj <- epi_inference[
  association == "partial_spearman_epithelial_and_library_adjusted"
]

shift_adj <- shift_summary[
  association == "partial_rank_residual_surface"
]

n_shift_p05 <- sum(
  shift_adj$toroidal_two_sided_p < 0.05,
  na.rm = TRUE
)

n_shift_below_null97 <- sum(
  shift_adj$observed_rho >
    shift_adj$toroidal_null_q975,
  na.rm = TRUE
)

interpretation <- paste0(
  "After conservative exclusion of both duplicated processed payloads, the ",
  "previous adjusted specimen-level coupling remains the primary spatial ",
  "summary (n=12; median rho=",
  sprintf("%.3f", primary_dup$median_rho),
  "; exact two-sided sign-flip P=",
  formatC(primary_dup$exact_two_sided_p, format = "g", digits = 3),
  "). In the epithelial-enriched proxy subset, the adjusted median rho was ",
  sprintf("%.3f", epi_adj$median_rho),
  " (exact two-sided P=",
  formatC(epi_adj$exact_two_sided_p, format = "g", digits = 3),
  "). Toroidal-shift nulls were evaluated within each conservative specimen; ",
  n_shift_p05,
  "/12 had adjusted spatial-alignment P<0.05 and ",
  n_shift_below_null97,
  "/12 had observed adjusted rho above the specimen-specific 97.5th null ",
  "percentile. These data are same-study orthogonal spatial support for local ",
  "proliferation coupling, not an independent validation cohort and not proof ",
  "of malignant-cell-specific metabolic activity."
)

writeLines(
  c(
    "status=PASS_E5_SPATIAL_ROBUSTNESS_UPGRADE",
    "script_version=E5_V1",
    "dataset=GSE223501_same_study_SlideSeq2",
    "primary_duplicate_handling=exclude_PA056_and_KRAS_11",
    "duplicate_sensitivity_scenarios=exclude_both;retain_PA056_only;retain_KRAS11_only;retain_both_descriptive",
    "epithelial_enriched_proxy=within_specimen_top_quartile_epithelial_score",
    "malignant_specific_bins_claimed=FALSE",
    "spatial_null=toroidal_shift_of_proliferation_surface",
    paste0(
      "spatial_null_min_overlap_fraction=",
      MIN_OVERLAP_FRACTION
    ),
    paste0(
      "spatial_null_min_overlap_grids=",
      MIN_OVERLAP_GRIDS
    ),
    "grid_width=400_deposited_coordinate_units",
    "micron_conversion_asserted=FALSE",
    paste0(
      "adjusted_toroidal_p_lt_0.05_specimens=",
      n_shift_p05,
      "/12"
    ),
    paste0(
      "adjusted_observed_above_null_97.5pct_specimens=",
      n_shift_below_null97,
      "/12"
    ),
    paste0(
      "reviewer_facing_interpretation=",
      interpretation
    ),
    "independent_spatial_cohort_claimed=FALSE",
    "primary_spatial_endpoint_replaced=FALSE",
    "gene_sets_changed=FALSE"
  ),
  file.path(
    OUT,
    "STEP_E5_COMPLETE.txt"
  )
)

cat(
  "\n============================================================\n",
  "STEP E5: SPATIAL ROBUSTNESS UPGRADE\n",
  "============================================================\n",
  "Duplicate-payload sensitivity:\n",
  sep = ""
)

print(
  dup_stats
)

cat(
  "\nEpithelial-enriched proxy inference:\n"
)

print(
  epi_inference
)

cat(
  "\nToroidal-shift adjusted summary:\n"
)

print(
  shift_adj[
    ,
    .(
      sample,
      cohort,
      eligible_grids,
      admissible_nonzero_shifts,
      observed_rho,
      toroidal_null_q025,
      toroidal_null_q975,
      toroidal_two_sided_p
    )
  ]
)

cat(
  "\nRepresentative-map selection:\n"
)

print(
  rep_rule
)

cat(
  "\nInterpretation:\n",
  interpretation,
  "\n\nStatus: PASS_E5_SPATIAL_ROBUSTNESS_UPGRADE\n\n",
  "Key output files:\n",
  file.path(OUT, "STEP_E5_01_duplicate_payload_sensitivity.tsv"), "\n",
  file.path(OUT, "STEP_E5_03_epithelial_enriched_proxy_inference.tsv"), "\n",
  file.path(OUT, "STEP_E5_04_toroidal_shift_null_summary.tsv"), "\n",
  file.path(OUT, "STEP_E5_06_coordinate_scale_transparency.tsv"), "\n",
  file.path(OUT, "STEP_E5_07_representative_map_selection.tsv"), "\n",
  file.path(OUT, "STEP_E5_COMPLETE.txt"), "\n",
  "And visually inspect:\n",
  atlas_pdf, "\n",
  rep_pdf, "\n",
  sep = ""
)
