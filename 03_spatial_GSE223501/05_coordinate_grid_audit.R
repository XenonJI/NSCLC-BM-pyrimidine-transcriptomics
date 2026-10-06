# Coordinate-only grid audit
# Evaluates candidate spatial grid widths using coordinates and bead counts only. Expression
# correlations are not used in this step; the primary grid scale is selected from spatial coverage.

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
  200,
  250,
  300,
  400,
  500,
  600
)

BIN_THRESHOLDS <- c(
  50L,
  100L,
  250L,
  500L
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

banner("STEP C5A：空间网格聚合设计审计（坐标-only）")

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

sample_width_list <- list()
occupancy_list <- list()

k <- 0L
m <- 0L

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

  if (!file.exists(target_file)) {
    stop(
      "缺少STEP C3v2 target RDS：\n",
      target_file
    )
  }

  d <- readRDS(
    target_file
  )

  required <- c(
    "bead_key",
    "x",
    "y"
  )

  missing <- setdiff(
    required,
    names(d)
  )

  if (length(missing) > 0L) {
    stop(
      sample_name,
      " target RDS缺少：",
      paste(
        missing,
        collapse = ", "
      )
    )
  }

  if (
    anyNA(d$x) ||
      anyNA(d$y)
  ) {
    stop(
      sample_name,
      " 的x/y存在NA。"
    )
  }

  x_min <- min(d$x)
  x_max <- max(d$x)
  y_min <- min(d$y)
  y_max <- max(d$y)

  cat(
    "\n[",
    i,
    "/14] ",
    sample_name,
    " | beads=",
    nrow(d),
    " | x-range=",
    round(x_min, 1),
    "–",
    round(x_max, 1),
    " | y-range=",
    round(y_min, 1),
    "–",
    round(y_max, 1),
    "\n",
    sep = ""
  )

  for (w in GRID_WIDTHS) {

    # Per-sample origin anchored at that sample's observed minimum.
    # The width itself is fixed across all samples.
    gx <- floor(
      (d$x - x_min) /
        w
    )

    gy <- floor(
      (d$y - y_min) /
        w
    )

    occ <- data.table(
      gx = gx,
      gy = gy
    )[
      ,
      .(
        n_beads = .N
      ),
      by = .(
        gx,
        gy
      )
    ]

    n_bins <- nrow(
      occ
    )

    total_beads <- sum(
      occ$n_beads
    )

    if (
      total_beads !=
        nrow(d)
    ) {
      stop(
        "网格计数未覆盖全部beads：",
        sample_name,
        " width=",
        w
      )
    }

    k <- k + 1L

    row <- data.table(
      sample = sample_name,
      cohort = cohort,
      grid_width = w,
      total_beads = total_beads,
      nonempty_bins = n_bins,
      min_beads_per_bin =
        min(
          occ$n_beads
        ),
      q10_beads_per_bin =
        as.numeric(
          quantile(
            occ$n_beads,
            0.10,
            names = FALSE
          )
        ),
      q25_beads_per_bin =
        as.numeric(
          quantile(
            occ$n_beads,
            0.25,
            names = FALSE
          )
        ),
      median_beads_per_bin =
        median(
          occ$n_beads
        ),
      mean_beads_per_bin =
        mean(
          occ$n_beads
        ),
      q75_beads_per_bin =
        as.numeric(
          quantile(
            occ$n_beads,
            0.75,
            names = FALSE
          )
        ),
      q90_beads_per_bin =
        as.numeric(
          quantile(
            occ$n_beads,
            0.90,
            names = FALSE
          )
        ),
      max_beads_per_bin =
        max(
          occ$n_beads
        )
    )

    for (thr in BIN_THRESHOLDS) {

      keep <- occ$n_beads >=
        thr

      row[
        ,
        (paste0(
          "bins_ge_",
          thr
        )) := sum(
          keep
        )
      ]

      row[
        ,
        (paste0(
          "fraction_bins_ge_",
          thr
        )) := mean(
          keep
        )
      ]

      row[
        ,
        (paste0(
          "fraction_beads_in_bins_ge_",
          thr
        )) := sum(
          occ$n_beads[
            keep
          ]
        ) /
          total_beads
      ]
    }

    sample_width_list[[k]] <- row

    # Save full occupancy distribution for later deterministic review.
    m <- m + 1L

    occupancy_list[[m]] <- occ[
      ,
      `:=`(
        sample = sample_name,
        cohort = cohort,
        grid_width = w
      )
    ][
      ,
      .(
        sample,
        cohort,
        grid_width,
        gx,
        gy,
        n_beads
      )
    ]

    cat(
      "  width=",
      w,
      " | bins=",
      n_bins,
      " | median beads/bin=",
      round(
        median(
          occ$n_beads
        ),
        1
      ),
      " | bins>=250=",
      sum(
        occ$n_beads >=
          250L
      ),
      " | beads retained @>=250=",
      round(
        100 *
          sum(
            occ$n_beads[
              occ$n_beads >=
                250L
            ]
          ) /
          total_beads,
        1
      ),
      "%\n",
      sep = ""
    )
  }

  rm(
    d,
    gx,
    gy,
    occ
  )

  gc(
    verbose = FALSE
  )
}

sample_width <- rbindlist(
  sample_width_list,
  fill = TRUE
)

occupancy <- rbindlist(
  occupancy_list,
  fill = TRUE
)

setorder(
  sample_width,
  grid_width,
  cohort,
  sample
)

# ------------------------------------------------------------
# Across-sample summaries by width
# ------------------------------------------------------------

width_summary <- sample_width[
  ,
  .(
    samples = .N,
    median_nonempty_bins =
      median(
        nonempty_bins
      ),
    min_nonempty_bins =
      min(
        nonempty_bins
      ),
    max_nonempty_bins =
      max(
        nonempty_bins
      ),
    median_of_sample_median_beads_per_bin =
      median(
        median_beads_per_bin
      ),
    min_sample_median_beads_per_bin =
      min(
        median_beads_per_bin
      ),
    max_sample_median_beads_per_bin =
      max(
        median_beads_per_bin
      ),
    median_bins_ge_50 =
      median(
        bins_ge_50
      ),
    min_bins_ge_50 =
      min(
        bins_ge_50
      ),
    median_fraction_beads_in_bins_ge_50 =
      median(
        fraction_beads_in_bins_ge_50
      ),
    median_bins_ge_100 =
      median(
        bins_ge_100
      ),
    min_bins_ge_100 =
      min(
        bins_ge_100
      ),
    median_fraction_beads_in_bins_ge_100 =
      median(
        fraction_beads_in_bins_ge_100
      ),
    median_bins_ge_250 =
      median(
        bins_ge_250
      ),
    min_bins_ge_250 =
      min(
        bins_ge_250
      ),
    median_fraction_beads_in_bins_ge_250 =
      median(
        fraction_beads_in_bins_ge_250
      ),
    median_bins_ge_500 =
      median(
        bins_ge_500
      ),
    min_bins_ge_500 =
      min(
        bins_ge_500
      ),
    median_fraction_beads_in_bins_ge_500 =
      median(
        fraction_beads_in_bins_ge_500
      )
  ),
  by = grid_width
]

setorder(
  width_summary,
  grid_width
)

# ------------------------------------------------------------
# Cohort-stratified descriptive summary
# ------------------------------------------------------------

cohort_width_summary <- sample_width[
  ,
  .(
    samples = .N,
    median_nonempty_bins =
      median(
        nonempty_bins
      ),
    median_of_sample_median_beads_per_bin =
      median(
        median_beads_per_bin
      ),
    median_bins_ge_100 =
      median(
        bins_ge_100
      ),
    median_bins_ge_250 =
      median(
        bins_ge_250
      ),
    median_bins_ge_500 =
      median(
        bins_ge_500
      ),
    median_fraction_beads_in_bins_ge_250 =
      median(
        fraction_beads_in_bins_ge_250
      )
  ),
  by = .(
    grid_width,
    cohort
  )
]

setorder(
  cohort_width_summary,
  grid_width,
  cohort
)

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------

fwrite(
  sample_width,
  file.path(
    OUT_DIR,
    "STEP_C5A_01_grid_width_QC_by_sample.tsv"
  ),
  sep = "\t"
)

fwrite(
  width_summary,
  file.path(
    OUT_DIR,
    "STEP_C5A_02_grid_width_QC_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  cohort_width_summary,
  file.path(
    OUT_DIR,
    "STEP_C5A_03_grid_width_QC_by_cohort.tsv"
  ),
  sep = "\t"
)

fwrite(
  occupancy,
  file.path(
    OUT_DIR,
    "STEP_C5A_04_full_bin_occupancy.tsv"
  ),
  sep = "\t"
)

status <- if (
  nrow(sample_width) ==
    14L *
      length(
        GRID_WIDTHS
      ) &&
    all(
      sample_width$total_beads >
        0
    )
) {
  "PASS_COORDINATE_ONLY_GRID_AUDIT"
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
    "samples=14",
    paste0(
      "candidate_grid_widths=",
      paste(
        GRID_WIDTHS,
        collapse = ","
      )
    ),
    paste0(
      "bin_thresholds=",
      paste(
        BIN_THRESHOLDS,
        collapse = ","
      )
    ),
    "expression_values_used=FALSE",
    "module_scoring_performed=FALSE",
    "correlation_testing_performed=FALSE",
    "p_values_computed=FALSE",
    "purpose=lock spatial aggregation scale before viewing expression correlations"
  ),
  file.path(
    OUT_DIR,
    "STEP_C5A_COMPLETE.txt"
  )
)

banner("STEP C5A结果")

print(
  width_summary
)

cat(
  "\nStatus: ",
  status,
  "\n\n",
  "Key output files:\n",
  file.path(
    OUT_DIR,
    "STEP_C5A_01_grid_width_QC_by_sample.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C5A_02_grid_width_QC_summary.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C5A_03_grid_width_QC_by_cohort.tsv"
  ),
  "\n",
  sep = ""
)
