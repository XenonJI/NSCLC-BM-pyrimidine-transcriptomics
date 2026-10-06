# GSE223501 matrix and target-gene audit
# Confirms count-matrix orientation, bead/coordinate dimensions, and availability of the
# prespecified pyrimidine, proliferation, and epithelial marker genes.

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

OUT_DIR <- file.path(
  ROOT,
  "results/reviewer_defense"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

TARGET_GROUPS <- list(
  pyrimidine = c(
    "TYMS",
    "UMPS"
  ),
  proliferation_nonoverlap = c(
    "MKI67",
    "PCNA",
    "TOP2A",
    "UBE2C",
    "CENPF",
    "BIRC5",
    "CCNB1"
  ),
  epithelial_content = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "MUC1"
  )
)

TARGET_GENES <- unique(
  unlist(
    TARGET_GROUPS,
    use.names = FALSE
  )
)

banner <- function(x) {
  cat(
    "\n========================================\n",
    x,
    "\n========================================\n",
    sep = ""
  )
}

clean_csv_field <- function(x) {
  x <- trimws(
    as.character(x)
  )

  x <- sub(
    '^"',
    "",
    x
  )

  x <- sub(
    '"$',
    "",
    x
  )

  trimws(x)
}

split_csv_header <- function(x) {
  # GEO files use simple CSV headers; this parser is intentionally
  # lightweight because only field names are needed here.
  fields <- strsplit(
    x,
    ",",
    fixed = TRUE
  )[[1]]

  clean_csv_field(
    fields
  )
}

count_gz_data_lines <- function(
  path,
  skip_header = TRUE,
  chunk_size = 2000L
) {
  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  if (skip_header) {
    readLines(
      con,
      n = 1L,
      warn = FALSE
    )
  }

  n_total <- 0L

  repeat {
    x <- readLines(
      con,
      n = chunk_size,
      warn = FALSE
    )

    if (length(x) == 0L) {
      break
    }

    n_total <- n_total +
      length(x)
  }

  n_total
}

scan_first_column_for_targets <- function(
  path,
  targets,
  chunk_size = 20L
) {
  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  header <- readLines(
    con,
    n = 1L,
    warn = FALSE
  )

  target_upper <- toupper(
    targets
  )

  counts <- setNames(
    integer(
      length(targets)
    ),
    target_upper
  )

  first_labels <- character()
  n_rows <- 0L

  repeat {
    lines <- readLines(
      con,
      n = chunk_size,
      warn = FALSE
    )

    if (length(lines) == 0L) {
      break
    }

    n_rows <- n_rows +
      length(lines)

    first_field <- sub(
      ",.*$",
      "",
      lines
    )

    first_field <- clean_csv_field(
      first_field
    )

    if (length(first_labels) < 10L) {
      needed <- 10L -
        length(first_labels)

      first_labels <- c(
        first_labels,
        head(
          first_field,
          needed
        )
      )
    }

    first_upper <- toupper(
      first_field
    )

    tab <- table(
      first_upper[
        first_upper %in%
          target_upper
      ]
    )

    if (length(tab) > 0L) {
      idx <- match(
        names(tab),
        names(counts)
      )

      counts[idx] <- counts[idx] +
        as.integer(tab)
    }
  }

  list(
    target_counts = counts,
    n_data_rows = n_rows,
    first_labels = first_labels
  )
}

banner("STEP C2：矩阵方向 + 目标基因覆盖审计")

# ------------------------------------------------------------
# 1. 输入检查
# ------------------------------------------------------------

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

required_columns <- c(
  "gsm",
  "sample_raw",
  "official_paper_id",
  "cohort",
  "counts_file",
  "coords_file",
  "pair_complete",
  "official_metadata_match"
)

missing_columns <- setdiff(
  required_columns,
  names(inv)
)

if (length(missing_columns) > 0L) {
  stop(
    "STEP C1 inventory缺少字段：",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

if (
  nrow(inv) != 14L ||
    !all(inv$pair_complete) ||
    !all(inv$official_metadata_match)
) {
  stop(
    "STEP C1没有达到14个完整配对并全部匹配官方metadata。"
  )
}

# ------------------------------------------------------------
# 2. 每个样本低内存审计
# ------------------------------------------------------------

format_results <- list()
coverage_results <- list()

for (i in seq_len(nrow(inv))) {

  sample_name <- inv$official_paper_id[i]
  counts_path <- inv$counts_file[i]
  coords_path <- inv$coords_file[i]

  cat(
    "\n[",
    i,
    "/",
    nrow(inv),
    "] ",
    sample_name,
    "\n",
    sep = ""
  )

  if (!file.exists(counts_path)) {
    stop(
      "counts文件不存在：\n",
      counts_path
    )
  }

  if (!file.exists(coords_path)) {
    stop(
      "coords文件不存在：\n",
      coords_path
    )
  }

  # ----- counts header -----
  con_counts <- gzfile(
    counts_path,
    open = "rt"
  )

  counts_header_line <- readLines(
    con_counts,
    n = 1L,
    warn = FALSE
  )

  close(
    con_counts
  )

  if (length(counts_header_line) != 1L) {
    stop(
      "无法读取counts header：",
      counts_path
    )
  }

  counts_header <- split_csv_header(
    counts_header_line
  )

  header_upper <- toupper(
    counts_header
  )

  header_target_counts <- vapply(
    TARGET_GENES,
    function(g) {
      sum(
        header_upper ==
          toupper(g)
      )
    },
    integer(1)
  )

  n_targets_in_header <- sum(
    header_target_counts > 0L
  )

  # ----- coords header and row count -----
  con_coords <- gzfile(
    coords_path,
    open = "rt"
  )

  coords_header_line <- readLines(
    con_coords,
    n = 1L,
    warn = FALSE
  )

  close(
    con_coords
  )

  coords_header <- split_csv_header(
    coords_header_line
  )

  n_coords_rows <- count_gz_data_lines(
    coords_path,
    skip_header = TRUE,
    chunk_size = 5000L
  )

  # ----- orientation -----
  if (n_targets_in_header >= 2L) {

    orientation <- "beads_x_genes"

    n_counts_data_rows <- count_gz_data_lines(
      counts_path,
      skip_header = TRUE,
      # A row may be long when genes are columns.
      chunk_size = 20L
    )

    estimated_beads <- n_counts_data_rows

    target_counts <- setNames(
      header_target_counts,
      TARGET_GENES
    )

    first_labels <- NA_character_

  } else {

    # Likely genes are rows. Scan only the first field of each CSV row.
    scan_result <- scan_first_column_for_targets(
      counts_path,
      targets = TARGET_GENES,
      chunk_size = 20L
    )

    n_counts_data_rows <- scan_result$n_data_rows

    target_counts <- scan_result$target_counts

    names(target_counts) <- TARGET_GENES

    first_labels <- paste(
      scan_result$first_labels,
      collapse = " | "
    )

    n_targets_first_col <- sum(
      target_counts > 0L
    )

    if (n_targets_first_col >= 2L) {
      orientation <- "genes_x_beads"

      # first column is gene identifier; remaining header fields are beads.
      estimated_beads <- max(
        length(counts_header) - 1L,
        0L
      )
    } else {
      orientation <- "UNRESOLVED"
      estimated_beads <- NA_integer_
    }
  }

  bead_count_match <- (
    !is.na(estimated_beads) &&
      estimated_beads ==
        n_coords_rows
  )

  format_results[[i]] <- data.table(
    gsm = inv$gsm[i],
    sample = sample_name,
    cohort = inv$cohort[i],
    counts_header_fields =
      length(counts_header),
    counts_data_rows =
      n_counts_data_rows,
    coords_header_fields =
      length(coords_header),
    coords_data_rows =
      n_coords_rows,
    orientation =
      orientation,
    estimated_beads_from_counts =
      estimated_beads,
    coords_beads =
      n_coords_rows,
    bead_count_match =
      bead_count_match,
    targets_found_in_counts_header =
      n_targets_in_header,
    first_10_first_column_labels =
      first_labels
  )

  coverage_results[[i]] <- data.table(
    gsm = inv$gsm[i],
    sample = sample_name,
    cohort = inv$cohort[i],
    gene = TARGET_GENES,
    group = vapply(
      TARGET_GENES,
      function(g) {
        names(
          TARGET_GROUPS
        )[
          vapply(
            TARGET_GROUPS,
            function(z) {
              g %in% z
            },
            logical(1)
          )
        ][1]
      },
      character(1)
    ),
    occurrence_count =
      as.integer(
        target_counts[
          TARGET_GENES
        ]
      ),
    present =
      as.integer(
        target_counts[
          TARGET_GENES
        ]
      ) > 0L,
    duplicated_identifier =
      as.integer(
        target_counts[
          TARGET_GENES
        ]
      ) > 1L
  )

  cat(
    "  orientation: ",
    orientation,
    "\n",
    "  beads(counts): ",
    estimated_beads,
    "\n",
    "  beads(coords): ",
    n_coords_rows,
    "\n",
    "  bead count match: ",
    bead_count_match,
    "\n",
    "  target genes present: ",
    sum(
      target_counts > 0L
    ),
    "/",
    length(TARGET_GENES),
    "\n",
    sep = ""
  )

  rm(
    counts_header,
    header_upper,
    header_target_counts,
    target_counts
  )

  gc(
    verbose = FALSE
  )
}

format_audit <- rbindlist(
  format_results,
  fill = TRUE
)

coverage <- rbindlist(
  coverage_results,
  fill = TRUE
)

# ------------------------------------------------------------
# 3. 跨样本覆盖汇总
# ------------------------------------------------------------

coverage_summary <- coverage[
  ,
  .(
    samples_present =
      sum(
        present
      ),
    samples_total =
      .N,
    all_14_present =
      all(
        present
      ),
    any_duplicate_identifier =
      any(
        duplicated_identifier
      )
  ),
  by = .(
    group,
    gene
  )
]

setorder(
  coverage_summary,
  group,
  gene
)

group_summary <- coverage[
  ,
  .(
    genes_in_group =
      uniqueN(
        gene
      ),
    sample_gene_pairs =
      .N,
    present_pairs =
      sum(
        present
      ),
    complete_all_samples =
      all(
        present
      )
  ),
  by = group
]

# ------------------------------------------------------------
# 4. QC锁
# ------------------------------------------------------------

all_orientation_resolved <- all(
  format_audit$orientation !=
    "UNRESOLVED"
)

one_orientation_only <- (
  uniqueN(
    format_audit$orientation
  ) == 1L
)

all_beads_match <- all(
  format_audit$bead_count_match
)

all_pyrimidine_present <- coverage[
  group ==
    "pyrimidine",
  all(
    present
  )
]

all_proliferation_present <- coverage[
  group ==
    "proliferation_nonoverlap",
  all(
    present
  )
]

all_epithelial_present <- coverage[
  group ==
    "epithelial_content",
  all(
    present
  )
]

no_duplicate_target_ids <- !any(
  coverage$duplicated_identifier
)

qc_pass <- all(
  all_orientation_resolved,
  one_orientation_only,
  all_beads_match,
  all_pyrimidine_present,
  all_proliferation_present,
  all_epithelial_present,
  no_duplicate_target_ids
)

status <- if (qc_pass) {
  "PASS_MATRIX_FORMAT_AND_TARGET_COVERAGE"
} else {
  "REVIEW_REQUIRED"
}

# ------------------------------------------------------------
# 5. 输出
# ------------------------------------------------------------

fwrite(
  format_audit,
  file.path(
    OUT_DIR,
    "STEP_C2_01_matrix_format_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  coverage,
  file.path(
    OUT_DIR,
    "STEP_C2_02_target_gene_coverage_by_sample.tsv"
  ),
  sep = "\t"
)

fwrite(
  coverage_summary,
  file.path(
    OUT_DIR,
    "STEP_C2_03_target_gene_coverage_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  group_summary,
  file.path(
    OUT_DIR,
    "STEP_C2_04_target_group_coverage_summary.tsv"
  ),
  sep = "\t"
)

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
      "orientation=",
      paste(
        unique(
          format_audit$orientation
        ),
        collapse = ";"
      )
    ),
    paste0(
      "all_14_bead_counts_match_coords=",
      all_beads_match
    ),
    paste0(
      "TYMS_UMPS_present_all_samples=",
      all_pyrimidine_present
    ),
    paste0(
      "proliferation_genes_present_all_samples=",
      all_proliferation_present
    ),
    paste0(
      "epithelial_genes_present_all_samples=",
      all_epithelial_present
    ),
    paste0(
      "duplicate_target_gene_identifiers=",
      !no_duplicate_target_ids
    ),
    "full_expression_matrix_loaded=FALSE",
    "statistical_testing_performed=FALSE"
  ),
  file.path(
    OUT_DIR,
    "STEP_C2_COMPLETE.txt"
  )
)

banner("STEP C2结果")

print(
  format_audit[
    ,
    .(
      sample,
      cohort,
      orientation,
      estimated_beads_from_counts,
      coords_beads,
      bead_count_match
    )
  ]
)

cat(
  "\n目标基因跨14样本覆盖：\n"
)

print(
  coverage_summary
)

cat(
  "\nStatus: ",
  status,
  "\n",
  sep = ""
)

cat(
  "\nKey output files:\n",
  file.path(
    OUT_DIR,
    "STEP_C2_01_matrix_format_audit.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C2_03_target_gene_coverage_summary.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C2_04_target_group_coverage_summary.tsv"
  ),
  "\n",
  sep = ""
)
