# Extract fixed spatial target genes
# Streams the fixed target genes from GSE223501, harmonizes the documented terminal barcode
# suffix formats, requires complete one-to-one barcode matching, and saves compact count/coordinate objects.

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

C2_FORMAT_FILE <- file.path(
  ROOT,
  "results/reviewer_defense",
  "STEP_C2_01_matrix_format_audit.tsv"
)

OUT_DIR <- file.path(
  ROOT,
  "results/reviewer_defense"
)

TARGET_OUT_DIR <- file.path(
  ROOT,
  "data/processed/GSE223501_spatial_targets"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TARGET_OUT_DIR,
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

split_csv_line <- function(x) {
  clean_csv_field(
    strsplit(
      x,
      ",",
      fixed = TRUE
    )[[1]]
  )
}

first_csv_field <- function(x) {
  clean_csv_field(
    sub(
      ",.*$",
      "",
      x
    )
  )
}

safe_stub <- function(x) {
  gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x
  )
}

canon_bead <- function(x) {
  z <- toupper(
    clean_csv_field(x)
  )

  # Harmonize terminal library suffix only:
  # ABC.1 -> ABC_1
  # ABC-1 -> ABC_1
  z <- sub(
    "[.-]([0-9]+)$",
    "_\\1",
    z,
    perl = TRUE
  )

  z[
    !nzchar(z)
  ] <- NA_character_

  z
}

# ------------------------------------------------------------
# Stream target rows only.
# ------------------------------------------------------------

extract_target_rows <- function(
  counts_file,
  target_genes,
  chunk_lines = 20L
) {
  con <- gzfile(
    counts_file,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  header_line <- readLines(
    con,
    n = 1L,
    warn = FALSE
  )

  if (length(header_line) != 1L) {
    stop(
      "无法读取counts header：",
      counts_file
    )
  }

  header_fields <- split_csv_line(
    header_line
  )

  if (length(header_fields) < 2L) {
    stop(
      "counts header字段数异常：",
      counts_file
    )
  }

  bead_ids_raw <- header_fields[-1L]
  bead_keys <- canon_bead(
    bead_ids_raw
  )

  if (
    anyNA(bead_keys) ||
      any(!nzchar(bead_keys))
  ) {
    stop(
      "counts中存在无法标准化的bead barcode：",
      counts_file
    )
  }

  if (anyDuplicated(bead_keys)) {
    stop(
      "counts barcode标准化后出现重复key：",
      counts_file
    )
  }

  target_upper <- toupper(
    target_genes
  )

  found <- setNames(
    vector(
      "list",
      length(target_genes)
    ),
    target_upper
  )

  n_rows_scanned <- 0L

  repeat {
    lines <- readLines(
      con,
      n = chunk_lines,
      warn = FALSE
    )

    if (length(lines) == 0L) {
      break
    }

    n_rows_scanned <- n_rows_scanned +
      length(lines)

    genes_here <- toupper(
      first_csv_field(
        lines
      )
    )

    hit <- which(
      genes_here %in%
        target_upper
    )

    if (length(hit) == 0L) {
      next
    }

    for (h in hit) {
      gene <- genes_here[h]

      if (!is.null(found[[gene]])) {
        stop(
          "目标基因在counts中出现重复行：",
          gene,
          "\n文件：",
          counts_file
        )
      }

      fields <- split_csv_line(
        lines[h]
      )

      if (
        length(fields) !=
          length(bead_ids_raw) + 1L
      ) {
        stop(
          "目标基因 ",
          gene,
          " 的字段数与header不一致。\n文件：",
          counts_file
        )
      }

      values <- suppressWarnings(
        as.numeric(
          fields[-1L]
        )
      )

      if (anyNA(values)) {
        stop(
          "目标基因 ",
          gene,
          " 存在无法解析为数值的count。\n文件：",
          counts_file
        )
      }

      found[[gene]] <- values
    }

    if (
      all(
        !vapply(
          found,
          is.null,
          logical(1)
        )
      )
    ) {
      break
    }
  }

  missing <- names(found)[
    vapply(
      found,
      is.null,
      logical(1)
    )
  ]

  if (length(missing) > 0L) {
    stop(
      "目标基因未全部提取：",
      paste(
        missing,
        collapse = ", "
      ),
      "\n文件：",
      counts_file
    )
  }

  mat <- do.call(
    rbind,
    found
  )

  mat <- mat[
    toupper(target_genes),
    ,
    drop = FALSE
  ]

  rownames(mat) <- target_genes

  list(
    bead_ids_raw = bead_ids_raw,
    bead_keys = bead_keys,
    counts = mat,
    n_rows_scanned = n_rows_scanned
  )
}

# ------------------------------------------------------------
# Match coords by canonicalized barcode.
# ------------------------------------------------------------

match_coords_to_counts <- function(
  coords_file,
  counts_bead_ids_raw,
  counts_bead_keys
) {
  coords <- fread(
    coords_file,
    data.table = TRUE,
    showProgress = FALSE
  )

  if (
    nrow(coords) !=
      length(
        counts_bead_keys
      )
  ) {
    stop(
      "coords行数与counts bead数不一致：",
      coords_file
    )
  }

  # Identify barcode-like column by normalized overlap.
  overlap_table <- rbindlist(
    lapply(
      names(coords),
      function(nm) {
        raw <- as.character(
          coords[[nm]]
        )

        key <- canon_bead(
          raw
        )

        n_match <- sum(
          !is.na(key) &
            key %chin%
              counts_bead_keys
        )

        data.table(
          column = nm,
          class = paste(
            class(
              coords[[nm]]
            ),
            collapse = "/"
          ),
          raw_exact_overlap =
            sum(
              raw %chin%
                counts_bead_ids_raw
            ),
          canonical_overlap =
            n_match
        )
      }
    )
  )

  best_idx <- which.max(
    overlap_table$canonical_overlap
  )

  barcode_col <- overlap_table$column[
    best_idx
  ]

  canonical_overlap <- overlap_table$canonical_overlap[
    best_idx
  ]

  if (
    canonical_overlap !=
      length(
        counts_bead_keys
      )
  ) {
    stop(
      "barcode标准化后仍未达到100%匹配。\n",
      "最佳coords列：",
      barcode_col,
      "\n",
      "匹配数：",
      canonical_overlap,
      "/",
      length(
        counts_bead_keys
      ),
      "\n文件：",
      coords_file
    )
  }

  coords_ids_raw <- as.character(
    coords[[
      barcode_col
    ]]
  )

  coords_keys <- canon_bead(
    coords_ids_raw
  )

  if (anyNA(coords_keys)) {
    stop(
      "coords barcode标准化后存在NA：",
      coords_file
    )
  }

  if (anyDuplicated(coords_keys)) {
    stop(
      "coords barcode标准化后出现重复key：",
      coords_file
    )
  }

  set_identical <- setequal(
    counts_bead_keys,
    coords_keys
  )

  if (!set_identical) {
    stop(
      "counts与coords标准化barcode集合不一致：",
      coords_file
    )
  }

  original_order_identical_after_canon <- identical(
    counts_bead_keys,
    coords_keys
  )

  idx <- match(
    counts_bead_keys,
    coords_keys
  )

  if (anyNA(idx)) {
    stop(
      "标准化barcode match产生NA。"
    )
  }

  coords_aligned <- coords[
    idx
  ]

  coords_keys_aligned <- canon_bead(
    coords_aligned[[
      barcode_col
    ]]
  )

  if (!identical(
    counts_bead_keys,
    coords_keys_aligned
  )) {
    stop(
      "重排coords后barcode顺序仍不一致。"
    )
  }

  candidate_cols <- setdiff(
    names(coords_aligned),
    barcode_col
  )

  numeric_candidates <- candidate_cols[
    vapply(
      coords_aligned[
        ,
        ..candidate_cols
      ],
      is.numeric,
      logical(1)
    )
  ]

  if (length(numeric_candidates) < 2L) {
    stop(
      "无法识别两个数值坐标列：",
      coords_file
    )
  }

  # Prefer exact x/y names if available.
  x_col <- if (
    "x" %in%
      names(coords_aligned)
  ) {
    "x"
  } else {
    numeric_candidates[1L]
  }

  y_col <- if (
    "y" %in%
      names(coords_aligned)
  ) {
    "y"
  } else {
    setdiff(
      numeric_candidates,
      x_col
    )[1L]
  }

  list(
    overlap_table = overlap_table,
    barcode_col = barcode_col,
    x_col = x_col,
    y_col = y_col,
    raw_exact_overlap =
      overlap_table[
        column ==
          barcode_col,
        raw_exact_overlap
      ],
    canonical_overlap =
      canonical_overlap,
    original_order_identical_after_canon =
      original_order_identical_after_canon,
    aligned = data.table(
      bead_id_counts_raw =
        counts_bead_ids_raw,
      bead_id_coords_raw =
        as.character(
          coords_aligned[[
            barcode_col
          ]]
        ),
      bead_key =
        counts_bead_keys,
      x = as.numeric(
        coords_aligned[[
          x_col
        ]]
      ),
      y = as.numeric(
        coords_aligned[[
          y_col
        ]]
      )
    )
  )
}

banner("STEP C3 v2：规范化barcode后提取目标基因")

# ------------------------------------------------------------
# 1. Input checks
# ------------------------------------------------------------

if (!file.exists(INVENTORY_FILE)) {
  stop(
    "缺少STEP C1 inventory：\n",
    INVENTORY_FILE
  )
}

if (!file.exists(C2_FORMAT_FILE)) {
  stop(
    "缺少STEP C2 format audit：\n",
    C2_FORMAT_FILE
  )
}

inv <- fread(
  INVENTORY_FILE,
  data.table = TRUE,
  showProgress = FALSE
)

fmt <- fread(
  C2_FORMAT_FILE,
  data.table = TRUE,
  showProgress = FALSE
)

if (
  nrow(inv) != 14L ||
    nrow(fmt) != 14L
) {
  stop(
    "STEP C1/C2样本数不是14。"
  )
}

if (!all(
  fmt$orientation ==
    "genes_x_beads"
)) {
  stop(
    "STEP C2并非全部 genes_x_beads。"
  )
}

if (!all(
  fmt$bead_count_match
)) {
  stop(
    "STEP C2存在counts/coords bead数不一致。"
  )
}

# ------------------------------------------------------------
# 2. Process one sample at a time
# ------------------------------------------------------------

sample_audit_list <- vector(
  "list",
  nrow(inv)
)

gene_qc_list <- vector(
  "list",
  nrow(inv)
)

barcode_overlap_list <- vector(
  "list",
  nrow(inv)
)

for (i in seq_len(nrow(inv))) {

  sample_name <- inv$official_paper_id[i]
  cohort <- inv$cohort[i]
  counts_file <- inv$counts_file[i]
  coords_file <- inv$coords_file[i]

  cat(
    "\n[",
    i,
    "/",
    nrow(inv),
    "] ",
    sample_name,
    " (",
    cohort,
    ")\n",
    sep = ""
  )

  ext <- extract_target_rows(
    counts_file = counts_file,
    target_genes = TARGET_GENES,
    chunk_lines = 20L
  )

  coord_match <- match_coords_to_counts(
    coords_file =
      coords_file,
    counts_bead_ids_raw =
      ext$bead_ids_raw,
    counts_bead_keys =
      ext$bead_keys
  )

  barcode_overlap_list[[i]] <- copy(
    coord_match$overlap_table
  )[
    ,
    `:=`(
      sample = sample_name,
      cohort = cohort
    )
  ]

  target_mat <- ext$counts

  if (
    ncol(target_mat) !=
      nrow(
        coord_match$aligned
      )
  ) {
    stop(
      "目标矩阵与coords行数不一致：",
      sample_name
    )
  }

  target_dt <- as.data.table(
    t(
      target_mat
    )
  )

  setnames(
    target_dt,
    TARGET_GENES
  )

  compact <- cbind(
    data.table(
      gsm = inv$gsm[i],
      sample = sample_name,
      cohort = cohort
    ),
    coord_match$aligned,
    target_dt
  )

  rds_file <- file.path(
    TARGET_OUT_DIR,
    paste0(
      safe_stub(
        sample_name
      ),
      "_target_counts_coords.rds"
    )
  )

  saveRDS(
    compact,
    rds_file,
    compress = "xz"
  )

  gene_qc <- rbindlist(
    lapply(
      TARGET_GENES,
      function(g) {
        x <- compact[[g]]

        data.table(
          gsm = inv$gsm[i],
          sample = sample_name,
          cohort = cohort,
          gene = g,
          group = names(
            TARGET_GROUPS
          )[
            vapply(
              TARGET_GROUPS,
              function(z) {
                g %in% z
              },
              logical(1)
            )
          ][1],
          n_beads = length(x),
          total_count = sum(x),
          mean_count = mean(x),
          median_count = median(x),
          max_count = max(x),
          nonzero_beads = sum(
            x > 0
          ),
          nonzero_fraction = mean(
            x > 0
          ),
          integer_like = all(
            abs(
              x -
                round(x)
            ) <
              1e-8
          )
        )
      }
    )
  )

  gene_qc_list[[i]] <- gene_qc

  sample_audit_list[[i]] <- data.table(
    gsm = inv$gsm[i],
    sample = sample_name,
    cohort = cohort,
    n_beads = ncol(
      target_mat
    ),
    target_genes_extracted =
      nrow(
        target_mat
      ),
    counts_rows_scanned =
      ext$n_rows_scanned,
    coords_barcode_column =
      coord_match$barcode_col,
    coords_x_column =
      coord_match$x_col,
    coords_y_column =
      coord_match$y_col,
    raw_exact_barcode_overlap =
      coord_match$raw_exact_overlap,
    canonical_barcode_overlap =
      coord_match$canonical_overlap,
    canonical_overlap_fraction =
      coord_match$canonical_overlap /
        ncol(
          target_mat
        ),
    original_order_identical_after_canon =
      coord_match$
        original_order_identical_after_canon,
    canonical_barcode_set_match = TRUE,
    aligned_order_match = TRUE,
    all_target_values_nonnegative =
      all(
        target_mat >= 0
      ),
    all_target_values_integer_like =
      all(
        abs(
          target_mat -
            round(
              target_mat
            )
        ) <
          1e-8
      ),
    output_rds = rds_file,
    output_rds_MB = round(
      file.info(
        rds_file
      )$size /
        1024^2,
      3
    )
  )

  cat(
    "  raw exact barcode overlap: ",
    coord_match$raw_exact_overlap,
    "/",
    ncol(
      target_mat
    ),
    "\n",
    "  canonical barcode overlap: ",
    coord_match$canonical_overlap,
    "/",
    ncol(
      target_mat
    ),
    "\n",
    "  canonical set match: TRUE\n",
    "  aligned order match: TRUE\n",
    "  14 target genes extracted: TRUE\n",
    "  saved: ",
    basename(
      rds_file
    ),
    "\n",
    sep = ""
  )

  rm(
    ext,
    coord_match,
    target_mat,
    target_dt,
    compact,
    gene_qc
  )

  gc(
    verbose = FALSE
  )
}

sample_audit <- rbindlist(
  sample_audit_list,
  fill = TRUE
)

gene_qc <- rbindlist(
  gene_qc_list,
  fill = TRUE
)

barcode_overlap_audit <- rbindlist(
  barcode_overlap_list,
  fill = TRUE
)

setorder(
  sample_audit,
  cohort,
  sample
)

setorder(
  gene_qc,
  group,
  gene,
  cohort,
  sample
)

# ------------------------------------------------------------
# 3. Summary QC
# ------------------------------------------------------------

gene_summary <- gene_qc[
  ,
  .(
    samples = .N,
    total_beads =
      sum(
        n_beads
      ),
    total_count_all_samples =
      sum(
        total_count
      ),
    median_nonzero_fraction =
      median(
        nonzero_fraction
      ),
    min_nonzero_fraction =
      min(
        nonzero_fraction
      ),
    max_nonzero_fraction =
      max(
        nonzero_fraction
      ),
    all_samples_integer_like =
      all(
        integer_like
      )
  ),
  by = .(
    group,
    gene
  )
]

qc_pass <- (
  nrow(sample_audit) == 14L &&
    all(
      sample_audit$target_genes_extracted ==
        length(TARGET_GENES)
    ) &&
    all(
      sample_audit$canonical_barcode_set_match
    ) &&
    all(
      sample_audit$aligned_order_match
    ) &&
    all(
      sample_audit$canonical_overlap_fraction ==
        1
    ) &&
    all(
      sample_audit$all_target_values_nonnegative
    )
)

status <- if (qc_pass) {
  "PASS_TARGET_EXTRACTION_AND_CANONICAL_BARCODE_ALIGNMENT"
} else {
  "REVIEW_REQUIRED"
}

# ------------------------------------------------------------
# 4. Outputs
# ------------------------------------------------------------

fwrite(
  sample_audit,
  file.path(
    OUT_DIR,
    "STEP_C3v2_01_target_extraction_sample_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_qc,
  file.path(
    OUT_DIR,
    "STEP_C3v2_02_target_gene_raw_count_QC_by_sample.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_summary,
  file.path(
    OUT_DIR,
    "STEP_C3v2_03_target_gene_raw_count_QC_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  barcode_overlap_audit,
  file.path(
    OUT_DIR,
    "STEP_C3v2_04_barcode_overlap_audit.tsv"
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
    "barcode_rule=terminal .<integer> and -<integer> normalized to _<integer>",
    "row_order_used_as_primary_matching_method=FALSE",
    "barcode_matching_method=canonicalized exact one-to-one key matching",
    paste0(
      "all_14_canonical_barcode_overlap_fraction_1=",
      all(
        sample_audit$canonical_overlap_fraction ==
          1
      )
    ),
    paste0(
      "all_14_canonical_sets_match=",
      all(
        sample_audit$canonical_barcode_set_match
      )
    ),
    paste0(
      "all_14_aligned_orders_match=",
      all(
        sample_audit$aligned_order_match
      )
    ),
    "full_expression_matrix_loaded=FALSE",
    "statistical_testing_performed=FALSE"
  ),
  file.path(
    OUT_DIR,
    "STEP_C3v2_COMPLETE.txt"
  )
)

banner("STEP C3 v2结果")

print(
  sample_audit[
    ,
    .(
      sample,
      cohort,
      n_beads,
      coords_barcode_column,
      raw_exact_barcode_overlap,
      canonical_barcode_overlap,
      canonical_overlap_fraction,
      original_order_identical_after_canon,
      all_target_values_integer_like
    )
  ]
)

cat(
  "\n目标基因raw-count QC汇总：\n"
)

print(
  gene_summary
)

cat(
  "\nStatus: ",
  status,
  "\n\n",
  "Key output files:\n",
  file.path(
    OUT_DIR,
    "STEP_C3v2_01_target_extraction_sample_audit.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C3v2_02_target_gene_raw_count_QC_by_sample.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C3v2_03_target_gene_raw_count_QC_summary.tsv"
  ),
  "\n",
  sep = ""
)
