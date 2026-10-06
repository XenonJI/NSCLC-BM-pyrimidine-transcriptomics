# Compute Slide-seq library depth
# Streams all-gene counts for each spatial specimen to calculate bead-level library size and
# detected-gene counts without loading the complete expression matrix into memory.

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
  TARGET_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  DEPTH_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

CHUNK_LINES <- 50L

TARGET_GENES <- c(
  "TYMS",
  "UMPS",
  "MKI67",
  "PCNA",
  "TOP2A",
  "UBE2C",
  "CENPF",
  "BIRC5",
  "CCNB1",
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "MUC1"
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
  clean_csv_field(
    strsplit(
      x,
      ",",
      fixed = TRUE
    )[[1]]
  )
}

canon_bead <- function(x) {
  z <- toupper(
    clean_csv_field(x)
  )

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

validate_existing_depth <- function(
  depth_file,
  target
) {
  if (!file.exists(depth_file)) {
    return(FALSE)
  }

  x <- tryCatch(
    readRDS(
      depth_file
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(x)) {
    return(FALSE)
  }

  required <- c(
    "bead_key",
    "library_size",
    "detected_genes"
  )

  if (!all(
    required %in%
      names(x)
  )) {
    return(FALSE)
  }

  if (
    nrow(x) !=
      nrow(target)
  ) {
    return(FALSE)
  }

  if (!identical(
    as.character(
      x$bead_key
    ),
    as.character(
      target$bead_key
    )
  )) {
    return(FALSE)
  }

  if (
    anyNA(
      x$library_size
    ) ||
      any(
        x$library_size < 0
      )
  ) {
    return(FALSE)
  }

  if (
    anyNA(
      x$detected_genes
    ) ||
      any(
        x$detected_genes < 0
      )
  ) {
    return(FALSE)
  }

  target_sum <- rowSums(
    as.matrix(
      target[
        ,
        ..TARGET_GENES
      ]
    )
  )

  if (any(
    target_sum >
      x$library_size
  )) {
    return(FALSE)
  }

  TRUE
}

compute_depth_streaming <- function(
  counts_file,
  target,
  chunk_lines = 50L,
  progress_prefix = ""
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

  header_fields <- split_csv_header(
    header_line
  )

  bead_key <- canon_bead(
    header_fields[-1L]
  )

  n_beads <- length(
    bead_key
  )

  if (
    nrow(target) !=
      n_beads
  ) {
    stop(
      "target RDS行数与counts bead数不一致：",
      counts_file
    )
  }

  if (!identical(
    bead_key,
    as.character(
      target$bead_key
    )
  )) {
    stop(
      "target RDS bead_key与counts header顺序不一致：",
      counts_file
    )
  }

  library_size <- numeric(
    n_beads
  )

  detected_genes <- integer(
    n_beads
  )

  genes_processed <- 0L
  chunks_processed <- 0L

  all_values_nonnegative <- TRUE
  all_values_integer_like <- TRUE

  start_time <- Sys.time()

  repeat {

    lines <- readLines(
      con,
      n = chunk_lines,
      warn = FALSE
    )

    if (length(lines) == 0L) {
      break
    }

    chunks_processed <- chunks_processed + 1L
    n_lines <- length(lines)

    numeric_text <- sub(
      "^[^,]*,",
      "",
      lines
    )

    values <- scan(
      text = paste(
        numeric_text,
        collapse = ","
      ),
      what = double(),
      sep = ",",
      quiet = TRUE
    )

    expected_n <- as.double(
      n_lines
    ) *
      as.double(
        n_beads
      )

    if (
      length(values) !=
        expected_n
    ) {
      stop(
        "数值解析长度异常：",
        counts_file,
        "\nchunk=",
        chunks_processed,
        "\n实际=",
        length(values),
        "\n预期=",
        expected_n
      )
    }

    if (anyNA(values)) {
      stop(
        "全基因count解析出现NA：",
        counts_file,
        "\nchunk=",
        chunks_processed
      )
    }

    if (any(
      values < 0
    )) {
      all_values_nonnegative <- FALSE
    }

    if (any(
      abs(
        values -
          round(values)
      ) >
        1e-8
    )) {
      all_values_integer_like <- FALSE
    }

    dim(values) <- c(
      n_beads,
      n_lines
    )

    library_size <- library_size +
      rowSums(
        values
      )

    detected_genes <- detected_genes +
      rowSums(
        values > 0
      )

    genes_processed <- genes_processed +
      n_lines

    if (
      chunks_processed %% 50L ==
        0L
    ) {
      elapsed_min <- as.numeric(
        difftime(
          Sys.time(),
          start_time,
          units = "mins"
        )
      )

      cat(
        progress_prefix,
        "chunks=",
        chunks_processed,
        " | genes=",
        genes_processed,
        " | elapsed=",
        round(
          elapsed_min,
          2
        ),
        " min\n",
        sep = ""
      )
    }

    rm(
      lines,
      numeric_text,
      values
    )

    gc(
      verbose = FALSE
    )
  }

  elapsed_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )

  target_sum <- rowSums(
    as.matrix(
      target[
        ,
        ..TARGET_GENES
      ]
    )
  )

  target_sum_not_exceed_library <- all(
    target_sum <=
      library_size
  )

  if (!target_sum_not_exceed_library) {
    stop(
      "14个目标基因count和出现大于全基因library size：",
      counts_file
    )
  }

  list(
    depth = data.table(
      bead_key =
        bead_key,
      library_size =
        library_size,
      detected_genes =
        detected_genes
    ),
    audit = list(
      n_beads =
        n_beads,
      genes_processed =
        genes_processed,
      chunks_processed =
        chunks_processed,
      elapsed_seconds =
        elapsed_seconds,
      all_values_nonnegative =
        all_values_nonnegative,
      all_values_integer_like =
        all_values_integer_like,
      target_sum_not_exceed_library =
        target_sum_not_exceed_library
    )
  )
}

banner("STEP C4B：14样本library-depth流式计算")

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

sample_audit_list <- vector(
  "list",
  nrow(inv)
)

global_start <- Sys.time()

for (i in seq_len(nrow(inv))) {

  sample_name <- inv$official_paper_id[i]
  cohort <- inv$cohort[i]
  counts_file <- inv$counts_file[i]

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
    "/",
    nrow(inv),
    "] ",
    sample_name,
    " (",
    cohort,
    ")\n",
    sep = ""
  )

  if (!file.exists(target_file)) {
    stop(
      "缺少STEP C3v2 target RDS：\n",
      target_file
    )
  }

  if (!file.exists(counts_file)) {
    stop(
      "counts文件不存在：\n",
      counts_file
    )
  }

  target <- readRDS(
    target_file
  )

  reused <- validate_existing_depth(
    depth_file = depth_file,
    target = target
  )

  if (reused) {

    depth <- readRDS(
      depth_file
    )

    cat(
      "  已存在合法depth RDS，直接复用。\n"
    )

    sample_audit_list[[i]] <- data.table(
      gsm = inv$gsm[i],
      sample = sample_name,
      cohort = cohort,
      n_beads = nrow(
        depth
      ),
      genes_processed = NA_integer_,
      chunks_processed = NA_integer_,
      chunk_lines = CHUNK_LINES,
      elapsed_seconds = 0,
      elapsed_minutes = 0,
      reused_existing = TRUE,
      min_library_size =
        min(
          depth$library_size
        ),
      median_library_size =
        median(
          depth$library_size
        ),
      mean_library_size =
        mean(
          depth$library_size
        ),
      max_library_size =
        max(
          depth$library_size
        ),
      zero_library_beads =
        sum(
          depth$library_size == 0
        ),
      nonzero_library_fraction =
        mean(
          depth$library_size > 0
        ),
      median_detected_genes =
        median(
          depth$detected_genes
        ),
      mean_detected_genes =
        mean(
          depth$detected_genes
        ),
      all_values_nonnegative = TRUE,
      all_values_integer_like = TRUE,
      target_sum_not_exceed_library = TRUE,
      output_rds = depth_file,
      output_rds_MB = round(
        file.info(
          depth_file
        )$size /
          1024^2,
        3
      )
    )

    rm(
      target,
      depth
    )

    gc(
      verbose = FALSE
    )

    next
  }

  res <- compute_depth_streaming(
    counts_file = counts_file,
    target = target,
    chunk_lines = CHUNK_LINES,
    progress_prefix = "  "
  )

  saveRDS(
    res$depth,
    depth_file,
    compress = "xz"
  )

  a <- res$audit

  sample_audit_list[[i]] <- data.table(
    gsm = inv$gsm[i],
    sample = sample_name,
    cohort = cohort,
    n_beads =
      a$n_beads,
    genes_processed =
      a$genes_processed,
    chunks_processed =
      a$chunks_processed,
    chunk_lines =
      CHUNK_LINES,
    elapsed_seconds =
      a$elapsed_seconds,
    elapsed_minutes =
      a$elapsed_seconds /
        60,
    reused_existing = FALSE,
    min_library_size =
      min(
        res$depth$library_size
      ),
    median_library_size =
      median(
        res$depth$library_size
      ),
    mean_library_size =
      mean(
        res$depth$library_size
      ),
    max_library_size =
      max(
        res$depth$library_size
      ),
    zero_library_beads =
      sum(
        res$depth$library_size == 0
      ),
    nonzero_library_fraction =
      mean(
        res$depth$library_size > 0
      ),
    median_detected_genes =
      median(
        res$depth$detected_genes
      ),
    mean_detected_genes =
      mean(
        res$depth$detected_genes
      ),
    all_values_nonnegative =
      a$all_values_nonnegative,
    all_values_integer_like =
      a$all_values_integer_like,
    target_sum_not_exceed_library =
      a$target_sum_not_exceed_library,
    output_rds = depth_file,
    output_rds_MB = round(
      file.info(
        depth_file
      )$size /
        1024^2,
      3
    )
  )

  cat(
    "  完成：",
    round(
      a$elapsed_seconds /
        60,
      2
    ),
    " min | median library=",
    round(
      median(
        res$depth$library_size
      ),
      1
    ),
    "\n",
    sep = ""
  )

  rm(
    target,
    res,
    a
  )

  gc(
    verbose = FALSE
  )
}

sample_audit <- rbindlist(
  sample_audit_list,
  fill = TRUE
)

setorder(
  sample_audit,
  cohort,
  sample
)

global_elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    global_start,
    units = "mins"
  )
)

# ------------------------------------------------------------
# Cohort-level descriptive QC only
# ------------------------------------------------------------

cohort_qc <- sample_audit[
  ,
  .(
    n_samples = .N,
    total_beads =
      sum(
        n_beads
      ),
    median_of_sample_median_library =
      median(
        median_library_size
      ),
    min_sample_median_library =
      min(
        median_library_size
      ),
    max_sample_median_library =
      max(
        median_library_size
      ),
    median_of_sample_median_detected_genes =
      median(
        median_detected_genes
      ),
    samples_with_zero_library_beads =
      sum(
        zero_library_beads > 0
      )
  ),
  by = cohort
]

qc_pass <- (
  nrow(sample_audit) == 14L &&
    all(
      sample_audit$n_beads > 0
    ) &&
    all(
      sample_audit$nonzero_library_fraction == 1
    ) &&
    all(
      sample_audit$all_values_nonnegative
    ) &&
    all(
      sample_audit$all_values_integer_like
    ) &&
    all(
      sample_audit$target_sum_not_exceed_library
    )
)

status <- if (qc_pass) {
  "PASS_ALL14_LIBRARY_DEPTH_STREAMING"
} else {
  "REVIEW_REQUIRED"
}

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------

fwrite(
  sample_audit,
  file.path(
    OUT_DIR,
    "STEP_C4B_01_all14_library_depth_sample_QC.tsv"
  ),
  sep = "\t"
)

fwrite(
  cohort_qc,
  file.path(
    OUT_DIR,
    "STEP_C4B_02_library_depth_cohort_descriptive_QC.tsv"
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
    "samples=14",
    paste0(
      "brain_metastasis_samples=",
      sum(
        sample_audit$cohort ==
          "Brain_metastasis"
      )
    ),
    paste0(
      "primary_tumor_samples=",
      sum(
        sample_audit$cohort ==
          "Primary_tumor"
      )
    ),
    paste0(
      "total_beads=",
      sum(
        sample_audit$n_beads
      )
    ),
    paste0(
      "reused_existing_samples=",
      sum(
        sample_audit$reused_existing
      )
    ),
    paste0(
      "newly_computed_samples=",
      sum(
        !sample_audit$reused_existing
      )
    ),
    paste0(
      "wall_clock_elapsed_minutes=",
      round(
        global_elapsed_minutes,
        4
      )
    ),
    paste0(
      "all_nonzero_library_fraction_1=",
      all(
        sample_audit$nonzero_library_fraction ==
          1
      )
    ),
    paste0(
      "all_values_integer_like=",
      all(
        sample_audit$all_values_integer_like
      )
    ),
    "normalization_or_scoring_performed=FALSE",
    "statistical_testing_performed=FALSE"
  ),
  file.path(
    OUT_DIR,
    "STEP_C4B_COMPLETE.txt"
  )
)

banner("STEP C4B结果")

print(
  sample_audit[
    ,
    .(
      sample,
      cohort,
      n_beads,
      reused_existing,
      elapsed_minutes,
      median_library_size,
      median_detected_genes,
      zero_library_beads
    )
  ]
)

cat(
  "\nCohort descriptive QC：\n"
)

print(
  cohort_qc
)

cat(
  "\nStatus: ",
  status,
  "\n",
  "Elapsed wall time: ",
  round(
    global_elapsed_minutes,
    2
  ),
  " minutes\n\n",
  "Key output files:\n",
  file.path(
    OUT_DIR,
    "STEP_C4B_01_all14_library_depth_sample_QC.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C4B_02_library_depth_cohort_descriptive_QC.tsv"
  ),
  "\n",
  sep = ""
)
