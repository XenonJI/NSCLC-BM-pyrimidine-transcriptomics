# GSE223499 fixed-module extraction
# Streams the fixed target genes from per-sample count files and constructs the initial
# specimen-level five-gene summaries. These are intermediate inputs to the corrected
# patient-aware analysis below. KRAS_17 is repaired with script 02, then this script is rerun.
# Set ROOT to the local project directory before running.

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  timeout = 7200,
  scipen = 999
)

ROOT <- "D:/LUAD_LM_PYRIMIDINE"

RAW_DIR <- file.path(
  ROOT,
  "data/raw/GSE223499"
)

COUNTS_DIR <- file.path(
  RAW_DIR,
  "counts_csv"
)

META_FILE <- file.path(
  RAW_DIR,
  "GSE223499_sn_integrated_data.csv.gz"
)

TABLE_DIR <- file.path(
  ROOT,
  "results/tables/GSE223499"
)

FIGURE_DIR <- file.path(
  ROOT,
  "results/figures/GSE223499"
)

LOCK_DIR <- file.path(
  ROOT,
  "results/locks"
)

PROCESSED_DIR <- file.path(
  ROOT,
  "data/processed"
)

CHECKPOINT_FILE <- file.path(
  PROCESSED_DIR,
  "GSE223499_STEP07C_sample_gene_counts_checkpoint.rds"
)

for (current_dir in c(
  TABLE_DIR,
  FIGURE_DIR,
  LOCK_DIR,
  PROCESSED_DIR
)) {
  dir.create(
    current_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

required_packages <- c(
  "data.table",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nInstall once with: install.packages(c(",
    paste(
      sprintf('"%s"', missing_packages),
      collapse = ", "
    ),
    "))"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

banner <- function(text) {
  cat(
    "\n========================================\n",
    text,
    "\n========================================\n",
    sep = ""
  )
}

strip_quotes <- function(x) {
  result <- trimws(
    as.character(x)
  )

  result <- sub(
    '^"',
    "",
    result
  )

  result <- sub(
    '"$',
    "",
    result
  )

  result
}

canonical_barcode <- function(x) {
  result <- strip_quotes(x)

  result <- sub(
    "^X(?=[ACGTN]{10,})",
    "",
    result,
    perl = TRUE
  )

  result <- sub(
    "\\.([0-9]+)$",
    "-\\1",
    result,
    perl = TRUE
  )

  extracted <- regmatches(
    result,
    regexpr(
      "[ACGTN]{10,}-[0-9]+$",
      result,
      perl = TRUE
    )
  )

  has_extract <- nzchar(
    extracted
  )

  result[
    has_extract
  ] <- extracted[
    has_extract
  ]

  result
}

barcode_core_sequence <- function(x) {
  result <- canonical_barcode(x)

  extracted <- regmatches(
    result,
    regexpr(
      "[ACGTN]{10,}",
      result,
      perl = TRUE
    )
  )

  extracted[
    !nzchar(extracted)
  ] <- NA_character_

  extracted
}

is_gzip_file <- function(path) {
  if (
    !file.exists(path) ||
      file.info(path)$size < 2
  ) {
    return(FALSE)
  }

  connection <- file(
    path,
    open = "rb"
  )

  on.exit(
    close(connection),
    add = TRUE
  )

  header <- readBin(
    connection,
    what = "raw",
    n = 2L
  )

  identical(
    as.integer(header),
    c(31L, 139L)
  )
}

safe_zscore <- function(x) {
  x <- as.numeric(x)

  current_sd <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(current_sd) ||
      current_sd <= 0
  ) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) /
    current_sd
}

hedges_g_two_group <- function(
  brain,
  primary
) {
  brain <- brain[
    is.finite(brain)
  ]

  primary <- primary[
    is.finite(primary)
  ]

  n1 <- length(brain)
  n0 <- length(primary)

  if (
    n1 < 2L ||
      n0 < 2L
  ) {
    return(NA_real_)
  }

  pooled_variance <- (
    (n1 - 1) *
      stats::var(brain) +
      (n0 - 1) *
        stats::var(primary)
  ) /
    (
      n1 +
        n0 -
        2
    )

  if (
    !is.finite(pooled_variance) ||
      pooled_variance <= 0
  ) {
    return(NA_real_)
  }

  d <- (
    mean(brain) -
      mean(primary)
  ) /
    sqrt(
      pooled_variance
    )

  df <- n1 +
    n0 -
    2

  correction <- 1 -
    3 /
      (
        4 *
          df -
          1
      )

  correction *
    d
}

one_sided_wilcox <- function(
  brain,
  primary
) {
  result <- tryCatch(
    stats::wilcox.test(
      x = brain,
      y = primary,
      alternative = "greater",
      exact = FALSE,
      correct = FALSE
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(result)) {
    return(NA_real_)
  }

  as.numeric(
    result$p.value
  )
}

one_sided_lm_p <- function(
  estimate,
  two_sided_p
) {
  if (
    !is.finite(estimate) ||
      !is.finite(two_sided_p)
  ) {
    return(NA_real_)
  }

  if (estimate > 0) {
    two_sided_p / 2
  } else {
    1 -
      two_sided_p / 2
  }
}

derive_batch <- function(sample) {
  fcase(
    grepl(
      "^PA",
      sample
    ),
    "PA",
    grepl(
      "^STK",
      sample
    ),
    "STK",
    grepl(
      "^KRAS",
      sample
    ),
    "KRAS",
    grepl(
      "^N",
      sample
    ),
    "N",
    default = "Other"
  )
}

extract_raw_barcode <- function(
  cell_id,
  sample
) {
  prefix <- paste0(
    sample,
    "_"
  )

  starts <- startsWith(
    cell_id,
    prefix
  )

  result <- cell_id

  result[starts] <- substring(
    cell_id[starts],
    nchar(
      prefix[starts]
    ) +
      1L
  )

  result
}

# ============================================================
# 阶段0：固定42份文件清单并检查下载完整性
# ============================================================

banner(
  "阶段0：检查42份counts文件"
)

manifest <- data.table(
  gsm = c(
    "GSM6957590", "GSM6957591", "GSM6957592",
    "GSM6957593", "GSM6957594", "GSM6957595",
    "GSM6957596", "GSM6957597", "GSM6957598",
    "GSM6957599", "GSM6957600", "GSM6957601",
    "GSM6957602", "GSM6957603", "GSM6957604",
    "GSM6957605", "GSM6957606", "GSM6957607",
    "GSM6957608", "GSM6957609", "GSM6957610",
    "GSM6957611", "GSM6957612",
    "GSM6957613", "GSM6957614", "GSM6957615",
    "GSM6957616", "GSM6957617", "GSM6957618",
    "GSM6957619", "GSM6957620", "GSM6957621",
    "GSM6957622", "GSM6957623",
    "GSM6957625", "GSM6957626", "GSM6957627",
    "GSM6957628",
    "GSM6957629", "GSM6957630", "GSM6957631",
    "GSM6957632"
  ),
  sample = c(
    "PA001", "PA004", "PA005", "PA019", "PA025",
    "PA034", "PA042", "PA043", "PA048", "PA054",
    "PA056", "PA060", "PA067", "PA068", "PA070",
    "PA072", "PA076", "PA080", "PA104", "PA125",
    "PA141",
    "N254", "N586",
    "STK_1", "STK_3", "STK_5dot1", "STK_5dot2",
    "STK_14", "STK_18",
    "STK_2", "STK_15", "STK_20", "STK_21",
    "STK_22dot2",
    "KRAS_6", "KRAS_7", "KRAS_8", "KRAS_17",
    "KRAS_10", "KRAS_11", "KRAS_12", "KRAS_13"
  )
)

manifest[
  ,
  filename := paste0(
    gsm,
    "_NSCLC_",
    sample,
    "_sn_counts.csv.gz"
  )
]

manifest[
  ,
  path := file.path(
    COUNTS_DIR,
    filename
  )
]

manifest[
  ,
  exists := file.exists(
    path
  )
]

manifest[
  ,
  bytes := fifelse(
    exists,
    as.numeric(
      file.info(path)$size
    ),
    NA_real_
  )
]

manifest[
  ,
  gzip_valid := vapply(
    path,
    is_gzip_file,
    logical(1)
  )
]

manifest[
  ,
  valid :=
    exists &
      bytes >=
        100 *
          1024 &
      gzip_valid
]

fwrite(
  manifest,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_counts_file_audit.tsv"
  ),
  sep = "\t"
)

cat(
  "有效文件：",
  manifest[
    valid == TRUE,
    .N
  ],
  "/42。\n",
  sep = ""
)

if (
  manifest[
    valid != TRUE,
    .N
  ] >
    0L
) {
  cat(
    "\n缺失或无效文件：\n"
  )

  print(
    manifest[
      valid != TRUE,
      .(
        gsm,
        sample,
        filename,
        exists,
        bytes,
        gzip_valid
      )
    ]
  )

  stop(
    "42份counts尚未全部通过检查。"
  )
}

if (!file.exists(META_FILE)) {
  stop(
    "找不到整合注释：\n",
    META_FILE
  )
}

manifest_signature <- paste(
  manifest$filename,
  manifest$bytes,
  collapse = "|"
)

# ============================================================
# 阶段1：重建固定作者肿瘤细胞队列
# ============================================================

banner(
  "阶段1：重建固定患者与作者肿瘤细胞"
)

metadata <- fread(
  META_FILE,
  showProgress = TRUE,
  data.table = TRUE,
  check.names = FALSE
)

required_columns <- c(
  "V1",
  "orig.ident",
  "PRIMARY vs BRAIN_METS vs CHEST_WALL_MET",
  "patient",
  "patient_long",
  "predicted_doublets",
  "nCount_RNA",
  "tumor_nontumor_major"
)

missing_columns <- setdiff(
  required_columns,
  names(metadata)
)

if (length(missing_columns) > 0L) {
  stop(
    "作者整合注释缺少字段：",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

setnames(
  metadata,
  old = c(
    "V1",
    "orig.ident",
    "PRIMARY vs BRAIN_METS vs CHEST_WALL_MET",
    "tumor_nontumor_major"
  ),
  new = c(
    "cell_id",
    "sample",
    "site_author",
    "tumor_major"
  )
)

metadata[
  ,
  `:=`(
    cell_id =
      trimws(
        as.character(cell_id)
      ),
    sample =
      trimws(
        as.character(sample)
      ),
    patient =
      trimws(
        as.character(patient)
      ),
    patient_long =
      trimws(
        as.character(patient_long)
      ),
    site_author =
      trimws(
        as.character(site_author)
      ),
    tumor_major =
      trimws(
        as.character(tumor_major)
      )
  )
]

metadata[
  patient %in%
    c(
      "",
      "NA",
      "<NA>"
    ),
  patient := NA_character_
]

metadata[
  ,
  patient_resolved := fifelse(
    is.na(patient),
    paste0(
      "SAMPLE_",
      sample
    ),
    patient
  )
]

metadata[
  ,
  cohort := fcase(
    site_author == "BRAIN_METS",
    "Brain_metastasis",
    site_author == "PRIMARY",
    "Primary_tumor",
    site_author == "CHEST_WALL_MET",
    "Chest_wall_metastasis",
    default =
      "Unresolved"
  )
]

metadata[
  ,
  batch := derive_batch(
    sample
  )
]

doublet_flag <- tolower(
  trimws(
    as.character(
      metadata$predicted_doublets
    )
  )
)

metadata[
  ,
  is_doublet :=
    doublet_flag %in%
      c(
        "true",
        "t",
        "1",
        "yes"
      )
]

metadata[
  ,
  analysis_nucleus :=
    cohort %in%
      c(
        "Brain_metastasis",
        "Primary_tumor"
      ) &
      tumor_major == "Tumor" &
      !is_doublet
]

metadata[
  ,
  raw_barcode :=
    extract_raw_barcode(
      cell_id,
      sample
    )
]

metadata[
  ,
  canonical_barcode :=
    canonical_barcode(
      raw_barcode
    )
]

sample_audit <- metadata[
  sample %in%
    manifest$sample,
  .(
    total_integrated_nuclei =
      .N,
    author_tumor_nuclei =
      sum(
        analysis_nucleus
      ),
    selected_library_size =
      sum(
        as.numeric(
          nCount_RNA[
            analysis_nucleus
          ]
        ),
        na.rm = TRUE
      )
  ),
  by = .(
    sample,
    patient =
      patient_resolved,
    cohort,
    batch
  )
]

sample_audit[
  ,
  primary_eligible :=
    author_tumor_nuclei >=
      50L
]

setorder(
  sample_audit,
  cohort,
  patient
)

fwrite(
  sample_audit,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_fixed_sample_audit.tsv"
  ),
  sep = "\t"
)

print(
  sample_audit
)

if (
  sample_audit[
    cohort == "Brain_metastasis",
    .N
  ] != 31L ||
    sample_audit[
      cohort == "Primary_tumor",
      .N
    ] != 11L
) {
  stop(
    "队列数量不是31例脑转移和11例原发肿瘤。"
  )
}

if (
  sample_audit[
    primary_eligible == FALSE,
    .N
  ] != 1L ||
    sample_audit[
      primary_eligible == FALSE,
      sample
    ] !=
      "STK_20"
) {
  stop(
    "低于50个作者肿瘤细胞核的样本不是唯一的STK_20。"
  )
}

# ============================================================
# 阶段2：在读取表达结果前冻结最终统计细节
# ============================================================

banner(
  "阶段2：冻结STEP07-C统计细节"
)

analysis_lock <- file.path(
  LOCK_DIR,
  "STEP_07C_GSE223499_ANALYSIS_LOCK.txt"
)

if (!file.exists(analysis_lock)) {
  writeLines(
    c(
      paste0(
        "created_at=",
        format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        )
      ),
      "expression_target_rows_scanned_at_lock_creation=FALSE",
      "fixed_module_genes=DHFR,DHODH,SHMT1,TYMS,UMPS",
      "primary_analysis=31_brain_vs_10_primary_with_at_least_50_author_tumor_nuclei",
      "all_sample_sensitivity=31_brain_vs_11_primary_including_STK_20",
      "primary_unit=resolved_patient",
      "normalization=log2_CPM_plus_1_using_sum_nCount_RNA_of_fixed_author_tumor_nuclei",
      "gene_standardization=within_each_analysis_population",
      "module_score=mean_z_score_of_the_five_fixed_genes",
      "primary_test=one_sided_Wilcoxon_rank_sum_without_continuity_correction",
      "primary_direction=Brain_metastasis_greater_than_Primary_tumor",
      "primary_p_threshold=0.05",
      "primary_Hedges_g_threshold=0.40",
      "direction_support=at_least_65_percent_of_brain_patients_above_primary_median",
      "cell_cycle_adjustment=residual_from_module_score_on_fixed_S_and_G2M_scores",
      "strict_validation_requires=primary_gate_pass_and_adjusted_one_sided_p_at_most_0.05",
      "batch_sensitivity=KRAS_and_STK_overlap_only_linear_model_with_batch_S_and_G2M",
      "batch_sensitivity_is_supportive_not_part_of_primary_gate",
      "individual_gene_tests=exploratory_with_BH_FDR",
      "barcode_primary_mapping=canonical_10x_barcode_with_suffix",
      "barcode_fallback=ignore_library_suffix_only_when_core_sequences_are_unique_and_100_percent_matched",
      "outcome_driven_exclusion=not_allowed",
      "threshold_or_gene_change_after_results=not_allowed",
      paste0(
        "manifest_signature=",
        manifest_signature
      )
    ),
    analysis_lock
  )
}

cat(
  "统计锁：\n",
  analysis_lock,
  "\n",
  sep = ""
)

# ============================================================
# 阶段3：固定基因
# ============================================================

banner(
  "阶段3：固定5基因和细胞周期基因"
)

module_genes <- c(
  "DHFR",
  "DHODH",
  "SHMT1",
  "TYMS",
  "UMPS"
)

s_genes <- c(
  "MCM5", "PCNA", "TYMS", "FEN1", "MCM2",
  "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
  "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU",
  "HELLS", "RFC2", "POLR1B", "NASP", "RAD51AP1",
  "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
  "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
  "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1",
  "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1",
  "CHAF1B", "BRIP1", "E2F8"
)

g2m_genes <- c(
  "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5",
  "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
  "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3",
  "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
  "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B",
  "GTSE1", "KIF20B", "HJURP", "CDCA3", "JPT1",
  "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
  "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2",
  "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
  "LBR", "CKAP5", "CENPE", "CTCF", "NEK2",
  "G2E3", "GAS2L3", "CBX5", "CENPA"
)

target_genes <- unique(
  toupper(
    c(
      module_genes,
      s_genes,
      g2m_genes
    )
  )
)

cat(
  "计划提取：",
  length(target_genes),
  "个唯一基因。\n",
  sep = ""
)

# ============================================================
# 阶段4：初始化或读取断点检查点
# ============================================================

barcode_amendment_file <- file.path(
  LOCK_DIR,
  "STEP_07C_GSE223499_BARCODE_MAPPING_AMENDMENT.txt"
)

if (!file.exists(barcode_amendment_file)) {
  writeLines(
    c(
      paste0(
        "created_at=",
        format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        )
      ),
      "reason=KRAS_17_full_barcode_suffix_mismatch_detected_before_gene_counts_were_extracted_for_that_sample",
      "biological_hypothesis_changed=FALSE",
      "module_genes_changed=FALSE",
      "statistical_gate_changed=FALSE",
      "sample_exclusion_changed=FALSE",
      "primary_mapping=canonical_10x_barcode_with_suffix",
      "technical_fallback=ignore_suffix_only_when_core_sequences_are_unique_in_both_metadata_and_counts_and_match_100_percent",
      "fallback_is_applied_before_reading_target_gene_rows_for_the_affected_sample"
    ),
    barcode_amendment_file
  )
}

banner(
  "阶段4：初始化逐样本低内存检查点"
)

checkpoint <- list(
  version =
    "STEP07C_v1",
  manifest_signature =
    manifest_signature,
  processed =
    list()
)

if (file.exists(CHECKPOINT_FILE)) {
  existing_checkpoint <- tryCatch(
    readRDS(
      CHECKPOINT_FILE
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    is.list(existing_checkpoint) &&
      identical(
        existing_checkpoint$version,
        "STEP07C_v1"
      ) &&
      identical(
        existing_checkpoint$manifest_signature,
        manifest_signature
      ) &&
      is.list(
        existing_checkpoint$processed
      )
  ) {
    checkpoint <-
      existing_checkpoint

    cat(
      "检测到有效检查点；已处理样本：",
      length(
        checkpoint$processed
      ),
      "。\n",
      sep = ""
    )
  } else {
    cat(
      "旧检查点不匹配，将从头扫描。\n"
    )
  }
}

# ============================================================
# 阶段5：逐样本流式提取固定基因
# ============================================================

banner(
  "阶段5：逐样本流式扫描42份counts"
)

process_one_sample <- function(
  sample_name,
  counts_path
) {
  sample_meta <- metadata[
    sample == sample_name
  ]

  target_meta <- sample_meta[
    analysis_nucleus == TRUE
  ]

  if (nrow(target_meta) < 1L) {
    stop(
      sample_name,
      "没有固定作者肿瘤细胞核。"
    )
  }

  if (
    anyDuplicated(
      target_meta$canonical_barcode
    )
  ) {
    stop(
      sample_name,
      "的目标完整barcode存在重复。"
    )
  }

  connection <- gzfile(
    counts_path,
    open = "rt"
  )

  header_line <- readLines(
    connection,
    n = 1L,
    warn = FALSE
  )

  if (length(header_line) != 1L) {
    close(connection)

    stop(
      sample_name,
      "无法读取counts表头。"
    )
  }

  header_fields <- strsplit(
    header_line,
    split = ",",
    fixed = TRUE
  )[[1]]

  header_fields <- strip_quotes(
    header_fields
  )

  header_barcodes <-
    canonical_barcode(
      header_fields[-1]
    )

  if (
    anyDuplicated(
      header_barcodes
    )
  ) {
    close(connection)

    stop(
      sample_name,
      "的counts完整表头标准化后存在重复barcode。"
    )
  }

  selected_zero_based <- match(
    target_meta$canonical_barcode,
    header_barcodes
  )

  matched_n <- sum(
    !is.na(
      selected_zero_based
    )
  )

  match_fraction <- matched_n /
    nrow(
      target_meta
    )

  match_strategy <-
    "canonical_barcode_with_suffix"

  full_suffix_match_fraction <-
    match_fraction

  if (match_fraction < 0.98) {
    target_core <-
      barcode_core_sequence(
        target_meta$canonical_barcode
      )

    header_core <-
      barcode_core_sequence(
        header_barcodes
      )

    diagnostic <- data.table(
      sample =
        sample_name,
      metadata_target_cells =
        nrow(
          target_meta
        ),
      counts_header_cells =
        length(
          header_barcodes
        ),
      full_suffix_matches =
        matched_n,
      full_suffix_match_fraction =
        full_suffix_match_fraction,
      target_core_missing =
        sum(
          is.na(
            target_core
          )
        ),
      header_core_missing =
        sum(
          is.na(
            header_core
          )
        ),
      target_core_duplicates =
        anyDuplicated(
          target_core
        ),
      header_core_duplicates =
        anyDuplicated(
          header_core
        )
    )

    fwrite(
      diagnostic,
      file.path(
        TABLE_DIR,
        paste0(
          "GSE223499_STEP07C_",
          sample_name,
          "_barcode_fallback_diagnostic.tsv"
        )
      ),
      sep = "\t"
    )

    if (
      any(
        is.na(
          target_core
        )
      ) ||
        any(
          is.na(
            header_core
          )
        )
    ) {
      close(connection)

      stop(
        sample_name,
        "无法从全部barcode中提取核心碱基序列。"
      )
    }

    if (
      anyDuplicated(
        target_core
      ) ||
        anyDuplicated(
          header_core
        )
    ) {
      close(connection)

      stop(
        sample_name,
        "忽略10x后缀后出现重复核心barcode，不能安全映射。"
      )
    }

    core_selected_zero_based <- match(
      target_core,
      header_core
    )

    core_matched_n <- sum(
      !is.na(
        core_selected_zero_based
      )
    )

    core_match_fraction <-
      core_matched_n /
        nrow(
          target_meta
        )

    suffix_audit <- data.table(
      metadata_canonical_barcode =
        target_meta$canonical_barcode,
      metadata_core =
        target_core,
      counts_column_index =
        core_selected_zero_based,
      counts_canonical_barcode =
        fifelse(
          is.na(
            core_selected_zero_based
          ),
          NA_character_,
          header_barcodes[
            core_selected_zero_based
          ]
        ),
      core_matched =
        !is.na(
          core_selected_zero_based
        )
    )

    fwrite(
      suffix_audit,
      file.path(
        TABLE_DIR,
        paste0(
          "GSE223499_STEP07C_",
          sample_name,
          "_suffix_ignored_mapping.tsv.gz"
        )
      ),
      sep = "\t",
      compress = "gzip"
    )

    if (
      core_match_fraction <
        1
    ) {
      close(connection)

      stop(
        sample_name,
        "忽略后缀后的核心barcode仍未达到100%匹配：",
        round(
          100 *
            core_match_fraction,
          3
        ),
        "%"
      )
    }

    selected_zero_based <-
      core_selected_zero_based

    matched_n <-
      core_matched_n

    match_fraction <-
      core_match_fraction

    match_strategy <-
      "unique_core_sequence_suffix_ignored"

    cat(
      sample_name,
      "：完整后缀匹配率为",
      round(
        100 *
          full_suffix_match_fraction,
        3
      ),
      "%；核心barcode唯一且100%匹配，采用后缀忽略映射。\n",
      sep = ""
    )
  }

  selected_field_positions <-
    selected_zero_based +
      1L

  gene_counts <- setNames(
    rep(
      0,
      length(
        target_genes
      )
    ),
    target_genes
  )

  gene_rows_found <- setNames(
    rep(
      0L,
      length(
        target_genes
      )
    ),
    target_genes
  )

  scanned_rows <- 0L

  repeat {
    lines <- readLines(
      connection,
      n = 500L,
      warn = FALSE
    )

    if (length(lines) == 0L) {
      break
    }

    scanned_rows <- scanned_rows +
      length(lines)

    first_fields <- sub(
      ",.*$",
      "",
      lines
    )

    first_fields <- toupper(
      strip_quotes(
        first_fields
      )
    )

    target_index <- which(
      first_fields %in%
        target_genes
    )

    if (length(target_index) > 0L) {
      for (line_index in target_index) {
        current_gene <-
          first_fields[line_index]

        fields <- strsplit(
          lines[line_index],
          split = ",",
          fixed = TRUE
        )[[1]]

        if (
          max(
            selected_field_positions
          ) >
            length(fields)
        ) {
          close(connection)

          stop(
            sample_name,
            "的基因行列数不足：",
            current_gene
          )
        }

        selected_values <- suppressWarnings(
          as.numeric(
            strip_quotes(
              fields[
                selected_field_positions
              ]
            )
          )
        )

        if (
          any(
            !is.finite(
              selected_values
            )
          )
        ) {
          close(connection)

          stop(
            sample_name,
            "基因",
            current_gene,
            "出现不可解析计数。"
          )
        }

        gene_counts[current_gene] <-
          gene_counts[current_gene] +
            sum(
              selected_values
            )

        gene_rows_found[current_gene] <-
          gene_rows_found[current_gene] +
            1L
      }
    }
  }

  close(connection)

  sample_info <- sample_audit[
    sample == sample_name
  ]

  list(
    sample =
      sample_name,
    patient =
      sample_info$patient,
    cohort =
      sample_info$cohort,
    batch =
      sample_info$batch,
    author_tumor_nuclei =
      sample_info$author_tumor_nuclei,
    selected_library_size =
      sample_info$selected_library_size,
    primary_eligible =
      sample_info$primary_eligible,
    header_cell_columns =
      length(
        header_barcodes
      ),
    matched_tumor_barcodes =
      matched_n,
    barcode_match_fraction =
      match_fraction,
    barcode_match_strategy =
      match_strategy,
    full_suffix_match_fraction =
      full_suffix_match_fraction,
    scanned_gene_rows =
      scanned_rows,
    gene_counts =
      gene_counts,
    gene_rows_found =
      gene_rows_found
  )
}

for (
  sample_index in
    seq_len(
      nrow(
        manifest
      )
    )
) {
  current_sample <-
    manifest$sample[
      sample_index
    ]

  current_path <-
    manifest$path[
      sample_index
    ]

  if (
    current_sample %in%
      names(
        checkpoint$processed
      )
  ) {
    cat(
      "[",
      sample_index,
      "/42] 跳过检查点样本：",
      current_sample,
      "\n",
      sep = ""
    )

    next
  }

  cat(
    "\n[",
    sample_index,
    "/42] 正在处理：",
    current_sample,
    "\n",
    sep = ""
  )

  current_result <- process_one_sample(
    sample_name =
      current_sample,
    counts_path =
      current_path
  )

  checkpoint$processed[[current_sample]] <- current_result

  saveRDS(
    checkpoint,
    CHECKPOINT_FILE
  )

  cat(
    "完成：目标细胞",
    current_result$matched_tumor_barcodes,
    "；barcode匹配率",
    round(
      100 *
        current_result$barcode_match_fraction,
      2
    ),
    "%；策略=",
    current_result$barcode_match_strategy,
    "；扫描基因行",
    format(
      current_result$scanned_gene_rows,
      big.mark = ","
    ),
    "。\n",
    sep = ""
  )

  gc(
    verbose = FALSE
  )
}

if (
  length(
    checkpoint$processed
  ) != 42L
) {
  stop(
    "检查点中不是42个处理完成样本。"
  )
}

# ============================================================
# 阶段6：合并样本并形成患者级伪批量表达
# ============================================================

banner(
  "阶段6：患者级伪批量表达"
)

sample_result_rows <- rbindlist(
  lapply(
    checkpoint$processed,
    function(current_result) {
      data.table(
        sample =
          current_result$sample,
        patient =
          current_result$patient,
        cohort =
          current_result$cohort,
        batch =
          current_result$batch,
        author_tumor_nuclei =
          current_result$author_tumor_nuclei,
        selected_library_size =
          current_result$selected_library_size,
        primary_eligible =
          current_result$primary_eligible,
        header_cell_columns =
          current_result$header_cell_columns,
        matched_tumor_barcodes =
          current_result$matched_tumor_barcodes,
        barcode_match_fraction =
          current_result$barcode_match_fraction,
        barcode_match_strategy =
          if (
            is.null(
              current_result$barcode_match_strategy
            )
          ) {
            "canonical_barcode_with_suffix"
          } else {
            current_result$barcode_match_strategy
          },
        full_suffix_match_fraction =
          if (
            is.null(
              current_result$full_suffix_match_fraction
            )
          ) {
            current_result$barcode_match_fraction
          } else {
            current_result$full_suffix_match_fraction
          },
        scanned_gene_rows =
          current_result$scanned_gene_rows
      )
    }
  ),
  fill = TRUE
)

sample_gene_counts <- rbindlist(
  lapply(
    checkpoint$processed,
    function(current_result) {
      data.table(
        sample =
          current_result$sample,
        patient =
          current_result$patient,
        cohort =
          current_result$cohort,
        batch =
          current_result$batch,
        primary_eligible =
          current_result$primary_eligible,
        gene =
          names(
            current_result$gene_counts
          ),
        raw_count =
          as.numeric(
            current_result$gene_counts
          ),
        rows_found =
          as.integer(
            current_result$gene_rows_found
          ),
        selected_library_size =
          current_result$selected_library_size,
        author_tumor_nuclei =
          current_result$author_tumor_nuclei
      )
    }
  ),
  fill = TRUE
)

fwrite(
  sample_result_rows,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_sample_processing_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  sample_gene_counts,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_sample_target_gene_counts.tsv.gz"
  ),
  sep = "\t",
  compress = "gzip"
)

patient_gene_counts <- sample_gene_counts[
  ,
  .(
    raw_count =
      sum(
        raw_count
      ),
    rows_found =
      sum(
        rows_found
      ),
    selected_library_size =
      sum(
        selected_library_size
      ),
    author_tumor_nuclei =
      sum(
        author_tumor_nuclei
      ),
    sample_count =
      uniqueN(sample),
    samples =
      paste(
        sort(
          unique(sample)
        ),
        collapse = ";"
      ),
    primary_eligible =
      all(
        primary_eligible
      )
  ),
  by = .(
    patient,
    cohort,
    batch,
    gene
  )
]

patient_gene_counts[
  ,
  log2_CPM_plus1 :=
    log2(
      raw_count /
        selected_library_size *
        1000000 +
        1
    )
]

fwrite(
  patient_gene_counts,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_patient_gene_expression_long.tsv"
  ),
  sep = "\t"
)

gene_coverage <- patient_gene_counts[
  ,
  .(
    patients =
      uniqueN(
        patient
      ),
    patients_with_row_found =
      uniqueN(
        patient[
          rows_found >
            0L
        ]
      ),
    total_raw_count =
      sum(
        raw_count
      ),
    patients_with_nonzero_count =
      uniqueN(
        patient[
          raw_count >
            0
        ]
      )
  ),
  by = gene
][
  order(
    -patients_with_row_found,
    gene
  )
]

fwrite(
  gene_coverage,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_target_gene_coverage.tsv"
  ),
  sep = "\t"
)

module_coverage <- gene_coverage[
  gene %in%
    module_genes
]

cat(
  "\n固定5基因覆盖：\n"
)

print(
  module_coverage
)

if (
  module_coverage[
    patients_with_row_found >
      0L,
    uniqueN(gene)
  ] != 5L
) {
  stop(
    "固定5基因中至少一个在全部样本均未检测到。"
  )
}

patient_info <- unique(
  patient_gene_counts[
    ,
    .(
      patient,
      cohort,
      batch,
      author_tumor_nuclei,
      primary_eligible,
      selected_library_size,
      sample_count,
      samples
    )
  ]
)

expression_wide <- dcast(
  patient_gene_counts,
  patient +
    cohort +
    batch +
    author_tumor_nuclei +
    primary_eligible +
    selected_library_size +
    sample_count +
    samples ~
    gene,
  value.var =
    "log2_CPM_plus1"
)

# ============================================================
# 阶段7：计算主要和敏感性评分
# ============================================================

banner(
  "阶段7：计算患者级固定模块评分"
)

calculate_scores <- function(
  wide_data,
  analysis_name
) {
  result <- copy(
    wide_data
  )

  module_matrix <- as.matrix(
    result[
      ,
      ..module_genes
    ]
  )

  module_z <- apply(
    module_matrix,
    2L,
    safe_zscore
  )

  if (is.null(dim(module_z))) {
    module_z <- matrix(
      module_z,
      ncol =
        length(
          module_genes
        )
    )
  }

  colnames(
    module_z
  ) <- module_genes

  result[
    ,
    module_score :=
      rowMeans(
        module_z,
        na.rm = TRUE
      )
  ]

  detected_s <- intersect(
    s_genes,
    names(result)
  )

  detected_g2m <- intersect(
    g2m_genes,
    names(result)
  )

  variable_s <- detected_s[
    vapply(
      result[
        ,
        ..detected_s
      ],
      function(x) {
        stats::sd(
          x,
          na.rm = TRUE
        ) >
          0
      },
      logical(1)
    )
  ]

  variable_g2m <- detected_g2m[
    vapply(
      result[
        ,
        ..detected_g2m
      ],
      function(x) {
        stats::sd(
          x,
          na.rm = TRUE
        ) >
          0
      },
      logical(1)
    )
  ]

  if (
    length(variable_s) <
      10L ||
      length(variable_g2m) <
        10L
  ) {
    stop(
      analysis_name,
      "可用细胞周期基因过少。"
    )
  }

  s_matrix <- as.matrix(
    result[
      ,
      ..variable_s
    ]
  )

  g2m_matrix <- as.matrix(
    result[
      ,
      ..variable_g2m
    ]
  )

  s_z <- apply(
    s_matrix,
    2L,
    safe_zscore
  )

  g2m_z <- apply(
    g2m_matrix,
    2L,
    safe_zscore
  )

  result[
    ,
    s_score :=
      rowMeans(
        s_z,
        na.rm = TRUE
      )
  ]

  result[
    ,
    g2m_score :=
      rowMeans(
        g2m_z,
        na.rm = TRUE
      )
  ]

  adjustment_model <- stats::lm(
    module_score ~
      s_score +
      g2m_score,
    data = result
  )

  result[
    ,
    adjusted_module_score :=
      stats::residuals(
        adjustment_model
      )
  ]

  result[
    ,
    analysis :=
      analysis_name
  ]

  attr(
    result,
    "s_genes_used"
  ) <- variable_s

  attr(
    result,
    "g2m_genes_used"
  ) <- variable_g2m

  result
}

main_wide <- expression_wide[
  primary_eligible == TRUE
]

all_wide <- copy(
  expression_wide
)

if (
  main_wide[
    cohort == "Brain_metastasis",
    .N
  ] != 31L ||
    main_wide[
      cohort == "Primary_tumor",
      .N
    ] != 10L
) {
  stop(
    "主分析不是31例脑转移对10例原发。"
  )
}

main_scores <- calculate_scores(
  main_wide,
  "Primary_31_vs_10"
)

all_scores <- calculate_scores(
  all_wide,
  "Sensitivity_31_vs_11"
)

patient_scores <- rbindlist(
  list(
    main_scores,
    all_scores
  ),
  fill = TRUE
)

score_columns <- c(
  "analysis",
  "patient",
  "samples",
  "cohort",
  "batch",
  "author_tumor_nuclei",
  "primary_eligible",
  "selected_library_size",
  "module_score",
  "s_score",
  "g2m_score",
  "adjusted_module_score"
)

fwrite(
  patient_scores[
    ,
    ..score_columns
  ],
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_patient_level_scores.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 阶段8：预注册验证闸门
# ============================================================

banner(
  "阶段8：正式验证闸门"
)

calculate_gate <- function(
  score_data,
  endpoint,
  score_column
) {
  brain <- score_data[
    cohort ==
      "Brain_metastasis",
    get(
      score_column
    )
  ]

  primary <- score_data[
    cohort ==
      "Primary_tumor",
    get(
      score_column
    )
  ]

  primary_median <- stats::median(
    primary,
    na.rm = TRUE
  )

  direction_count <- sum(
    brain >
      primary_median,
    na.rm = TRUE
  )

  direction_required <- ceiling(
    0.65 *
      length(
        brain
      )
  )

  p_value <- one_sided_wilcox(
    brain,
    primary
  )

  effect <- hedges_g_two_group(
    brain,
    primary
  )

  median_difference <- stats::median(
    brain,
    na.rm = TRUE
  ) -
    primary_median

  data.table(
    analysis =
      unique(
        score_data$analysis
      ),
    endpoint =
      endpoint,
    brain_patients =
      length(
        brain
      ),
    primary_patients =
      length(
        primary
      ),
    brain_above_primary_median =
      direction_count,
    direction_required =
      direction_required,
    brain_median =
      stats::median(
        brain,
        na.rm = TRUE
      ),
    primary_median =
      primary_median,
    median_difference =
      median_difference,
    one_sided_wilcox_p =
      p_value,
    hedges_g =
      effect,
    patient_count_pass =
      length(brain) %in%
        c(
          31L
        ) &&
        length(primary) %in%
          c(
            10L,
            11L
          ),
    direction_pass =
      median_difference >
        0 &&
        direction_count >=
          direction_required,
    p_pass =
      is.finite(
        p_value
      ) &&
        p_value <=
          0.05,
    effect_pass =
      is.finite(
        effect
      ) &&
        effect >=
          0.40
  )[
    ,
    endpoint_pass :=
      patient_count_pass &
        direction_pass &
        p_pass &
        effect_pass
  ]
}

main_raw_gate <- calculate_gate(
  main_scores,
  "Fixed_5_gene_module",
  "module_score"
)

main_adjusted_gate <- calculate_gate(
  main_scores,
  "Cell_cycle_adjusted_5_gene_module",
  "adjusted_module_score"
)

all_raw_gate <- calculate_gate(
  all_scores,
  "Fixed_5_gene_module",
  "module_score"
)

all_adjusted_gate <- calculate_gate(
  all_scores,
  "Cell_cycle_adjusted_5_gene_module",
  "adjusted_module_score"
)

gate_summary <- rbindlist(
  list(
    main_raw_gate,
    main_adjusted_gate,
    all_raw_gate,
    all_adjusted_gate
  ),
  fill = TRUE
)

strict_status <- if (
  isTRUE(
    main_raw_gate$endpoint_pass
  ) &&
    isTRUE(
      main_adjusted_gate$direction_pass
    ) &&
    isTRUE(
      main_adjusted_gate$p_pass
    )
) {
  "GSE223499_STRICT_VALIDATION_PASS"
} else if (
  isTRUE(
    main_raw_gate$endpoint_pass
  )
) {
  "GSE223499_PRIMARY_PASS_ADJUSTED_NOT_PASS"
} else {
  "GSE223499_VALIDATION_FAIL"
}

gate_summary[
  ,
  final_validation_status :=
    strict_status
]

print(
  gate_summary
)

cat(
  "\n正式状态：",
  strict_status,
  "\n",
  sep = ""
)

fwrite(
  gate_summary,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_validation_gate_summary.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 阶段9：批次重叠敏感性分析
# ============================================================

banner(
  "阶段9：KRAS/STK重叠批次敏感性分析"
)

overlap_wide <- expression_wide[
  primary_eligible == TRUE &
    batch %in%
      c(
        "KRAS",
        "STK"
      )
]

overlap_scores <- calculate_scores(
  overlap_wide,
  "Batch_overlap_KRAS_STK"
)

overlap_scores[
  ,
  brain_indicator :=
    as.integer(
      cohort ==
        "Brain_metastasis"
    )
]

batch_model <- stats::lm(
  module_score ~
    brain_indicator +
    factor(batch) +
    s_score +
    g2m_score,
  data =
    overlap_scores
)

batch_coef <- summary(
  batch_model
)$coefficients

brain_row <- batch_coef[
  "brain_indicator",
  ,
  drop = FALSE
]

batch_sensitivity <- data.table(
  analysis =
    "KRAS_STK_overlap_batch_adjusted",
  brain_patients =
    overlap_scores[
      cohort ==
        "Brain_metastasis",
      .N
    ],
  primary_patients =
    overlap_scores[
      cohort ==
        "Primary_tumor",
      .N
    ],
  coefficient =
    brain_row[
      1,
      "Estimate"
    ],
  standard_error =
    brain_row[
      1,
      "Std. Error"
    ],
  t_value =
    brain_row[
      1,
      "t value"
    ],
  two_sided_p =
    brain_row[
      1,
      "Pr(>|t|)"
    ]
)

batch_sensitivity[
  ,
  one_sided_positive_p :=
    one_sided_lm_p(
      coefficient,
      two_sided_p
    )
]

batch_sensitivity[
  ,
  positive_direction :=
    coefficient >
      0
]

print(
  batch_sensitivity
)

fwrite(
  batch_sensitivity,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_batch_overlap_sensitivity.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 阶段10：逐基因探索性审计
# ============================================================

banner(
  "阶段10：固定5基因逐基因审计"
)

gene_test_rows <- lapply(
  module_genes,
  function(current_gene) {
    brain <- main_scores[
      cohort ==
        "Brain_metastasis",
      get(
        current_gene
      )
    ]

    primary <- main_scores[
      cohort ==
        "Primary_tumor",
      get(
        current_gene
      )
    ]

    current_p <- one_sided_wilcox(
      brain,
      primary
    )

    data.table(
      gene =
        current_gene,
      brain_median =
        stats::median(
          brain,
          na.rm = TRUE
        ),
      primary_median =
        stats::median(
          primary,
          na.rm = TRUE
        ),
      median_difference =
        stats::median(
          brain,
          na.rm = TRUE
        ) -
          stats::median(
            primary,
            na.rm = TRUE
          ),
      one_sided_wilcox_p =
        current_p,
      hedges_g =
        hedges_g_two_group(
          brain,
          primary
        ),
      positive_direction =
        stats::median(
          brain,
          na.rm = TRUE
        ) >
          stats::median(
            primary,
            na.rm = TRUE
          )
    )
  }
)

gene_tests <- rbindlist(
  gene_test_rows
)

gene_tests[
  ,
  BH_FDR :=
    stats::p.adjust(
      one_sided_wilcox_p,
      method = "BH"
    )
]

setorder(
  gene_tests,
  one_sided_wilcox_p
)

print(
  gene_tests
)

fwrite(
  gene_tests,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_fixed5_gene_audit.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 阶段11：留一患者稳定性
# ============================================================

banner(
  "阶段11：主评分留一患者稳定性"
)

leave_one_rows <- lapply(
  main_scores$patient,
  function(omitted_patient) {
    current <- main_scores[
      patient !=
        omitted_patient
    ]

    brain <- current[
      cohort ==
        "Brain_metastasis",
      module_score
    ]

    primary <- current[
      cohort ==
        "Primary_tumor",
      module_score
    ]

    data.table(
      omitted_patient =
        omitted_patient,
      omitted_cohort =
        main_scores[
          patient ==
            omitted_patient,
          cohort
        ][1],
      brain_median =
        stats::median(
          brain,
          na.rm = TRUE
        ),
      primary_median =
        stats::median(
          primary,
          na.rm = TRUE
        ),
      direction_preserved =
        stats::median(
          brain,
          na.rm = TRUE
        ) >
          stats::median(
            primary,
            na.rm = TRUE
          ),
      one_sided_wilcox_p =
        one_sided_wilcox(
          brain,
          primary
        ),
      hedges_g =
        hedges_g_two_group(
          brain,
          primary
        )
    )
  }
)

leave_one_out <- rbindlist(
  leave_one_rows
)

fwrite(
  leave_one_out,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_leave_one_patient_out.tsv"
  ),
  sep = "\t"
)

cat(
  "留一后方向保持：",
  leave_one_out[
    direction_preserved == TRUE,
    .N
  ],
  "/",
  nrow(
    leave_one_out
  ),
  "。\n",
  sep = ""
)

# ============================================================
# 阶段12：患者级图形
# ============================================================

banner(
  "阶段12：生成患者级图形"
)

plot_data <- copy(
  main_scores
)

plot_data[
  ,
  cohort_label := factor(
    cohort,
    levels = c(
      "Primary_tumor",
      "Brain_metastasis"
    ),
    labels = c(
      "Primary LUAD",
      "Brain metastasis"
    )
  )
]

main_plot <- ggplot(
  plot_data,
  aes(
    x = cohort_label,
    y = module_score
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.2
  ) +
  labs(
    x = NULL,
    y = "Fixed five-gene module score",
    title = "GSE223499 patient-level validation",
    subtitle =
      "Author-defined tumor nuclei; patient as statistical unit"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "GSE223499_STEP07C_primary_module_score.png"
  ),
  plot = main_plot,
  width = 6.5,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "GSE223499_STEP07C_primary_module_score.pdf"
  ),
  plot = main_plot,
  width = 6.5,
  height = 5
)

adjusted_plot <- ggplot(
  plot_data,
  aes(
    x = cohort_label,
    y = adjusted_module_score
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.2
  ) +
  labs(
    x = NULL,
    y = "Cell-cycle-adjusted module score",
    title = "GSE223499 cell-cycle sensitivity"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "GSE223499_STEP07C_adjusted_module_score.png"
  ),
  plot = adjusted_plot,
  width = 6.5,
  height = 5,
  dpi = 300
)

# ============================================================
# 阶段13：完成锁
# ============================================================

banner(
  "STEP07-C全部完成"
)

completion_file <- file.path(
  LOCK_DIR,
  "STEP_07C_GSE223499_COMPLETE.txt"
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
      "validation_status=",
      strict_status
    ),
    "counts_files_valid=42/42",
    "full_expression_matrix_loaded=FALSE",
    "Seurat_object_created=FALSE",
    "primary_analysis=31_brain_vs_10_primary",
    "all_sample_sensitivity=31_brain_vs_11_primary",
    "primary_unit=resolved_patient",
    "fixed_module_genes=DHFR,DHODH,SHMT1,TYMS,UMPS",
    paste0(
      "main_one_sided_p=",
      main_raw_gate$one_sided_wilcox_p
    ),
    paste0(
      "main_hedges_g=",
      main_raw_gate$hedges_g
    ),
    paste0(
      "main_direction_count=",
      main_raw_gate$brain_above_primary_median,
      "/",
      main_raw_gate$brain_patients
    ),
    paste0(
      "adjusted_one_sided_p=",
      main_adjusted_gate$one_sided_wilcox_p
    ),
    paste0(
      "batch_overlap_coefficient=",
      batch_sensitivity$coefficient
    ),
    paste0(
      "batch_overlap_one_sided_p=",
      batch_sensitivity$one_sided_positive_p
    ),
    "outcome_driven_exclusion=none",
    "post_result_threshold_change=not_allowed",
    "barcode_suffix_fallback=unique_core_sequence_only_and_100_percent_match_required"
  ),
  completion_file
)

cat(
  "Status: ",
  strict_status,
  "\n",
  sep = ""
)

cat(
  "\nKey output files:\n"
)

cat(
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_validation_gate_summary.tsv"
  ),
  "\n",
  sep = ""
)

cat(
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_patient_level_scores.tsv"
  ),
  "\n",
  sep = ""
)

cat(
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_batch_overlap_sensitivity.tsv"
  ),
  "\n",
  sep = ""
)

cat(
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_fixed5_gene_audit.tsv"
  ),
  "\n",
  sep = ""
)

cat(
  completion_file,
  "\n",
  sep = ""
)
