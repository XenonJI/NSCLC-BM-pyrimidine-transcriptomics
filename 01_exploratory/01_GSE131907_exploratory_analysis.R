# GSE131907 exploratory analysis
# Low-memory extraction and aggregation of the fixed pyrimidine module and cell-cycle genes.
# This script is retained as the historical exploratory resource; it does not reconstruct
# the complete candidate-selection search space described in the manuscript.
# Set ROOT to the local project directory before running.

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  timeout = 7200
)

set.seed(20260806)

ROOT <- "D:/LUAD_LM_PYRIMIDINE"

RAW_DIR <- file.path(
  ROOT,
  "data/raw/GSE131907"
)

ANNOTATION_FILE <- file.path(
  RAW_DIR,
  "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
)

NORMALIZED_TEXT_FILE <- file.path(
  RAW_DIR,
  "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz"
)

TABLE_DIR <- file.path(
  ROOT,
  "results/tables"
)

FIGURE_DIR <- file.path(
  ROOT,
  "results/figures"
)

LOCK_DIR <- file.path(
  ROOT,
  "results/locks"
)

PROCESSED_DIR <- file.path(
  ROOT,
  "data/processed"
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  LOCK_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  PROCESSED_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

required_packages <- c(
  "data.table",
  "ggplot2",
  "digest"
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
    "缺少R包：",
    paste(
      missing_packages,
      collapse = ", "
    )
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

  x
}

safe_zscore <- function(x) {
  x <- as.numeric(x)

  current_sd <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(current_sd) ||
      current_sd == 0
  ) {
    return(
      rep(
        0,
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
  ) / current_sd
}

exact_rank_sum_greater <- function(
  brain_values,
  primary_values
) {
  brain_values <- brain_values[
    is.finite(brain_values)
  ]

  primary_values <- primary_values[
    is.finite(primary_values)
  ]

  n_brain <- length(
    brain_values
  )

  n_primary <- length(
    primary_values
  )

  if (
    n_brain < 2L ||
      n_primary < 2L
  ) {
    return(NA_real_)
  }

  all_values <- c(
    brain_values,
    primary_values
  )

  ranks <- rank(
    all_values,
    ties.method = "average"
  )

  observed <- sum(
    ranks[
      seq_len(n_brain)
    ]
  )

  combinations <- combn(
    length(all_values),
    n_brain
  )

  permutation_sums <- colSums(
    matrix(
      ranks[combinations],
      nrow = n_brain
    )
  )

  mean(
    permutation_sums >=
      observed - 1e-12
  )
}

hedges_g_independent <- function(
  brain_values,
  primary_values
) {
  brain_values <- brain_values[
    is.finite(brain_values)
  ]

  primary_values <- primary_values[
    is.finite(primary_values)
  ]

  n1 <- length(
    brain_values
  )

  n0 <- length(
    primary_values
  )

  if (
    n1 < 2L ||
      n0 < 2L
  ) {
    return(NA_real_)
  }

  pooled_sd <- sqrt(
    (
      (n1 - 1) *
        stats::var(brain_values) +
        (n0 - 1) *
        stats::var(primary_values)
    ) /
      (
        n1 +
          n0 -
          2
      )
  )

  if (
    !is.finite(pooled_sd) ||
      pooled_sd == 0
  ) {
    return(NA_real_)
  }

  cohen_d <- (
    mean(brain_values) -
      mean(primary_values)
  ) / pooled_sd

  correction <- 1 -
    3 / (
      4 *
        (
          n1 +
            n0
        ) -
        9
    )

  correction *
    cohen_d
}

calculate_endpoint <- function(
  endpoint_name,
  brain_values,
  primary_values
) {
  brain_values <- brain_values[
    is.finite(brain_values)
  ]

  primary_values <- primary_values[
    is.finite(primary_values)
  ]

  primary_median <- median(
    primary_values,
    na.rm = TRUE
  )

  data.table(
    endpoint = endpoint_name,
    brain_patients = length(
      brain_values
    ),
    primary_patients = length(
      primary_values
    ),
    brain_above_primary_median = sum(
      brain_values >
        primary_median
    ),
    brain_median = median(
      brain_values,
      na.rm = TRUE
    ),
    primary_median = primary_median,
    median_difference = median(
      brain_values,
      na.rm = TRUE
    ) -
      primary_median,
    exact_one_sided_rank_p =
      exact_rank_sum_greater(
        brain_values,
        primary_values
      ),
    hedges_g =
      hedges_g_independent(
        brain_values,
        primary_values
      )
  )
}

# ============================================================
# 阶段0：安全检查
# ============================================================

banner(
  "阶段0：16GB低内存安全检查"
)

dangerous_directory <- file.path(
  PROCESSED_DIR,
  "GSE131907_serialization_layers"
)

if (dir.exists(dangerous_directory)) {
  dangerous_files <- list.files(
    dangerous_directory,
    recursive = TRUE,
    full.names = TRUE
  )

  dangerous_size <- if (
    length(dangerous_files) > 0L
  ) {
    sum(
      file.info(
        dangerous_files
      )$size,
      na.rm = TRUE
    )
  } else {
    0
  }

  if (
    is.finite(dangerous_size) &&
      dangerous_size >
        100 *
          1024^2
  ) {
    stop(
      "检测到旧V5生成的大型解压目录：\n",
      dangerous_directory,
      "\n大小约：",
      round(
        dangerous_size /
          1024^3,
        2
      ),
      " GB。\n",
      "请关闭RStudio后删除整个目录，再运行本脚本。"
    )
  }
}

if (!file.exists(ANNOTATION_FILE)) {
  stop(
    "缺少注释文件：\n",
    ANNOTATION_FILE
  )
}

if (!file.exists(NORMALIZED_TEXT_FILE)) {
  stop(
    "缺少低内存分析所需的官方标准化文本矩阵：\n",
    NORMALIZED_TEXT_FILE,
    "\n请下载GEO中的：",
    "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
    "\n不要解压。"
  )
}

cat(
  "本脚本不会读取raw_UMI_matrix.rds.gz。\n"
)

cat(
  "本脚本不会在C盘或D盘生成完整解压矩阵。\n"
)

cat(
  "表达矩阵将通过gzfile逐行流式读取。\n"
)

# ============================================================
# 阶段1：冻结低内存方案修订
# ============================================================

banner(
  "阶段1：冻结低内存验证方案"
)

amendment_file <- file.path(
  LOCK_DIR,
  "STEP_06A_LOW_MEMORY_AMENDMENT_02.txt"
)

outcome_files <- c(
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_patient_level_scores.tsv"
  ),
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_validation_gate_summary.tsv"
  )
)

if (
  !file.exists(amendment_file) &&
    any(
      file.exists(
        outcome_files
      )
    )
) {
  stop(
    "已存在低内存表达结果，但不存在对应修订锁；",
    "为避免结果后修改，分析停止。"
  )
}

if (!file.exists(amendment_file)) {
  amendment_text <- c(
    paste0(
      "amended_at=",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    "reason=full RDS expands beyond 16GB hardware capacity",
    "amended_before_normalized_expression_outcome=TRUE",
    "data_source=official normalized log2TPM text matrix",
    "reading_strategy=stream gzip line by line without full decompression",
    "primary_unit=patient",
    "fixed_cohorts=10 brain metastasis patients and 11 primary LUAD patients",
    "fixed_malignant_definition=primary tS1/tS2/tS3 plus brain Malignant cells",
    "fixed_main_genes=CAD,CMPK1,CPS1,DHFR,DHFRP1,DHODH,MTOR,SHMT1,TYMS,UMPS",
    "fixed_core_genes=CAD,CMPK1,DHFR,DHODH,SHMT1,TYMS,UMPS",
    "patient_expression=mean official log2TPM across malignant cells per gene",
    "pathway_score=mean gene-wise z score across fixed patients",
    "cell_cycle_sensitivity=patient-level residual adjustment for S and G2M scores",
    "main_direction_rule=at least 8 of 10 brain patients above primary median",
    "main_p_rule=exact one-sided rank-sum p <= 0.05",
    "main_effect_rule=Hedges g >= 0.50",
    "support_rule=core and cell-cycle-adjusted endpoints must also pass",
    "all_fixed_patients_retained=TRUE",
    "low_cell_count_samples=retained and transparently flagged",
    "outcome_driven_patient_exclusion=not allowed",
    "gene_or_threshold_change_after_result=not allowed"
  )

  writeLines(
    amendment_text,
    amendment_file
  )
}

cat(
  "低内存修订锁：\n",
  amendment_file,
  "\n",
  sep = ""
)

# ============================================================
# 阶段2：读取注释并固定目标细胞
# ============================================================

banner(
  "阶段2：读取作者注释并固定目标恶性细胞"
)

annotation <- fread(
  ANNOTATION_FILE,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  showProgress = FALSE
)

required_annotation_columns <- c(
  "Index",
  "Sample",
  "Sample_Origin",
  "Cell_type",
  "Cell_subtype"
)

missing_annotation_columns <-
  setdiff(
    required_annotation_columns,
    names(annotation)
  )

if (
  length(
    missing_annotation_columns
  ) >
    0L
) {
  stop(
    "注释文件缺少字段：",
    paste(
      missing_annotation_columns,
      collapse = ", "
    )
  )
}

annotation[
  ,
  cell_id :=
    strip_quotes(
      Index
    )
]

annotation[
  ,
  sample :=
    trimws(
      as.character(
        Sample
      )
    )
]

annotation[
  ,
  origin :=
    trimws(
      as.character(
        Sample_Origin
      )
    )
]

annotation[
  ,
  broad_type :=
    trimws(
      as.character(
        Cell_type
      )
    )
]

annotation[
  ,
  subtype :=
    trimws(
      as.character(
        Cell_subtype
      )
    )
]

annotation[
  ,
  cohort := fifelse(
    origin == "mBrain",
    "Brain_metastasis",
    fifelse(
      origin == "tLung",
      "Primary_LUAD",
      "Other"
    )
  )
]

annotation[
  ,
  malignant := (
    broad_type ==
      "Epithelial cells"
  ) &
    (
      (
        cohort ==
          "Brain_metastasis" &
          subtype ==
            "Malignant cells"
      ) |
        (
          cohort ==
            "Primary_LUAD" &
            subtype %in%
              c(
                "tS1",
                "tS2",
                "tS3"
              )
        )
    )
]

target_annotation <- annotation[
  malignant &
    cohort %in%
      c(
        "Brain_metastasis",
        "Primary_LUAD"
      ),
  .(
    cell_id,
    patient = sample,
    cohort,
    subtype
  )
]

if (
  anyNA(
    target_annotation$cell_id
  ) ||
    anyDuplicated(
      target_annotation$cell_id
    )
) {
  stop(
    "目标细胞ID存在缺失或重复。"
  )
}

patient_cell_counts <-
  target_annotation[
    ,
    .N,
    by = .(
      patient,
      cohort
    )
  ]

setorder(
  patient_cell_counts,
  cohort,
  patient
)

patient_cell_counts[
  ,
  low_cell_flag :=
    N < 50L
]

print(
  patient_cell_counts
)

fwrite(
  patient_cell_counts,
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_target_cell_counts.tsv"
  ),
  sep = "\t"
)

brain_patients <- sort(
  unique(
    target_annotation[
      cohort ==
        "Brain_metastasis",
      patient
    ]
  )
)

primary_patients <- sort(
  unique(
    target_annotation[
      cohort ==
        "Primary_LUAD",
      patient
    ]
  )
)

if (
  length(brain_patients) != 10L ||
    length(primary_patients) != 11L
) {
  stop(
    "固定患者数异常：脑转移",
    length(brain_patients),
    "，原发",
    length(primary_patients),
    "。"
  )
}

all_patients <- c(
  primary_patients,
  brain_patients
)

patient_cohort <- c(
  setNames(
    rep(
      "Primary_LUAD",
      length(
        primary_patients
      )
    ),
    primary_patients
  ),
  setNames(
    rep(
      "Brain_metastasis",
      length(
        brain_patients
      )
    ),
    brain_patients
  )
)

cat(
  "固定目标细胞：",
  format(
    nrow(
      target_annotation
    ),
    big.mark = ","
  ),
  "个；患者：",
  length(
    all_patients
  ),
  "名。\n",
  sep = ""
)

# ============================================================
# 阶段3：固定目标基因
# ============================================================

banner(
  "阶段3：固定嘧啶与细胞周期基因"
)

main_genes <- c(
  "CAD",
  "CMPK1",
  "CPS1",
  "DHFR",
  "DHFRP1",
  "DHODH",
  "MTOR",
  "SHMT1",
  "TYMS",
  "UMPS"
)

core_genes <- c(
  "CAD",
  "CMPK1",
  "DHFR",
  "DHODH",
  "SHMT1",
  "TYMS",
  "UMPS"
)

load_fixed_cell_cycle_genes <- function() {
  # 首选：通过data()读取Seurat附带的数据集对象。
  # 数据集对象不保证存在于包namespace中，因此不使用
  # get(..., envir = asNamespace("Seurat"))。
  if (
    requireNamespace(
      "Seurat",
      quietly = TRUE
    )
  ) {
    cc_environment <- new.env(
      parent = emptyenv()
    )

    load_success <- tryCatch(
      {
        suppressWarnings(
          utils::data(
            list =
              "cc.genes.updated.2019",
            package = "Seurat",
            envir = cc_environment
          )
        )

        exists(
          "cc.genes.updated.2019",
          envir = cc_environment,
          inherits = FALSE
        )
      },
      error = function(e) {
        FALSE
      }
    )

    if (isTRUE(load_success)) {
      cc_object <- get(
        "cc.genes.updated.2019",
        envir = cc_environment,
        inherits = FALSE
      )

      if (
        is.list(cc_object) &&
          all(
            c(
              "s.genes",
              "g2m.genes"
            ) %in%
              names(cc_object)
          )
      ) {
        return(
          list(
            s.genes = unique(
              toupper(
                cc_object$s.genes
              )
            ),
            g2m.genes = unique(
              toupper(
                cc_object$g2m.genes
              )
            ),
            source =
              "Seurat::cc.genes.updated.2019 loaded via utils::data"
          )
        )
      }
    }
  }

  # 回退：冻结的Seurat 2019人类细胞周期基因列表。
  # 这样不依赖本地Seurat版本，也不会联网下载。
  fixed_s_genes <- c(
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

  fixed_g2m_genes <- c(
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

  list(
    s.genes = fixed_s_genes,
    g2m.genes = fixed_g2m_genes,
    source =
      "script-embedded frozen Seurat 2019 human cell-cycle list"
  )
}

cc_genes <- load_fixed_cell_cycle_genes()

s_genes <- cc_genes$s.genes
g2m_genes <- cc_genes$g2m.genes

cat(
  "细胞周期基因来源：",
  cc_genes$source,
  "\n",
  sep = ""
)

cat(
  "固定S期基因：",
  length(s_genes),
  "个；固定G2/M期基因：",
  length(g2m_genes),
  "个。\n",
  sep = ""
)

writeLines(
  c(
    paste0(
      "source=",
      cc_genes$source
    ),
    paste0(
      "S_genes=",
      paste(
        s_genes,
        collapse = ","
      )
    ),
    paste0(
      "G2M_genes=",
      paste(
        g2m_genes,
        collapse = ","
      )
    )
  ),
  file.path(
    LOCK_DIR,
    "STEP_06A_LOWMEM_CELL_CYCLE_GENE_LOCK.txt"
  )
)

genes_to_extract <- unique(
  c(
    main_genes,
    s_genes,
    g2m_genes
  )
)

cat(
  "计划提取基因数：",
  length(
    genes_to_extract
  ),
  "\n",
  sep = ""
)

# ============================================================
# 阶段4：流式读取标准化矩阵
# ============================================================

banner(
  "阶段4：流式扫描2.9GB压缩矩阵"
)

checkpoint_file <- file.path(
  PROCESSED_DIR,
  "GSE131907_lowmem_patient_gene_means.rds"
)

if (file.exists(checkpoint_file)) {
  cat(
    "检测到患者×基因检查点，直接复用：\n",
    checkpoint_file,
    "\n",
    sep = ""
  )

  checkpoint <- readRDS(
    checkpoint_file
  )

  patient_gene_means <-
    checkpoint$patient_gene_means

  detected_genes <-
    checkpoint$detected_genes
} else {
  connection <- gzfile(
    NORMALIZED_TEXT_FILE,
    open = "rt"
  )

  # 连接在流式读取完成后显式关闭。
  # 不在脚本顶层使用 on.exit()，以避免顶层环境错误。

  header_line <- readLines(
    connection,
    n = 1L,
    warn = FALSE
  )

  if (
    length(
      header_line
    ) != 1L
  ) {
    stop(
      "无法读取标准化矩阵表头。"
    )
  }

  header_fields <- strsplit(
    header_line,
    split = "\t",
    fixed = TRUE
  )[[1]]

  header_fields <- strip_quotes(
    header_fields
  )

  if (
    length(
      header_fields
    ) <
      100000L
  ) {
    stop(
      "矩阵表头列数异常：",
      length(
        header_fields
      )
    )
  }

  selected_positions <- match(
    target_annotation$cell_id,
    header_fields
  )

  alignment_fraction <- mean(
    !is.na(
      selected_positions
    )
  )

  cat(
    "目标细胞与矩阵表头匹配率：",
    round(
      100 *
        alignment_fraction,
      3
    ),
    "%\n",
    sep = ""
  )

  if (
    alignment_fraction <
      0.98
  ) {
    stop(
      "目标细胞与标准化矩阵匹配率低于98%。"
    )
  }

  target_annotation <- target_annotation[
    !is.na(
      selected_positions
    )
  ]

  selected_positions <- selected_positions[
    !is.na(
      selected_positions
    )
  ]

  selected_patient <- factor(
    target_annotation$patient,
    levels = all_patients
  )

  patient_gene_list <- list()

  detected_genes <- character()

  row_number <- 0L

  repeat {
    current_line <- readLines(
      connection,
      n = 1L,
      warn = FALSE
    )

    if (
      length(
        current_line
      ) == 0L
    ) {
      break
    }

    row_number <- row_number + 1L

    first_tab <- regexpr(
      "\t",
      current_line,
      fixed = TRUE
    )[1]

    if (
      is.na(
        first_tab
      ) ||
        first_tab < 2L
    ) {
      next
    }

    current_gene <- substr(
      current_line,
      1L,
      first_tab - 1L
    )

    current_gene <- toupper(
      strip_quotes(
        current_gene
      )
    )

    if (
      current_gene %in%
        genes_to_extract &&
        !current_gene %in%
          detected_genes
    ) {
      fields <- strsplit(
        current_line,
        split = "\t",
        fixed = TRUE
      )[[1]]

      if (
        max(
          selected_positions
        ) >
          length(
            fields
          )
      ) {
        stop(
          "基因 ",
          current_gene,
          " 的字段数少于表头列数。"
        )
      }

      selected_values <- suppressWarnings(
        as.numeric(
          fields[
            selected_positions
          ]
        )
      )

      if (
        anyNA(
          selected_values
        )
      ) {
        stop(
          "基因 ",
          current_gene,
          " 的目标表达值包含无法解析的字符。"
        )
      }

      sums <- rowsum(
        selected_values,
        group = selected_patient,
        reorder = TRUE
      )

      sum_by_patient <- setNames(
        as.numeric(
          sums[
            ,
            1
          ]
        ),
        rownames(
          sums
        )
      )

      count_by_patient <- table(
        selected_patient
      )

      patient_means <-
        sum_by_patient[
          all_patients
        ] /
          as.numeric(
            count_by_patient[
              all_patients
            ]
          )

      names(
        patient_means
      ) <- all_patients

      if (
        anyNA(patient_means) ||
          any(
            !is.finite(
              patient_means
            )
          )
      ) {
        stop(
          "基因 ",
          current_gene,
          " 的患者级表达汇总出现缺失或非有限值。"
        )
      }

      patient_gene_list[[current_gene]] <- patient_means

      detected_genes <- c(
        detected_genes,
        current_gene
      )

      cat(
        "已提取：",
        current_gene,
        "（",
        length(
          detected_genes
        ),
        "/",
        length(
          genes_to_extract
        ),
        "）\n",
        sep = ""
      )

      rm(
        fields,
        selected_values,
        sums,
        patient_means
      )

      if (
        length(
          detected_genes
        ) %%
          10L ==
          0L
      ) {
        gc()
      }
    }

    if (
      row_number %%
        2000L ==
        0L
    ) {
      cat(
        "已扫描基因行：",
        format(
          row_number,
          big.mark = ","
        ),
        "；已找到目标基因：",
        length(
          detected_genes
        ),
        "\n",
        sep = ""
      )
    }
  }

  close(
    connection
  )

  if (
    length(
      patient_gene_list
    ) == 0L
  ) {
    stop(
      "标准化矩阵中没有找到任何固定目标基因。"
    )
  }

  patient_gene_means <- do.call(
    cbind,
    patient_gene_list
  )

  rownames(
    patient_gene_means
  ) <- all_patients

  saveRDS(
    list(
      patient_gene_means =
        patient_gene_means,
      detected_genes =
        detected_genes,
      generated_at =
        Sys.time()
    ),
    checkpoint_file,
    compress = TRUE
  )

  cat(
    "低内存检查点已保存：\n",
    checkpoint_file,
    "\n",
    sep = ""
  )

  rm(
    patient_gene_list,
    header_line,
    header_fields,
    selected_positions
  )

  gc()
}

# ============================================================
# 阶段5：患者层面评分
# ============================================================

banner(
  "阶段5：患者层面通路和细胞周期评分"
)

detected_main <- intersect(
  main_genes,
  colnames(
    patient_gene_means
  )
)

detected_core <- intersect(
  core_genes,
  colnames(
    patient_gene_means
  )
)

detected_s <- intersect(
  s_genes,
  colnames(
    patient_gene_means
  )
)

detected_g2m <- intersect(
  g2m_genes,
  colnames(
    patient_gene_means
  )
)

gene_detection <- data.table(
  gene = genes_to_extract,
  category = fifelse(
    genes_to_extract %in%
      main_genes,
    "pyrimidine_main",
    fifelse(
      genes_to_extract %in%
        s_genes,
      "S_phase",
      "G2M_phase"
    )
  ),
  detected =
    genes_to_extract %in%
      colnames(
        patient_gene_means
      )
)

fwrite(
  gene_detection,
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_gene_detection.tsv"
  ),
  sep = "\t"
)

cat(
  "主基因检测：",
  length(
    detected_main
  ),
  "/10；核心基因：",
  length(
    detected_core
  ),
  "/7；S期：",
  length(
    detected_s
  ),
  "；G2M期：",
  length(
    detected_g2m
  ),
  "。\n",
  sep = ""
)

if (
  length(
    detected_main
  ) <
    9L
) {
  stop(
    "主基因检测少于9/10。"
  )
}

if (
  length(
    detected_core
  ) <
    7L
) {
  stop(
    "核心基因检测少于7/7。"
  )
}

if (
  length(
    detected_s
  ) <
    20L ||
    length(
      detected_g2m
    ) <
      20L
) {
  stop(
    "细胞周期基因检测不足，不能进行校正。"
  )
}

zscore_matrix <- apply(
  patient_gene_means,
  2,
  safe_zscore
)

if (
  is.null(
    dim(
      zscore_matrix
    )
  )
) {
  stop(
    "患者×基因标准化矩阵维度异常。"
  )
}

rownames(
  zscore_matrix
) <- rownames(
  patient_gene_means
)

main_score <- rowMeans(
  zscore_matrix[
    ,
    detected_main,
    drop = FALSE
  ]
)

core_score <- rowMeans(
  zscore_matrix[
    ,
    detected_core,
    drop = FALSE
  ]
)

s_score <- rowMeans(
  zscore_matrix[
    ,
    detected_s,
    drop = FALSE
  ]
)

g2m_score <- rowMeans(
  zscore_matrix[
    ,
    detected_g2m,
    drop = FALSE
  ]
)

adjustment_model <- lm(
  main_score ~
    s_score +
    g2m_score
)

adjusted_main_score <- as.numeric(
  residuals(
    adjustment_model
  )
)

names(
  adjusted_main_score
) <- names(
  main_score
)

patient_scores <- data.table(
  patient = all_patients,
  cohort = unname(
    patient_cohort[
      all_patients
    ]
  ),
  malignant_cells = patient_cell_counts$N[
    match(
      all_patients,
      patient_cell_counts$patient
    )
  ],
  low_cell_flag = patient_cell_counts$low_cell_flag[
    match(
      all_patients,
      patient_cell_counts$patient
    )
  ],
  main_score = as.numeric(
    main_score[
      all_patients
    ]
  ),
  core_score = as.numeric(
    core_score[
      all_patients
    ]
  ),
  s_score = as.numeric(
    s_score[
      all_patients
    ]
  ),
  g2m_score = as.numeric(
    g2m_score[
      all_patients
    ]
  ),
  adjusted_main_score = as.numeric(
    adjusted_main_score[
      all_patients
    ]
  )
)

setorder(
  patient_scores,
  cohort,
  patient
)

print(
  patient_scores
)

fwrite(
  patient_scores,
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_patient_level_scores.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 阶段6：验证闸门
# ============================================================

banner(
  "阶段6：独立验证阳性闸门"
)

brain <- patient_scores[
  cohort ==
    "Brain_metastasis"
]

primary <- patient_scores[
  cohort ==
    "Primary_LUAD"
]

main_endpoint <- calculate_endpoint(
  "Main_fixed_pyrimidine_score",
  brain$main_score,
  primary$main_score
)

core_endpoint <- calculate_endpoint(
  "Core_7_gene_sensitivity",
  brain$core_score,
  primary$core_score
)

adjusted_endpoint <- calculate_endpoint(
  "Cell_cycle_adjusted_main_score",
  brain$adjusted_main_score,
  primary$adjusted_main_score
)

gate <- rbindlist(
  list(
    main_endpoint,
    core_endpoint,
    adjusted_endpoint
  )
)

gate[
  ,
  patient_count_pass :=
    length(brain_patients) == 10L &
      length(primary_patients) == 11L
]

gate[
  ,
  direction_pass :=
    brain_above_primary_median >=
      8L
]

gate[
  ,
  p_pass :=
    is.finite(
      exact_one_sided_rank_p
    ) &
      exact_one_sided_rank_p <=
        0.05
]

gate[
  ,
  effect_pass :=
    is.finite(
      hedges_g
    ) &
      hedges_g >=
        0.50
]

gate[
  ,
  endpoint_pass :=
    patient_count_pass &
      direction_pass &
      p_pass &
      effect_pass
]

main_pass <- gate[
  endpoint ==
    "Main_fixed_pyrimidine_score",
  endpoint_pass
]

core_pass <- gate[
  endpoint ==
    "Core_7_gene_sensitivity",
  endpoint_pass
]

adjusted_pass <- gate[
  endpoint ==
    "Cell_cycle_adjusted_main_score",
  endpoint_pass
]

validation_status <- if (
  isTRUE(
    main_pass
  ) &&
    isTRUE(
      core_pass
    ) &&
    isTRUE(
      adjusted_pass
    )
) {
  "GSE131907_LOWMEM_VALIDATION_PASS"
} else {
  "GSE131907_LOWMEM_VALIDATION_FAIL"
}

gate[
  ,
  validation_status :=
    validation_status
]

print(
  gate
)

cat(
  "\n最终状态：",
  validation_status,
  "\n",
  sep = ""
)

fwrite(
  gate,
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_validation_gate_summary.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 阶段7：基因层面方向审计
# ============================================================

banner(
  "阶段7：固定主基因逐基因方向审计"
)

gene_direction_rows <- lapply(
  detected_main,
  function(gene) {
    values <- patient_gene_means[
      ,
      gene
    ]

    data.table(
      gene = gene,
      brain_median = median(
        values[
          brain_patients
        ],
        na.rm = TRUE
      ),
      primary_median = median(
        values[
          primary_patients
        ],
        na.rm = TRUE
      ),
      brain_minus_primary =
        median(
          values[
            brain_patients
          ],
          na.rm = TRUE
        ) -
        median(
          values[
            primary_patients
          ],
          na.rm = TRUE
        )
    )
  }
)

gene_direction <- rbindlist(
  gene_direction_rows
)

gene_direction[
  ,
  positive_direction :=
    brain_minus_primary >
      0
]

print(
  gene_direction
)

fwrite(
  gene_direction,
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_gene_direction_audit.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 阶段8：留一患者稳定性
# ============================================================

banner(
  "阶段8：留一患者稳定性检查"
)

leave_one_rows <- list()

for (
  omitted_patient in
    patient_scores$patient
) {
  current <- patient_scores[
    patient !=
      omitted_patient
  ]

  current_brain <- current[
    cohort ==
      "Brain_metastasis",
    main_score
  ]

  current_primary <- current[
    cohort ==
      "Primary_LUAD",
    main_score
  ]

  approximate_p <- suppressWarnings(
    wilcox.test(
      current_brain,
      current_primary,
      alternative = "greater",
      exact = FALSE,
      correct = FALSE
    )$p.value
  )

  leave_one_rows[[omitted_patient]] <- data.table(
    omitted_patient =
      omitted_patient,
    omitted_cohort =
      patient_scores[
        patient ==
          omitted_patient,
        cohort
      ][1],
    brain_median =
      median(
        current_brain
      ),
    primary_median =
      median(
        current_primary
      ),
    direction_preserved =
      median(
        current_brain
      ) >
        median(
          current_primary
        ),
    approximate_one_sided_p =
      approximate_p,
    hedges_g =
      hedges_g_independent(
        current_brain,
        current_primary
      )
  )
}

leave_one <- rbindlist(
  leave_one_rows
)

print(
  leave_one
)

fwrite(
  leave_one,
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_leave_one_patient_out.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 阶段9：绘图
# ============================================================

banner(
  "阶段9：生成患者级图形"
)

patient_scores[
  ,
  cohort := factor(
    cohort,
    levels = c(
      "Primary_LUAD",
      "Brain_metastasis"
    )
  )
]

p_main <- ggplot(
  patient_scores,
  aes(
    x = cohort,
    y = main_score
  )
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.8
  ) +
  geom_text(
    aes(
      label = patient
    ),
    nudge_x = 0.18,
    check_overlap = TRUE,
    size = 2.7
  ) +
  theme_classic(
    base_size = 13
  ) +
  labs(
    title =
      "GSE131907 low-memory independent validation",
    subtitle =
      "Fixed de novo pyrimidine biosynthesis score",
    x = NULL,
    y =
      "Patient-level gene-wise z-score mean"
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE131907_lowmem_primary_vs_brain_main.png"
  ),
  p_main,
  width = 8,
  height = 5.5,
  dpi = 300
)

p_adjusted <- ggplot(
  patient_scores,
  aes(
    x = cohort,
    y = adjusted_main_score
  )
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.8
  ) +
  geom_text(
    aes(
      label = patient
    ),
    nudge_x = 0.18,
    check_overlap = TRUE,
    size = 2.7
  ) +
  theme_classic(
    base_size = 13
  ) +
  labs(
    title =
      "Cell-cycle-adjusted sensitivity analysis",
    x = NULL,
    y =
      "Patient-level residual pyrimidine score"
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE131907_lowmem_cell_cycle_adjusted.png"
  ),
  p_adjusted,
  width = 8,
  height = 5.5,
  dpi = 300
)

# ============================================================
# 阶段10：完成锁
# ============================================================

completion_file <- file.path(
  LOCK_DIR,
  "STEP_06A_LOWMEM_COMPLETE.txt"
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
      validation_status
    ),
    "method=streamed official normalized log2TPM text",
    "full_matrix_loaded=FALSE",
    "full_matrix_decompressed_to_disk=FALSE",
    paste0(
      "main_genes_detected=",
      length(
        detected_main
      ),
      "/10"
    ),
    paste0(
      "core_genes_detected=",
      length(
        detected_core
      ),
      "/7"
    ),
    paste0(
      "brain_patients=",
      nrow(
        brain
      )
    ),
    paste0(
      "primary_patients=",
      nrow(
        primary
      )
    ),
    "primary_unit=patient",
    "outcome_driven_patient_exclusion=none",
    "threshold_change_after_result=not_allowed"
  ),
  completion_file
)

banner(
  "STEP06-A低内存版全部完成"
)

cat(
  "患者级结果：\n",
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_patient_level_scores.tsv"
  ),
  "\n",
  sep = ""
)

cat(
  "验证闸门：\n",
  file.path(
    TABLE_DIR,
    "GSE131907_lowmem_validation_gate_summary.tsv"
  ),
  "\n",
  sep = ""
)

cat(
  "完成锁：\n",
  completion_file,
  "\n",
  sep = ""
)
