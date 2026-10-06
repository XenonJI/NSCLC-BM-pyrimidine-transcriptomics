# GSE223499 fixed robustness analyses
# Generates the historical robustness outputs used to assemble the primary analysis table.
# The final patient mapping and inferential corrections are applied in downstream scripts.
# The fixed gene set and prespecified sample thresholds are unchanged.

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

ROOT <- "D:/LUAD_LM_PYRIMIDINE"

INPUT_FILE <- file.path(
  ROOT,
  "results/tables/GSE223499",
  "GSE223499_STEP07C_patient_gene_expression_long.tsv"
)

STEP07_GATE_FILE <- file.path(
  ROOT,
  "results/tables/GSE223499",
  "GSE223499_STEP07C_validation_gate_summary.tsv"
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

for (current_dir in c(
  TABLE_DIR,
  FIGURE_DIR,
  LOCK_DIR
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

  degrees_freedom <- n1 +
    n0 -
    2

  correction <- 1 -
    3 /
      (
        4 *
          degrees_freedom -
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

one_sided_p_from_lm <- function(
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
    two_sided_p /
      2
  } else {
    1 -
      two_sided_p /
        2
  }
}

extract_lm_term <- function(
  model,
  term_name
) {
  coefficient_table <- summary(
    model
  )$coefficients

  if (
    !term_name %in%
      rownames(
        coefficient_table
      )
  ) {
    return(
      data.table(
        coefficient =
          NA_real_,
        standard_error =
          NA_real_,
        t_value =
          NA_real_,
        two_sided_p =
          NA_real_,
        one_sided_positive_p =
          NA_real_,
        positive_direction =
          NA
      )
    )
  }

  estimate <- coefficient_table[
    term_name,
    "Estimate"
  ]

  two_sided_p <- coefficient_table[
    term_name,
    "Pr(>|t|)"
  ]

  data.table(
    coefficient =
      as.numeric(
        estimate
      ),
    standard_error =
      as.numeric(
        coefficient_table[
          term_name,
          "Std. Error"
        ]
      ),
    t_value =
      as.numeric(
        coefficient_table[
          term_name,
          "t value"
        ]
      ),
    two_sided_p =
      as.numeric(
        two_sided_p
      ),
    one_sided_positive_p =
      one_sided_p_from_lm(
        estimate =
          estimate,
        two_sided_p =
          two_sided_p
      ),
    positive_direction =
      estimate >
        0
  )
}

module_genes <- c(
  "DHFR",
  "DHODH",
  "SHMT1",
  "TYMS",
  "UMPS"
)

module_genes_without_TYMS <- setdiff(
  module_genes,
  "TYMS"
)

s_genes_original <- c(
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

g2m_genes_original <- c(
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

s_genes_nonoverlap <- setdiff(
  s_genes_original,
  module_genes
)

g2m_genes_nonoverlap <- setdiff(
  g2m_genes_original,
  module_genes
)

fixed_cell_thresholds <- c(
  50L,
  200L,
  500L,
  1000L
)

# ============================================================
# 阶段0：输入审计及分析锁
# ============================================================

banner(
  "阶段0：输入审计及STEP08分析锁"
)

if (!file.exists(INPUT_FILE)) {
  stop(
    "找不到STEP07患者级基因表达长表：\n",
    INPUT_FILE
  )
}

if (!file.exists(STEP07_GATE_FILE)) {
  stop(
    "找不到STEP07验证闸门文件：\n",
    STEP07_GATE_FILE
  )
}

step07_gate <- fread(
  STEP07_GATE_FILE,
  data.table = TRUE,
  showProgress = FALSE
)

expected_primary <- step07_gate[
  analysis ==
    "Primary_31_vs_10" &
    endpoint ==
      "Fixed_5_gene_module"
]

if (
  nrow(expected_primary) != 1L ||
    !isTRUE(
      expected_primary$endpoint_pass
    )
) {
  stop(
    "STEP07主要5基因模块并非唯一阳性记录，停止STEP08。"
  )
}

analysis_lock_file <- file.path(
  LOCK_DIR,
  "STEP_08_GSE223499_ROBUSTNESS_ANALYSIS_LOCK.txt"
)

if (!file.exists(analysis_lock_file)) {
  writeLines(
    c(
      paste0(
        "created_at=",
        format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        )
      ),
      "STEP07_results_are_immutable=TRUE",
      "fixed_primary_module=DHFR,DHODH,SHMT1,TYMS,UMPS",
      "primary_STEP08_question=does_the_five_gene_module_remain_positive_after_cell_cycle_adjustment_without_gene_overlap",
      "nonoverlap_rule=remove_all_five_module_genes_from_S_and_G2M_sets",
      "primary_adjusted_model=module_score_brain_indicator_S_nonoverlap_G2M_nonoverlap",
      "primary_adjusted_success=brain_coefficient_positive_and_one_sided_p_at_most_0.05",
      "residualized_rank_test=supportive_only",
      "KRAS17_exclusion=fixed_sensitivity_analysis",
      "cell_thresholds=50,200,500,1000",
      "batch_analyses=KRAS_only,STK_only,KRAS_STK_overlap,full_batch_adjusted",
      "four_gene_without_TYMS=post_hoc_sensitivity_only_not_primary_validation",
      "no_new_gene_selection=TRUE",
      "no_outcome_driven_threshold_change=TRUE",
      "stop_after_STEP08=TRUE"
    ),
    analysis_lock_file
  )
}

cat(
  "STEP08分析锁：\n",
  analysis_lock_file,
  "\n",
  sep = ""
)

# ============================================================
# 阶段1：读取患者级基因表达长表
# ============================================================

banner(
  "阶段1：读取患者级基因表达"
)

expression_long <- fread(
  INPUT_FILE,
  data.table = TRUE,
  showProgress = TRUE
)

required_columns <- c(
  "patient",
  "cohort",
  "batch",
  "gene",
  "log2_CPM_plus1",
  "author_tumor_nuclei",
  "primary_eligible",
  "samples",
  "selected_library_size"
)

missing_columns <- setdiff(
  required_columns,
  names(
    expression_long
  )
)

if (length(missing_columns) > 0L) {
  stop(
    "患者级表达长表缺少字段：",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

expression_long[
  ,
  `:=`(
    patient =
      as.character(
        patient
      ),
    cohort =
      as.character(
        cohort
      ),
    batch =
      as.character(
        batch
      ),
    gene =
      toupper(
        as.character(
          gene
        )
      ),
    samples =
      as.character(
        samples
      ),
    log2_CPM_plus1 =
      as.numeric(
        log2_CPM_plus1
      ),
    author_tumor_nuclei =
      as.integer(
        author_tumor_nuclei
      ),
    selected_library_size =
      as.numeric(
        selected_library_size
      )
  )
]

expression_long <- expression_long[
  cohort %in%
    c(
      "Brain_metastasis",
      "Primary_tumor"
    )
]

if (
  uniqueN(
    expression_long$patient
  ) != 42L
) {
  stop(
    "患者级表达长表不是42个比较患者。"
  )
}

patient_info <- unique(
  expression_long[
    ,
    .(
      patient,
      samples,
      cohort,
      batch,
      author_tumor_nuclei,
      primary_eligible,
      selected_library_size
    )
  ]
)

if (
  patient_info[
    cohort ==
      "Brain_metastasis",
    .N
  ] != 31L ||
    patient_info[
      cohort ==
        "Primary_tumor",
      .N
    ] != 11L
) {
  stop(
    "患者数量不是31例脑转移和11例原发。"
  )
}

expression_wide <- dcast(
  expression_long,
  patient +
    samples +
    cohort +
    batch +
    author_tumor_nuclei +
    primary_eligible +
    selected_library_size ~
    gene,
  value.var =
    "log2_CPM_plus1"
)

missing_module_genes <- setdiff(
  module_genes,
  names(
    expression_wide
  )
)

if (length(missing_module_genes) > 0L) {
  stop(
    "缺少固定模块基因：",
    paste(
      missing_module_genes,
      collapse = ", "
    )
  )
}

# ============================================================
# 通用评分和统计函数
# ============================================================

calculate_patient_scores <- function(
  wide_data,
  module_gene_set,
  analysis_name
) {
  result <- copy(
    wide_data
  )

  missing_current_module <- setdiff(
    module_gene_set,
    names(
      result
    )
  )

  if (length(missing_current_module) > 0L) {
    stop(
      analysis_name,
      "缺少模块基因：",
      paste(
        missing_current_module,
        collapse = ", "
      )
    )
  }

  module_matrix <- as.matrix(
    result[
      ,
      ..module_gene_set
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
          module_gene_set
        )
    )
  }

  colnames(
    module_z
  ) <- module_gene_set

  result[
    ,
    module_score :=
      rowMeans(
        module_z,
        na.rm = TRUE
      )
  ]

  available_s <- intersect(
    s_genes_nonoverlap,
    names(
      result
    )
  )

  available_g2m <- intersect(
    g2m_genes_nonoverlap,
    names(
      result
    )
  )

  variable_s <- available_s[
    vapply(
      result[
        ,
        ..available_s
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

  variable_g2m <- available_g2m[
    vapply(
      result[
        ,
        ..available_g2m
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
      "无重叠细胞周期基因数量不足。"
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
    s_score_nonoverlap :=
      rowMeans(
        s_z,
        na.rm = TRUE
      )
  ]

  result[
    ,
    g2m_score_nonoverlap :=
      rowMeans(
        g2m_z,
        na.rm = TRUE
      )
  ]

  result[
    ,
    brain_indicator :=
      as.integer(
        cohort ==
          "Brain_metastasis"
      )
  ]

  residual_model <- stats::lm(
    module_score ~
      s_score_nonoverlap +
      g2m_score_nonoverlap,
    data =
      result
  )

  result[
    ,
    residualized_module_score :=
      stats::residuals(
        residual_model
      )
  ]

  result[
    ,
    analysis :=
      analysis_name
  ]

  attr(
    result,
    "module_genes"
  ) <- module_gene_set

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

calculate_unadjusted_gate <- function(
  score_data,
  endpoint_name
) {
  brain <- score_data[
    cohort ==
      "Brain_metastasis",
    module_score
  ]

  primary <- score_data[
    cohort ==
      "Primary_tumor",
    module_score
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
      endpoint_name,
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
      direction_pass &
        p_pass &
        effect_pass
  ]
}

calculate_adjusted_results <- function(
  score_data,
  analysis_label,
  include_batch = FALSE
) {
  if (include_batch) {
    model <- stats::lm(
      module_score ~
        brain_indicator +
        factor(batch) +
        s_score_nonoverlap +
        g2m_score_nonoverlap,
      data =
        score_data
    )
  } else {
    model <- stats::lm(
      module_score ~
        brain_indicator +
        s_score_nonoverlap +
        g2m_score_nonoverlap,
      data =
        score_data
    )
  }

  coefficient_result <- extract_lm_term(
    model,
    "brain_indicator"
  )

  brain_residual <- score_data[
    cohort ==
      "Brain_metastasis",
    residualized_module_score
  ]

  primary_residual <- score_data[
    cohort ==
      "Primary_tumor",
    residualized_module_score
  ]

  residual_primary_median <- stats::median(
    primary_residual,
    na.rm = TRUE
  )

  residual_direction_count <- sum(
    brain_residual >
      residual_primary_median,
    na.rm = TRUE
  )

  residual_direction_required <- ceiling(
    0.65 *
      length(
        brain_residual
      )
  )

  cbind(
    data.table(
      analysis =
        analysis_label,
      brain_patients =
        score_data[
          cohort ==
            "Brain_metastasis",
          .N
        ],
      primary_patients =
        score_data[
          cohort ==
            "Primary_tumor",
          .N
        ],
      batches =
        paste(
          sort(
            unique(
              score_data$batch
            )
          ),
          collapse = ","
        )
    ),
    coefficient_result,
    data.table(
      adjusted_primary_pass =
        isTRUE(
          coefficient_result$positive_direction
        ) &&
          is.finite(
            coefficient_result$one_sided_positive_p
          ) &&
          coefficient_result$one_sided_positive_p <=
            0.05,
      residual_brain_median =
        stats::median(
          brain_residual,
          na.rm = TRUE
        ),
      residual_primary_median =
        residual_primary_median,
      residual_median_difference =
        stats::median(
          brain_residual,
          na.rm = TRUE
        ) -
          residual_primary_median,
      residual_brain_above_primary_median =
        residual_direction_count,
      residual_direction_required =
        residual_direction_required,
      residual_one_sided_wilcox_p =
        one_sided_wilcox(
          brain_residual,
          primary_residual
        ),
      residual_hedges_g =
        hedges_g_two_group(
          brain_residual,
          primary_residual
        ),
      residual_direction_pass =
        stats::median(
          brain_residual,
          na.rm = TRUE
        ) >
          residual_primary_median &&
          residual_direction_count >=
            residual_direction_required
    )
  )
}

run_analysis_set <- function(
  current_wide,
  module_gene_set,
  analysis_name,
  include_batch_model = FALSE
) {
  scores <- calculate_patient_scores(
    wide_data =
      current_wide,
    module_gene_set =
      module_gene_set,
    analysis_name =
      analysis_name
  )

  unadjusted <- calculate_unadjusted_gate(
    score_data =
      scores,
    endpoint_name =
      paste0(
        length(
          module_gene_set
        ),
        "_gene_module"
      )
  )

  adjusted <- calculate_adjusted_results(
    score_data =
      scores,
    analysis_label =
      analysis_name,
    include_batch =
      include_batch_model
  )

  list(
    scores =
      scores,
    unadjusted =
      unadjusted,
    adjusted =
      adjusted
  )
}

# ============================================================
# 阶段2：主要无重叠细胞周期校正
# ============================================================

banner(
  "阶段2：主要无重叠细胞周期校正"
)

main_wide <- expression_wide[
  author_tumor_nuclei >=
    50L
]

if (
  main_wide[
    cohort ==
      "Brain_metastasis",
    .N
  ] != 31L ||
    main_wide[
      cohort ==
        "Primary_tumor",
      .N
    ] != 10L
) {
  stop(
    "STEP08主要分析不是31例脑转移对10例原发。"
  )
}

main_result <- run_analysis_set(
  current_wide =
    main_wide,
  module_gene_set =
    module_genes,
  analysis_name =
    "STEP08_primary_31_vs_10_nooverlap_cycle",
  include_batch_model =
    FALSE
)

main_batch_adjusted_result <-
  calculate_adjusted_results(
    score_data =
      main_result$scores,
    analysis_label =
      "STEP08_primary_31_vs_10_full_batch_cycle_adjusted",
    include_batch =
      TRUE
  )

primary_summary <- cbind(
  main_result$unadjusted,
  main_result$adjusted[
    ,
    !c(
      "analysis",
      "brain_patients",
      "primary_patients"
    )
  ]
)

fwrite(
  primary_summary,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_primary_nooverlap_adjustment.tsv"
  ),
  sep = "\t"
)

fwrite(
  main_batch_adjusted_result,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_full_batch_cycle_adjusted.tsv"
  ),
  sep = "\t"
)

print(
  primary_summary
)

cat(
  "\n全样本批次+无重叠细胞周期模型：\n"
)

print(
  main_batch_adjusted_result
)

# ============================================================
# 阶段3：31对11全样本敏感性分析
# ============================================================

banner(
  "阶段3：31对11全样本无重叠校正"
)

all_result <- run_analysis_set(
  current_wide =
    expression_wide,
  module_gene_set =
    module_genes,
  analysis_name =
    "STEP08_sensitivity_31_vs_11_nooverlap_cycle",
  include_batch_model =
    FALSE
)

# ============================================================
# 阶段4：删除KRAS_17固定敏感性分析
# ============================================================

banner(
  "阶段4：删除KRAS_17敏感性分析"
)

without_kras17_wide <- main_wide[
  samples !=
    "KRAS_17"
]

if (
  without_kras17_wide[
    cohort ==
      "Brain_metastasis",
    .N
  ] != 30L ||
    without_kras17_wide[
      cohort ==
        "Primary_tumor",
      .N
    ] != 10L
) {
  stop(
    "删除KRAS_17后不是30例脑转移对10例原发。"
  )
}

without_kras17_result <- run_analysis_set(
  current_wide =
    without_kras17_wide,
  module_gene_set =
    module_genes,
  analysis_name =
    "STEP08_without_KRAS17_30_vs_10",
  include_batch_model =
    FALSE
)

kras17_exclusion_summary <- cbind(
  without_kras17_result$unadjusted,
  without_kras17_result$adjusted[
    ,
    !c(
      "analysis",
      "brain_patients",
      "primary_patients"
    )
  ]
)

fwrite(
  kras17_exclusion_summary,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_KRAS17_exclusion.tsv"
  ),
  sep = "\t"
)

print(
  kras17_exclusion_summary
)

# ============================================================
# 阶段5：固定细胞数阈值
# ============================================================

banner(
  "阶段5：固定细胞数阈值敏感性分析"
)

threshold_rows <- lapply(
  fixed_cell_thresholds,
  function(current_threshold) {
    current_wide <- expression_wide[
      author_tumor_nuclei >=
        current_threshold
    ]

    current_name <- paste0(
      "STEP08_cell_threshold_",
      current_threshold
    )

    if (
      current_wide[
        cohort ==
          "Brain_metastasis",
        .N
      ] <
        3L ||
        current_wide[
          cohort ==
            "Primary_tumor",
          .N
        ] <
          3L
    ) {
      return(
        data.table(
          analysis =
            current_name,
          threshold =
            current_threshold,
          brain_patients =
            current_wide[
              cohort ==
                "Brain_metastasis",
              .N
            ],
          primary_patients =
            current_wide[
              cohort ==
                "Primary_tumor",
              .N
            ],
          evaluable =
            FALSE
        )
      )
    }

    current_result <- run_analysis_set(
      current_wide =
        current_wide,
      module_gene_set =
        module_genes,
      analysis_name =
        current_name,
      include_batch_model =
        FALSE
    )

    cbind(
      data.table(
        threshold =
          current_threshold,
        evaluable =
          TRUE
      ),
      current_result$unadjusted,
      current_result$adjusted[
        ,
        !c(
          "analysis",
          "brain_patients",
          "primary_patients"
        )
      ]
    )
  }
)

threshold_summary <- rbindlist(
  threshold_rows,
  fill = TRUE
)

fwrite(
  threshold_summary,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_cell_count_thresholds.tsv"
  ),
  sep = "\t"
)

print(
  threshold_summary
)

# ============================================================
# 阶段6：批次内和重叠批次分析
# ============================================================

banner(
  "阶段6：批次内和重叠批次分析"
)

batch_rows <- list()

for (current_batch in c(
  "KRAS",
  "STK"
)) {
  current_wide <- main_wide[
    batch ==
      current_batch
  ]

  if (
    current_wide[
      cohort ==
        "Brain_metastasis",
      .N
    ] >=
      3L &&
      current_wide[
        cohort ==
          "Primary_tumor",
        .N
      ] >=
        3L
  ) {
    current_result <- run_analysis_set(
      current_wide =
        current_wide,
      module_gene_set =
        module_genes,
      analysis_name =
        paste0(
          "STEP08_",
          current_batch,
          "_within_batch"
        ),
      include_batch_model =
        FALSE
    )

    batch_rows[[current_batch]] <- cbind(
      current_result$unadjusted,
      current_result$adjusted[
        ,
        !c(
          "analysis",
          "brain_patients",
          "primary_patients"
        )
      ]
    )
  }
}

overlap_wide <- main_wide[
  batch %in%
    c(
      "KRAS",
      "STK"
    )
]

overlap_result <- run_analysis_set(
  current_wide =
    overlap_wide,
  module_gene_set =
    module_genes,
  analysis_name =
    "STEP08_KRAS_STK_overlap",
  include_batch_model =
    TRUE
)

batch_rows[["KRAS_STK_overlap"]] <- cbind(
  overlap_result$unadjusted,
  overlap_result$adjusted[
    ,
    !c(
      "analysis",
      "brain_patients",
      "primary_patients"
    )
  ]
)

batch_summary <- rbindlist(
  batch_rows,
  fill = TRUE
)

fwrite(
  batch_summary,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_batch_stratified.tsv"
  ),
  sep = "\t"
)

print(
  batch_summary
)

# ============================================================
# 阶段7：4基因无TYMS事后敏感性
# ============================================================

banner(
  "阶段7：4基因无TYMS事后敏感性"
)

four_gene_result <- run_analysis_set(
  current_wide =
    main_wide,
  module_gene_set =
    module_genes_without_TYMS,
  analysis_name =
    "STEP08_posthoc_4gene_without_TYMS",
  include_batch_model =
    FALSE
)

four_gene_summary <- cbind(
  four_gene_result$unadjusted,
  four_gene_result$adjusted[
    ,
    !c(
      "analysis",
      "brain_patients",
      "primary_patients"
    )
  ]
)

four_gene_summary[
  ,
  interpretation_rule :=
    "post_hoc_sensitivity_only_not_primary_validation"
]

fwrite(
  four_gene_summary,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_4gene_without_TYMS_sensitivity.tsv"
  ),
  sep = "\t"
)

print(
  four_gene_summary
)

# ============================================================
# 阶段8：相关性与综合结果
# ============================================================

banner(
  "阶段8：相关性与综合结果"
)

correlation_summary <- data.table(
  analysis =
    "STEP08_primary_31_vs_10",
  module_vs_S_nonoverlap =
    stats::cor(
      main_result$scores$module_score,
      main_result$scores$s_score_nonoverlap,
      use = "complete.obs",
      method = "pearson"
    ),
  module_vs_G2M_nonoverlap =
    stats::cor(
      main_result$scores$module_score,
      main_result$scores$g2m_score_nonoverlap,
      use = "complete.obs",
      method = "pearson"
    ),
  S_vs_G2M_nonoverlap =
    stats::cor(
      main_result$scores$s_score_nonoverlap,
      main_result$scores$g2m_score_nonoverlap,
      use = "complete.obs",
      method = "pearson"
    )
)

fwrite(
  correlation_summary,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_score_correlations.tsv"
  ),
  sep = "\t"
)

patient_score_columns <- c(
  "analysis",
  "patient",
  "samples",
  "cohort",
  "batch",
  "author_tumor_nuclei",
  "selected_library_size",
  "module_score",
  "s_score_nonoverlap",
  "g2m_score_nonoverlap",
  "residualized_module_score"
)

all_score_outputs <- rbindlist(
  list(
    main_result$scores[
      ,
      ..patient_score_columns
    ],
    all_result$scores[
      ,
      ..patient_score_columns
    ],
    without_kras17_result$scores[
      ,
      ..patient_score_columns
    ],
    overlap_result$scores[
      ,
      ..patient_score_columns
    ],
    four_gene_result$scores[
      ,
      ..patient_score_columns
    ]
  ),
  fill = TRUE
)

fwrite(
  all_score_outputs,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_patient_scores.tsv"
  ),
  sep = "\t"
)

overall_summary <- rbindlist(
  list(
    cbind(
      main_result$unadjusted,
      main_result$adjusted[
        ,
        !c(
          "analysis",
          "brain_patients",
          "primary_patients"
        )
      ]
    ),
    cbind(
      all_result$unadjusted,
      all_result$adjusted[
        ,
        !c(
          "analysis",
          "brain_patients",
          "primary_patients"
        )
      ]
    ),
    kras17_exclusion_summary,
    batch_summary,
    four_gene_summary
  ),
  fill = TRUE
)

fwrite(
  overall_summary,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_all_robustness_results.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 阶段9：固定解释状态
# ============================================================

banner(
  "阶段9：固定解释状态"
)

main_adjusted_pass <- isTRUE(
  main_result$adjusted$adjusted_primary_pass
)

kras17_raw_pass <- isTRUE(
  without_kras17_result$unadjusted$endpoint_pass
)

kras17_adjusted_positive <- isTRUE(
  without_kras17_result$adjusted$positive_direction
)

threshold_direction_consistency <- all(
  threshold_summary[
    evaluable == TRUE,
    median_difference >
      0
  ],
  na.rm = TRUE
)

if (
  main_adjusted_pass &&
    kras17_raw_pass &&
    kras17_adjusted_positive
) {
  final_status <-
    "STEP08_ROBUST_POSITIVE_INDEPENDENT_OF_NONOVERLAP_CELL_CYCLE"
} else if (
  isTRUE(
    main_result$unadjusted$endpoint_pass
  ) &&
    kras17_raw_pass &&
    threshold_direction_consistency
) {
  final_status <-
    "STEP08_ROBUST_PRIMARY_POSITIVE_WITH_PROLIFERATION_DEPENDENCE"
} else {
  final_status <-
    "STEP08_PRIMARY_POSITIVE_NOT_ROBUST_ACROSS_FIXED_SENSITIVITIES"
}

status_table <- data.table(
  final_status =
    final_status,
  STEP07_primary_positive =
    isTRUE(
      expected_primary$endpoint_pass
    ),
  STEP08_primary_unadjusted_pass =
    isTRUE(
      main_result$unadjusted$endpoint_pass
    ),
  STEP08_nonoverlap_adjusted_coefficient =
    main_result$adjusted$coefficient,
  STEP08_nonoverlap_adjusted_one_sided_p =
    main_result$adjusted$one_sided_positive_p,
  STEP08_nonoverlap_adjusted_pass =
    main_adjusted_pass,
  KRAS17_excluded_unadjusted_pass =
    kras17_raw_pass,
  KRAS17_excluded_adjusted_positive =
    kras17_adjusted_positive,
  threshold_direction_consistency =
    threshold_direction_consistency,
  four_gene_analysis_role =
    "post_hoc_sensitivity_only"
)

fwrite(
  status_table,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_final_status.tsv"
  ),
  sep = "\t"
)

print(
  status_table
)

# ============================================================
# 阶段10：图形
# ============================================================

banner(
  "阶段10：生成STEP08图形"
)

plot_data <- copy(
  main_result$scores
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

raw_plot <- ggplot(
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
    y = "Fixed five-gene score",
    title = "STEP08 unadjusted module score"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE223499_STEP08_unadjusted_module.png"
  ),
  raw_plot,
  width = 6.5,
  height = 5,
  dpi = 300
)

residual_plot <- ggplot(
  plot_data,
  aes(
    x = cohort_label,
    y = residualized_module_score
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
    y = "Residualized five-gene score",
    title = "STEP08 non-overlapping cell-cycle sensitivity"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE223499_STEP08_nonoverlap_residualized_module.png"
  ),
  residual_plot,
  width = 6.5,
  height = 5,
  dpi = 300
)

# ============================================================
# 阶段11：完成锁
# ============================================================

banner(
  "STEP08全部完成"
)

completion_file <- file.path(
  LOCK_DIR,
  "STEP_08_GSE223499_COMPLETE.txt"
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
      "final_status=",
      final_status
    ),
    "STEP07_results_modified=FALSE",
    "new_gene_selection_performed=FALSE",
    "fixed_module=DHFR,DHODH,SHMT1,TYMS,UMPS",
    "S_gene_overlap_removed=TYMS",
    "G2M_gene_overlap_removed=none",
    paste0(
      "primary_unadjusted_p=",
      main_result$unadjusted$one_sided_wilcox_p
    ),
    paste0(
      "primary_unadjusted_hedges_g=",
      main_result$unadjusted$hedges_g
    ),
    paste0(
      "nonoverlap_adjusted_coefficient=",
      main_result$adjusted$coefficient
    ),
    paste0(
      "nonoverlap_adjusted_one_sided_p=",
      main_result$adjusted$one_sided_positive_p
    ),
    paste0(
      "nonoverlap_residual_wilcox_p=",
      main_result$adjusted$residual_one_sided_wilcox_p
    ),
    paste0(
      "without_KRAS17_unadjusted_p=",
      without_kras17_result$unadjusted$one_sided_wilcox_p
    ),
    paste0(
      "without_KRAS17_adjusted_one_sided_p=",
      without_kras17_result$adjusted$one_sided_positive_p
    ),
    paste0(
      "full_batch_adjusted_one_sided_p=",
      main_batch_adjusted_result$one_sided_positive_p
    ),
    paste0(
      "threshold_direction_consistency=",
      threshold_direction_consistency
    ),
    "four_gene_without_TYMS_role=post_hoc_sensitivity_only",
    "post_result_model_change=not_allowed",
    "analysis_stops_after_STEP08=TRUE"
  ),
  completion_file
)

cat(
  "Status: ",
  final_status,
  "\n",
  sep = ""
)

cat(
  "\nKey output files:\n",
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_primary_nooverlap_adjustment.tsv"
  ),
  "\n",
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_KRAS17_exclusion.tsv"
  ),
  "\n",
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_cell_count_thresholds.tsv"
  ),
  "\n",
  file.path(
    TABLE_DIR,
    "GSE223499_STEP08_batch_stratified.tsv"
  ),
  "\n",
  completion_file,
  "\n",
  sep = ""
)
