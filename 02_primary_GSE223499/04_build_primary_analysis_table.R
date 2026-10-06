# Build the primary GSE223499 analysis table
# Extracts the prespecified primary strata from the earlier GSE223499 outputs and merges
# them into the 41-specimen analysis table used by the final patient-aware analysis.
# This step performs data assembly and consistency checks only.

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/LUAD_LM_PYRIMIDINE"

INPUT_DIR <- file.path(
  ROOT,
  "results/tables/GSE223499"
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

FILE_07C <- file.path(
  INPUT_DIR,
  "GSE223499_STEP07C_patient_level_scores.tsv"
)

FILE_08 <- file.path(
  INPUT_DIR,
  "GSE223499_STEP08_patient_scores.tsv"
)

PRIMARY_ANALYSIS_07C <- "Primary_31_vs_10"

PRIMARY_ANALYSIS_08 <-
  "STEP08_primary_31_vs_10_nooverlap_cycle"

cat("\n========================================\n")
cat("STEP A3：锁定31 BM vs 10 Primary患者级主分析表\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# 1. 输入检查
# ------------------------------------------------------------

required_files <- c(
  FILE_07C,
  FILE_08
)

if (!all(file.exists(required_files))) {
  stop(
    "缺少输入文件：\n",
    paste(
      required_files[
        !file.exists(required_files)
      ],
      collapse = "\n"
    )
  )
}

dt07 <- fread(
  FILE_07C,
  data.table = TRUE,
  showProgress = FALSE
)

dt08 <- fread(
  FILE_08,
  data.table = TRUE,
  showProgress = FALSE
)

required_07 <- c(
  "analysis",
  "patient",
  "samples",
  "cohort",
  "batch",
  "author_tumor_nuclei",
  "selected_library_size",
  "module_score",
  "s_score",
  "g2m_score",
  "adjusted_module_score"
)

required_08 <- c(
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

missing_07 <- setdiff(
  required_07,
  names(dt07)
)

missing_08 <- setdiff(
  required_08,
  names(dt08)
)

if (length(missing_07) > 0L) {
  stop(
    "STEP07C缺少字段：",
    paste(
      missing_07,
      collapse = ", "
    )
  )
}

if (length(missing_08) > 0L) {
  stop(
    "STEP08缺少字段：",
    paste(
      missing_08,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# 2. 审计两个文件中的所有 analysis strata
# ------------------------------------------------------------

summarize_strata <- function(dt, source_name) {
  dt[
    ,
    .(
      rows = .N,
      unique_patients = uniqueN(patient),
      n_brain = sum(
        cohort == "Brain_metastasis",
        na.rm = TRUE
      ),
      n_primary = sum(
        cohort == "Primary_tumor",
        na.rm = TRUE
      ),
      batches = paste(
        sort(
          unique(
            batch[
              !is.na(batch)
            ]
          )
        ),
        collapse = ";"
      )
    ),
    by = analysis
  ][
    ,
    source_file := source_name
  ][]
}

strata07 <- summarize_strata(
  dt07,
  basename(FILE_07C)
)

strata08 <- summarize_strata(
  dt08,
  basename(FILE_08)
)

strata_audit <- rbindlist(
  list(
    strata07,
    strata08
  ),
  fill = TRUE
)

setcolorder(
  strata_audit,
  c(
    "source_file",
    "analysis",
    "rows",
    "unique_patients",
    "n_brain",
    "n_primary",
    "batches"
  )
)

cat("STEP07C / STEP08 analysis strata审计：\n")
print(strata_audit)

fwrite(
  strata_audit,
  file.path(
    OUT_DIR,
    "STEP_A3_analysis_strata_audit.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# 3. 严格锁定预定主分析层
# ------------------------------------------------------------

if (
  !PRIMARY_ANALYSIS_07C %in%
    dt07$analysis
) {
  stop(
    "STEP07C中没有精确找到：",
    PRIMARY_ANALYSIS_07C,
    "\n现有analysis：\n",
    paste(
      unique(dt07$analysis),
      collapse = "\n"
    )
  )
}

if (
  !PRIMARY_ANALYSIS_08 %in%
    dt08$analysis
) {
  stop(
    "STEP08中没有精确找到：",
    PRIMARY_ANALYSIS_08,
    "\n现有analysis：\n",
    paste(
      unique(dt08$analysis),
      collapse = "\n"
    )
  )
}

primary07 <- dt07[
  analysis ==
    PRIMARY_ANALYSIS_07C
]

primary08 <- dt08[
  analysis ==
    PRIMARY_ANALYSIS_08
]

# ------------------------------------------------------------
# 4. 主分析样本数和唯一患者严格审计
# ------------------------------------------------------------

audit_primary_subset <- function(
  dt,
  label
) {
  out <- data.table(
    source = label,
    rows = nrow(dt),
    unique_patients =
      uniqueN(dt$patient),
    duplicate_patient_rows =
      nrow(dt) -
      uniqueN(dt$patient),
    n_brain = sum(
      dt$cohort ==
        "Brain_metastasis",
      na.rm = TRUE
    ),
    n_primary = sum(
      dt$cohort ==
        "Primary_tumor",
      na.rm = TRUE
    ),
    unique_batches =
      uniqueN(dt$batch)
  )

  out
}

primary_subset_audit <- rbindlist(
  list(
    audit_primary_subset(
      primary07,
      "STEP07C_Primary_31_vs_10"
    ),
    audit_primary_subset(
      primary08,
      "STEP08_primary_31_vs_10_nooverlap_cycle"
    )
  )
)

print(primary_subset_audit)

fwrite(
  primary_subset_audit,
  file.path(
    OUT_DIR,
    "STEP_A3_primary_subset_audit.tsv"
  ),
  sep = "\t"
)

for (nm in c("primary07", "primary08")) {
  current <- get(nm)

  if (nrow(current) != 41L) {
    stop(
      nm,
      " 主分析行数不是41，而是：",
      nrow(current)
    )
  }

  if (
    uniqueN(
      current$patient
    ) != 41L
  ) {
    stop(
      nm,
      " 主分析并非41名唯一患者。"
    )
  }

  if (
    anyDuplicated(
      current$patient
    )
  ) {
    stop(
      nm,
      " 存在重复患者ID。"
    )
  }

  if (
    sum(
      current$cohort ==
        "Brain_metastasis"
    ) != 31L
  ) {
    stop(
      nm,
      " 脑转移患者数不是31。"
    )
  }

  if (
    sum(
      current$cohort ==
        "Primary_tumor"
    ) != 10L
  ) {
    stop(
      nm,
      " 原发肿瘤患者数不是10。"
    )
  }
}

# ------------------------------------------------------------
# 5. 检查 STEP07C 和 STEP08 是否对应同一批41名患者
# ------------------------------------------------------------

patients07 <- sort(
  primary07$patient
)

patients08 <- sort(
  primary08$patient
)

if (!identical(
  patients07,
  patients08
)) {
  only07 <- setdiff(
    patients07,
    patients08
  )

  only08 <- setdiff(
    patients08,
    patients07
  )

  stop(
    "STEP07C与STEP08主分析患者集合不一致。\n",
    "仅STEP07C：",
    paste(
      only07,
      collapse = ", "
    ),
    "\n仅STEP08：",
    paste(
      only08,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# 6. 合并非重叠细胞周期信息
# ------------------------------------------------------------

merge08 <- primary08[
  ,
  .(
    patient,
    cohort_08 = cohort,
    batch_08 = batch,
    samples_08 = samples,
    module_score_08 =
      module_score,
    s_score_nonoverlap =
      s_score_nonoverlap,
    g2m_score_nonoverlap =
      g2m_score_nonoverlap,
    residualized_module_score =
      residualized_module_score
  )
]

locked <- merge(
  primary07[
    ,
    .(
      patient,
      samples,
      cohort,
      batch,
      author_tumor_nuclei,
      primary_eligible,
      selected_library_size,
      module_score,
      s_score_original =
        s_score,
      g2m_score_original =
        g2m_score,
      adjusted_module_score_original =
        adjusted_module_score
    )
  ],
  merge08,
  by = "patient",
  all = FALSE,
  sort = FALSE
)

if (nrow(locked) != 41L) {
  stop(
    "合并后不是41名患者：",
    nrow(locked)
  )
}

# ------------------------------------------------------------
# 7. 检查两版本元数据和module score是否完全一致
# ------------------------------------------------------------

locked[
  ,
  cohort_match :=
    cohort ==
      cohort_08
]

locked[
  ,
  batch_match :=
    batch ==
      batch_08
]

locked[
  ,
  sample_match :=
    as.character(samples) ==
      as.character(samples_08)
]

locked[
  ,
  module_score_abs_diff :=
    abs(
      module_score -
        module_score_08
    )
]

if (!all(
  locked$cohort_match
)) {
  stop(
    "STEP07C与STEP08存在cohort不一致。"
  )
}

if (!all(
  locked$batch_match
)) {
  stop(
    "STEP07C与STEP08存在batch不一致。"
  )
}

if (!all(
  locked$sample_match
)) {
  stop(
    "STEP07C与STEP08存在sample不一致。"
  )
}

if (
  max(
    locked$module_score_abs_diff,
    na.rm = TRUE
  ) >
    1e-10
) {
  stop(
    "STEP07C与STEP08的module_score不一致；最大绝对差：",
    max(
      locked$module_score_abs_diff,
      na.rm = TRUE
    )
  )
}

# ------------------------------------------------------------
# 8. 形成后续 robustness 唯一锁定输入表
# ------------------------------------------------------------

locked_final <- locked[
  ,
  .(
    patient,
    samples,
    cohort,
    batch,
    author_tumor_nuclei,
    primary_eligible,
    selected_library_size,
    module_score,
    s_score_original,
    g2m_score_original,
    adjusted_module_score_original,
    s_score_nonoverlap,
    g2m_score_nonoverlap,
    residualized_module_score
  )
]

setorder(
  locked_final,
  batch,
  cohort,
  patient
)

# ------------------------------------------------------------
# 9. Batch结构审计
# ------------------------------------------------------------

batch_structure <- locked_final[
  ,
  .(
    n_total = .N,
    n_brain = sum(
      cohort ==
        "Brain_metastasis"
    ),
    n_primary = sum(
      cohort ==
        "Primary_tumor"
    ),
    median_module_score =
      median(
        module_score,
        na.rm = TRUE
      )
  ),
  by = batch
]

batch_structure[
  ,
  overlapping_batch :=
    n_brain > 0L &
      n_primary > 0L
]

setorder(
  batch_structure,
  batch
)

cat(
  "\n固定主分析的Batch结构：\n"
)

print(
  batch_structure
)

# ------------------------------------------------------------
# 10. 保存
# ------------------------------------------------------------

locked_file <- file.path(
  OUT_DIR,
  "STEP_A3_locked_primary_31BM_10Primary_patient_table.tsv"
)

batch_file <- file.path(
  OUT_DIR,
  "STEP_A3_batch_structure.tsv"
)

fwrite(
  locked_final,
  locked_file,
  sep = "\t"
)

fwrite(
  batch_structure,
  batch_file,
  sep = "\t"
)

lock_file <- file.path(
  OUT_DIR,
  "STEP_A3_INPUT_LOCK.txt"
)

writeLines(
  c(
    paste0(
      "created_at=",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "STEP07C_analysis=",
      PRIMARY_ANALYSIS_07C
    ),
    paste0(
      "STEP08_analysis=",
      PRIMARY_ANALYSIS_08
    ),
    "inferential_unit=patient",
    "n_brain_metastasis=31",
    "n_primary_tumor=10",
    "primary_module_score_source=STEP07C_Primary_31_vs_10",
    "nonoverlap_cellcycle_source=STEP08_primary_31_vs_10_nooverlap_cycle",
    "patient_exclusion_after_result=none",
    "primary_endpoint_replacement=FALSE",
    "purpose=reviewer_defense_batch_sensitivity_only"
  ),
  lock_file
)

cat(
  "\n========================================\n",
  "STEP A3 complete\n",
  "========================================\n",
  "Locked patient-level table:\n",
  locked_file,
  "\n\nBatch structure:\n",
  batch_file,
  "\n\nInput lock:\n",
  lock_file,
  "\n\nReview the locked table and batch structure before downstream analyses.\n",
  "Key audit files: STEP_A3_primary_subset_audit.tsv and ",
  "STEP_A3_batch_structure.tsv.\n",
  sep = ""
)
