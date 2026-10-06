# Match the locked GSE223499 analysis set to official sample metadata
# Uses the locked specimen table and the parent-study Supplementary Table 1, with
# punctuation-insensitive sample-ID canonicalization where required. No histology is inferred
# from expression data and no analysis endpoint is changed.

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

ROOT <- "D:/LUAD_LM_PYRIMIDINE"

SUPP_FILE <- file.path(
  ROOT,
  "data/reference/GSE223499/Supplementary_Table_1.xlsx"
)

LOCKED_FILE <- file.path(
  ROOT,
  "results/reviewer_defense/STEP_A3_locked_primary_31BM_10Primary_patient_table.tsv"
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

banner <- function(x) {
  cat(
    "\n========================================\n",
    x,
    "\n========================================\n",
    sep = ""
  )
}

canon_id <- function(x) {
  z <- toupper(
    trimws(
      as.character(x)
    )
  )

  # Harmonize textual "dot" and punctuation.
  z <- gsub(
    "DOT",
    ".",
    z,
    fixed = TRUE
  )

  # Remove common technical prefixes only at the beginning.
  z <- sub(
    "^NSCLC[_-]*",
    "",
    z
  )

  z <- sub(
    "^NSCL[_-]*",
    "",
    z
  )

  z <- sub(
    "^SAMPLE[_-]*",
    "",
    z
  )

  # Ignore punctuation/spacing differences.
  z <- gsub(
    "[^A-Z0-9]+",
    "",
    z
  )

  z[
    !nzchar(z)
  ] <- NA_character_

  z
}

normalize_site <- function(x) {
  z <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA_character_,
    length(z)
  )

  out[
    grepl(
      "BRAIN",
      z
    )
  ] <- "Brain_metastasis"

  out[
    grepl(
      "PRIMARY",
      z
    )
  ] <- "Primary_tumor"

  out
}

banner("STEP B2 v2：稳健Paper ID匹配")

# ------------------------------------------------------------
# 1. Input checks
# ------------------------------------------------------------

if (!file.exists(SUPP_FILE)) {
  stop(
    "未找到官方Supplementary Table 1：\n",
    SUPP_FILE
  )
}

if (!file.exists(LOCKED_FILE)) {
  stop(
    "未找到STEP A3锁定患者表：\n",
    LOCKED_FILE
  )
}

# ------------------------------------------------------------
# 2. Read official table with the verified header position
# ------------------------------------------------------------

supp <- as.data.table(
  read_excel(
    SUPP_FILE,
    sheet = "Sheet1",
    skip = 2,
    .name_repair = "unique"
  )
)

names(supp) <- trimws(
  names(supp)
)

required_supp <- c(
  "Paper ID",
  "Sample type(s)"
)

missing_supp <- setdiff(
  required_supp,
  names(supp)
)

if (length(missing_supp) > 0L) {
  stop(
    "官方表缺少预期列：",
    paste(
      missing_supp,
      collapse = ", "
    )
  )
}

histology_cols <- grep(
  "hist|pathol|diagnos|histologic|histology|cancer type|tumou?r type|subtype",
  names(supp),
  ignore.case = TRUE,
  value = TRUE
)

official <- data.table(
  official_paper_id =
    trimws(
      as.character(
        supp[["Paper ID"]]
      )
    ),
  official_sample_type =
    trimws(
      as.character(
        supp[["Sample type(s)"]]
      )
    )
)

official <- official[
  !is.na(official_paper_id) &
    nzchar(official_paper_id)
]

official[
  ,
  match_key :=
    canon_id(
      official_paper_id
    )
]

official[
  ,
  official_site :=
    normalize_site(
      official_sample_type
    )
]

# ------------------------------------------------------------
# 3. Audit official-key uniqueness
# ------------------------------------------------------------

official_key_audit <- official[
  ,
  .(
    n_rows = .N,
    n_raw_ids =
      uniqueN(
        official_paper_id
      ),
    raw_ids =
      paste(
        unique(
          official_paper_id
        ),
        collapse = ";"
      ),
    n_site_labels =
      uniqueN(
        official_site[
          !is.na(
            official_site
          )
        ]
      ),
    site_labels =
      paste(
        unique(
          official_site[
            !is.na(
              official_site
            )
          ]
        ),
        collapse = ";"
      )
  ),
  by = match_key
]

conflicting_official <- official_key_audit[
  n_raw_ids > 1L |
    n_site_labels > 1L
]

fwrite(
  official_key_audit,
  file.path(
    OUT_DIR,
    "STEP_B2v2_01_official_ID_key_audit.tsv"
  ),
  sep = "\t"
)

if (nrow(conflicting_official) > 0L) {
  fwrite(
    conflicting_official,
    file.path(
      OUT_DIR,
      "STEP_B2v2_OFFICIAL_KEY_CONFLICTS.tsv"
    ),
    sep = "\t"
  )

  stop(
    "官方Paper ID标准化后出现冲突；已输出 STEP_B2v2_OFFICIAL_KEY_CONFLICTS.tsv。"
  )
}

official_unique <- official[
  ,
  .(
    official_paper_id =
      official_paper_id[1],
    official_sample_type =
      official_sample_type[1],
    official_site =
      official_site[1]
  ),
  by = match_key
]

# ------------------------------------------------------------
# 4. Read locked 41-patient table
# ------------------------------------------------------------

locked <- fread(
  LOCKED_FILE,
  data.table = TRUE,
  showProgress = FALSE
)

required_locked <- c(
  "patient",
  "samples",
  "cohort",
  "batch"
)

missing_locked <- setdiff(
  required_locked,
  names(locked)
)

if (length(missing_locked) > 0L) {
  stop(
    "锁定患者表缺少列：",
    paste(
      missing_locked,
      collapse = ", "
    )
  )
}

if (
  nrow(locked) != 41L ||
    uniqueN(
      locked$patient
    ) != 41L
) {
  stop(
    "锁定患者表必须为41名唯一患者。"
  )
}

# PRIMARY mapping key = samples.
locked[
  ,
  sample_match_key :=
    canon_id(
      samples
    )
]

# Patient-derived key is diagnostic only.
locked[
  ,
  patient_audit_key :=
    canon_id(
      patient
    )
]

if (anyNA(
  locked$sample_match_key
)) {
  fwrite(
    locked[
      is.na(
        sample_match_key
      )
    ],
    file.path(
      OUT_DIR,
      "STEP_B2v2_UNRESOLVED_SAMPLE_KEYS.tsv"
    ),
    sep = "\t"
  )

  stop(
    "存在无法从samples字段标准化的ID；已输出诊断表。"
  )
}

if (anyDuplicated(
  locked$sample_match_key
)) {
  dup <- locked[
    duplicated(
      sample_match_key
    ) |
      duplicated(
        sample_match_key,
        fromLast = TRUE
      )
  ]

  fwrite(
    dup,
    file.path(
      OUT_DIR,
      "STEP_B2v2_DUPLICATE_SAMPLE_KEYS.tsv"
    ),
    sep = "\t"
  )

  stop(
    "samples标准化后出现重复ID；已输出诊断表。"
  )
}

# Audit cases where patient technical name points to a different canonical ID.
locked[
  ,
  patient_vs_sample_key_same :=
    !is.na(
      patient_audit_key
    ) &
      patient_audit_key ==
        sample_match_key
]

id_mapping_audit <- locked[
  ,
  .(
    patient,
    samples,
    cohort,
    batch,
    patient_audit_key,
    sample_match_key,
    patient_vs_sample_key_same
  )
]

fwrite(
  id_mapping_audit,
  file.path(
    OUT_DIR,
    "STEP_B2v2_02_locked_ID_mapping_audit.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# 5. Deterministic match using sample_match_key only
# ------------------------------------------------------------

matched <- merge(
  locked,
  official_unique,
  by.x = "sample_match_key",
  by.y = "match_key",
  all.x = TRUE,
  sort = FALSE
)

matched[
  ,
  matched_to_official :=
    !is.na(
      official_paper_id
    )
]

matched[
  ,
  site_comparable :=
    !is.na(
      official_site
    )
]

matched[
  ,
  site_match :=
    site_comparable &
      cohort ==
        official_site
]

setorder(
  matched,
  cohort,
  batch,
  samples
)

fwrite(
  matched,
  file.path(
    OUT_DIR,
    "STEP_B2v2_03_locked41_to_official_metadata_match.tsv"
  ),
  sep = "\t"
)

unmatched <- matched[
  matched_to_official == FALSE
]

if (nrow(unmatched) > 0L) {
  fwrite(
    unmatched,
    file.path(
      OUT_DIR,
      "STEP_B2v2_UNMATCHED_LOCKED_SAMPLES.tsv"
    ),
    sep = "\t"
  )
}

site_mismatch <- matched[
  matched_to_official == TRUE &
    site_comparable == TRUE &
    site_match == FALSE
]

if (nrow(site_mismatch) > 0L) {
  fwrite(
    site_mismatch,
    file.path(
      OUT_DIR,
      "STEP_B2v2_SITE_MISMATCHES.tsv"
    ),
    sep = "\t"
  )
}

# ------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------

summary_match <- data.table(
  locked_patients =
    nrow(matched),
  matched_to_official =
    sum(
      matched$matched_to_official
    ),
  unmatched =
    sum(
      !matched$matched_to_official
    ),
  site_comparable =
    sum(
      matched$site_comparable
    ),
  site_matches =
    sum(
      matched$site_match,
      na.rm = TRUE
    ),
  site_mismatches =
    nrow(
      site_mismatch
    ),
  histology_candidate_columns =
    length(
      histology_cols
    )
)

print(
  summary_match
)

fwrite(
  summary_match,
  file.path(
    OUT_DIR,
    "STEP_B2v2_04_match_summary.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# 7. Histology-field status
# ------------------------------------------------------------

if (length(histology_cols) == 0L) {
  histology_status <- c(
    "histology_column_present=FALSE",
    "supplementary_table1_can_verify_sample_identity_and_site=TRUE",
    "supplementary_table1_can_verify_patient_level_LUAD_histology=FALSE",
    "marker_based_histology_inference=PROHIBITED",
    "safe_scope_pending_external_official_text_verification=NSCLC",
    "next_step=verify LUAD scope from original article text and/or other official supplementary metadata"
  )
} else {
  histology_status <- c(
    "histology_column_present=TRUE",
    paste0(
      "histology_columns=",
      paste(
        histology_cols,
        collapse = ";"
      )
    ),
    "next_step=patient-level extraction of official histology field"
  )
}

writeLines(
  histology_status,
  file.path(
    OUT_DIR,
    "STEP_B2v2_05_histology_field_status.txt"
  )
)

# ------------------------------------------------------------
# 8. Completion status
# ------------------------------------------------------------

overall_status <- if (
  summary_match$matched_to_official == 41L &&
    summary_match$site_mismatches == 0L
) {
  "OFFICIAL_SAMPLE_AND_SITE_MATCH_COMPLETE"
} else {
  "MATCH_REQUIRES_REVIEW"
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
    "matching_primary_key=locked samples -> official Paper ID",
    "patient_field_used_as_primary_key=FALSE",
    "punctuation_insensitive_matching=TRUE",
    paste0(
      "matched_to_official=",
      summary_match$matched_to_official
    ),
    paste0(
      "site_matches=",
      summary_match$site_matches
    ),
    paste0(
      "site_mismatches=",
      summary_match$site_mismatches
    ),
    paste0(
      "histology_candidate_columns=",
      summary_match$histology_candidate_columns
    ),
    paste0(
      "status=",
      overall_status
    ),
    "primary_endpoint_changed=FALSE"
  ),
  file.path(
    OUT_DIR,
    "STEP_B2v2_COMPLETE.txt"
  )
)

cat(
  "\n========================================\n",
  "STEP B2 v2 complete\n",
  "========================================\n",
  "Status: ",
  overall_status,
  "\n\n",
  "Key output files:\n",
  file.path(
    OUT_DIR,
    "STEP_B2v2_02_locked_ID_mapping_audit.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_B2v2_03_locked41_to_official_metadata_match.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_B2v2_04_match_summary.tsv"
  ),
  "\n\nReview the console summary together with these files.\n",
  sep = ""
)
