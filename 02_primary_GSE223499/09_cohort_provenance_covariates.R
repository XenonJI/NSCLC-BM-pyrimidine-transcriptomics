# Cohort provenance and covariate audit
# Reconciles the processed GSE223499 cohort with the source metadata, documents the matched
# patient relationships and exclusions, and assembles the technical/clinical covariates used
# in downstream sensitivity models.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

needed <- c("data.table", "readxl")

missing <- needed[
  !vapply(
    needed,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0L) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing, collapse = ", "),
      ". Install once manually, then rerun."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

ROOT <- "D:/LUAD_LM_PYRIMIDINE"
OUT <- file.path(ROOT, "results", "major_revision")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

OFFICIAL_XLSX <- file.path(
  ROOT,
  "data",
  "reference",
  "GSE223499",
  "Supplementary_Table_1.xlsx"
)

S1_CANDIDATES <- c(
  file.path(
    OUT,
    "STEP_E1_10_rebuilt_Supplementary_Table_S1_master_provenance.tsv"
  ),
  file.path(
    ROOT,
    "results",
    "final_figures_tables",
    "supplementary_source_tables",
    "Supplementary_Table_S1_FINAL.tsv"
  )
)

S1_FILE <- S1_CANDIDATES[
  file.exists(
    S1_CANDIDATES
  )
][1]

if (!file.exists(OFFICIAL_XLSX)) {
  stop(
    "Official Supplementary Table 1 not found:\n",
    OFFICIAL_XLSX,
    call. = FALSE
  )
}

if (is.na(S1_FILE) || !nzchar(S1_FILE)) {
  stop(
    "Could not find a current 43-specimen provenance table.",
    call. = FALSE
  )
}

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
norm_id <- function(x) {
  z <- toupper(
    trimws(
      as.character(x)
    )
  )

  z <- gsub(
    "DOT",
    "_",
    z,
    fixed = TRUE
  )

  gsub(
    "[^A-Z0-9]+",
    "",
    z
  )
}

clean_colname <- function(x) {
  z <- tolower(
    trimws(
      as.character(x)
    )
  )

  z <- gsub(
    "[\r\n]+",
    " ",
    z
  )

  z <- gsub(
    "[^a-z0-9]+",
    "_",
    z
  )

  z <- gsub(
    "^_+|_+$",
    "",
    z
  )

  make.unique(
    z,
    sep = "_dup"
  )
}

nonempty_fraction <- function(x) {
  z <- trimws(
    as.character(x)
  )

  mean(
    !is.na(z) &
      nzchar(z)
  )
}

candidate_cols <- function(nms, patterns) {
  hit <- rep(FALSE, length(nms))

  for (pp in patterns) {
    hit <- hit |
      grepl(
        pp,
        nms,
        ignore.case = TRUE,
        perl = TRUE
      )
  }

  nms[hit]
}

first_existing <- function(nms, candidates) {
  z <- candidates[
    candidates %in% nms
  ]

  if (length(z) == 0L) {
    NA_character_
  } else {
    z[1]
  }
}

# ----------------------------------------------------------------------
# 1. Read current processed-site provenance
# ----------------------------------------------------------------------
s1 <- fread(
  S1_FILE,
  data.table = TRUE,
  showProgress = FALSE
)

# Accept either old or E1-rebuilt naming.
sample_col <- first_existing(
  names(s1),
  c("sample", "samples")
)

tumor_n_col <- first_existing(
  names(s1),
  c("author_tumor_cells", "author_tumor_nuclei")
)

site_col <- first_existing(
  names(s1),
  c("site_author", "site_analysis")
)

if (anyNA(
  c(
    sample_col,
    tumor_n_col,
    site_col
  )
)) {
  stop(
    "Current provenance source lacks required sample/site/tumor-count fields.",
    call. = FALSE
  )
}

base <- s1[
  ,
  .(
    sample =
      as.character(
        get(sample_col)
      ),
    author_tumor_nuclei =
      as.numeric(
        get(tumor_n_col)
      ),
    site_source =
      as.character(
        get(site_col)
      )
  )
]

# Normalize site labels from either old or rebuilt S1.
base[
  ,
  site_author := fifelse(
    grepl(
      "brain",
      tolower(site_source)
    ),
    "BRAIN_METS",
    fifelse(
      grepl(
        "primary",
        tolower(site_source)
      ),
      "PRIMARY",
      fifelse(
        grepl(
          "chest",
          tolower(site_source)
        ),
        "CHEST_WALL_MET",
        site_source
      )
    )
  )
]

base[
  ,
  sample_key := norm_id(sample)
]

if (uniqueN(base$sample) != 43L) {
  warning(
    "Expected 43 source specimens but found ",
    uniqueN(base$sample),
    ". Audit will continue."
  )
}

# ----------------------------------------------------------------------
# 2. Auto-detect the real header row in Supplementary Table 1
# ----------------------------------------------------------------------
sheets <- excel_sheets(
  OFFICIAL_XLSX
)

header_candidates <- list()
raw_sheets <- list()

header_keywords <- c(
  "sample",
  "patient",
  "subject",
  "sex",
  "gender",
  "age",
  "site",
  "brain",
  "primary",
  "histology",
  "mutation",
  "kras",
  "stk11",
  "10x",
  "chemistry",
  "library"
)

for (sh in sheets) {
  raw <- as.data.table(
    read_excel(
      OFFICIAL_XLSX,
      sheet = sh,
      col_names = FALSE,
      .name_repair = "minimal"
    )
  )

  raw_sheets[[sh]] <- raw

  nscan <- min(
    20L,
    nrow(raw)
  )

  for (rr in seq_len(nscan)) {
    vals <- tolower(
      trimws(
        as.character(
          unlist(
            raw[rr],
            use.names = FALSE
          )
        )
      )
    )

    vals <- vals[
      !is.na(vals) &
        nzchar(vals)
    ]

    if (length(vals) == 0L) next

    keyword_hits <- sum(
      vapply(
        header_keywords,
        function(k) {
          any(
            grepl(
              k,
              vals,
              fixed = TRUE
            )
          )
        },
        logical(1)
      )
    )

    sample_id_hits <- sum(
      norm_id(vals) %chin%
        base$sample_key
    )

    header_candidates[[length(header_candidates) + 1L]] <- data.table(
      sheet = sh,
      row = rr,
      keyword_hits = keyword_hits,
      sample_id_hits = sample_id_hits,
      nonempty_cells = length(vals),
      preview = paste(
        head(vals, 12L),
        collapse = " | "
      )
    )
  }
}

header_audit <- rbindlist(
  header_candidates,
  fill = TRUE
)

setorder(
  header_audit,
  -keyword_hits,
  sample_id_hits,
  row
)

fwrite(
  header_audit,
  file.path(
    OUT,
    "STEP_E4_01_official_header_row_audit.tsv"
  ),
  sep = "\t"
)

if (nrow(header_audit) == 0L) {
  stop(
    "Could not audit Supplementary Table 1 header rows.",
    call. = FALSE
  )
}

best_header <- header_audit[1]

# Require a reasonably plausible header row.
if (best_header$keyword_hits < 2L) {
  stop(
    "Could not confidently identify the official table header row. ",
    "Inspect STEP_E4_01_official_header_row_audit.tsv.",
    call. = FALSE
  )
}

official_sheet <- best_header$sheet
header_row <- best_header$row

raw <- raw_sheets[[official_sheet]]

header_vals <- as.character(
  unlist(
    raw[header_row],
    use.names = FALSE
  )
)

header_vals[
  is.na(header_vals) |
    !nzchar(
      trimws(
        header_vals
      )
    )
] <- paste0(
  "unnamed_",
  which(
    is.na(header_vals) |
      !nzchar(
        trimws(
          header_vals
        )
      )
  )
)

colnames_raw <- clean_colname(
  header_vals
)

if (header_row >= nrow(raw)) {
  stop("Detected header row is the final workbook row.", call. = FALSE)
}

official <- copy(
  raw[
    (header_row + 1L):nrow(raw)
  ]
)

setnames(
  official,
  colnames_raw
)

# Drop completely empty rows.
row_nonempty <- apply(
  official,
  1,
  function(z) {
    any(
      !is.na(z) &
        nzchar(
          trimws(
            as.character(z)
          )
        )
    )
  }
)

official <- official[
  row_nonempty
]

# ----------------------------------------------------------------------
# 3. Identify official sample-ID column by maximum exact overlap
# ----------------------------------------------------------------------
id_overlap <- rbindlist(
  lapply(
    names(official),
    function(cn) {
      vals <- norm_id(
        official[[cn]]
      )

      data.table(
        column = cn,
        overlap_with_43_ids = sum(
          unique(vals) %chin%
            base$sample_key
        ),
        nonempty_fraction = nonempty_fraction(
          official[[cn]]
        ),
        unique_nonempty = uniqueN(
          trimws(
            as.character(
              official[[cn]]
            )
          )[
            !is.na(
              official[[cn]]
            ) &
              nzchar(
                trimws(
                  as.character(
                    official[[cn]]
                  )
                )
              )
          ]
        )
      )
    }
  )
)

setorder(
  id_overlap,
  -overlap_with_43_ids,
  -nonempty_fraction
)

fwrite(
  id_overlap,
  file.path(
    OUT,
    "STEP_E4_02_official_sample_id_column_audit.tsv"
  ),
  sep = "\t"
)

best_id <- id_overlap[1]

if (best_id$overlap_with_43_ids < 40L) {
  stop(
    "Official table could not be matched to at least 40/43 processed sample IDs. ",
    "Inspect STEP_E4_02_official_sample_id_column_audit.tsv.",
    call. = FALSE
  )
}

official_sample_col <- best_id$column

official[
  ,
  official_sample_id :=
    as.character(
      get(
        official_sample_col
      )
    )
]

official[
  ,
  sample_key := norm_id(
    official_sample_id
  )
]

official <- official[
  sample_key %chin%
    base$sample_key
]

official <- unique(
  official,
  by = "sample_key"
)

master <- merge(
  base,
  official,
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)

# ----------------------------------------------------------------------
# 4. Corrected patient identity + analysis cohort
# ----------------------------------------------------------------------
# Prefer official patient field when present; otherwise use sample identity
# and explicitly correct the two matched pairs.
patient_candidates <- candidate_cols(
  names(master),
  c(
    "^patient$",
    "patient_id",
    "^subject$",
    "subject_id",
    "^donor$",
    "case_id"
  )
)

official_patient_col <- if (
  length(patient_candidates) > 0L
) {
  patient_candidates[1]
} else {
  NA_character_
}

master[
  ,
  corrected_patient_id := sample
]

if (!is.na(official_patient_col)) {
  z <- trimws(
    as.character(
      master[[official_patient_col]]
    )
  )

  use <- !is.na(z) &
    nzchar(z)

  master[
    use,
    corrected_patient_id := z[use]
  ]
}

# Source-study known matched pairs override missing/inconsistent patient fields.
master[
  sample %chin% c("PA060", "N254"),
  corrected_patient_id := "MATCHED_PATIENT_PA060_N254"
]

master[
  sample %chin% c("PA067", "N586"),
  corrected_patient_id := "MATCHED_PATIENT_PA067_N586"
]

master[
  ,
  matched_pair_id := fifelse(
    sample %chin% c("PA060", "N254"),
    "PA060_N254",
    fifelse(
      sample %chin% c("PA067", "N586"),
      "PA067_N586",
      ""
    )
  )
]

master[
  ,
  site_analysis := fifelse(
    site_author == "BRAIN_METS",
    "Brain_metastasis",
    fifelse(
      site_author == "PRIMARY",
      "Primary_tumor",
      "Excluded_non_target_site"
    )
  )
]

master[
  ,
  threshold_pass_50 := author_tumor_nuclei >= 50L
]

master[
  ,
  included_primary_analysis :=
    site_author %chin% c("BRAIN_METS", "PRIMARY") &
    threshold_pass_50
]

master[
  ,
  exclusion_reason := fifelse(
    site_author == "CHEST_WALL_MET",
    "Excluded from PT-vs-BM target-site comparison because the author-provided processed site field is CHEST_WALL_MET",
    fifelse(
      !threshold_pass_50,
      "Excluded from the locked primary >=50 malignant-nucleus analysis",
      ""
    )
  )
]

master[
  ,
  sample_prefix_stratum := fifelse(
    grepl("^KRAS", sample, ignore.case = TRUE),
    "KRAS",
    fifelse(
      grepl("^STK", sample, ignore.case = TRUE),
      "STK",
      fifelse(
        grepl("^PA", sample, ignore.case = TRUE),
        "PA",
        fifelse(
          grepl("^N", sample, ignore.case = TRUE),
          "N",
          "Other"
        )
      )
    )
  )
]

master[
  ,
  prefix_stratum_interpretation :=
    "Sample-prefix/processing stratum; not assumed to be a pure technical sequencing batch unless independently documented"
]

master[
  ,
  malignant_annotation_interpretation :=
    "Malignant/tumor identity was inherited from the parent study's curated annotation and was not independently re-called from CNV in this reanalysis"
]

# ----------------------------------------------------------------------
# 5. Inventory clinical / biological / technical fields
# ----------------------------------------------------------------------
nms <- names(master)

field_sets <- list(
  age = candidate_cols(
    nms,
    c("^age$", "age_", "_age")
  ),
  sex = candidate_cols(
    nms,
    c("^sex$", "^gender$", "sex_", "gender_")
  ),
  histology = candidate_cols(
    nms,
    c("histolog", "subtype", "adenocarc", "squamous")
  ),
  genotype_driver = candidate_cols(
    nms,
    c(
      "kras",
      "stk11",
      "egfr",
      "tp53",
      "alk",
      "ros1",
      "braf",
      "met",
      "ret",
      "driver",
      "mutation",
      "genotype"
    )
  ),
  chemistry_library = candidate_cols(
    nms,
    c(
      "10x",
      "chem",
      "library",
      "protocol",
      "kit",
      "platform",
      "sequenc",
      "lane",
      "run",
      "prep"
    )
  ),
  timing_treatment = candidate_cols(
    nms,
    c(
      "time",
      "date",
      "treatment",
      "therapy",
      "naive",
      "collection"
    )
  )
)

field_inventory <- rbindlist(
  lapply(
    names(field_sets),
    function(category) {
      cols <- field_sets[[category]]

      if (length(cols) == 0L) {
        return(
          data.table(
            category = category,
            column = NA_character_,
            nonempty_fraction_43 = NA_real_,
            unique_nonempty_43 = NA_integer_,
            nonempty_fraction_primary_analysis = NA_real_,
            unique_nonempty_primary_analysis = NA_integer_
          )
        )
      }

      rbindlist(
        lapply(
          cols,
          function(cn) {
            v43 <- master[[cn]]
            v41 <- master[
              included_primary_analysis == TRUE,
              get(cn)
            ]

            data.table(
              category = category,
              column = cn,
              nonempty_fraction_43 = nonempty_fraction(v43),
              unique_nonempty_43 = uniqueN(
                trimws(
                  as.character(v43)
                )[
                  !is.na(v43) &
                    nzchar(
                      trimws(
                        as.character(v43)
                      )
                    )
                ]
              ),
              nonempty_fraction_primary_analysis = nonempty_fraction(v41),
              unique_nonempty_primary_analysis = uniqueN(
                trimws(
                  as.character(v41)
                )[
                  !is.na(v41) &
                    nzchar(
                      trimws(
                        as.character(v41)
                      )
                    )
                ]
              )
            )
          }
        )
      )
    }
  ),
  fill = TRUE
)

fwrite(
  field_inventory,
  file.path(
    OUT,
    "STEP_E4_03_clinical_technical_field_inventory.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 6. Explicit N561 provenance record
# ----------------------------------------------------------------------
n561 <- master[
  sample == "N561"
]

if (nrow(n561) != 1L) {
  stop(
    "N561 was not uniquely identified in the 43-specimen reconstruction.",
    call. = FALSE
  )
}

fwrite(
  n561,
  file.path(
    OUT,
    "STEP_E4_04_N561_PROVENANCE_AUDIT.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 7. Cohort flow
# ----------------------------------------------------------------------
flow <- data.table(
  stage = c(
    "Retrieved processed specimens",
    "Processed site = BRAIN_METS",
    "Processed site = PRIMARY",
    "Processed site = CHEST_WALL_MET",
    "Target BM/PT specimens before >=50 threshold",
    "Locked primary BM specimens >=50",
    "Locked primary PT specimens >=50",
    "Locked primary specimens total",
    "Locked primary unique patients after matched-pair correction"
  ),
  n = c(
    uniqueN(master$sample),
    sum(master$site_author == "BRAIN_METS"),
    sum(master$site_author == "PRIMARY"),
    sum(master$site_author == "CHEST_WALL_MET"),
    sum(master$site_author %chin% c("BRAIN_METS", "PRIMARY")),
    sum(
      master$site_author == "BRAIN_METS" &
        master$threshold_pass_50
    ),
    sum(
      master$site_author == "PRIMARY" &
        master$threshold_pass_50
    ),
    sum(master$included_primary_analysis),
    uniqueN(
      master[
        included_primary_analysis == TRUE,
        corrected_patient_id
      ]
    )
  )
)

fwrite(
  flow,
  file.path(
    OUT,
    "STEP_E4_05_COHORT_FLOW.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 8. reporting master provenance table
# ----------------------------------------------------------------------
front <- c(
  "sample",
  "corrected_patient_id",
  "matched_pair_id",
  "site_author",
  "site_analysis",
  "author_tumor_nuclei",
  "threshold_pass_50",
  "included_primary_analysis",
  "exclusion_reason",
  "sample_prefix_stratum",
  "prefix_stratum_interpretation",
  "official_sample_id",
  "malignant_annotation_interpretation"
)

front <- front[
  front %in%
    names(master)
]

setcolorder(
  master,
  c(
    front,
    setdiff(
      names(master),
      front
    )
  )
)

fwrite(
  master,
  file.path(
    OUT,
    "STEP_E4_06_MASTER_SAMPLE_PROVENANCE_WITH_OFFICIAL_METADATA.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 9. Decide whether a genuine documented technical covariate is available
# ----------------------------------------------------------------------
technical_available <- field_inventory[
  category == "chemistry_library" &
    !is.na(column) &
    nonempty_fraction_primary_analysis >= 0.80 &
    unique_nonempty_primary_analysis >= 2L
]

clinical_available <- field_inventory[
  category %chin% c("age", "sex", "genotype_driver") &
    !is.na(column) &
    nonempty_fraction_primary_analysis >= 0.80
]

status <- "PASS_E4_PROVENANCE_COVARIATE_AUDIT"

recommendation <- if (nrow(technical_available) > 0L) {
  paste0(
    "At least one documented chemistry/library/platform field has >=80% coverage ",
    "and variation in the locked primary cohort. Use STEP E4 outputs to define ",
    "a genuinely technical covariate sensitivity model; do not equate KRAS/STK/PA/N ",
    "prefixes with sequencing batch."
  )
} else {
  paste0(
    "No well-covered varying chemistry/library/platform field was automatically ",
    "identified. Retain KRAS/STK/PA/N only as sample-prefix/processing strata and ",
    "describe residual technical/genomic confounding transparently rather than ",
    "claiming full technical-batch adjustment."
  )
}

writeLines(
  c(
    paste0("status=", status),
    paste0("official_sheet=", official_sheet),
    paste0("detected_header_row=", header_row),
    paste0("official_sample_column=", official_sample_col),
    paste0("official_sample_overlap=", best_id$overlap_with_43_ids),
    paste0("processed_site_BM=", sum(master$site_author == "BRAIN_METS")),
    paste0("processed_site_PRIMARY=", sum(master$site_author == "PRIMARY")),
    paste0("processed_site_CHEST_WALL_MET=", sum(master$site_author == "CHEST_WALL_MET")),
    paste0(
      "primary_analysis_unique_patients_corrected=",
      uniqueN(
        master[
          included_primary_analysis == TRUE,
          corrected_patient_id
        ]
      )
    ),
    paste0(
      "N561_processed_site=",
      n561$site_author[1]
    ),
    paste0(
      "documented_technical_fields_eligible_for_followup=",
      paste(
        technical_available$column,
        collapse = ";"
      )
    ),
    paste0(
      "clinical_fields_with_ge80pct_coverage=",
      paste(
        unique(clinical_available$column),
        collapse = ";"
      )
    ),
    paste0("recommendation=", recommendation),
    "module_score_outcomes_inspected_for_covariate_selection=FALSE",
    "biological_model_fitted=FALSE",
    "primary_endpoint_changed=FALSE"
  ),
  file.path(
    OUT,
    "STEP_E4_COMPLETE.txt"
  )
)

cat(
  "\n============================================================\n",
  "STEP E4: OFFICIAL COHORT PROVENANCE + COVARIATE AUDIT\n",
  "============================================================\n",
  "Official workbook sheet: ", official_sheet, "\n",
  "Detected true header row: ", header_row, "\n",
  "Official sample-ID column: ", official_sample_col, "\n",
  "Matched processed sample IDs: ", best_id$overlap_with_43_ids, "/43\n\n",
  "Cohort flow:\n",
  sep = ""
)

print(flow)

cat(
  "\nN561 processed-site record:\n"
)

print(
  n561[
    ,
    intersect(
      c(
        "sample",
        "site_author",
        "site_analysis",
        "official_sample_id",
        names(n561)
      ),
      names(n561)
    ),
    with = FALSE
  ]
)

cat(
  "\nCandidate clinical/technical fields:\n"
)

print(field_inventory)

cat(
  "\nRecommendation:\n",
  recommendation,
  "\n\nStatus: ",
  status,
  "\n\nKey output files:\n",
  file.path(OUT, "STEP_E4_03_clinical_technical_field_inventory.tsv"), "\n",
  file.path(OUT, "STEP_E4_04_N561_PROVENANCE_AUDIT.tsv"), "\n",
  file.path(OUT, "STEP_E4_05_COHORT_FLOW.tsv"), "\n",
  file.path(OUT, "STEP_E4_06_MASTER_SAMPLE_PROVENANCE_WITH_OFFICIAL_METADATA.tsv"), "\n",
  file.path(OUT, "STEP_E4_COMPLETE.txt"), "\n",
  sep = ""
)
