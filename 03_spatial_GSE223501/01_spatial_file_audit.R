# GSE223501 spatial file audit
# Verifies the Slide-seq count/coordinate file inventory and source-study sample metadata
# before any spatial expression analysis is performed.

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

ROOT <- "D:/LUAD_LM_PYRIMIDINE"

SPATIAL_DIR <- file.path(
  ROOT,
  "data/raw/GSE223501"
)

TAR_FILE <- file.path(
  SPATIAL_DIR,
  "GSE223501_RAW.tar"
)

SUPP_FILE <- file.path(
  ROOT,
  "data/reference/GSE223499/Supplementary_Table_1.xlsx"
)

OUT_DIR <- file.path(
  ROOT,
  "results/reviewer_defense"
)

dir.create(
  SPATIAL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
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

  z <- gsub(
    "DOT",
    ".",
    z,
    fixed = TRUE
  )

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

parse_spatial_filename <- function(path) {
  b <- basename(path)

  gsm <- sub(
    "^(GSM[0-9]+).*",
    "\\1",
    b
  )

  sample <- b

  sample <- sub(
    "^GSM[0-9]+_NSCLC_",
    "",
    sample,
    ignore.case = TRUE
  )

  sample <- sub(
    "_slideseq_(counts|coords)\\.csv\\.gz$",
    "",
    sample,
    ignore.case = TRUE
  )

  data.table(
    gsm = gsm,
    sample_raw = sample,
    sample_key = canon_id(sample)
  )
}

banner("STEP C1：GSE223501 Slide-seq 文件与样本审计")

# ------------------------------------------------------------
# 1. 检查Supplementary Table 1
# ------------------------------------------------------------

if (!file.exists(SUPP_FILE)) {
  stop(
    "缺少官方Supplementary Table 1：\n",
    SUPP_FILE
  )
}

# ------------------------------------------------------------
# 2. 检查/解压RAW.tar
# ------------------------------------------------------------

count_files_before <- list.files(
  SPATIAL_DIR,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "_slideseq_counts\\.csv\\.gz$",
  ignore.case = TRUE
)

coord_files_before <- list.files(
  SPATIAL_DIR,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "_slideseq_coords\\.csv\\.gz$",
  ignore.case = TRUE
)

if (
  length(count_files_before) == 0L &&
    length(coord_files_before) == 0L
) {
  if (!file.exists(TAR_FILE)) {
    stop(
      paste0(
        "尚未找到 GSE223501_RAW.tar，也没有找到已解压的 Slide-seq 文件。\n\n",
        "请人工下载 GSE223501_RAW.tar 并保存为：\n",
        TAR_FILE
      )
    )
  }

  cat(
    "尚未发现已解压文件，现在自动解压：\n",
    TAR_FILE,
    "\n\n",
    sep = ""
  )

  utils::untar(
    TAR_FILE,
    exdir = SPATIAL_DIR
  )

  gc()
} else {
  cat(
    "已经发现解压后的Slide-seq文件，跳过再次解压。\n"
  )
}

# ------------------------------------------------------------
# 3. 文件清单
# ------------------------------------------------------------

count_files <- list.files(
  SPATIAL_DIR,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "_slideseq_counts\\.csv\\.gz$",
  ignore.case = TRUE
)

coord_files <- list.files(
  SPATIAL_DIR,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "_slideseq_coords\\.csv\\.gz$",
  ignore.case = TRUE
)

cat(
  "\ncounts files: ",
  length(count_files),
  "\ncoords files: ",
  length(coord_files),
  "\n",
  sep = ""
)

if (
  length(count_files) == 0L ||
    length(coord_files) == 0L
) {
  stop(
    "没有同时找到 counts 与 coords 文件。"
  )
}

# ------------------------------------------------------------
# 4. 解析文件名
# ------------------------------------------------------------

count_info <- rbindlist(
  lapply(
    count_files,
    function(f) {
      x <- parse_spatial_filename(f)

      x[
        ,
        `:=`(
          counts_file = f,
          counts_size_MB = round(
            file.info(f)$size /
              1024^2,
            3
          )
        )
      ]

      x
    }
  )
)

coord_info <- rbindlist(
  lapply(
    coord_files,
    function(f) {
      x <- parse_spatial_filename(f)

      x[
        ,
        `:=`(
          coords_file = f,
          coords_size_MB = round(
            file.info(f)$size /
              1024^2,
            3
          )
        )
      ]

      x
    }
  )
)

if (anyDuplicated(
  count_info$sample_key
)) {
  stop(
    "counts文件中出现重复sample key。"
  )
}

if (anyDuplicated(
  coord_info$sample_key
)) {
  stop(
    "coords文件中出现重复sample key。"
  )
}

paired <- merge(
  count_info,
  coord_info[
    ,
    .(
      sample_key,
      gsm_coords = gsm,
      sample_raw_coords = sample_raw,
      coords_file,
      coords_size_MB
    )
  ],
  by = "sample_key",
  all = TRUE,
  sort = FALSE
)

paired[
  ,
  counts_present :=
    !is.na(
      counts_file
    )
]

paired[
  ,
  coords_present :=
    !is.na(
      coords_file
    )
]

paired[
  ,
  pair_complete :=
    counts_present &
      coords_present
]

paired[
  ,
  gsm_match :=
    !is.na(gsm) &
      !is.na(gsm_coords) &
      gsm ==
        gsm_coords
]

# ------------------------------------------------------------
# 5. 读取官方Supplementary Table 1样本信息
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

required_cols <- c(
  "Paper ID",
  "Sample type(s)"
)

missing_cols <- setdiff(
  required_cols,
  names(supp)
)

if (length(missing_cols) > 0L) {
  stop(
    "Supplementary Table 1缺少字段：",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

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
  sample_key :=
    canon_id(
      official_paper_id
    )
]

official[
  ,
  cohort :=
    normalize_site(
      official_sample_type
    )
]

official_unique <- official[
  ,
  .(
    official_paper_id =
      official_paper_id[1],
    official_sample_type =
      official_sample_type[1],
    cohort =
      {
        x <- unique(
          cohort[
            !is.na(cohort)
          ]
        )

        if (length(x) == 0L) {
          NA_character_
        } else {
          x[1]
        }
      }
  ),
  by = sample_key
]

paired <- merge(
  paired,
  official_unique,
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)

paired[
  ,
  official_metadata_match :=
    !is.na(
      official_paper_id
    )
]

setorder(
  paired,
  cohort,
  official_paper_id
)

# ------------------------------------------------------------
# 6. 汇总
# ------------------------------------------------------------

summary_dt <- data.table(
  n_count_files =
    length(count_files),
  n_coord_files =
    length(coord_files),
  n_unique_spatial_samples =
    nrow(paired),
  n_complete_pairs =
    sum(
      paired$pair_complete
    ),
  n_gsm_matches =
    sum(
      paired$gsm_match,
      na.rm = TRUE
    ),
  n_official_metadata_matches =
    sum(
      paired$official_metadata_match
    ),
  n_brain_metastasis =
    sum(
      paired$cohort ==
        "Brain_metastasis",
      na.rm = TRUE
    ),
  n_primary_tumor =
    sum(
      paired$cohort ==
        "Primary_tumor",
      na.rm = TRUE
    )
)

print(
  paired[
    ,
    .(
      gsm,
      sample_raw,
      official_paper_id,
      cohort,
      counts_size_MB,
      coords_size_MB,
      pair_complete,
      official_metadata_match
    )
  ]
)

cat(
  "\nSTEP C1 summary：\n"
)

print(
  summary_dt
)

# ------------------------------------------------------------
# 7. 严格QC状态
# ------------------------------------------------------------

qc_pass <- (
  summary_dt$n_count_files == 14L &&
    summary_dt$n_coord_files == 14L &&
    summary_dt$n_unique_spatial_samples == 14L &&
    summary_dt$n_complete_pairs == 14L &&
    summary_dt$n_gsm_matches == 14L &&
    summary_dt$n_official_metadata_matches == 14L
)

status <- if (
  isTRUE(
    qc_pass
  )
) {
  "PASS_14_SPATIAL_SAMPLES_COMPLETE"
} else {
  "REVIEW_REQUIRED"
}

# ------------------------------------------------------------
# 8. 输出
# ------------------------------------------------------------

fwrite(
  paired,
  file.path(
    OUT_DIR,
    "STEP_C1_01_GSE223501_spatial_file_inventory.tsv"
  ),
  sep = "\t"
)

fwrite(
  summary_dt,
  file.path(
    OUT_DIR,
    "STEP_C1_02_GSE223501_spatial_inventory_summary.tsv"
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
      "count_files=",
      summary_dt$n_count_files
    ),
    paste0(
      "coord_files=",
      summary_dt$n_coord_files
    ),
    paste0(
      "complete_pairs=",
      summary_dt$n_complete_pairs
    ),
    paste0(
      "official_metadata_matches=",
      summary_dt$n_official_metadata_matches
    ),
    paste0(
      "brain_metastasis_samples=",
      summary_dt$n_brain_metastasis
    ),
    paste0(
      "primary_tumor_samples=",
      summary_dt$n_primary_tumor
    ),
    "expression_matrix_loaded=FALSE",
    "statistical_testing_performed=FALSE"
  ),
  file.path(
    OUT_DIR,
    "STEP_C1_COMPLETE.txt"
  )
)

cat(
  "\n========================================\n",
  "STEP C1 complete\n",
  "========================================\n",
  "Status: ",
  status,
  "\n\n",
  "Key output files:\n",
  file.path(
    OUT_DIR,
    "STEP_C1_01_GSE223501_spatial_file_inventory.tsv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "STEP_C1_02_GSE223501_spatial_inventory_summary.tsv"
  ),
  "\n",
  sep = ""
)
