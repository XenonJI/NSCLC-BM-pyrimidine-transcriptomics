# KRAS_17 Matrix Market reconstruction
# Reconstructs target-gene counts for KRAS_17 from the official Matrix Market files using
# the 1,000 exact unique barcode matches and updates the checkpoint used by script 01.
# No biological endpoint or gene set is changed.

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

ROOT <- "D:/LUAD_LM_PYRIMIDINE"

RAW_DIR <- file.path(
  ROOT,
  "data/raw/GSE223499"
)

DATA_DIR <- file.path(
  RAW_DIR,
  "counts_csv"
)

META_FILE <- file.path(
  RAW_DIR,
  "GSE223499_sn_integrated_data.csv.gz"
)

BARCODES_FILE <- file.path(
  DATA_DIR,
  "GSM6957628_NSCLC_KRAS_17_sn_barcodes.tsv.gz"
)

FEATURES_FILE <- file.path(
  DATA_DIR,
  "GSM6957628_NSCLC_KRAS_17_sn_features.tsv.gz"
)

MATRIX_FILE <- file.path(
  DATA_DIR,
  "GSM6957628_NSCLC_KRAS_17_sn_matrix.mtx.gz"
)

PROCESSED_DIR <- file.path(
  ROOT,
  "data/processed"
)

CHECKPOINT_FILE <- file.path(
  PROCESSED_DIR,
  "GSE223499_STEP07C_sample_gene_counts_checkpoint.rds"
)

TABLE_DIR <- file.path(
  ROOT,
  "results/tables/GSE223499"
)

LOCK_DIR <- file.path(
  ROOT,
  "results/locks"
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  LOCK_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("缺少data.table包。")
}

suppressPackageStartupMessages({
  library(data.table)
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

banner(
  "阶段0：检查MTX三件套和已有检查点"
)

required_files <- c(
  META_FILE,
  BARCODES_FILE,
  FEATURES_FILE,
  MATRIX_FILE,
  CHECKPOINT_FILE
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files) > 0L) {
  stop(
    "缺少文件：\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

gzip_files <- c(
  META_FILE,
  BARCODES_FILE,
  FEATURES_FILE,
  MATRIX_FILE
)

invalid_gzip <- gzip_files[
  !vapply(
    gzip_files,
    is_gzip_file,
    logical(1)
  )
]

if (length(invalid_gzip) > 0L) {
  stop(
    "以下文件不是有效gzip：\n",
    paste(
      invalid_gzip,
      collapse = "\n"
    )
  )
}

checkpoint <- readRDS(
  CHECKPOINT_FILE
)

if (
  !is.list(checkpoint) ||
    !identical(
      checkpoint$version,
      "STEP07C_v1"
    ) ||
    !is.list(
      checkpoint$processed
    )
) {
  stop(
    "STEP07-C检查点格式不正确。"
  )
}

cat(
  "现有检查点样本数：",
  length(
    checkpoint$processed
  ),
  "。\n",
  sep = ""
)

if (
  "KRAS_17" %in%
    names(
      checkpoint$processed
    )
) {
  stop(
    "检查点中已经存在KRAS_17。为避免覆盖，请先核查。"
  )
}

banner(
  "阶段1：读取KRAS_17整合注释"
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
    "整合注释缺少字段：",
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

kras17_meta <- metadata[
  sample == "KRAS_17" &
    site_author == "BRAIN_METS"
]

if (nrow(kras17_meta) != 1504L) {
  stop(
    "整合注释中的KRAS_17不是1504个细胞，而是：",
    nrow(kras17_meta)
  )
}

kras17_meta[
  ,
  raw_barcode := sub(
    "^KRAS_17_",
    "",
    cell_id
  )
]

kras17_meta[
  ,
  canonical_barcode :=
    canonical_barcode(
      raw_barcode
    )
]

kras17_meta[
  ,
  selected_tumor :=
    tumor_major == "Tumor" &
      !is_doublet
]

if (
  kras17_meta[
    selected_tumor == TRUE,
    .N
  ] != 1040L
) {
  stop(
    "KRAS_17固定作者肿瘤细胞不是1040个。"
  )
}

if (
  anyDuplicated(
    kras17_meta$canonical_barcode
  )
) {
  stop(
    "KRAS_17整合注释barcode存在重复。"
  )
}

cat(
  "KRAS_17总细胞：1504；固定作者肿瘤细胞：1040。\n"
)

banner(
  "阶段2：读取官方barcodes和features"
)

barcode_connection <- gzfile(
  BARCODES_FILE,
  open = "rt"
)

official_barcodes <- tryCatch(
  readLines(
    barcode_connection,
    warn = FALSE
  ),
  finally = {
    close(
      barcode_connection
    )
  }
)

official_barcodes <-
  canonical_barcode(
    official_barcodes
  )

official_barcode_n <- length(
  official_barcodes
)

if (official_barcode_n < 1504L) {
  stop(
    "官方KRAS_17 barcode数少于整合注释的1504个：",
    official_barcode_n
  )
}

if (
  anyDuplicated(
    official_barcodes
  )
) {
  stop(
    "官方KRAS_17 barcode存在重复。"
  )
}

all_match <- sum(
  kras17_meta$canonical_barcode %in%
    official_barcodes
)

metadata_unmatched <- kras17_meta[
  !canonical_barcode %in%
    official_barcodes
]

metadata_matched <- kras17_meta[
  canonical_barcode %in%
    official_barcodes
]

selected_meta_all <- kras17_meta[
  selected_tumor == TRUE
]

selected_meta <- selected_meta_all[
  canonical_barcode %in%
    official_barcodes
]

unmatched_selected_meta <- selected_meta_all[
  !canonical_barcode %in%
    official_barcodes
]

matched_tumor_n <- nrow(
  selected_meta
)

metadata_tumor_n <- nrow(
  selected_meta_all
)

unmatched_tumor_n <- nrow(
  unmatched_selected_meta
)

if (
  all_match != 1464L ||
    nrow(
      metadata_unmatched
    ) != 40L
) {
  stop(
    "KRAS_17当前barcode结构与已审计结果不一致：匹配",
    all_match,
    "/1504；未匹配",
    nrow(
      metadata_unmatched
    ),
    "。请停止并复核输入文件。"
  )
}

if (
  matched_tumor_n != 1000L ||
    metadata_tumor_n != 1040L ||
    unmatched_tumor_n != 40L
) {
  stop(
    "KRAS_17肿瘤细胞覆盖与已审计结果不一致：",
    matched_tumor_n,
    "/",
    metadata_tumor_n,
    "；未匹配",
    unmatched_tumor_n,
    "。"
  )
}

if (
  anyDuplicated(
    metadata_matched$canonical_barcode
  ) ||
    anyDuplicated(
      official_barcodes
    )
) {
  stop(
    "精确匹配集合存在重复barcode，不能继续。"
  )
}

fwrite(
  metadata_unmatched[
    ,
    .(
      cell_id,
      raw_barcode,
      canonical_barcode,
      selected_tumor,
      nCount_RNA
    )
  ],
  file.path(
    TABLE_DIR,
    "GSE223499_KRAS17_MTX_unmatched_40_cells_final_audit.tsv"
  ),
  sep = "\t"
)

extra_official_barcodes <- setdiff(
  official_barcodes,
  kras17_meta$canonical_barcode
)

if (length(extra_official_barcodes) > 0L) {
  fwrite(
    data.table(
      official_extra_barcode =
        extra_official_barcodes
    ),
    file.path(
      TABLE_DIR,
      "GSE223499_KRAS17_MTX_extra_official_barcodes.tsv"
    ),
    sep = "\t"
  )
}

selected_columns <- match(
  selected_meta$canonical_barcode,
  official_barcodes
)

if (
  length(
    selected_columns
  ) != 1000L ||
    any(
      is.na(
        selected_columns
      )
    ) ||
    anyDuplicated(
      selected_columns
    )
) {
  stop(
    "无法建立1000个精确匹配肿瘤细胞的唯一MTX列映射。"
  )
}

cat(
  "官方MTX barcode：",
  official_barcode_n,
  "；整合注释精确匹配：",
  all_match,
  "/1504；肿瘤细胞精确匹配：",
  matched_tumor_n,
  "/",
  metadata_tumor_n,
  "；排除未匹配肿瘤细胞：",
  unmatched_tumor_n,
  "。\n",
  sep = ""
)

features_connection <- gzfile(
  FEATURES_FILE,
  open = "rt"
)

features <- tryCatch(
  read.delim(
    features_connection,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  finally = {
    close(
      features_connection
    )
  }
)

if (ncol(features) < 2L) {
  stop(
    "features文件少于两列。"
  )
}

feature_symbols <- toupper(
  trimws(
    as.character(
      features[[2]]
    )
  )
)

target_feature_indices <- which(
  feature_symbols %in%
    target_genes
)

if (length(target_feature_indices) < 1L) {
  stop(
    "features中没有找到固定基因。"
  )
}

feature_gene_lookup <- setNames(
  feature_symbols[
    target_feature_indices
  ],
  as.character(
    target_feature_indices
  )
)

gene_rows_found <- table(
  factor(
    feature_symbols[
      target_feature_indices
    ],
    levels =
      target_genes
  )
)

gene_rows_found <- setNames(
  as.integer(
    gene_rows_found
  ),
  target_genes
)

cat(
  "features总行数：",
  nrow(features),
  "；固定目标feature行：",
  length(
    target_feature_indices
  ),
  "；MTX列映射仅使用1000个精确匹配的作者肿瘤细胞。\n",
  sep = ""
)

banner(
  "阶段3：流式扫描Matrix Market坐标"
)

matrix_connection <- gzfile(
  MATRIX_FILE,
  open = "rt"
)

first_line <- readLines(
  matrix_connection,
  n = 1L,
  warn = FALSE
)

if (
  length(first_line) != 1L ||
    !grepl(
      "^%%MatrixMarket",
      first_line
    )
) {
  close(matrix_connection)

  stop(
    "MTX缺少MatrixMarket表头。"
  )
}

repeat {
  dimension_line <- readLines(
    matrix_connection,
    n = 1L,
    warn = FALSE
  )

  if (length(dimension_line) == 0L) {
    close(matrix_connection)

    stop(
      "MTX没有维度行。"
    )
  }

  if (
    !startsWith(
      dimension_line,
      "%"
    )
  ) {
    break
  }
}

dimensions <- as.integer(
  strsplit(
    trimws(
      dimension_line
    ),
    "\\s+"
  )[[1]]
)

if (length(dimensions) != 3L) {
  close(matrix_connection)

  stop(
    "MTX维度行不是三个整数。"
  )
}

feature_n <- dimensions[1]
barcode_n <- dimensions[2]
nonzero_n <- dimensions[3]

if (
  feature_n != nrow(features) ||
    barcode_n != length(
      official_barcodes
    )
) {
  close(matrix_connection)

  stop(
    "MTX维度与features/barcodes不一致：",
    paste(
      dimensions,
      collapse = " x "
    )
  )
}

cat(
  "MTX维度：",
  feature_n,
  " features × ",
  barcode_n,
  " cells；非零项 ",
  format(
    nonzero_n,
    big.mark = ","
  ),
  "。\n",
  sep = ""
)

target_feature_set <- as.integer(
  target_feature_indices
)

selected_column_set <- as.integer(
  selected_columns
)

gene_counts <- setNames(
  rep(
    0,
    length(
      target_genes
    )
  ),
  target_genes
)

entries_read <- 0L
target_entries_used <- 0L

repeat {
  current_lines <- readLines(
    matrix_connection,
    n = 100000L,
    warn = FALSE
  )

  if (length(current_lines) == 0L) {
    break
  }

  block <- fread(
    text =
      paste(
        current_lines,
        collapse = "\n"
      ),
    header = FALSE,
    col.names = c(
      "feature_index",
      "cell_index",
      "count"
    ),
    showProgress = FALSE
  )

  entries_read <- entries_read +
    nrow(block)

  hit <- block[
    feature_index %in%
      target_feature_set &
      cell_index %in%
        selected_column_set
  ]

  if (nrow(hit) > 0L) {
    hit[
      ,
      gene :=
        feature_gene_lookup[
          as.character(
            feature_index
          )
        ]
    ]

    hit_sum <- hit[
      ,
      .(
        raw_count =
          sum(
            as.numeric(
              count
            )
          )
      ),
      by = gene
    ]

    gene_counts[
      hit_sum$gene
    ] <- gene_counts[
      hit_sum$gene
    ] +
      hit_sum$raw_count

    target_entries_used <-
      target_entries_used +
        nrow(hit)
  }

  rm(
    block,
    hit
  )

  gc(
    verbose = FALSE
  )
}

close(matrix_connection)

if (entries_read != nonzero_n) {
  stop(
    "实际读取的MTX非零项数与维度行不一致：",
    entries_read,
    " vs ",
    nonzero_n
  )
}

cat(
  "MTX扫描完成；使用固定目标非零项：",
  format(
    target_entries_used,
    big.mark = ","
  ),
  "。\n",
  sep = ""
)

module_detected <- module_genes[
  gene_rows_found[
    module_genes
  ] >
    0L
]

if (length(module_detected) != 5L) {
  stop(
    "固定5基因没有全部出现在features中：",
    paste(
      module_detected,
      collapse = ", "
    )
  )
}

banner(
  "阶段4：写入KRAS_17检查点"
)

patient_value <- unique(
  kras17_meta$patient
)

patient_value <- patient_value[
  !is.na(
    patient_value
  ) &
    nzchar(
      patient_value
    )
]

if (length(patient_value) != 1L) {
  stop(
    "KRAS_17作者patient字段不是唯一值。"
  )
}

selected_library_size <- sum(
  as.numeric(
    selected_meta$nCount_RNA
  ),
  na.rm = TRUE
)

kras17_result <- list(
  sample =
    "KRAS_17",
  patient =
    patient_value,
  cohort =
    "Brain_metastasis",
  batch =
    "KRAS",
  author_tumor_nuclei =
    1000L,
  metadata_author_tumor_nuclei =
    1040L,
  unmatched_author_tumor_nuclei =
    40L,
  selected_library_size =
    selected_library_size,
  primary_eligible =
    TRUE,
  header_cell_columns =
    official_barcode_n,
  matched_tumor_barcodes =
    1000L,
  barcode_match_fraction =
    1000 /
      1040,
  barcode_match_strategy =
    "official_MTX_exact_subset_1000_of_1040",
  full_suffix_match_fraction =
    1000 /
      1040,
  scanned_gene_rows =
    feature_n,
  gene_counts =
    gene_counts,
  gene_rows_found =
    gene_rows_found,
  expression_source =
    "GSM6957628_matrix_mtx"
)

timestamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

backup_file <- paste0(
  CHECKPOINT_FILE,
  ".before_KRAS17_MTX_",
  timestamp,
  ".rds"
)

if (
  !file.copy(
    CHECKPOINT_FILE,
    backup_file,
    overwrite = FALSE
  )
) {
  stop(
    "无法备份原检查点。"
  )
}

checkpoint$processed[["KRAS_17"]] <-
  kras17_result

saveRDS(
  checkpoint,
  CHECKPOINT_FILE
)

verification <- readRDS(
  CHECKPOINT_FILE
)

if (
  !"KRAS_17" %in%
    names(
      verification$processed
    ) ||
    verification$processed[["KRAS_17"]]$matched_tumor_barcodes !=
      1000L
) {
  stop(
    "写入后的检查点复核失败。原检查点备份位于：",
    backup_file
  )
}

audit <- data.table(
  sample =
    "KRAS_17",
  expression_source =
    "official_MTX_exact_subset",
  metadata_cells =
    1504L,
  official_MTX_barcodes =
    official_barcode_n,
  exact_metadata_matches =
    all_match,
  metadata_match_fraction =
    all_match /
      1504,
  metadata_author_tumor_cells =
    metadata_tumor_n,
  exact_matched_author_tumor_cells =
    matched_tumor_n,
  unmatched_author_tumor_cells =
    unmatched_tumor_n,
  tumor_match_fraction =
    matched_tumor_n /
      metadata_tumor_n,
  selected_MTX_columns =
    length(
      selected_columns
    ),
  official_extra_barcodes =
    length(
      extra_official_barcodes
    ),
  features =
    feature_n,
  MTX_nonzero_entries =
    nonzero_n,
  target_nonzero_entries =
    target_entries_used,
  matched_cells_library_size =
    selected_library_size,
  fixed_module_genes_detected =
    paste(
      module_detected,
      collapse = ","
    ),
  checkpoint_samples_after_patch =
    length(
      verification$processed
    ),
  checkpoint_backup =
    backup_file
)

fwrite(
  audit,
  file.path(
    TABLE_DIR,
    "GSE223499_STEP07C_KRAS17_MTX_checkpoint_patch_audit.tsv"
  ),
  sep = "\t"
)

amendment_file <- file.path(
  LOCK_DIR,
  "STEP_07C_GSE223499_KRAS17_MTX_EXACT_SUBSET_AMENDMENT.txt"
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
    "reason=deposited_integrated_metadata_and_official_KRAS17_MTX_have_a_40_tumor_cell_version_difference",
    "decision_time=before_reading_KRAS17_expression_matrix_entries",
    "KRAS17_expression_source=official_GSM6957628_matrix_mtx_barcodes_features",
    "KRAS17_metadata_cells=1504",
    paste0(
      "KRAS17_official_MTX_barcodes=",
      official_barcode_n
    ),
    "KRAS17_exact_metadata_barcode_match=1464/1504",
    "KRAS17_metadata_author_tumor_cells=1040",
    "KRAS17_exact_matched_author_tumor_cells=1000/1040",
    "KRAS17_unmatched_author_tumor_cells=40",
    "KRAS17_rule=use_only_exact_unique_barcode_matches_for_both_gene_counts_and_nCount_RNA_denominator",
    "KRAS17_unmatched_cells_used=FALSE",
    "KRAS17_official_extra_cells_used=FALSE",
    "KRAS17_patient_retained=TRUE",
    "KRAS17_leave_one_patient_out_sensitivity_required=TRUE",
    "module_genes_changed=FALSE",
    "cell_cycle_genes_changed=FALSE",
    "statistical_gate_changed=FALSE",
    "outcome_values_inspected_before_amendment=FALSE",
    "checkpoint_backup_created=TRUE",
    paste0(
      "checkpoint_backup=",
      backup_file
    )
  ),
  amendment_file
)

banner(
  "KRAS_17精确匹配1000细胞检查点补丁完成"
)

print(
  audit
)

cat(
  "\n现在重新运行：\n",
  'source("D:/LUAD_LM_PYRIMIDINE/scripts/run_07C_GSE223499_lowmem_5gene_validation_v2_resume.R")',
  "\n",
  sep = ""
)
