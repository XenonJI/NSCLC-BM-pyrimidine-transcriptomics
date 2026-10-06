# Cell-state composition and G0/G1-like sensitivity
# Reconstructs marker-defined cycling states from fixed non-overlapping S/G2M marker sets,
# quantifies the cycling-like fraction, and evaluates the fixed five-gene score within the
# G0/G1-like malignant-nucleus compartment.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/LUAD_LM_PYRIMIDINE"
OUT <- file.path(ROOT, "results", "major_revision")
CACHE <- file.path(OUT, "E2B_target_gene_cache")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

META <- file.path(
  ROOT,
  "data",
  "raw",
  "GSE223499",
  "GSE223499_sn_integrated_data.csv.gz"
)

COUNTS_DIR <- file.path(
  ROOT,
  "data",
  "raw",
  "GSE223499",
  "counts_csv"
)

ZREF_FILE <- file.path(
  OUT,
  "STEP_E1_05_fixed_zscore_reference_parameters.tsv"
)

if (!file.exists(META)) stop("Missing metadata: ", META)
if (!dir.exists(COUNTS_DIR)) stop("Missing counts directory: ", COUNTS_DIR)
if (!file.exists(ZREF_FILE)) {
  stop(
    "Missing E1 fixed z-reference:\n",
    ZREF_FILE,
    "\nRun STEP E1 first."
  )
}

# ----------------------------------------------------------------------
# 1. Fixed gene sets
# ----------------------------------------------------------------------
module_genes <- c(
  "DHFR",
  "DHODH",
  "SHMT1",
  "TYMS",
  "UMPS"
)

s_genes_original <- c(
  "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2",
  "MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2",
  "RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7",
  "POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1",
  "TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1",
  "CHAF1B","BRIP1","E2F8"
)

g2m_genes_original <- c(
  "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80",
  "CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A",
  "SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E",
  "TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK",
  "CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8",
  "ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5",
  "CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"
)

# Non-overlap is fixed before this post hoc analysis.
s_genes <- setdiff(
  s_genes_original,
  module_genes
)

g2m_genes <- setdiff(
  g2m_genes_original,
  module_genes
)

if (!identical(
  intersect(s_genes_original, module_genes),
  "TYMS"
)) {
  stop("Unexpected S-phase/module overlap.")
}

if (length(
  intersect(g2m_genes_original, module_genes)
) != 0L) {
  stop("Unexpected G2M/module overlap.")
}

cycle_genes <- unique(c(s_genes, g2m_genes))
target_genes <- unique(c(module_genes, cycle_genes))

# ----------------------------------------------------------------------
# 2. Helpers
# ----------------------------------------------------------------------
strip_quotes <- function(x) {
  x <- sub('^"', "", x)
  x <- sub('"$', "", x)
  x
}

clean_barcode <- function(x) {
  trimws(
    strip_quotes(
      as.character(x)
    )
  )
}

barcode_core16 <- function(x) {
  z <- toupper(
    clean_barcode(
      x
    )
  )

  m <- regexpr(
    "[ACGT]{16}",
    z,
    perl = TRUE
  )

  out <- rep(
    NA_character_,
    length(z)
  )

  ok <- m > 0L

  out[ok] <- substring(
    z[ok],
    m[ok],
    m[ok] + 15L
  )

  out
}

parse_sample_from_count_filename <- function(x) {
  z <- basename(x)
  z <- sub("^GSM[0-9]+_NSCLC_", "", z, ignore.case = TRUE)
  z <- sub("_sn_counts\\.csv\\.gz$", "", z, ignore.case = TRUE)
  z
}

hedges_g_ci <- function(x_bm, x_pt) {
  x_bm <- x_bm[is.finite(x_bm)]
  x_pt <- x_pt[is.finite(x_pt)]

  n1 <- length(x_bm)
  n0 <- length(x_pt)
  df <- n1 + n0 - 2L

  sp <- sqrt(
    (
      (n1 - 1L) * var(x_bm) +
      (n0 - 1L) * var(x_pt)
    ) / df
  )

  d <- (mean(x_bm) - mean(x_pt)) / sp
  J <- 1 - 3 / (4 * df - 1)
  g <- J * d

  var_d <- (
    (n1 + n0) / (n1 * n0)
  ) + (
    d^2 / (2 * df)
  )

  se_g <- sqrt(J^2 * var_d)

  data.table(
    hedges_g = g,
    hedges_g_se = se_g,
    hedges_g_95CI_low = g - qnorm(0.975) * se_g,
    hedges_g_95CI_high = g + qnorm(0.975) * se_g
  )
}

rank_biserial_two_sample <- function(x_bm, x_pt) {
  n1 <- length(x_bm)
  n0 <- length(x_pt)
  ranks <- rank(c(x_bm, x_pt), ties.method = "average")
  R1 <- sum(ranks[seq_len(n1)])
  U1 <- R1 - n1 * (n1 + 1) / 2
  2 * U1 / (n1 * n0) - 1
}

summarize_two_group <- function(
  dt,
  score_col,
  analysis_name
) {
  bm <- dt[
    cohort == "Brain_metastasis",
    get(score_col)
  ]
  pt <- dt[
    cohort == "Primary_tumor",
    get(score_col)
  ]

  bm <- bm[is.finite(bm)]
  pt <- pt[is.finite(pt)]

  if (length(bm) < 2L || length(pt) < 2L) {
    stop("Too few observations for: ", analysis_name)
  }

  w <- suppressWarnings(
    wilcox.test(
      bm,
      pt,
      alternative = "two.sided",
      exact = FALSE,
      correct = FALSE,
      conf.int = TRUE,
      conf.level = 0.95
    )
  )

  g <- hedges_g_ci(bm, pt)

  data.table(
    analysis = analysis_name,
    outcome = score_col,
    n_BM = length(bm),
    n_PT = length(pt),
    BM_median = median(bm),
    PT_median = median(pt),
    median_difference = median(bm) - median(pt),
    two_sided_wilcoxon_p = unname(w$p.value),
    hodges_lehmann_shift = if (!is.null(w$estimate)) {
      unname(w$estimate)
    } else {
      NA_real_
    },
    hodges_lehmann_95CI_low = if (!is.null(w$conf.int)) {
      unname(w$conf.int[1])
    } else {
      NA_real_
    },
    hodges_lehmann_95CI_high = if (!is.null(w$conf.int)) {
      unname(w$conf.int[2])
    } else {
      NA_real_
    },
    rank_biserial = rank_biserial_two_sample(bm, pt),
    hedges_g = g$hedges_g,
    hedges_g_95CI_low = g$hedges_g_95CI_low,
    hedges_g_95CI_high = g$hedges_g_95CI_high
  )
}

cluster_robust_cr1 <- function(
  fit,
  cluster,
  term_name
) {
  X <- model.matrix(fit)
  e <- residuals(fit)
  cluster <- as.character(cluster)

  if (nrow(X) != length(cluster)) {
    stop("Cluster vector length does not match model rows.")
  }

  keep <- complete.cases(X, e, cluster)
  X <- X[keep, , drop = FALSE]
  e <- e[keep]
  cluster <- cluster[keep]

  N <- nrow(X)
  K <- ncol(X)
  clusters <- unique(cluster)
  G <- length(clusters)

  if (G <= K) stop("Too few clusters.")

  bread <- solve(crossprod(X))
  meat <- matrix(0, nrow = K, ncol = K)

  for (cl in clusters) {
    idx <- cluster == cl
    Xg <- X[idx, , drop = FALSE]
    eg <- e[idx]
    sg <- crossprod(Xg, eg)
    meat <- meat + sg %*% t(sg)
  }

  correction <- (G / (G - 1)) * ((N - 1) / (N - K))
  vc <- correction * bread %*% meat %*% bread

  b <- coef(fit)
  se <- sqrt(diag(vc))

  est <- unname(b[term_name])
  se0 <- unname(se[term_name])
  df <- G - 1
  t0 <- est / se0
  p2 <- 2 * pt(-abs(t0), df = df)
  crit <- qt(0.975, df = df)

  data.table(
    term = term_name,
    estimate = est,
    CR1_SE = se0,
    df_clusters = df,
    t_CR1 = t0,
    two_sided_p_CR1 = p2,
    CI95_low_CR1 = est - crit * se0,
    CI95_high_CR1 = est + crit * se0,
    n_observations = N,
    n_patient_clusters = G
  )
}

extract_target_rows_from_csv <- function(
  path,
  target_genes,
  selected_columns,
  expected_n_columns
) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)

  header <- readLines(con, n = 1L, warn = FALSE)
  if (length(header) != 1L) {
    stop("Cannot read count header: ", path)
  }

  found <- vector("list", length(target_genes))
  names(found) <- target_genes

  repeat {
    lines <- readLines(con, n = 2000L, warn = FALSE)
    if (length(lines) == 0L) break

    comma_pos <- regexpr(",", lines, fixed = TRUE)
    valid <- comma_pos > 1L

    gene_raw <- rep(NA_character_, length(lines))
    gene_raw[valid] <- substring(
      lines[valid],
      1L,
      comma_pos[valid] - 1L
    )
    gene_raw <- toupper(strip_quotes(gene_raw))

    hit_idx <- which(
      !is.na(gene_raw) &
      gene_raw %chin% target_genes
    )

    if (length(hit_idx) == 0L) next

    for (ii in hit_idx) {
      gg <- gene_raw[ii]
      fields <- strsplit(
        lines[ii],
        ",",
        fixed = TRUE
      )[[1]]

      if (length(fields) != expected_n_columns + 1L) {
        stop(
          "Unexpected field count for gene ",
          gg,
          " in ",
          basename(path),
          ": expected ",
          expected_n_columns + 1L,
          " but observed ",
          length(fields)
        )
      }

      vals <- suppressWarnings(
        as.numeric(
          strip_quotes(fields[-1L])
        )
      )

      if (anyNA(vals)) {
        stop(
          "Non-numeric count encountered for ",
          gg,
          " in ",
          basename(path)
        )
      }

      vsel <- vals[selected_columns]

      if (is.null(found[[gg]])) {
        found[[gg]] <- vsel
      } else {
        found[[gg]] <- found[[gg]] + vsel
      }
    }
  }

  nsel <- length(selected_columns)

  mat <- matrix(
    0,
    nrow = length(target_genes),
    ncol = nsel,
    dimnames = list(target_genes, NULL)
  )

  for (gg in target_genes) {
    if (!is.null(found[[gg]])) {
      mat[gg, ] <- found[[gg]]
    }
  }

  mat
}

# ----------------------------------------------------------------------
# 3. True GSE223499 metadata
# ----------------------------------------------------------------------
md <- fread(
  META,
  showProgress = TRUE
)

required_md <- c(
  "V1",
  "orig.ident",
  "patient",
  "PRIMARY vs BRAIN_METS vs CHEST_WALL_MET",
  "tumor_nontumor_major",
  "predicted_doublets",
  "nCount_RNA"
)

missing_md <- setdiff(required_md, names(md))
if (length(missing_md) > 0L) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_md, collapse = ", ")
  )
}

setnames(
  md,
  c(
    "V1",
    "orig.ident",
    "PRIMARY vs BRAIN_METS vs CHEST_WALL_MET"
  ),
  c(
    "barcode",
    "sample",
    "site_author"
  )
)

md[
  ,
  barcode := clean_barcode(barcode)
]

md[
  ,
  site_author := as.character(site_author)
]

md[
  ,
  cohort := fifelse(
    site_author == "BRAIN_METS",
    "Brain_metastasis",
    fifelse(
      site_author == "PRIMARY",
      "Primary_tumor",
      NA_character_
    )
  )
]

md[
  ,
  predicted_doublets_logical := tolower(
    as.character(predicted_doublets)
  ) %chin% c("true", "t", "1", "yes")
]

eligible <- md[
  !is.na(cohort) &
  tumor_nontumor_major == "Tumor" &
  predicted_doublets_logical == FALSE
]

eligible[
  ,
  corrected_patient_id := sample
]

eligible[
  sample %chin% c("PA060", "N254"),
  corrected_patient_id := "MATCHED_PATIENT_PA060_N254"
]

eligible[
  sample %chin% c("PA067", "N586"),
  corrected_patient_id := "MATCHED_PATIENT_PA067_N586"
]

# Locked >=50 specimen set; STK_20 is therefore excluded.
sample_n <- eligible[
  ,
  .(
    malignant_nuclei = .N,
    cohort = cohort[1],
    corrected_patient_id = corrected_patient_id[1]
  ),
  by = sample
]

analysis_samples <- sample_n[
  malignant_nuclei >= 50L
]

if (
  sum(analysis_samples$cohort == "Brain_metastasis") != 31L ||
  sum(analysis_samples$cohort == "Primary_tumor") != 10L
) {
  stop(
    "Expected locked >=50 set of 31 BM / 10 PT but observed ",
    sum(analysis_samples$cohort == "Brain_metastasis"),
    " BM / ",
    sum(analysis_samples$cohort == "Primary_tumor"),
    " PT."
  )
}

# ----------------------------------------------------------------------
# 4. Count-file map + strict barcode harmonization preflight
# ----------------------------------------------------------------------
# The integrated metadata uses a global barcode namespace, whereas the
# per-sample CSV headers can use a sample-local suffix convention. Therefore:
#   1) exact cleaned barcode matching is attempted first;
#   2) if exact matching is incomplete, a suffix-insensitive 16-nt core
#      fallback is permitted ONLY when the 16-nt cores are unique within
#      both the malignant metadata subset and the raw count columns, and
#      100% of malignant metadata nuclei map to a raw count column.
# This mirrors the previously documented barcode-harmonization principle
# and prevents ambiguous many-to-one joins.

count_files <- list.files(
  COUNTS_DIR,
  pattern = "_sn_counts\\.csv\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

count_map <- data.table(
  sample = vapply(
    count_files,
    parse_sample_from_count_filename,
    character(1)
  ),
  count_file = count_files
)

if (anyDuplicated(count_map$sample)) {
  stop("Duplicate raw count file mapping detected.")
}

analysis_samples <- merge(
  analysis_samples,
  count_map,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(analysis_samples$count_file)) {
  stop(
    "Missing count file for: ",
    paste(
      analysis_samples[is.na(count_file), sample],
      collapse = ", "
    )
  )
}

join_audit <- vector(
  "list",
  nrow(analysis_samples)
)

for (i in seq_len(nrow(analysis_samples))) {
  ss <- analysis_samples$sample[i]
  ff <- analysis_samples$count_file[i]

  con <- gzfile(
    ff,
    open = "rt"
  )

  header <- readLines(
    con,
    n = 1L,
    warn = FALSE
  )

  close(
    con
  )

  fields <- strsplit(
    header,
    ",",
    fixed = TRUE
  )[[1]]

  count_barcodes <- clean_barcode(
    fields[-1L]
  )

  md_barcodes <- eligible[
    sample == ss,
    barcode
  ]

  exact_idx <- match(
    md_barcodes,
    count_barcodes
  )

  exact_n <- sum(
    !is.na(
      exact_idx
    )
  )

  exact_fraction <- exact_n /
    length(
      md_barcodes
    )

  md_core <- barcode_core16(
    md_barcodes
  )

  count_core <- barcode_core16(
    count_barcodes
  )

  md_core_valid <- all(
    !is.na(
      md_core
    )
  )

  count_core_valid <- all(
    !is.na(
      count_core
    )
  )

  md_core_unique <- md_core_valid &&
    !anyDuplicated(
      md_core
    )

  count_core_unique <- count_core_valid &&
    !anyDuplicated(
      count_core
    )

  core_idx <- if (
    md_core_unique &&
    count_core_unique
  ) {
    match(
      md_core,
      count_core
    )
  } else {
    rep(
      NA_integer_,
      length(
        md_core
      )
    )
  }

  core_n <- sum(
    !is.na(
      core_idx
    )
  )

  core_fraction <- core_n /
    length(
      md_barcodes
    )

  join_method <- if (
    exact_fraction ==
      1
  ) {
    "exact_cleaned_barcode"
  } else if (
    md_core_unique &&
    count_core_unique &&
    core_fraction ==
      1
  ) {
    "unique_core16_100pct_fallback"
  } else {
    "unresolved"
  }

  selected_match_n <- if (
    join_method ==
      "exact_cleaned_barcode"
  ) {
    exact_n
  } else if (
    join_method ==
      "unique_core16_100pct_fallback"
  ) {
    core_n
  } else {
    max(
      exact_n,
      core_n
    )
  }

  selected_fraction <- selected_match_n /
    length(
      md_barcodes
    )

  join_audit[[i]] <- data.table(
    sample = ss,
    cohort = analysis_samples$cohort[i],
    malignant_nuclei_metadata = length(
      md_barcodes
    ),
    count_columns = length(
      count_barcodes
    ),
    exact_cleaned_matches = exact_n,
    exact_cleaned_fraction = exact_fraction,
    metadata_core16_unique = md_core_unique,
    count_core16_unique = count_core_unique,
    core16_matches = core_n,
    core16_match_fraction = core_fraction,
    join_method = join_method,
    selected_match_n = selected_match_n,
    selected_match_fraction = selected_fraction,
    example_metadata_barcode = md_barcodes[1],
    example_metadata_core16 = md_core[1],
    example_count_barcode = count_barcodes[1],
    example_count_core16 = count_core[1],
    count_file = ff
  )
}

join_audit <- rbindlist(
  join_audit
)

# KRAS_17 remains excluded a priori from E2B because its deposited CSV
# discordance was documented before this post hoc analysis.
join_audit[
  ,
  E2B_included :=
    sample != "KRAS_17" &
    join_method != "unresolved" &
    selected_match_fraction == 1
]

join_audit[
  ,
  exclusion_reason := fifelse(
    sample == "KRAS_17",
    "Documented pre-expression metadata/count-file discordance; E2B excludes the deposited CSV",
    fifelse(
      join_method == "unresolved",
      "Neither exact barcode nor unique 16-nt core achieved 100% malignant-nucleus mapping",
      ""
    )
  )
]

fwrite(
  join_audit,
  file.path(
    OUT,
    "STEP_E2B_01_cell_level_barcode_join_audit.tsv"
  ),
  sep = "\t"
)

if (
  any(
    join_audit[
      sample != "KRAS_17",
      join_method
    ] ==
      "unresolved"
  )
) {
  bad <- join_audit[
    sample != "KRAS_17" &
    join_method == "unresolved",
    sample
  ]

  stop(
    "Strict barcode harmonization remains unresolved for non-KRAS_17 specimen(s): ",
    paste(
      bad,
      collapse = ", "
    ),
    "\nInspect STEP_E2B_01_cell_level_barcode_join_audit.tsv."
  )
}

e2b_samples <- analysis_samples[
  sample %chin%
    join_audit[E2B_included == TRUE, sample]
]

e2b_samples <- merge(
  e2b_samples,
  join_audit[
    E2B_included == TRUE,
    .(
      sample,
      join_method
    )
  ],
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

# Expected: 30 BM + 10 PT after KRAS_17 technical exclusion.
if (
  sum(e2b_samples$cohort == "Brain_metastasis") != 30L ||
  sum(e2b_samples$cohort == "Primary_tumor") != 10L
) {
  stop(
    "Expected E2B set of 30 BM / 10 PT after KRAS_17 exclusion."
  )
}

# ----------------------------------------------------------------------
# 5. Stream only fixed target genes and cache compact matrices
# ----------------------------------------------------------------------
cycle_accum <- data.table(
  gene = cycle_genes,
  n = 0,
  sum_x = 0,
  sum_x2 = 0
)

availability <- vector("list", nrow(e2b_samples))

for (i in seq_len(nrow(e2b_samples))) {
  ss <- e2b_samples$sample[i]
  ff <- e2b_samples$count_file[i]

  cat(
    "[E2B extraction ",
    i,
    "/",
    nrow(e2b_samples),
    "] ",
    ss,
    "\n",
    sep = ""
  )

  md_ss <- eligible[
    sample == ss
  ]

  con <- gzfile(ff, open = "rt")
  header <- readLines(con, n = 1L, warn = FALSE)
  close(con)

  fields <- strsplit(
    header,
    ",",
    fixed = TRUE
  )[[1]]

  count_barcodes <- clean_barcode(
    fields[-1L]
  )

  method_ss <- e2b_samples[
    sample == ss,
    join_method
  ][1]

  if (
    method_ss ==
      "exact_cleaned_barcode"
  ) {
    idx <- match(
      md_ss$barcode,
      count_barcodes
    )
  } else if (
    method_ss ==
      "unique_core16_100pct_fallback"
  ) {
    md_core <- barcode_core16(
      md_ss$barcode
    )

    count_core <- barcode_core16(
      count_barcodes
    )

    if (
      anyNA(
        md_core
      ) ||
      anyNA(
        count_core
      ) ||
      anyDuplicated(
        md_core
      ) ||
      anyDuplicated(
        count_core
      )
    ) {
      stop(
        "Core16 uniqueness changed unexpectedly after preflight for ",
        ss
      )
    }

    idx <- match(
      md_core,
      count_core
    )
  } else {
    stop(
      "Unexpected join method for ",
      ss,
      ": ",
      method_ss
    )
  }

  if (
    anyNA(
      idx
    )
  ) {
    stop(
      "Unexpected unmatched malignant barcode after strict preflight for ",
      ss
    )
  }

  counts_target <- extract_target_rows_from_csv(
    path = ff,
    target_genes = target_genes,
    selected_columns = idx,
    expected_n_columns = length(count_barcodes)
  )

  lib <- as.numeric(md_ss$nCount_RNA)

  if (any(!is.finite(lib)) || any(lib <= 0)) {
    stop("Invalid nCount_RNA in ", ss)
  }

  lognorm <- log1p(
    sweep(
      counts_target,
      2,
      lib / 1e4,
      "/"
    )
  )

  # Accumulate global cell-level mean/SD reference for cycle markers.
  for (gg in cycle_genes) {
    x <- lognorm[gg, ]
    jj <- match(gg, cycle_accum$gene)

    cycle_accum$n[jj] <- cycle_accum$n[jj] + length(x)
    cycle_accum$sum_x[jj] <- cycle_accum$sum_x[jj] + sum(x)
    cycle_accum$sum_x2[jj] <- cycle_accum$sum_x2[jj] + sum(x^2)
  }

  availability[[i]] <- data.table(
    sample = ss,
    malignant_nuclei = ncol(counts_target),
    detected_module_genes = sum(rowSums(counts_target[module_genes, , drop = FALSE]) > 0),
    detected_S_genes = sum(rowSums(counts_target[s_genes, , drop = FALSE]) > 0),
    detected_G2M_genes = sum(rowSums(counts_target[g2m_genes, , drop = FALSE]) > 0)
  )

  saveRDS(
    list(
      sample = ss,
      cohort = e2b_samples$cohort[i],
      corrected_patient_id = e2b_samples$corrected_patient_id[i],
      barcode = md_ss$barcode,
      libsize = lib,
      counts = counts_target
    ),
    file.path(
      CACHE,
      paste0(ss, "_E2B_targets.rds")
    ),
    compress = TRUE
  )

  rm(
    counts_target,
    lognorm,
    md_ss,
    lib
  )
  gc(verbose = FALSE)
}

availability <- rbindlist(availability)

fwrite(
  availability,
  file.path(
    OUT,
    "STEP_E2B_02_target_gene_availability.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 6. Global fixed cell-level reference for non-overlapping S/G2M genes
# ----------------------------------------------------------------------
cycle_accum[
  ,
  mean_log1p10k := sum_x / n
]

cycle_accum[
  ,
  variance_log1p10k := pmax(
    (
      sum_x2 -
      (sum_x^2 / n)
    ) / pmax(n - 1, 1),
    0
  )
]

cycle_accum[
  ,
  sd_log1p10k := sqrt(variance_log1p10k)
]

cycle_accum[
  ,
  usable := is.finite(sd_log1p10k) &
    sd_log1p10k > 1e-8
]

fwrite(
  cycle_accum,
  file.path(
    OUT,
    "STEP_E2B_03_fixed_cellcycle_gene_reference.tsv"
  ),
  sep = "\t"
)

usable_s <- intersect(
  s_genes,
  cycle_accum[usable == TRUE, gene]
)

usable_g2m <- intersect(
  g2m_genes,
  cycle_accum[usable == TRUE, gene]
)

if (length(usable_s) < 20L) {
  stop(
    "Too few usable S-phase markers: ",
    length(usable_s)
  )
}

if (length(usable_g2m) < 20L) {
  stop(
    "Too few usable G2M markers: ",
    length(usable_g2m)
  )
}

# ----------------------------------------------------------------------
# 7. Fixed E1 module reference
# ----------------------------------------------------------------------
zref <- fread(
  ZREF_FILE,
  showProgress = FALSE
)

required_ref <- c(
  "gene",
  "reference_mean",
  "reference_sd"
)

if (!all(required_ref %in% names(zref))) {
  stop("E1 z-reference file has unexpected columns.")
}

zref[
  ,
  gene := toupper(as.character(gene))
]

if (!all(module_genes %chin% zref$gene)) {
  stop("E1 z-reference is missing one or more module genes.")
}

# ----------------------------------------------------------------------
# 8. Reconstruct marker-based phase state + G0/G1-like pseudobulk
# ----------------------------------------------------------------------
sample_results <- list()
gene_results <- list()

for (i in seq_len(nrow(e2b_samples))) {
  ss <- e2b_samples$sample[i]

  obj <- readRDS(
    file.path(
      CACHE,
      paste0(ss, "_E2B_targets.rds")
    )
  )

  counts_target <- obj$counts
  lib <- obj$libsize

  lognorm <- log1p(
    sweep(
      counts_target,
      2,
      lib / 1e4,
      "/"
    )
  )

  cycle_ref <- cycle_accum[
    match(rownames(lognorm), gene)
  ]

  zcycle <- matrix(
    NA_real_,
    nrow = nrow(lognorm),
    ncol = ncol(lognorm),
    dimnames = dimnames(lognorm)
  )

  cyc_rows <- which(
    rownames(lognorm) %chin%
      cycle_accum[usable == TRUE, gene]
  )

  for (rr in cyc_rows) {
    gg <- rownames(lognorm)[rr]
    refrow <- cycle_accum[gene == gg]

    zcycle[rr, ] <- (
      lognorm[rr, ] -
      refrow$mean_log1p10k
    ) /
      refrow$sd_log1p10k
  }

  S_score <- colMeans(
    zcycle[
      usable_s,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  G2M_score <- colMeans(
    zcycle[
      usable_g2m,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  phase_like <- fifelse(
    S_score <= 0 &
      G2M_score <= 0,
    "G0G1_like",
    fifelse(
      S_score >= G2M_score,
      "S_like",
      "G2M_like"
    )
  )

  g1 <- phase_like == "G0G1_like"

  n_total <- length(g1)
  n_g1 <- sum(g1)
  n_cycle <- n_total - n_g1

  sample_results[[i]] <- data.table(
    sample = ss,
    cohort = obj$cohort,
    corrected_patient_id = obj$corrected_patient_id,
    malignant_nuclei = n_total,
    G0G1_like_nuclei = n_g1,
    cycling_like_nuclei = n_cycle,
    G0G1_like_fraction = n_g1 / n_total,
    cycling_like_fraction = n_cycle / n_total,
    S_like_fraction = mean(phase_like == "S_like"),
    G2M_like_fraction = mean(phase_like == "G2M_like"),
    mean_S_score = mean(S_score),
    mean_G2M_score = mean(G2M_score)
  )

  # G0/G1-like pseudobulk only if >= 20 nuclei; >=50 is primary E2B
  # reporting threshold, while >=20 is retained as a descriptive sensitivity.
  if (n_g1 > 0L) {
    denom <- sum(lib[g1])

    for (gg in module_genes) {
      raw_sum <- sum(
        counts_target[
          gg,
          g1
        ]
      )

      expr <- log2(
        raw_sum /
          denom *
          1e6 +
          1
      )

      ref <- zref[gene == gg]

      gene_results[[length(gene_results) + 1L]] <- data.table(
        sample = ss,
        cohort = obj$cohort,
        corrected_patient_id = obj$corrected_patient_id,
        G0G1_like_nuclei = n_g1,
        gene = gg,
        raw_sum = raw_sum,
        G0G1_library_size = denom,
        log2_CPM_plus1 = expr,
        z_fixed_E1_reference = (
          expr -
          ref$reference_mean
        ) /
          ref$reference_sd
      )
    }
  }

  rm(
    obj,
    counts_target,
    lognorm,
    zcycle,
    S_score,
    G2M_score,
    phase_like,
    g1
  )
  gc(verbose = FALSE)
}

sample_results <- rbindlist(sample_results)
gene_results <- rbindlist(gene_results)

fwrite(
  sample_results,
  file.path(
    OUT,
    "STEP_E2B_04_specimen_cellcycle_composition.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_results,
  file.path(
    OUT,
    "STEP_E2B_05_G0G1_like_gene_pseudobulk.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 9. G0/G1-like module score, same fixed E1 gene reference
# ----------------------------------------------------------------------
g1_module <- gene_results[
  ,
  .(
    cohort = cohort[1],
    corrected_patient_id = corrected_patient_id[1],
    G0G1_like_nuclei = G0G1_like_nuclei[1],
    module_score_G0G1_fixed = mean(z_fixed_E1_reference)
  ),
  by = sample
]

fwrite(
  g1_module,
  file.path(
    OUT,
    "STEP_E2B_06_G0G1_like_module_score.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 10. Statistical reporting
# ----------------------------------------------------------------------
# A) Cycling composition
cycling_stats <- summarize_two_group(
  sample_results,
  score_col = "cycling_like_fraction",
  analysis_name = "Reviewer-requested malignant cycling-like fraction"
)

# B) G0/G1-like module, >=50 G0/G1-like nuclei
g1_50 <- g1_module[
  G0G1_like_nuclei >= 50L
]

if (
  sum(g1_50$cohort == "Brain_metastasis") < 2L ||
  sum(g1_50$cohort == "Primary_tumor") < 2L
) {
  stop(
    "Too few specimens with >=50 G0/G1-like malignant nuclei."
  )
}

g1_50_stats <- summarize_two_group(
  g1_50,
  score_col = "module_score_G0G1_fixed",
  analysis_name = "G0/G1-like malignant-cell pseudobulk, >=50 nuclei"
)

# C) G0/G1-like module, >=20 descriptive sensitivity
g1_20 <- g1_module[
  G0G1_like_nuclei >= 20L
]

g1_20_stats <- summarize_two_group(
  g1_20,
  score_col = "module_score_G0G1_fixed",
  analysis_name = "G0/G1-like malignant-cell pseudobulk, >=20 nuclei sensitivity"
)

# D) corrected patient-cluster CR1 for G0/G1-like >=50
g1_50[
  ,
  site_BM := as.integer(
    cohort == "Brain_metastasis"
  )
]

fit_g1_cluster <- lm(
  module_score_G0G1_fixed ~ site_BM,
  data = g1_50
)

g1_cluster <- cluster_robust_cr1(
  fit = fit_g1_cluster,
  cluster = g1_50$corrected_patient_id,
  term_name = "site_BM"
)

# E) independence-preserving subset: remove both matched patient pairs
pair_samples <- c(
  "PA060",
  "N254",
  "PA067",
  "N586"
)

g1_independent <- g1_50[
  !sample %chin% pair_samples
]

g1_independent_stats <- summarize_two_group(
  g1_independent,
  score_col = "module_score_G0G1_fixed",
  analysis_name = "G0/G1-like >=50, independence-preserving exclusion of both matched pairs"
)

all_module_stats <- rbindlist(
  list(
    cycling_stats,
    g1_50_stats,
    g1_20_stats,
    g1_independent_stats
  ),
  fill = TRUE
)

fwrite(
  all_module_stats,
  file.path(
    OUT,
    "STEP_E2B_07_cellcycle_and_G0G1_module_statistics.tsv"
  ),
  sep = "\t"
)

fwrite(
  g1_cluster,
  file.path(
    OUT,
    "STEP_E2B_08_G0G1_patient_cluster_CR1.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 11. G0/G1-like gene-level decomposition, >=50 nuclei
# ----------------------------------------------------------------------
g1_gene_50 <- gene_results[
  G0G1_like_nuclei >= 50L
]

gene_stats <- rbindlist(
  lapply(
    module_genes,
    function(gg) {
      d <- g1_gene_50[
        gene == gg
      ]

      out <- summarize_two_group(
        d,
        score_col = "log2_CPM_plus1",
        analysis_name = paste0(
          "G0/G1-like gene decomposition: ",
          gg
        )
      )

      out[
        ,
        gene := gg
      ]

      out
    }
  ),
  fill = TRUE
)

gene_stats[
  ,
  BH_FDR_across_5_genes := p.adjust(
    two_sided_wilcoxon_p,
    method = "BH"
  )
]

setcolorder(
  gene_stats,
  c(
    "gene",
    setdiff(names(gene_stats), "gene")
  )
)

fwrite(
  gene_stats,
  file.path(
    OUT,
    "STEP_E2B_09_G0G1_gene_level_statistics.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 12. reporting interpretation guard
# ----------------------------------------------------------------------
g1_p <- g1_50_stats$two_sided_wilcoxon_p
g1_eff <- g1_50_stats$median_difference
cycle_p <- cycling_stats$two_sided_wilcoxon_p
cycle_eff <- cycling_stats$median_difference

interpretation <- if (
  is.finite(g1_p) &&
  g1_p < 0.05 &&
  g1_eff > 0
) {
  paste0(
    "The five-gene signal remains detectable within the marker-defined ",
    "G0/G1-like malignant compartment. This argues against a purely ",
    "compositional explanation, although the analysis remains post hoc and ",
    "does not establish proliferation-independent metabolic activation."
  )
} else {
  paste0(
    "The five-gene site difference is not retained within the marker-defined ",
    "G0/G1-like malignant compartment. This is consistent with the main ",
    "cell-cycle-adjustment result and supports the interpretation that much ",
    "of the observed signal reflects proliferation-related cell-state ",
    "composition/transcriptional variation."
  )
}

writeLines(
  c(
    paste0("status=PASS_E2B_G0G1_SENSITIVITY_COMPLETE"),
    "script_version=E2B_V2",
    "analysis_type=POST_HOC_REVIEWER_REQUESTED_SENSITIVITY",
    paste0("exact_barcode_samples=", sum(join_audit$join_method == "exact_cleaned_barcode")),
    paste0("core16_fallback_samples=", sum(join_audit$join_method == "unique_core16_100pct_fallback")),
    "author_phase_available=FALSE",
    "phase_method=marker_based_G0G1_like_not_author_phase",
    "cellcycle_reference=global_z_of_log1p_1e4_normalized_nonoverlapping_Seurat2019_markers",
    "G0G1_like_rule=S_score_le_0_AND_G2M_score_le_0",
    "KRAS17_cell_level_E2B_included=FALSE",
    paste0(
      "E2B_specimens_BM=",
      sum(sample_results$cohort == "Brain_metastasis")
    ),
    paste0(
      "E2B_specimens_PT=",
      sum(sample_results$cohort == "Primary_tumor")
    ),
    paste0(
      "G0G1_ge50_BM=",
      sum(g1_50$cohort == "Brain_metastasis")
    ),
    paste0(
      "G0G1_ge50_PT=",
      sum(g1_50$cohort == "Primary_tumor")
    ),
    paste0(
      "cycling_fraction_two_sided_p=",
      signif(cycle_p, 10)
    ),
    paste0(
      "cycling_fraction_median_difference_BM_minus_PT=",
      signif(cycle_eff, 10)
    ),
    paste0(
      "G0G1_module_two_sided_p=",
      signif(g1_p, 10)
    ),
    paste0(
      "G0G1_module_median_difference_BM_minus_PT=",
      signif(g1_eff, 10)
    ),
    paste0(
      "reviewer_facing_interpretation=",
      interpretation
    ),
    "primary_endpoint_replaced=FALSE",
    "gene_set_changed=FALSE",
    "result_driven_exclusion=FALSE"
  ),
  file.path(
    OUT,
    "STEP_E2B_COMPLETE.txt"
  )
)

cat(
  "\n============================================================\n",
  "STEP E2B V2: CELL-CYCLE COMPOSITION + G0/G1-LIKE SENSITIVITY\n",
  "============================================================\n",
  "E2B technically concordant specimens: ",
  sum(sample_results$cohort == "Brain_metastasis"),
  " BM / ",
  sum(sample_results$cohort == "Primary_tumor"),
  " PT\n",
  "G0/G1-like >=50 specimens: ",
  sum(g1_50$cohort == "Brain_metastasis"),
  " BM / ",
  sum(g1_50$cohort == "Primary_tumor"),
  " PT\n\n",
  "Cycling-like fraction comparison:\n",
  sep = ""
)

print(cycling_stats)

cat(
  "\nG0/G1-like module comparison (>=50 nuclei):\n"
)

print(g1_50_stats)

cat(
  "\nCorrected patient-cluster CR1 model:\n"
)

print(g1_cluster)

cat(
  "\nIndependence-preserving G0/G1-like sensitivity:\n"
)

print(g1_independent_stats)

cat(
  "\nGene-level G0/G1-like decomposition:\n"
)

print(
  gene_stats[
    ,
    .(
      gene,
      median_difference,
      two_sided_wilcoxon_p,
      BH_FDR_across_5_genes,
      hedges_g,
      hedges_g_95CI_low,
      hedges_g_95CI_high
    )
  ]
)

cat(
  "\nReviewer-facing interpretation:\n",
  interpretation,
  "\n\nStatus: PASS_E2B_G0G1_SENSITIVITY_COMPLETE\n\n",
  "Key output files:\n",
  file.path(OUT, "STEP_E2B_01_cell_level_barcode_join_audit.tsv"), "\n",
  file.path(OUT, "STEP_E2B_04_specimen_cellcycle_composition.tsv"), "\n",
  file.path(OUT, "STEP_E2B_07_cellcycle_and_G0G1_module_statistics.tsv"), "\n",
  file.path(OUT, "STEP_E2B_08_G0G1_patient_cluster_CR1.tsv"), "\n",
  file.path(OUT, "STEP_E2B_09_G0G1_gene_level_statistics.tsv"), "\n",
  file.path(OUT, "STEP_E2B_COMPLETE.txt"), "\n",
  sep = ""
)
