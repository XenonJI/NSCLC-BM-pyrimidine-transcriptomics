# Spatial sample-identity audit
# Compares deposited spatial count/coordinate payloads to identify duplicated processed content,
# including the PA056/KRAS_11 duplicate pair used in the conservative spatial sensitivity analysis.

rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))

ROOT <- "D:/LUAD_LM_PYRIMIDINE"
INV <- file.path(ROOT, "results/reviewer_defense", "STEP_C1_01_GSE223501_spatial_file_inventory.tsv")
TARGET_DIR <- file.path(ROOT, "data/processed/GSE223501_spatial_targets")
DEPTH_DIR <- file.path(ROOT, "data/processed/GSE223501_spatial_library_depth")
OUT <- file.path(ROOT, "results/reviewer_defense")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

GENES <- c("TYMS","UMPS","MKI67","PCNA","TOP2A","UBE2C","CENPF","BIRC5","CCNB1",
           "EPCAM","KRT8","KRT18","KRT19","MUC1")

safe_stub <- function(x) gsub("[^A-Za-z0-9._-]+","_",x)

gz_text_identical <- function(f1, f2, chunk=5L) {
  c1 <- gzfile(f1, "rt"); c2 <- gzfile(f2, "rt")
  on.exit(close(c1), add=TRUE); on.exit(close(c2), add=TRUE)
  seen <- 0L
  repeat {
    a <- readLines(c1, n=chunk, warn=FALSE)
    b <- readLines(c2, n=chunk, warn=FALSE)
    if (length(a) != length(b)) return(list(identical=FALSE, first_difference_line=seen+1L))
    if (length(a) == 0L) return(list(identical=TRUE, first_difference_line=NA_integer_))
    z <- which(a != b)
    if (length(z)) return(list(identical=FALSE, first_difference_line=seen+z[1]))
    seen <- seen + length(a)
  }
}

cat("\nSTEP C6A: sample identity / duplicate-content audit\n")

if (!file.exists(INV)) stop("Missing inventory: ", INV)
inv <- fread(INV)
if (nrow(inv) != 14L) stop("Inventory does not contain 14 samples.")

# Raw compressed-file fingerprints
fp <- copy(inv)
fp[, counts_md5 := unname(tools::md5sum(counts_file))]
fp[, coords_md5 := unname(tools::md5sum(coords_file))]
fp[, counts_bytes := file.info(counts_file)$size]
fp[, coords_bytes := file.info(coords_file)$size]

# Load only compact objects
compact <- vector("list", nrow(inv)); names(compact) <- inv$official_paper_id
depths <- vector("list", nrow(inv)); names(depths) <- inv$official_paper_id

for (i in seq_len(nrow(inv))) {
  s <- inv$official_paper_id[i]
  tf <- file.path(TARGET_DIR, paste0(safe_stub(s), "_target_counts_coords.rds"))
  df <- file.path(DEPTH_DIR, paste0(safe_stub(s), "_library_depth.rds"))
  if (!file.exists(tf)) stop("Missing target RDS: ", tf)
  if (!file.exists(df)) stop("Missing depth RDS: ", df)
  d <- readRDS(tf); dep <- readRDS(df)
  need <- c("bead_key","x","y",GENES)
  miss <- setdiff(need, names(d))
  if (length(miss)) stop(s, " target RDS missing: ", paste(miss, collapse=", "))
  compact[[s]] <- d[, ..need]
  depths[[s]] <- dep[, .(bead_key, library_size, detected_genes)]
}

pairs <- combn(inv$official_paper_id, 2, simplify=FALSE)

pair_audit <- rbindlist(lapply(pairs, function(p) {
  a <- compact[[p[1]]]; b <- compact[[p[2]]]
  da <- depths[[p[1]]]; db <- depths[[p[2]]]
  same_n <- nrow(a) == nrow(b)
  bead_same <- same_n && identical(a$bead_key, b$bead_key)
  coord_same <- same_n && identical(a$x,b$x) && identical(a$y,b$y)
  target_same <- same_n && all(vapply(GENES, function(g) identical(a[[g]], b[[g]]), logical(1)))
  depth_same <- nrow(da)==nrow(db) &&
                identical(da$bead_key,db$bead_key) &&
                identical(da$library_size,db$library_size) &&
                identical(da$detected_genes,db$detected_genes)
  i1 <- match(p[1], inv$official_paper_id); i2 <- match(p[2], inv$official_paper_id)
  data.table(
    sample_1=p[1], cohort_1=inv$cohort[i1],
    sample_2=p[2], cohort_2=inv$cohort[i2],
    same_n_beads=same_n,
    bead_key_identical=bead_same,
    coords_identical=coord_same,
    all_14_target_counts_identical=target_same,
    compact_content_identical=(bead_same && coord_same && target_same),
    library_depth_content_identical=depth_same,
    compressed_counts_md5_identical=(fp$counts_md5[i1] == fp$counts_md5[i2]),
    compressed_coords_md5_identical=(fp$coords_md5[i1] == fp$coords_md5[i2])
  )
}))

# Any suspicious pair, plus PA056 vs KRAS_11 explicitly
flag <- pair_audit[
  compact_content_identical | library_depth_content_identical |
  compressed_counts_md5_identical | compressed_coords_md5_identical |
  ((sample_1=="PA056" & sample_2=="KRAS_11") | (sample_1=="KRAS_11" & sample_2=="PA056"))
]

raw_compare <- data.table()
if (nrow(flag)) {
  tmp <- vector("list", nrow(flag))
  for (i in seq_len(nrow(flag))) {
    s1 <- flag$sample_1[i]; s2 <- flag$sample_2[i]
    i1 <- match(s1, inv$official_paper_id); i2 <- match(s2, inv$official_paper_id)
    cat("\nRaw comparison: ", s1, " vs ", s2, "\n", sep="")
    cc <- gz_text_identical(inv$coords_file[i1], inv$coords_file[i2], chunk=1000L)
    rc <- gz_text_identical(inv$counts_file[i1], inv$counts_file[i2], chunk=5L)
    cat("  coords identical: ", cc$identical, "\n", sep="")
    cat("  counts identical: ", rc$identical, "\n", sep="")
    tmp[[i]] <- data.table(
      sample_1=s1, sample_2=s2,
      raw_coords_decompressed_identical=cc$identical,
      raw_coords_first_difference_line=cc$first_difference_line,
      raw_counts_decompressed_identical=rc$identical,
      raw_counts_first_difference_line=rc$first_difference_line
    )
  }
  raw_compare <- rbindlist(tmp, fill=TRUE)
}

exact_compact_dup <- pair_audit[compact_content_identical & library_depth_content_identical]
exact_raw_dup <- if (nrow(raw_compare)) raw_compare[
  raw_coords_decompressed_identical & raw_counts_decompressed_identical
] else data.table()

status <- if (nrow(exact_compact_dup) || nrow(exact_raw_dup)) {
  "CRITICAL_DUPLICATE_REVIEW_REQUIRED"
} else {
  "PASS_NO_DUPLICATE_SPECIMEN_CONTENT"
}

fwrite(fp, file.path(OUT, "STEP_C6A_01_raw_file_fingerprints.tsv"), sep="\t")
fwrite(pair_audit, file.path(OUT, "STEP_C6A_02_pairwise_compact_identity_audit.tsv"), sep="\t")
fwrite(raw_compare, file.path(OUT, "STEP_C6A_03_flagged_raw_decompressed_comparison.tsv"), sep="\t")

writeLines(c(
  paste0("status=", status),
  paste0("samples=", nrow(inv)),
  paste0("pairwise_comparisons=", nrow(pair_audit)),
  paste0("compact_exact_duplicate_pairs=", nrow(exact_compact_dup)),
  paste0("raw_exact_duplicate_pairs=", nrow(exact_raw_dup)),
  "biological_endpoint_testing_performed=FALSE",
  "correlation_testing_performed=FALSE",
  "p_values_computed=FALSE"
), file.path(OUT, "STEP_C6A_COMPLETE.txt"))

cat("\nPA056 vs KRAS_11:\n")
print(pair_audit[
  (sample_1=="PA056" & sample_2=="KRAS_11") |
  (sample_1=="KRAS_11" & sample_2=="PA056")
])
if (nrow(raw_compare)) {
  cat("\nFlagged raw comparisons:\n")
  print(raw_compare)
}
cat("\nStatus: ", status, "\n", sep="")
