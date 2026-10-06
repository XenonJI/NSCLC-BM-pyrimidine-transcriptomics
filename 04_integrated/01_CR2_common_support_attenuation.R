# Small-sample, common-support, and attenuation analyses
# Applies CR2/Satterthwaite inference, metadata-defined common-support restrictions,
# proliferation-conditioned models, formal attenuation tests, influence diagnostics, and
# multiplicity correction for the specimen-level spatial shift tests.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260818)

ROOT <- "D:/LUAD_LM_PYRIMIDINE"
OUT  <- file.path(ROOT, "results", "major_revision", "E8_critical_reviewer_defense")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# 0. Packages
# ----------------------------------------------------------------------
required <- c("data.table", "clubSandwich")
missing_required <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_required) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_required, collapse = ", "),
    "\nInstall once with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_required), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(clubSandwich)
})

HAS_WILD_BOOT <- requireNamespace("fwildclusterboot", quietly = TRUE)

if (HAS_WILD_BOOT && requireNamespace("dqrng", quietly = TRUE)) {
  dqrng::dqset.seed(20260818)
}

cat("\n============================================================\n")
cat("STEP E8: INTEGRATED ROBUSTNESS AND ATTENUATION ANALYSIS\n")
cat("============================================================\n")
cat("fwildclusterboot available: ", HAS_WILD_BOOT, "\n", sep = "")

# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------
find_exact_file <- function(basename_target) {
  roots <- c(
    file.path(ROOT, "results", "major_revision"),
    file.path(ROOT, "results", "reviewer_defense"),
    file.path(ROOT, "results"),
    ROOT
  )

  hits <- character(0)

  for (rr in roots) {
    if (!dir.exists(rr)) next
    hh <- list.files(
      rr,
      recursive = TRUE,
      full.names = TRUE,
      pattern = paste0("^", gsub("\\.", "\\\\.", basename_target), "$")
    )
    hits <- unique(c(hits, hh))
  }

  if (length(hits) == 0L) {
    stop("Required file not found anywhere under project: ", basename_target)
  }

  # Prefer major_revision, then reviewer_defense, then shortest path.
  ord <- order(
    !grepl("[/\\\\]major_revision[/\\\\]", hits),
    !grepl("[/\\\\]reviewer_defense[/\\\\]", hits),
    nchar(hits)
  )
  hits[ord][1]
}

safe_fread <- function(base) {
  f <- find_exact_file(base)
  cat("Input: ", f, "\n", sep = "")
  fread(f, showProgress = FALSE)
}

cr2_site <- function(dat, formula, term = "site_BM", label = "") {
  fit <- lm(formula, data = dat)

  tt <- clubSandwich::coef_test(
    fit,
    vcov = "CR2",
    cluster = dat$corrected_patient_id,
    test = "Satterthwaite",
    alternative = "two-sided"
  )

  tt <- as.data.frame(tt)

  coef_name_col <- grep(
    "^Coef\\.?$|^coefficient$|^term$",
    names(tt),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(coef_name_col) >= 1L) {
    tt$term_internal <- as.character(tt[[coef_name_col[1]]])
  } else {
    tt$term_internal <- rownames(tt)
  }

  rr <- tt[tt$term_internal == term, , drop = FALSE]

  if (nrow(rr) != 1L) {
    stop(
      "Could not uniquely extract term ", term,
      " from model: ", deparse(formula)
    )
  }

  pick_num <- function(patterns) {
    nms <- names(rr)
    for (pat in patterns) {
      hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
      if (length(hit) >= 1L) return(as.numeric(rr[[hit[1]]]))
    }
    NA_real_
  }

  beta <- unname(coef(fit)[term])
  se   <- pick_num(c("^SE$", "std"))
  df   <- pick_num(c("df_Satt", "d.f"))
  p    <- pick_num(c("p_Satt", "p-val", "p_value", "^p$"))

  crit <- if (is.finite(df)) qt(0.975, df = df) else qnorm(0.975)

  data.table(
    analysis = label,
    formula = paste(deparse(formula), collapse = ""),
    n_obs = nrow(dat),
    n_clusters = uniqueN(dat$corrected_patient_id),
    n_BM = sum(dat$site_BM == 1L),
    n_PT = sum(dat$site_BM == 0L),
    estimate_site_BM = beta,
    CR2_SE = se,
    Satterthwaite_df = df,
    two_sided_p_CR2 = p,
    CI95_low_CR2 = beta - crit * se,
    CI95_high_CR2 = beta + crit * se
  )
}

wild_boot_site <- function(dat, formula, term = "site_BM", label = "") {
  if (!HAS_WILD_BOOT) {
    return(data.table(
      analysis = label,
      status = "NOT_RUN_fwildclusterboot_not_installed",
      B = NA_integer_,
      p_wild_cluster = NA_real_,
      CI95_low_wild = NA_real_,
      CI95_high_wild = NA_real_
    ))
  }

  tryCatch({
    fit <- lm(formula, data = dat)

    # Current fwildclusterboot syntax supports lm objects with character
    # clustid / param. Reproducibility seed is set above.
    boot <- fwildclusterboot::boottest(
      fit,
      param = term,
      clustid = "corrected_patient_id",
      B = 9999,
      impose_null = TRUE
    )

    td <- tryCatch(
      generics::tidy(boot),
      error = function(e) NULL
    )

    if (!is.null(td) && nrow(td) >= 1L) {
      return(data.table(
        analysis = label,
        status = "PASS",
        B = 9999L,
        p_wild_cluster = as.numeric(td$p.value[1]),
        CI95_low_wild = if ("conf.low" %in% names(td)) as.numeric(td$conf.low[1]) else NA_real_,
        CI95_high_wild = if ("conf.high" %in% names(td)) as.numeric(td$conf.high[1]) else NA_real_
      ))
    }

    # Fallback if tidy method changes.
    pp <- tryCatch(fwildclusterboot::pval(boot), error = function(e) NA_real_)
    cc <- tryCatch(stats::confint(boot), error = function(e) c(NA_real_, NA_real_))

    data.table(
      analysis = label,
      status = "PASS_fallback_extraction",
      B = 9999L,
      p_wild_cluster = as.numeric(pp[1]),
      CI95_low_wild = as.numeric(cc[1]),
      CI95_high_wild = as.numeric(cc[2])
    )
  }, error = function(e) {
    data.table(
      analysis = label,
      status = paste0("SKIPPED_ERROR: ", conditionMessage(e)),
      B = 9999L,
      p_wild_cluster = NA_real_,
      CI95_low_wild = NA_real_,
      CI95_high_wild = NA_real_
    )
  })
}

smd_continuous <- function(x_bm, x_pt) {
  den <- sqrt((var(x_bm, na.rm = TRUE) + var(x_pt, na.rm = TRUE)) / 2)
  if (!is.finite(den) || den == 0) return(NA_real_)
  (mean(x_bm, na.rm = TRUE) - mean(x_pt, na.rm = TRUE)) / den
}

smd_binary <- function(x_bm, x_pt) {
  p1 <- mean(x_bm, na.rm = TRUE)
  p0 <- mean(x_pt, na.rm = TRUE)
  den <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  if (!is.finite(den) || den == 0) return(NA_real_)
  (p1 - p0) / den
}

fit_site_plain <- function(dat, outcome = "module_score") {
  if (nrow(dat) < 4L || uniqueN(dat$site_BM) < 2L) return(NA_real_)
  fit <- lm(reformulate("site_BM", outcome), data = dat)
  unname(coef(fit)["site_BM"])
}

# ----------------------------------------------------------------------
# 2. Load compact source-of-truth inputs
# ----------------------------------------------------------------------
score_src <- safe_fread("STEP_B2v2_03_locked41_to_official_metadata_match.tsv")
meta_src  <- safe_fread("STEP_E4_06_MASTER_SAMPLE_PROVENANCE_WITH_OFFICIAL_METADATA.tsv")
cycle_src <- safe_fread("STEP_E2B_04_specimen_cellcycle_composition.tsv")
g1_src    <- safe_fread("STEP_E2B_06_G0G1_like_module_score.tsv")
spatial_src <- safe_fread("STEP_E5_04_toroidal_shift_null_summary.tsv")

needed_score <- c("samples", "module_score")
needed_meta <- c(
  "sample", "corrected_patient_id", "site_analysis",
  "included_primary_analysis", "sample_prefix_stratum",
  "snrna_seq_chemistry", "age_at_resection_of_profiled_specimen", "sex"
)
needed_cycle <- c(
  "sample", "corrected_patient_id", "cohort",
  "malignant_nuclei", "G0G1_like_nuclei", "cycling_like_nuclei",
  "cycling_like_fraction"
)
needed_g1 <- c(
  "sample", "corrected_patient_id", "cohort",
  "G0G1_like_nuclei", "module_score_G0G1_fixed"
)

stopifnot(all(needed_score %in% names(score_src)))
stopifnot(all(needed_meta %in% names(meta_src)))
stopifnot(all(needed_cycle %in% names(cycle_src)))
stopifnot(all(needed_g1 %in% names(g1_src)))

score_src <- score_src[, .(
  sample = as.character(samples),
  module_score = as.numeric(module_score)
)]

main <- merge(
  meta_src[included_primary_analysis == TRUE, ..needed_meta],
  score_src,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

if (nrow(main) != 41L || any(!is.finite(main$module_score))) {
  stop("Main 41-specimen score merge failed.")
}

main[, site_BM := as.integer(site_analysis == "Brain_metastasis")]
main[, chemistry_V2 := as.integer(grepl("V2", snrna_seq_chemistry, fixed = TRUE))]
main[, sex_M := as.integer(toupper(sex) == "M")]
main[, age := as.numeric(age_at_resection_of_profiled_specimen)]
main[, age_c := age - mean(age, na.rm = TRUE)]
main[, prefix := as.character(sample_prefix_stratum)]

if (sum(main$site_BM) != 31L || sum(main$site_BM == 0L) != 10L) {
  stop("Unexpected primary cohort counts after merge.")
}

fwrite(main, file.path(OUT, "STEP_E8_01_analysis_dataset_41.tsv"), sep = "\t")

# ----------------------------------------------------------------------
# 3. Positivity / common-support audit
# ----------------------------------------------------------------------
chemistry_cross <- main[, .N, by = .(site_analysis, snrna_seq_chemistry)]
prefix_cross <- main[, .N, by = .(site_analysis, prefix)]

age_by_site <- main[, .(
  n = .N,
  mean_age = mean(age),
  sd_age = sd(age),
  median_age = median(age),
  min_age = min(age),
  max_age = max(age)
), by = site_analysis]

age_lo <- max(
  main[site_BM == 1L, min(age)],
  main[site_BM == 0L, min(age)]
)
age_hi <- min(
  main[site_BM == 1L, max(age)],
  main[site_BM == 0L, max(age)]
)

smd_table <- rbindlist(list(
  data.table(
    covariate = "Age (years)",
    type = "continuous",
    SMD_BM_minus_PT = smd_continuous(
      main[site_BM == 1L, age],
      main[site_BM == 0L, age]
    )
  ),
  data.table(
    covariate = "Chemistry: V2",
    type = "binary",
    SMD_BM_minus_PT = smd_binary(
      main[site_BM == 1L, chemistry_V2],
      main[site_BM == 0L, chemistry_V2]
    )
  ),
  data.table(
    covariate = "Sex: male",
    type = "binary",
    SMD_BM_minus_PT = smd_binary(
      main[site_BM == 1L, sex_M],
      main[site_BM == 0L, sex_M]
    )
  )
))

for (pp in sort(unique(main$prefix))) {
  tmp <- data.table(
    covariate = paste0("Prefix: ", pp),
    type = "binary",
    SMD_BM_minus_PT = smd_binary(
      as.integer(main[site_BM == 1L, prefix] == pp),
      as.integer(main[site_BM == 0L, prefix] == pp)
    )
  )
  smd_table <- rbind(smd_table, tmp, fill = TRUE)
}

fwrite(chemistry_cross, file.path(OUT, "STEP_E8_02_chemistry_by_site.tsv"), sep = "\t")
fwrite(prefix_cross, file.path(OUT, "STEP_E8_03_prefix_by_site.tsv"), sep = "\t")
fwrite(age_by_site, file.path(OUT, "STEP_E8_04_age_by_site.tsv"), sep = "\t")
fwrite(smd_table, file.path(OUT, "STEP_E8_05_covariate_standardized_mean_differences.tsv"), sep = "\t")

writeLines(
  c(
    paste0("empirical_age_common_support_low=", age_lo),
    paste0("empirical_age_common_support_high=", age_hi),
    "definition=intersection of observed age ranges in BM and PT; defined without outcome information"
  ),
  file.path(OUT, "STEP_E8_06_age_common_support_definition.txt")
)

# Prespecified/common-support subsets
sub_v2 <- main[chemistry_V2 == 1L]
sub_prefix <- main[prefix %chin% c("KRAS", "STK")]
sub_age <- main[age >= age_lo & age <= age_hi]

# ----------------------------------------------------------------------
# 4. Core site-identification models with CR2 Satterthwaite
# ----------------------------------------------------------------------
site_models <- rbindlist(list(
  cr2_site(
    main,
    module_score ~ site_BM,
    label = "All 41: site only"
  ),
  cr2_site(
    main,
    module_score ~ site_BM + chemistry_V2 + age_c + sex_M,
    label = "All 41: site + chemistry + age + sex"
  ),
  cr2_site(
    sub_v2,
    module_score ~ site_BM,
    label = "Common chemistry only: 5' V2"
  ),
  cr2_site(
    sub_prefix,
    module_score ~ site_BM + factor(prefix),
    label = "Common-prefix subset: KRAS/STK + prefix adjustment"
  ),
  cr2_site(
    sub_age,
    module_score ~ site_BM + age_c,
    label = paste0("Empirical age common support ", age_lo, "-", age_hi, " years + age")
  )
), fill = TRUE)

fwrite(site_models, file.path(OUT, "STEP_E8_07_site_identification_CR2_models.tsv"), sep = "\t")

# Optional wild-cluster bootstrap for the same central models
wild_models <- rbindlist(list(
  wild_boot_site(
    main,
    module_score ~ site_BM,
    label = "All 41: site only"
  ),
  wild_boot_site(
    main,
    module_score ~ site_BM + chemistry_V2 + age_c + sex_M,
    label = "All 41: site + chemistry + age + sex"
  ),
  wild_boot_site(
    sub_v2,
    module_score ~ site_BM,
    label = "Common chemistry only: 5' V2"
  ),
  wild_boot_site(
    sub_prefix,
    module_score ~ site_BM + factor(prefix),
    label = "Common-prefix subset: KRAS/STK + prefix adjustment"
  )
), fill = TRUE)

fwrite(wild_models, file.path(OUT, "STEP_E8_08_site_identification_wild_cluster_bootstrap.tsv"), sep = "\t")

# ----------------------------------------------------------------------
# 5. Influence analyses: patient / prefix / chemistry
# ----------------------------------------------------------------------
patients <- unique(main$corrected_patient_id)

loo_patient <- rbindlist(lapply(patients, function(pid) {
  dd <- main[corrected_patient_id != pid]
  data.table(
    omitted_patient = pid,
    n_obs = nrow(dd),
    n_clusters = uniqueN(dd$corrected_patient_id),
    estimate_site_BM = fit_site_plain(dd)
  )
}))

prefixes <- sort(unique(main$prefix))
loo_prefix <- rbindlist(lapply(prefixes, function(pp) {
  dd <- main[prefix != pp]
  data.table(
    omitted_prefix = pp,
    n_obs = nrow(dd),
    n_BM = sum(dd$site_BM == 1L),
    n_PT = sum(dd$site_BM == 0L),
    estimate_site_BM = fit_site_plain(dd)
  )
}))

chem_levels <- sort(unique(main$snrna_seq_chemistry))
leave_chem <- rbindlist(lapply(chem_levels, function(cc) {
  dd <- main[snrna_seq_chemistry != cc]
  data.table(
    omitted_chemistry = cc,
    retained_chemistry = paste(sort(unique(dd$snrna_seq_chemistry)), collapse = ";"),
    n_obs = nrow(dd),
    n_BM = sum(dd$site_BM == 1L),
    n_PT = sum(dd$site_BM == 0L),
    estimate_site_BM = fit_site_plain(dd)
  )
}))

fwrite(loo_patient, file.path(OUT, "STEP_E8_09_leave_one_patient_out.tsv"), sep = "\t")
fwrite(loo_prefix, file.path(OUT, "STEP_E8_10_leave_one_prefix_out.tsv"), sep = "\t")
fwrite(leave_chem, file.path(OUT, "STEP_E8_11_leave_one_chemistry_stratum_out.tsv"), sep = "\t")

# ----------------------------------------------------------------------
# 6. Formal proliferation-conditioning models
# ----------------------------------------------------------------------
cycle <- copy(cycle_src)
cycle[, sample := as.character(sample)]

cycle_dat <- merge(
  cycle,
  main[, .(
    sample, corrected_patient_id, module_score, site_BM,
    chemistry_V2, age_c, sex_M, prefix
  )],
  by = c("sample", "corrected_patient_id"),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(cycle_dat) != 40L || any(!is.finite(cycle_dat$module_score))) {
  stop("40-specimen cycling/main-score merge failed.")
}

# A. Does continuous cycling fraction account for the site-associated score?
prolif_score_models <- rbindlist(list(
  cr2_site(
    cycle_dat,
    module_score ~ site_BM,
    label = "40 concordant specimens: site only"
  ),
  cr2_site(
    cycle_dat,
    module_score ~ site_BM + cycling_like_fraction,
    label = "40 concordant: site + continuous cycling fraction"
  ),
  cr2_site(
    cycle_dat,
    module_score ~ site_BM + cycling_like_fraction + chemistry_V2 + age_c + sex_M,
    label = "40 concordant: site + cycling + chemistry + age + sex"
  )
), fill = TRUE)

fwrite(
  prolif_score_models,
  file.path(OUT, "STEP_E8_12_module_score_conditioned_on_cycling_CR2.tsv"),
  sep = "\t"
)

# B. Is the cycling fraction itself still associated with site after
# documented technical / clinical composition is considered?
# Each specimen remains one biological replicate.
cycle_dat[, cycling_logit := qlogis(
  (cycling_like_nuclei + 0.5) / (malignant_nuclei + 1)
)]

cycling_models_raw <- rbindlist(list(
  cr2_site(
    cycle_dat,
    cycling_like_fraction ~ site_BM,
    label = "Cycling fraction: site only"
  ),
  cr2_site(
    cycle_dat,
    cycling_like_fraction ~ site_BM + chemistry_V2,
    label = "Cycling fraction: site + chemistry"
  ),
  cr2_site(
    cycle_dat,
    cycling_like_fraction ~ site_BM + age_c,
    label = "Cycling fraction: site + age"
  ),
  cr2_site(
    cycle_dat,
    cycling_like_fraction ~ site_BM + chemistry_V2 + age_c + sex_M,
    label = "Cycling fraction: site + chemistry + age + sex"
  )
), fill = TRUE)

cycling_models_logit <- rbindlist(list(
  cr2_site(
    cycle_dat,
    cycling_logit ~ site_BM,
    label = "Empirical-logit cycling fraction: site only"
  ),
  cr2_site(
    cycle_dat,
    cycling_logit ~ site_BM + chemistry_V2,
    label = "Empirical-logit cycling fraction: site + chemistry"
  ),
  cr2_site(
    cycle_dat,
    cycling_logit ~ site_BM + chemistry_V2 + age_c + sex_M,
    label = "Empirical-logit cycling fraction: full covariates"
  )
), fill = TRUE)

# Common-support cycling analyses
cycle_v2 <- cycle_dat[chemistry_V2 == 1L]
cycle_prefix <- cycle_dat[prefix %chin% c("KRAS", "STK")]

cycling_common_support <- rbindlist(list(
  cr2_site(
    cycle_v2,
    cycling_like_fraction ~ site_BM,
    label = "Cycling fraction: V2-only common chemistry"
  ),
  cr2_site(
    cycle_prefix,
    cycling_like_fraction ~ site_BM + factor(prefix),
    label = "Cycling fraction: KRAS/STK common-prefix + prefix adjustment"
  )
), fill = TRUE)

fwrite(cycling_models_raw, file.path(OUT, "STEP_E8_13_cycling_fraction_covariate_CR2.tsv"), sep = "\t")
fwrite(cycling_models_logit, file.path(OUT, "STEP_E8_14_cycling_fraction_empirical_logit_CR2.tsv"), sep = "\t")
fwrite(cycling_common_support, file.path(OUT, "STEP_E8_15_cycling_fraction_common_support_CR2.tsv"), sep = "\t")

# ----------------------------------------------------------------------
# 7. Formal attenuation: all malignant vs marker-defined G0/G1-like
# ----------------------------------------------------------------------
g1 <- copy(g1_src)
g1[, sample := as.character(sample)]
g1 <- g1[G0G1_like_nuclei >= 50L]

paired <- merge(
  cycle_dat[, .(
    sample, corrected_patient_id, site_BM, module_score,
    chemistry_V2, age_c, sex_M
  )],
  g1[, .(
    sample, module_score_G0G1_fixed
  )],
  by = "sample",
  all = FALSE,
  sort = FALSE
)

if (nrow(paired) < 35L) {
  stop("Unexpectedly small matched all-vs-G0G1 specimen set.")
}

# 7A. Stacked paired compartment interaction:
# interaction = change in site coefficient in G0/G1-like relative to all malignant.
stacked <- rbindlist(list(
  paired[, .(
    sample, corrected_patient_id, site_BM,
    chemistry_V2, age_c, sex_M,
    state_G0G1 = 0L,
    score = module_score
  )],
  paired[, .(
    sample, corrected_patient_id, site_BM,
    chemistry_V2, age_c, sex_M,
    state_G0G1 = 1L,
    score = module_score_G0G1_fixed
  )]
))

fit_interaction <- lm(
  score ~ site_BM * state_G0G1,
  data = stacked
)

interaction_test <- clubSandwich::coef_test(
  fit_interaction,
  vcov = "CR2",
  cluster = stacked$corrected_patient_id,
  test = "Satterthwaite",
  alternative = "two-sided"
)

interaction_df <- as.data.frame(interaction_test)

interaction_coef_col <- grep(
  "^Coef\\.?$|^coefficient$|^term$",
  names(interaction_df),
  ignore.case = TRUE,
  value = TRUE
)

if (length(interaction_coef_col) >= 1L) {
  interaction_df$term <- as.character(interaction_df[[interaction_coef_col[1]]])
} else {
  interaction_df$term <- rownames(interaction_df)
}

interaction_out <- as.data.table(interaction_df)

fwrite(
  interaction_out,
  file.path(OUT, "STEP_E8_16_all_vs_G0G1_site_by_state_interaction_CR2.tsv"),
  sep = "\t"
)

# 7B. Cluster bootstrap of Delta beta = beta_all - beta_G0G1.
# Resampling unit = corrected patient cluster, preserving both specimens
# when a matched patient is sampled.
cluster_ids <- unique(paired$corrected_patient_id)
B_DELTA <- 10000L

beta_all_obs <- fit_site_plain(
  paired[, .(module_score, site_BM)]
)
beta_g1_obs <- {
  fit <- lm(module_score_G0G1_fixed ~ site_BM, data = paired)
  unname(coef(fit)["site_BM"])
}
delta_obs <- beta_all_obs - beta_g1_obs

boot_delta <- rep(NA_real_, B_DELTA)

for (b in seq_len(B_DELTA)) {
  draw <- sample(cluster_ids, length(cluster_ids), replace = TRUE)

  boot_rows <- rbindlist(lapply(seq_along(draw), function(j) {
    z <- paired[corrected_patient_id == draw[j]]
    z[, bootstrap_cluster := paste0("draw_", j)]
    z
  }))

  if (uniqueN(boot_rows$site_BM) < 2L) next

  ba <- unname(coef(lm(module_score ~ site_BM, data = boot_rows))["site_BM"])
  bg <- unname(coef(lm(module_score_G0G1_fixed ~ site_BM, data = boot_rows))["site_BM"])

  if (is.finite(ba) && is.finite(bg)) {
    boot_delta[b] <- ba - bg
  }
}

boot_delta_ok <- boot_delta[is.finite(boot_delta)]

if (length(boot_delta_ok) < 0.95 * B_DELTA) {
  warning("More than 5% of cluster bootstrap draws were invalid.")
}

ci_delta <- quantile(boot_delta_ok, probs = c(0.025, 0.975), names = FALSE)
p_delta <- 2 * min(
  (1 + sum(boot_delta_ok <= 0)) / (length(boot_delta_ok) + 1),
  (1 + sum(boot_delta_ok >= 0)) / (length(boot_delta_ok) + 1)
)
p_delta <- min(1, p_delta)

delta_summary <- data.table(
  n_specimens = nrow(paired),
  n_patient_clusters = uniqueN(paired$corrected_patient_id),
  beta_all_same_specimens = beta_all_obs,
  beta_G0G1_same_specimens = beta_g1_obs,
  delta_beta_all_minus_G0G1 = delta_obs,
  bootstrap_valid_draws = length(boot_delta_ok),
  bootstrap_CI95_low = ci_delta[1],
  bootstrap_CI95_high = ci_delta[2],
  bootstrap_two_sided_p_delta = p_delta
)

fwrite(
  delta_summary,
  file.path(OUT, "STEP_E8_17_delta_beta_cluster_bootstrap.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 8. Spatial toroidal-shift multiplicity
# ----------------------------------------------------------------------
sp_adj <- spatial_src[
  association == "partial_rank_residual_surface"
]

if (nrow(sp_adj) != 12L) {
  stop("Expected exactly 12 adjusted toroidal-shift specimen rows.")
}

sp_adj[, BH_FDR_12 := p.adjust(toroidal_two_sided_p, method = "BH")]
sp_adj[, survives_BH_0_05 := BH_FDR_12 < 0.05]

fwrite(
  sp_adj,
  file.path(OUT, "STEP_E8_18_spatial_toroidal_shift_BH_FDR.tsv"),
  sep = "\t"
)

spatial_summary <- sp_adj[, .(
  n_specimens = .N,
  n_nominal_P_lt_0_05 = sum(toroidal_two_sided_p < 0.05),
  n_BH_FDR_lt_0_05 = sum(BH_FDR_12 < 0.05),
  min_nominal_P = min(toroidal_two_sided_p),
  min_BH_FDR = min(BH_FDR_12),
  specimen_with_min_P = sample[which.min(toroidal_two_sided_p)]
)]

fwrite(
  spatial_summary,
  file.path(OUT, "STEP_E8_19_spatial_multiplicity_summary.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 9. Duplicate-threshold reporting audit
# ----------------------------------------------------------------------
sens <- safe_fread("STEP_E1_07_fixed_score_sensitivity_results.tsv")
duplicate_threshold_note <- data.table(
  item = c(
    "Primary >=50 vs minimum 200",
    "G0/G1 >=50 vs >=20"
  ),
  recommendation = c(
    "Do not present as two independent robustness checks if all values are identical. Keep >=50 as primary and state that all eligible specimens also exceeded 200 in this dataset, or move the redundant threshold to a footnote.",
    "Do not present >=20 as additional evidence when it yields exactly the same specimens/results as >=50. Keep >=50 as the reported analysis and note that the >=20 sensitivity was identical."
  )
)

fwrite(
  duplicate_threshold_note,
  file.path(OUT, "STEP_E8_20_duplicate_threshold_reporting_guardrail.tsv"),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 10. Machine-readable interpretation guardrails
# ----------------------------------------------------------------------
full_cr2 <- site_models[analysis == "All 41: site + chemistry + age + sex"]
cycling_cond <- prolif_score_models[
  analysis == "40 concordant: site + continuous cycling fraction"
]

interaction_row_name <- grep(
  "site_BM:state_G0G1",
  interaction_out$term,
  value = TRUE
)

interaction_p <- NA_real_
interaction_beta <- NA_real_

if (length(interaction_row_name) == 1L) {
  rr <- interaction_out[term == interaction_row_name]

  bcol <- grep("^beta$|estimate", names(rr), ignore.case = TRUE, value = TRUE)
  if (length(bcol) >= 1L) interaction_beta <- as.numeric(rr[[bcol[1]]])

  pcol <- grep("p_Satt|p-val|p_value", names(rr), ignore.case = TRUE, value = TRUE)
  if (length(pcol) >= 1L) interaction_p <- as.numeric(rr[[pcol[1]]])
}

guardrails <- c(
  "status=PASS_E8_INTEGRATED_ROBUSTNESS_COMPLETE",
  "primary_endpoint_replaced=FALSE",
  "gene_set_changed=FALSE",
  "result_driven_exclusion=FALSE",
  "two_sided_inference=TRUE",
  paste0("age_common_support=", age_lo, " to ", age_hi, " years"),
  paste0(
    "full_covariate_CR2_site_beta=",
    signif(full_cr2$estimate_site_BM, 6),
    "; p=",
    signif(full_cr2$two_sided_p_CR2, 6)
  ),
  paste0(
    "cycling_conditioned_site_beta=",
    signif(cycling_cond$estimate_site_BM, 6),
    "; p=",
    signif(cycling_cond$two_sided_p_CR2, 6)
  ),
  paste0(
    "all_vs_G0G1_interaction_beta=",
    signif(interaction_beta, 6),
    "; p=",
    signif(interaction_p, 6)
  ),
  paste0(
    "delta_beta_bootstrap=",
    signif(delta_summary$delta_beta_all_minus_G0G1, 6),
    "; CI=",
    signif(delta_summary$bootstrap_CI95_low, 6),
    " to ",
    signif(delta_summary$bootstrap_CI95_high, 6),
    "; p=",
    signif(delta_summary$bootstrap_two_sided_p_delta, 6)
  ),
  paste0(
    "spatial_nominal_significant=",
    spatial_summary$n_nominal_P_lt_0_05,
    "/12"
  ),
  paste0(
    "spatial_BH_significant=",
    spatial_summary$n_BH_FDR_lt_0_05,
    "/12; min_BH_FDR=",
    signif(spatial_summary$min_BH_FDR, 6)
  ),
  "wording_guardrail_site=Use 'site-associated difference' unless common-support and covariate analyses support stronger wording; do not claim an independent site effect when the full model is compatible with zero.",
  "wording_guardrail_G0G1=Do not infer absence of a residual effect from P>0.05. Report attenuation point estimate, CI, and the formal interaction/delta-beta test.",
  "wording_guardrail_spatial=Treat toroidal-shift inference as the local-alignment analysis; report BH correction across 12 specimen-specific tests.",
  "wording_guardrail_proliferation=Use 'proliferation-coupled' as a title-level claim only if the formal site-by-state / delta-beta evidence supports effect attenuation and the cycling-fraction association is not explained by documented composition."
)

writeLines(
  guardrails,
  file.path(OUT, "STEP_E8_COMPLETE.txt")
)

# ----------------------------------------------------------------------
# 11. Console summary
# ----------------------------------------------------------------------
cat("\n============================================================\n")
cat("STEP E8 COMPLETE\n")
cat("============================================================\n\n")

cat("CR2 site-identification models:\n")
print(site_models)

cat("\nModule score conditioned on continuous cycling fraction:\n")
print(prolif_score_models)

cat("\nCycling-fraction covariate models:\n")
print(cycling_models_raw)

cat("\nFormal delta-beta bootstrap:\n")
print(delta_summary)

cat("\nSpatial multiplicity:\n")
print(spatial_summary)

if (!HAS_WILD_BOOT) {
  cat(
    "\nNOTE: fwildclusterboot was not installed, so wild-cluster bootstrap was safely skipped.\n",
    "Optional install:\n",
    "install.packages('fwildclusterboot')\n",
    "Then rerun this same script.\n",
    sep = ""
  )
}

cat("\nStatus: PASS_E8_INTEGRATED_ROBUSTNESS_COMPLETE\n")
cat("\nKey output files:\n")
cat(file.path(OUT, "STEP_E8_07_site_identification_CR2_models.tsv"), "\n")
cat(file.path(OUT, "STEP_E8_12_module_score_conditioned_on_cycling_CR2.tsv"), "\n")
cat(file.path(OUT, "STEP_E8_13_cycling_fraction_covariate_CR2.tsv"), "\n")
cat(file.path(OUT, "STEP_E8_15_cycling_fraction_common_support_CR2.tsv"), "\n")
cat(file.path(OUT, "STEP_E8_16_all_vs_G0G1_site_by_state_interaction_CR2.tsv"), "\n")
cat(file.path(OUT, "STEP_E8_17_delta_beta_cluster_bootstrap.tsv"), "\n")
cat(file.path(OUT, "STEP_E8_18_spatial_toroidal_shift_BH_FDR.tsv"), "\n")
cat(file.path(OUT, "STEP_E8_19_spatial_multiplicity_summary.tsv"), "\n")
cat(file.path(OUT, "STEP_E8_COMPLETE.txt"), "\n")
