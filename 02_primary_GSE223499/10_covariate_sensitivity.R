# Technical and clinical covariate sensitivity
# Fits the prespecified specimen-level site models with corrected patient clustering and
# documented chemistry, age, and sex covariates. Covariates are selected from source metadata,
# not from outcome associations.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/LUAD_LM_PYRIMIDINE"
OUT <- file.path(ROOT, "results", "major_revision")
RD <- file.path(ROOT, "results", "reviewer_defense")

LOCKED_FILE <- file.path(
  RD,
  "STEP_A3_locked_primary_31BM_10Primary_patient_table.tsv"
)

MASTER_FILE <- file.path(
  OUT,
  "STEP_E4_06_MASTER_SAMPLE_PROVENANCE_WITH_OFFICIAL_METADATA.tsv"
)

for (f in c(LOCKED_FILE, MASTER_FILE)) {
  if (!file.exists(f)) {
    stop("Required file not found:\n", f, call. = FALSE)
  }
}

# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------
cluster_robust_cr1 <- function(
  fit,
  cluster,
  term_name = "site_BM"
) {
  mf <- model.frame(fit)
  X <- model.matrix(fit)
  e <- residuals(fit)

  if (length(cluster) != nrow(mf)) {
    stop(
      "Cluster vector length does not match model rows.",
      call. = FALSE
    )
  }

  keep <- complete.cases(
    X,
    e,
    cluster
  )

  X <- X[
    keep,
    ,
    drop = FALSE
  ]

  e <- e[
    keep
  ]

  cluster <- as.character(
    cluster[
      keep
    ]
  )

  N <- nrow(X)
  K <- ncol(X)
  clusters <- unique(cluster)
  G <- length(clusters)

  if (G <= K) {
    stop(
      "Too few clusters relative to model parameters.",
      call. = FALSE
    )
  }

  XtX <- crossprod(X)

  if (qr(XtX)$rank < K) {
    stop(
      "Model matrix is rank deficient.",
      call. = FALSE
    )
  }

  bread <- solve(XtX)

  meat <- matrix(
    0,
    nrow = K,
    ncol = K
  )

  for (cl in clusters) {
    idx <- cluster == cl
    Xg <- X[
      idx,
      ,
      drop = FALSE
    ]
    eg <- e[
      idx
    ]

    sg <- crossprod(
      Xg,
      eg
    )

    meat <- meat +
      sg %*%
        t(sg)
  }

  correction <- (
    G /
      (G - 1)
  ) *
    (
      (N - 1) /
        (N - K)
    )

  vc <- correction *
    bread %*%
      meat %*%
      bread

  beta <- coef(fit)
  se <- sqrt(diag(vc))

  if (!term_name %in% names(beta)) {
    stop(
      "Requested coefficient not found: ",
      term_name,
      call. = FALSE
    )
  }

  est <- unname(
    beta[
      term_name
    ]
  )

  se0 <- unname(
    se[
      term_name
    ]
  )

  df <- G - 1L
  t0 <- est / se0
  p2 <- 2 * pt(
    -abs(t0),
    df = df
  )

  crit <- qt(
    0.975,
    df = df
  )

  data.table(
    estimate_site_BM = est,
    CR1_SE = se0,
    t_CR1 = t0,
    df_clusters = df,
    two_sided_p_CR1 = p2,
    CI95_low_CR1 = est - crit * se0,
    CI95_high_CR1 = est + crit * se0,
    n_observations = N,
    n_patient_clusters = G
  )
}

fit_model <- function(
  data,
  model_id,
  model_label,
  formula_text
) {
  vars <- all.vars(
    as.formula(
      formula_text
    )
  )

  needed <- unique(
    c(
      vars,
      "corrected_patient_id"
    )
  )

  d <- data[
    complete.cases(
      data[
        ,
        ..needed
      ]
    )
  ]

  fit <- lm(
    as.formula(
      formula_text
    ),
    data = d
  )

  cr <- cluster_robust_cr1(
    fit = fit,
    cluster = d$corrected_patient_id,
    term_name = "site_BM"
  )

  data.table(
    model_id = model_id,
    model_label = model_label,
    formula = formula_text,
    cr
  )
}

mean_sd <- function(x) {
  x <- x[
    is.finite(x)
  ]

  paste0(
    sprintf("%.2f", mean(x)),
    " (",
    sprintf("%.2f", sd(x)),
    ")"
  )
}

median_iqr <- function(x) {
  x <- x[
    is.finite(x)
  ]

  q <- quantile(
    x,
    c(0.25, 0.75),
    na.rm = TRUE
  )

  paste0(
    sprintf("%.1f", median(x)),
    " [",
    sprintf("%.1f", q[1]),
    ", ",
    sprintf("%.1f", q[2]),
    "]"
  )
}

# ----------------------------------------------------------------------
# 2. Read and join locked outcome + corrected official metadata
# ----------------------------------------------------------------------
locked <- fread(
  LOCKED_FILE,
  showProgress = FALSE
)

master <- fread(
  MASTER_FILE,
  showProgress = FALSE
)

required_locked <- c(
  "samples",
  "cohort",
  "module_score"
)

required_master <- c(
  "sample",
  "corrected_patient_id",
  "included_primary_analysis",
  "snrna_seq_chemistry",
  "age_at_resection_of_profiled_specimen",
  "sex",
  "other_mutation_status",
  "sample_prefix_stratum"
)

if (!all(required_locked %in% names(locked))) {
  stop(
    "Locked table lacks required columns.",
    call. = FALSE
  )
}

if (!all(required_master %in% names(master))) {
  stop(
    "E4 master table lacks required columns.",
    call. = FALSE
  )
}

meta <- master[
  included_primary_analysis == TRUE,
  ..required_master
]

setnames(
  meta,
  "sample",
  "samples"
)

d <- merge(
  locked,
  meta,
  by = "samples",
  all.x = TRUE,
  sort = FALSE
)

if (
  nrow(d) != 41L ||
  anyNA(d$corrected_patient_id) ||
  anyNA(d$snrna_seq_chemistry) ||
  anyNA(d$age_at_resection_of_profiled_specimen) ||
  anyNA(d$sex)
) {
  stop(
    "Expected complete 41-specimen metadata for patient ID, chemistry, age, and sex.",
    call. = FALSE
  )
}

if (
  uniqueN(d$corrected_patient_id) != 39L
) {
  stop(
    "Expected 39 corrected patient clusters.",
    call. = FALSE
  )
}

d[
  ,
  site_BM := as.integer(
    cohort ==
      "Brain_metastasis"
  )
]

d[
  ,
  chemistry_V2 := as.integer(
    snrna_seq_chemistry ==
      "5' V2"
  )
]

d[
  ,
  sex_M := as.integer(
    toupper(
      as.character(sex)
    ) ==
      "M"
  )
]

d[
  ,
  age_c := as.numeric(
    age_at_resection_of_profiled_specimen
  ) -
    mean(
      as.numeric(
        age_at_resection_of_profiled_specimen
      )
    )
]

# ----------------------------------------------------------------------
# 3. Descriptive covariate distribution
# ----------------------------------------------------------------------
chem_tab <- d[
  ,
  .N,
  by = .(
    cohort,
    snrna_seq_chemistry
  )
]

chem_tab[
  ,
  fraction_within_site :=
    N /
      sum(N),
  by = cohort
]

fwrite(
  chem_tab,
  file.path(
    OUT,
    "STEP_E4B_01_chemistry_by_site.tsv"
  ),
  sep = "\t"
)

sex_tab <- d[
  ,
  .N,
  by = .(
    cohort,
    sex
  )
]

sex_tab[
  ,
  fraction_within_site :=
    N /
      sum(N),
  by = cohort
]

age_summary <- d[
  ,
  .(
    n = .N,
    age_mean_sd = mean_sd(
      as.numeric(
        age_at_resection_of_profiled_specimen
      )
    ),
    age_median_IQR = median_iqr(
      as.numeric(
        age_at_resection_of_profiled_specimen
      )
    ),
    age_min = min(
      as.numeric(
        age_at_resection_of_profiled_specimen
      )
    ),
    age_max = max(
      as.numeric(
        age_at_resection_of_profiled_specimen
      )
    )
  ),
  by = cohort
]

prefix_tab <- d[
  ,
  .N,
  by = .(
    cohort,
    sample_prefix_stratum
  )
]

# Broad mutation descriptors are descriptive only.
d[
  ,
  mutation_text :=
    toupper(
      trimws(
        as.character(
          other_mutation_status
        )
      )
    )
]

d[
  ,
  KRAS_status_descriptive := fifelse(
    is.na(mutation_text) |
      !nzchar(mutation_text),
    "Unknown/not stated",
    fifelse(
      grepl(
        "KRAS WT",
        mutation_text,
        fixed = TRUE
      ),
      "KRAS WT stated",
      fifelse(
        grepl(
          "KRAS",
          mutation_text,
          fixed = TRUE
        ),
        "KRAS alteration stated",
        "Other driver/alteration stated"
      )
    )
  )
]

kras_tab <- d[
  ,
  .N,
  by = .(
    cohort,
    KRAS_status_descriptive
  )
]

baseline_long <- rbindlist(
  list(
    chem_tab[
      ,
      .(
        domain = "Documented technical",
        variable = "snRNA-seq chemistry",
        cohort,
        level = snrna_seq_chemistry,
        n = N,
        fraction = fraction_within_site
      )
    ],
    sex_tab[
      ,
      .(
        domain = "Clinical",
        variable = "Sex",
        cohort,
        level = as.character(sex),
        n = N,
        fraction = fraction_within_site
      )
    ],
    prefix_tab[
      ,
      .(
        domain = "Processing/genomic stratum",
        variable = "Sample prefix stratum",
        cohort,
        level = sample_prefix_stratum,
        n = N,
        fraction = N / sum(N)
      ),
      by = cohort
    ],
    kras_tab[
      ,
      .(
        domain = "Genomic descriptive",
        variable = "KRAS status from granular mutation text",
        cohort,
        level = KRAS_status_descriptive,
        n = N,
        fraction = N / sum(N)
      ),
      by = cohort
    ]
  ),
  fill = TRUE
)

fwrite(
  baseline_long,
  file.path(
    OUT,
    "STEP_E4B_02_baseline_categorical_distributions.tsv"
  ),
  sep = "\t"
)

fwrite(
  age_summary,
  file.path(
    OUT,
    "STEP_E4B_03_age_by_site.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 4. Prespecified patient-cluster CR1 models
# ----------------------------------------------------------------------
models <- rbindlist(
  list(
    fit_model(
      d,
      model_id = "M0",
      model_label = "Site only; corrected patient clustering",
      formula_text =
        "module_score ~ site_BM"
    ),
    fit_model(
      d,
      model_id = "M1",
      model_label = "Site + documented snRNA-seq chemistry",
      formula_text =
        "module_score ~ site_BM + chemistry_V2"
    ),
    fit_model(
      d,
      model_id = "M2",
      model_label = "Site + age at resection",
      formula_text =
        "module_score ~ site_BM + age_c"
    ),
    fit_model(
      d,
      model_id = "M3",
      model_label = "Site + sex",
      formula_text =
        "module_score ~ site_BM + sex_M"
    ),
    fit_model(
      d,
      model_id = "M4",
      model_label = "Site + chemistry + age + sex",
      formula_text =
        "module_score ~ site_BM + chemistry_V2 + age_c + sex_M"
    )
  ),
  fill = TRUE
)

fwrite(
  models,
  file.path(
    OUT,
    "STEP_E4B_04_patient_cluster_CR1_covariate_models.tsv"
  ),
  sep = "\t"
)

# ----------------------------------------------------------------------
# 5. Collinearity / design transparency
# ----------------------------------------------------------------------
chem_by_site <- d[
  ,
  .N,
  by = .(
    cohort,
    snrna_seq_chemistry
  )
]

age_by_site_numeric <- d[
  ,
  .(
    n = .N,
    mean_age = mean(
      as.numeric(
        age_at_resection_of_profiled_specimen
      )
    ),
    sd_age = sd(
      as.numeric(
        age_at_resection_of_profiled_specimen
      )
    ),
    median_age = median(
      as.numeric(
        age_at_resection_of_profiled_specimen
      )
    )
  ),
  by = cohort
]

design_note <- c(
  "KRAS/STK/PA/N sample prefixes were not entered as a pure technical batch covariate.",
  "The documented technical covariate is snRNA-seq chemistry (5' V1.1 vs 5' V2).",
  "Age and sex are treated as biological/clinical covariates, not technical covariates.",
  "Granular mutation text is summarized descriptively only because of heterogeneous categories and limited sample size.",
  "The M4 model is a multivariable sensitivity model; attenuation after adjustment should not be interpreted as proof of causal mediation or absence of a site-associated phenotype."
)

writeLines(
  design_note,
  file.path(
    OUT,
    "STEP_E4B_05_MODEL_INTERPRETATION_NOTES.txt"
  )
)

# ----------------------------------------------------------------------
# 6. Status and reporting wording
# ----------------------------------------------------------------------
m0 <- models[model_id == "M0"]
m1 <- models[model_id == "M1"]
m2 <- models[model_id == "M2"]
m4 <- models[model_id == "M4"]

interpretation <- paste0(
  "The site coefficient was ",
  sprintf("%.3f", m0$estimate_site_BM),
  " in the patient-clustered site-only model (two-sided P=",
  formatC(m0$two_sided_p_CR1, format = "g", digits = 3),
  "). After adjustment for documented snRNA-seq chemistry it was ",
  sprintf("%.3f", m1$estimate_site_BM),
  " (P=",
  formatC(m1$two_sided_p_CR1, format = "g", digits = 3),
  "); after age adjustment it was ",
  sprintf("%.3f", m2$estimate_site_BM),
  " (P=",
  formatC(m2$two_sided_p_CR1, format = "g", digits = 3),
  "); and in the chemistry+age+sex sensitivity model it was ",
  sprintf("%.3f", m4$estimate_site_BM),
  " (P=",
  formatC(m4$two_sided_p_CR1, format = "g", digits = 3),
  "). These models should be interpreted as sensitivity analyses in a small, ",
  "imbalanced observational cohort rather than as causal decomposition."
)

writeLines(
  c(
    "status=PASS_E4B_TECHNICAL_CLINICAL_COVARIATE_SENSITIVITY",
    "script_version=E4B_V1",
    "outcome=locked_fixed_five_gene_module_score",
    "inference=manual_CR1_patient_clustered_two_sided_t",
    "documented_technical_covariate=snrna_seq_chemistry",
    "clinical_covariates=age_at_resection_of_profiled_specimen;sex",
    "sample_prefix_stratum_treated_as_pure_technical_batch=FALSE",
    "granular_mutation_text_in_multivariable_model=FALSE",
    paste0(
      "M0_site_only_beta=",
      signif(m0$estimate_site_BM, 10)
    ),
    paste0(
      "M0_site_only_p=",
      signif(m0$two_sided_p_CR1, 10)
    ),
    paste0(
      "M1_site_plus_chemistry_beta=",
      signif(m1$estimate_site_BM, 10)
    ),
    paste0(
      "M1_site_plus_chemistry_p=",
      signif(m1$two_sided_p_CR1, 10)
    ),
    paste0(
      "M2_site_plus_age_beta=",
      signif(m2$estimate_site_BM, 10)
    ),
    paste0(
      "M2_site_plus_age_p=",
      signif(m2$two_sided_p_CR1, 10)
    ),
    paste0(
      "M4_site_plus_chemistry_age_sex_beta=",
      signif(m4$estimate_site_BM, 10)
    ),
    paste0(
      "M4_site_plus_chemistry_age_sex_p=",
      signif(m4$two_sided_p_CR1, 10)
    ),
    paste0(
      "reviewer_facing_interpretation=",
      interpretation
    ),
    "covariates_selected_from_outcome_results=FALSE",
    "primary_endpoint_replaced=FALSE"
  ),
  file.path(
    OUT,
    "STEP_E4B_COMPLETE.txt"
  )
)

cat(
  "\n============================================================\n",
  "STEP E4B: TECHNICAL + CLINICAL COVARIATE SENSITIVITY\n",
  "============================================================\n",
  "snRNA-seq chemistry by site:\n",
  sep = ""
)

print(
  chem_by_site
)

cat(
  "\nAge by site:\n"
)

print(
  age_by_site_numeric
)

cat(
  "\nPatient-cluster CR1 site coefficients:\n"
)

print(
  models[
    ,
    .(
      model_id,
      model_label,
      estimate_site_BM,
      CR1_SE,
      CI95_low_CR1,
      CI95_high_CR1,
      two_sided_p_CR1,
      n_observations,
      n_patient_clusters
    )
  ]
)

cat(
  "\nInterpretation:\n",
  interpretation,
  "\n\nStatus: PASS_E4B_TECHNICAL_CLINICAL_COVARIATE_SENSITIVITY\n\n",
  "Key output files:\n",
  file.path(
    OUT,
    "STEP_E4B_01_chemistry_by_site.tsv"
  ),
  "\n",
  file.path(
    OUT,
    "STEP_E4B_03_age_by_site.tsv"
  ),
  "\n",
  file.path(
    OUT,
    "STEP_E4B_04_patient_cluster_CR1_covariate_models.tsv"
  ),
  "\n",
  file.path(
    OUT,
    "STEP_E4B_COMPLETE.txt"
  ),
  "\n",
  sep = ""
)
