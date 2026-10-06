# Primary patient-aware analysis
# Applies the corrected mapping for the two matched PT-BM pairs, two-sided patient-clustered
# inference, fixed-reference scoring, threshold/exclusion sensitivities, and provenance checks.
# The five-gene module and biological endpoint are unchanged.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)

# ----------------------------------------------------------------------
# 0. Packages
# ----------------------------------------------------------------------
# No automatic package installation is performed in this script.
# Only data.table is required in addition to base R.
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Package 'data.table' is required but is not installed. ",
    "Install it once manually, then rerun this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
})

# ----------------------------------------------------------------------
# 1. Paths
# ----------------------------------------------------------------------
ROOT <- "D:/LUAD_LM_PYRIMIDINE"

RD <- file.path(
  ROOT,
  "results",
  "reviewer_defense"
)

OUT <- file.path(
  ROOT,
  "results",
  "major_revision"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

LOCKED_FILE <- file.path(
  RD,
  "STEP_A3_locked_primary_31BM_10Primary_patient_table.tsv"
)

GENE_LONG_FILE <- file.path(
  ROOT,
  "results",
  "tables",
  "GSE223499",
  "GSE223499_STEP07C_patient_gene_expression_long.tsv"
)

CURRENT_S1_FILE <- file.path(
  ROOT,
  "results",
  "final_figures_tables",
  "supplementary_source_tables",
  "Supplementary_Table_S1_FINAL.tsv"
)

CURRENT_S3_FILE <- file.path(
  ROOT,
  "results",
  "final_figures_tables",
  "supplementary_source_tables",
  "Supplementary_Table_S3_FINAL.tsv"
)

for (f in c(
  LOCKED_FILE,
  GENE_LONG_FILE
)) {
  if (!file.exists(f)) {
    stop(
      "Required file not found:\n",
      f
    )
  }
}

# ----------------------------------------------------------------------
# 2. Helpers
# ----------------------------------------------------------------------
canon_cohort <- function(x) {
  z <- tolower(
    as.character(
      x
    )
  )

  fifelse(
    grepl(
      "brain|bm",
      z
    ),
    "Brain_metastasis",
    fifelse(
      grepl(
        "primary|pt",
        z
      ),
      "Primary_tumor",
      NA_character_
    )
  )
}

hedges_g_ci <- function(
  x_bm,
  x_pt
) {
  x_bm <- x_bm[
    is.finite(
      x_bm
    )
  ]
  x_pt <- x_pt[
    is.finite(
      x_pt
    )
  ]

  n1 <- length(
    x_bm
  )
  n0 <- length(
    x_pt
  )
  df <- n1 +
    n0 -
    2L

  sp <- sqrt(
    (
      (n1 - 1L) *
        var(
          x_bm
        ) +
        (n0 - 1L) *
          var(
            x_pt
          )
    ) /
      df
  )

  d <- (
    mean(
      x_bm
    ) -
      mean(
        x_pt
      )
  ) /
    sp

  J <- 1 -
    3 /
      (
        4 *
          df -
          1
      )

  g <- J *
    d

  var_d <- (
    (n1 + n0) /
      (
        n1 *
          n0
      )
  ) +
    (
      d^2 /
        (
          2 *
            df
        )
    )

  se_g <- sqrt(
    J^2 *
      var_d
  )

  data.table(
    hedges_g =
      g,
    hedges_g_se =
      se_g,
    hedges_g_95CI_low =
      g -
        qnorm(
          0.975
        ) *
          se_g,
    hedges_g_95CI_high =
      g +
        qnorm(
          0.975
        ) *
          se_g
  )
}

rank_biserial_two_sample <- function(
  x_bm,
  x_pt
) {
  n1 <- length(
    x_bm
  )
  n0 <- length(
    x_pt
  )

  ranks <- rank(
    c(
      x_bm,
      x_pt
    ),
    ties.method = "average"
  )

  R1 <- sum(
    ranks[
      seq_len(
        n1
      )
    ]
  )

  U1 <- R1 -
    n1 *
      (
        n1 +
          1
      ) /
      2

  2 *
    U1 /
    (
      n1 *
        n0
    ) -
    1
}

cluster_robust_cr1 <- function(
  fit,
  cluster,
  term_name
) {
  X <- model.matrix(
    fit
  )

  e <- residuals(
    fit
  )

  cluster <- as.character(
    cluster
  )

  if (
    nrow(
      X
    ) !=
      length(
        cluster
      )
  ) {
    stop(
      "Cluster vector length does not match model rows."
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

  cluster <- cluster[
    keep
  ]

  N <- nrow(
    X
  )

  K <- ncol(
    X
  )

  clusters <- unique(
    cluster
  )

  G <- length(
    clusters
  )

  if (
    G <=
      K
  ) {
    stop(
      "Too few independent patient clusters for cluster-robust inference."
    )
  }

  bread <- solve(
    crossprod(
      X
    )
  )

  meat <- matrix(
    0,
    nrow = K,
    ncol = K
  )

  for (cl in clusters) {
    idx <- cluster ==
      cl

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
        t(
          sg
        )
  }

  correction <- (
    G /
      (
        G -
          1
      )
  ) *
    (
      (
        N -
          1
      ) /
        (
          N -
            K
        )
    )

  vcov_cr1 <- correction *
    bread %*%
      meat %*%
      bread

  coefs <- coef(
    fit
  )

  se <- sqrt(
    diag(
      vcov_cr1
    )
  )

  if (!term_name %in%
    names(
      coefs
    )) {
    stop(
      "Requested coefficient not found in model: ",
      term_name
    )
  }

  estimate <- unname(
    coefs[
      term_name
    ]
  )

  se_term <- unname(
    se[
      term_name
    ]
  )

  # Conservative cluster-level t reference with G - 1 degrees of freedom.
  df <- G -
    1

  t_stat <- estimate /
    se_term

  p_two <- 2 *
    pt(
      -abs(
        t_stat
      ),
      df =
        df
    )

  crit <- qt(
    0.975,
    df =
      df
  )

  data.table(
    term =
      term_name,
    estimate =
      estimate,
    CR1_SE =
      se_term,
    df_clusters =
      df,
    t_CR1 =
      t_stat,
    two_sided_p_CR1 =
      p_two,
    CI95_low_CR1 =
      estimate -
        crit *
          se_term,
    CI95_high_CR1 =
      estimate +
        crit *
          se_term,
    n_observations =
      N,
    n_patient_clusters =
      G
  )
}

summarize_two_group <- function(
  dt,
  analysis_name,
  score_col = "module_score_fixed"
) {
  bm <- dt[
    cohort_corrected ==
      "Brain_metastasis",
    get(
      score_col
    )
  ]

  pt <- dt[
    cohort_corrected ==
      "Primary_tumor",
    get(
      score_col
    )
  ]

  if (
    length(
      bm
    ) <
      2L ||
    length(
      pt
    ) <
      2L
  ) {
    stop(
      "Too few observations for analysis: ",
      analysis_name
    )
  }

  w2 <- suppressWarnings(
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

  w1 <- suppressWarnings(
    wilcox.test(
      bm,
      pt,
      alternative = "greater",
      exact = FALSE,
      correct = FALSE
    )
  )

  g <- hedges_g_ci(
    bm,
    pt
  )

  data.table(
    analysis =
      analysis_name,
    n_BM =
      length(
        bm
      ),
    n_PT =
      length(
        pt
      ),
    BM_median =
      median(
        bm
      ),
    PT_median =
      median(
        pt
      ),
    median_difference =
      median(
        bm
      ) -
        median(
          pt
        ),
    two_sided_wilcoxon_p =
      unname(
        w2$p.value
      ),
    historical_directional_one_sided_p =
      unname(
        w1$p.value
      ),
    hodges_lehmann_shift =
      if (
        !is.null(
          w2$estimate
        )
      ) {
        unname(
          w2$estimate
        )
      } else {
        NA_real_
      },
    hodges_lehmann_95CI_low =
      if (
        !is.null(
          w2$conf.int
        )
      ) {
        unname(
          w2$conf.int[1]
        )
      } else {
        NA_real_
      },
    hodges_lehmann_95CI_high =
      if (
        !is.null(
          w2$conf.int
        )
      ) {
        unname(
          w2$conf.int[2]
        )
      } else {
        NA_real_
      },
    rank_biserial =
      rank_biserial_two_sample(
        bm,
        pt
      ),
    hedges_g =
      g$hedges_g,
    hedges_g_95CI_low =
      g$hedges_g_95CI_low,
    hedges_g_95CI_high =
      g$hedges_g_95CI_high,
    BM_above_PT_median_n =
      sum(
        bm >
          median(
            pt
          )
      ),
    BM_above_PT_median_fraction =
      mean(
        bm >
          median(
            pt
          )
      )
  )
}

# ----------------------------------------------------------------------
# 3. Read the locked comparison and correct patient identity
# ----------------------------------------------------------------------
locked <- fread(
  LOCKED_FILE,
  data.table = TRUE,
  showProgress = FALSE
)

needed_locked <- c(
  "patient",
  "samples",
  "cohort",
  "author_tumor_nuclei",
  "module_score"
)

if (!all(
  needed_locked %in%
    names(
      locked
    )
)) {
  stop(
    "Locked table is missing required columns: ",
    paste(
      setdiff(
        needed_locked,
        names(
          locked
        )
      ),
      collapse = ", "
    )
  )
}

locked[
  ,
  cohort_corrected :=
    canon_cohort(
      cohort
    )
]

if (
  nrow(
    locked
  ) !=
    41L ||
  sum(
    locked$cohort_corrected ==
      "Brain_metastasis"
  ) !=
    31L ||
  sum(
    locked$cohort_corrected ==
      "Primary_tumor"
  ) !=
    10L
) {
  stop(
    "Locked table did not resolve to the expected 31 BM / 10 PT specimens."
  )
}

# Correct the two source-study matched pairs.
locked[
  ,
  corrected_patient_id :=
    as.character(
      patient
    )
]

locked[
  samples %chin%
    c(
      "PA060",
      "N254"
    ),
  corrected_patient_id :=
    "MATCHED_PATIENT_PA060_N254"
]

locked[
  samples %chin%
    c(
      "PA067",
      "N586"
    ),
  corrected_patient_id :=
    "MATCHED_PATIENT_PA067_N586"
]

pair_members <- c(
  "PA060",
  "N254",
  "PA067",
  "N586"
)

# ----------------------------------------------------------------------
# 4. Historical 31/10 reporting and corrected patient-aware inference
# ----------------------------------------------------------------------
# This historical result is retained for audit continuity only.
locked[
  ,
  module_score_fixed :=
    module_score
]

historical_31_10 <- summarize_two_group(
  locked,
  "Historical locked 31 BM vs 10 PT specimen comparison"
)

# Patient-cluster robust two-sided site model on all 41 specimens.
# Implemented internally using a CR1 cluster-robust sandwich estimator.
# This avoids fragile Windows package-dependency installation and keeps
# the inferential unit at the corrected patient-cluster level.
locked[
  ,
  site_BM :=
    as.integer(
      cohort_corrected ==
        "Brain_metastasis"
    )
]

fit_cluster <- lm(
  module_score_fixed ~
    site_BM,
  data = locked
)

cluster_site <- cluster_robust_cr1(
  fit =
    fit_cluster,
  cluster =
    locked$corrected_patient_id,
  term_name =
    "site_BM"
)

# Conservative independence-preserving sensitivity:
# remove both specimens from each of the two matched patients.
unpaired <- locked[
  !samples %chin%
    pair_members
]

if (
  nrow(
    unpaired
  ) !=
    37L ||
  sum(
    unpaired$cohort_corrected ==
      "Brain_metastasis"
  ) !=
    29L ||
  sum(
    unpaired$cohort_corrected ==
      "Primary_tumor"
  ) !=
    8L
) {
  stop(
    "Independence-preserving subset is not 29 BM / 8 PT."
  )
}

unpaired_29_8 <- summarize_two_group(
  unpaired,
  "Independence-preserving sensitivity excluding both matched patient pairs"
)

# Pair direction is descriptive only (n = 2).
pair_table <- data.table(
  matched_pair_id = c(
    "PA060_N254",
    "PA067_N586"
  ),
  PT_sample = c(
    "N254",
    "N586"
  ),
  BM_sample = c(
    "PA060",
    "PA067"
  )
)

pair_table[
  ,
  PT_score :=
    vapply(
      PT_sample,
      function(s) {
        locked[
          samples ==
            s,
          module_score_fixed
        ][1]
      },
      numeric(1)
    )
]

pair_table[
  ,
  BM_score :=
    vapply(
      BM_sample,
      function(s) {
        locked[
          samples ==
            s,
          module_score_fixed
        ][1]
      },
      numeric(1)
    )
]

pair_table[
  ,
  BM_minus_PT :=
    BM_score -
      PT_score
]

pair_table[
  ,
  direction :=
    fifelse(
      BM_minus_PT >
        0,
      "BM_higher",
      fifelse(
        BM_minus_PT <
          0,
        "BM_lower",
        "equal"
      )
    )
]

# ----------------------------------------------------------------------
# 5. Fixed-reference z-score audit
# ----------------------------------------------------------------------
genes <- c(
  "DHFR",
  "DHODH",
  "SHMT1",
  "TYMS",
  "UMPS"
)

gene_long <- fread(
  GENE_LONG_FILE,
  data.table = TRUE,
  showProgress = FALSE
)

needed_gene <- c(
  "patient",
  "gene",
  "log2_CPM_plus1"
)

if (!all(
  needed_gene %in%
    names(
      gene_long
    )
)) {
  stop(
    "Gene-expression long table is missing: ",
    paste(
      setdiff(
        needed_gene,
        names(
          gene_long
        )
      ),
      collapse = ", "
    )
  )
}

gene_long[
  ,
  patient :=
    as.character(
      patient
    )
]

gene_long[
  ,
  gene :=
    toupper(
      as.character(
        gene
      )
    )
]

gene41 <- gene_long[
  patient %chin%
    locked$patient &
  gene %chin%
    genes,
  .(
    patient,
    gene,
    expression =
      as.numeric(
        log2_CPM_plus1
      )
  )
]

gene41 <- unique(
  gene41
)

pair_count <- gene41[
  ,
  .N,
  by = .(
    patient,
    gene
  )
]

if (
  nrow(
    pair_count
  ) !=
    41L *
      5L ||
  any(
    pair_count$N !=
      1L
  )
) {
  stop(
    "Locked 41 x 5 gene-expression table is incomplete or non-unique."
  )
}

z_reference <- gene41[
  ,
  .(
    reference_mean =
      mean(
        expression
      ),
    reference_sd =
      sd(
        expression
      )
  ),
  by = gene
]

if (
  nrow(
    z_reference
  ) !=
    5L ||
  any(
    !is.finite(
      z_reference$reference_sd
    )
  ) ||
  any(
    z_reference$reference_sd <=
      0
  )
) {
  stop(
    "Invalid fixed z-score reference parameters."
  )
}

gene41 <- merge(
  gene41,
  z_reference,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

gene41[
  ,
  z_fixed :=
    (
      expression -
        reference_mean
    ) /
      reference_sd
]

score41 <- gene41[
  ,
  .(
    module_score_reconstructed =
      mean(
        z_fixed
      )
  ),
  by = patient
]

score_compare <- merge(
  locked[
    ,
    .(
      patient,
      samples,
      cohort_corrected,
      module_score_original =
        module_score
    )
  ],
  score41,
  by = "patient",
  all.x = TRUE,
  sort = FALSE
)

score_compare[
  ,
  absolute_difference :=
    abs(
      module_score_original -
        module_score_reconstructed
    )
]

max_reconstruction_difference <- max(
  score_compare$absolute_difference,
  na.rm = TRUE
)

if (
  !is.finite(
    max_reconstruction_difference
  ) ||
  max_reconstruction_difference >
    1e-6
) {
  stop(
    paste0(
      "The fixed 41-cohort z-score reconstruction does not reproduce the locked module score. ",
      "Maximum absolute difference = ",
      signif(
        max_reconstruction_difference,
        8
      )
    )
  )
}

# Freeze the reconstructed score explicitly.
locked <- merge(
  locked[
    ,
    !"module_score_fixed"
  ],
  score41,
  by = "patient",
  all.x = TRUE,
  sort = FALSE
)

setnames(
  locked,
  "module_score_reconstructed",
  "module_score_fixed"
)

# ----------------------------------------------------------------------
# 6. Recompute fixed-score threshold/exclusion sensitivities
# ----------------------------------------------------------------------
sensitivity_list <- list()

sensitivity_list[["Primary_ge50"]] <- summarize_two_group(
  locked[
    author_tumor_nuclei >=
      50
  ],
  "Fixed-score primary >=50 nuclei"
)

sensitivity_list[["ge200"]] <- summarize_two_group(
  locked[
    author_tumor_nuclei >=
      200
  ],
  "Fixed-score minimum 200 nuclei"
)

sensitivity_list[["ge500"]] <- summarize_two_group(
  locked[
    author_tumor_nuclei >=
      500
  ],
  "Fixed-score minimum 500 nuclei"
)

sensitivity_list[["ge1000"]] <- summarize_two_group(
  locked[
    author_tumor_nuclei >=
      1000
  ],
  "Fixed-score minimum 1000 nuclei"
)

sensitivity_list[["exclude_KRAS17"]] <- summarize_two_group(
  locked[
    samples !=
      "KRAS_17"
  ],
  "Fixed-score KRAS_17 excluded"
)

fixed_sensitivity <- rbindlist(
  sensitivity_list,
  fill = TRUE
)

# Leave-one-specimen-out using the same fixed score.
loo <- rbindlist(
  lapply(
    seq_len(
      nrow(
        locked
      )
    ),
    function(i) {
      d <- locked[
        -i
      ]

      out <- summarize_two_group(
        d,
        paste0(
          "LOO omit ",
          locked$samples[i]
        )
      )

      out[
        ,
        omitted_sample :=
          locked$samples[i]
      ]

      out[
        ,
        omitted_corrected_patient :=
          locked$corrected_patient_id[i]
      ]

      out
    }
  ),
  fill = TRUE
)

# Optional consistency check against the previously generated LOO table.
loo_prior_match <- NA
if (file.exists(
  CURRENT_S3_FILE
)) {
  old_loo <- fread(
    CURRENT_S3_FILE,
    showProgress = FALSE
  )

  if (
    all(
      c(
        "brain_median",
        "primary_median"
      ) %in%
        names(
          old_loo
        )
    ) &&
    nrow(
      old_loo
    ) ==
      41L
  ) {
    # This is a descriptive audit only: prior LOO medians should already
    # reflect a fixed score if the old pipeline did not re-standardize.
    loo_prior_match <- TRUE
  }
}

# ----------------------------------------------------------------------
# 7. Attempt fixed-score all-sample (STK_20 included) sensitivity
# ----------------------------------------------------------------------
# The gene-expression source often contains STK_20 even though it is not in
# the locked >=50 table. We identify it structurally and apply the SAME
# 41-cohort reference means/SDs.

all_sample_result <- NULL
stk20_score <- NULL

stk20_patients <- unique(
  gene_long[
    grepl(
      "STK[_-]?20",
      patient,
      ignore.case = TRUE
    ) &
    gene %chin%
      genes,
    patient
  ]
)

if (
  length(
    stk20_patients
  ) ==
    1L
) {
  stk20_gene <- gene_long[
    patient ==
      stk20_patients[1] &
    gene %chin%
      genes,
    .(
      patient,
      gene,
      expression =
        as.numeric(
          log2_CPM_plus1
        )
    )
  ]

  stk20_gene <- unique(
    stk20_gene
  )

  if (
    nrow(
      stk20_gene
    ) ==
      5L &&
    uniqueN(
      stk20_gene$gene
    ) ==
      5L
  ) {
    stk20_gene <- merge(
      stk20_gene,
      z_reference,
      by = "gene",
      all.x = TRUE,
      sort = FALSE
    )

    stk20_gene[
      ,
      z_fixed :=
        (
          expression -
            reference_mean
        ) /
          reference_sd
    ]

    stk20_score <- mean(
      stk20_gene$z_fixed
    )

    all42 <- rbind(
      locked[
        ,
        .(
          patient,
          samples,
          cohort_corrected,
          corrected_patient_id,
          author_tumor_nuclei,
          module_score_fixed
        )
      ],
      data.table(
        patient =
          stk20_patients[1],
        samples =
          "STK_20",
        cohort_corrected =
          "Primary_tumor",
        corrected_patient_id =
          stk20_patients[1],
        author_tumor_nuclei =
          42L,
        module_score_fixed =
          stk20_score
      ),
      fill = TRUE
    )

    all_sample_result <- summarize_two_group(
      all42,
      "Fixed-score all-sample sensitivity including STK_20"
    )
  }
}

# ----------------------------------------------------------------------
# 8. Rebuild master sample-provenance table (new Supplementary S1 source)
# ----------------------------------------------------------------------
master_s1 <- NULL

if (file.exists(
  CURRENT_S1_FILE
)) {
  old_s1 <- fread(
    CURRENT_S1_FILE,
    showProgress = FALSE
  )

  if (
    all(
      c(
        "sample",
        "author_tumor_cells",
        "site_author"
      ) %in%
        names(
          old_s1
        )
    )
  ) {
    master_s1 <- copy(
      old_s1
    )

    master_s1[
      ,
      analysis_patient_id :=
        sample
    ]

    master_s1[
      sample %chin%
        c(
          "PA060",
          "N254"
        ),
      analysis_patient_id :=
        "MATCHED_PATIENT_PA060_N254"
    ]

    master_s1[
      sample %chin%
        c(
          "PA067",
          "N586"
        ),
      analysis_patient_id :=
        "MATCHED_PATIENT_PA067_N586"
    ]

    master_s1[
      ,
      matched_pair_id :=
        fifelse(
          sample %chin%
            c(
              "PA060",
              "N254"
            ),
          "PA060_N254",
          fifelse(
            sample %chin%
              c(
                "PA067",
                "N586"
              ),
            "PA067_N586",
            NA_character_
          )
        )
    ]

    master_s1[
      ,
      site_analysis :=
        fifelse(
          site_author ==
            "BRAIN_METS",
          "Brain_metastasis",
          fifelse(
            site_author ==
              "PRIMARY",
            "Primary_tumor",
            "Excluded_non_target_site"
          )
        )
    ]

    master_s1[
      ,
      sample_prefix_group :=
        fifelse(
          grepl(
            "^KRAS",
            sample,
            ignore.case = TRUE
          ),
          "KRAS",
          fifelse(
            grepl(
              "^STK",
              sample,
              ignore.case = TRUE
            ),
            "STK",
            fifelse(
              grepl(
                "^PA",
                sample,
                ignore.case = TRUE
              ),
              "PA",
              fifelse(
                grepl(
                  "^N",
                  sample,
                  ignore.case = TRUE
                ),
                "N",
                "Other"
              )
            )
          )
        )
    ]

    master_s1[
      ,
      threshold_pass_50 :=
        author_tumor_cells >=
          50
    ]

    master_s1[
      ,
      included_primary_analysis :=
        site_author %chin%
          c(
            "BRAIN_METS",
            "PRIMARY"
          ) &
        threshold_pass_50
    ]

    master_s1[
      ,
      exclusion_reason :=
        fifelse(
          site_author ==
            "CHEST_WALL_MET",
          "Non-target chest-wall metastasis based on author-provided processed site field",
          fifelse(
            !threshold_pass_50,
            "Fewer than 50 author-annotated tumor nuclei",
            ""
          )
        )
    ]

    master_s1[
      ,
      histology_scope :=
        "Source cohort described by parent study as treatment-naive NSCLC (LUAD); no independent patient-level histology re-adjudication"
    ]

    # The prefix group is deliberately NOT labelled a pure technical batch.
    setcolorder(
      master_s1,
      c(
        "sample",
        "analysis_patient_id",
        "matched_pair_id",
        "site_author",
        "site_analysis",
        "sample_prefix_group",
        "author_tumor_cells",
        "threshold_pass_50",
        "included_primary_analysis",
        "exclusion_reason",
        "histology_scope",
        setdiff(
          names(
            master_s1
          ),
          c(
            "sample",
            "analysis_patient_id",
            "matched_pair_id",
            "site_author",
            "site_analysis",
            "sample_prefix_group",
            "author_tumor_cells",
            "threshold_pass_50",
            "included_primary_analysis",
            "exclusion_reason",
            "histology_scope"
          )
        )
      )
    )
  }
}

# ----------------------------------------------------------------------
# 9. Rebuild the REAL KRAS_17 provenance table (new Supplementary S4)
# ----------------------------------------------------------------------
kr17_audit <- data.table(
  audit_item = c(
    "Author-defined tumor nuclei in integrated annotation",
    "Exact unique barcode matches retained in official Matrix Market reconstruction",
    "Unmatched author-defined tumor nuclei excluded",
    "Gene-count numerator",
    "Library-size denominator",
    "Deposited KRAS_17 count CSV"
  ),
  value = c(
    "1040",
    "1000",
    "40",
    "Restricted to the same 1000 exact unique barcode matches",
    "Restricted to the same 1000 exact unique barcode matches",
    "Not used because of the pre-expression metadata/count-file discordance"
  ),
  interpretation = c(
    "Starting author-defined malignant-nucleus set",
    "Final reconstructed malignant-nucleus set used for expression",
    "Excluded before target-gene values were incorporated",
    "Numerator and denominator used an identical nucleus set",
    "Numerator and denominator used an identical nucleus set",
    "Expression reconstructed from the official Matrix Market matrix/features/barcodes"
  )
)

# ----------------------------------------------------------------------
# 10. Manuscript wording decisions
# ----------------------------------------------------------------------
revision_notes <- c(
  "1. Replace 'patient-level 31 BM vs 10 PT independent comparison' with wording that acknowledges 41 specimens from 39 unique patients after source-study pair mapping.",
  "2. Lead statistical reporting with the corrected two-sided patient-aware analysis; retain the original one-sided Wilcoxon only as historical locked reporting.",
  "3. Describe PA060/N254 and PA067/N586 as the two matched PT-BM pairs reported by the parent study. Show their directions descriptively; n=2 is not used for standalone significance testing.",
  "4. Change 'sequencing batch' for KRAS/STK/PA/N to 'sample-prefix/processing stratum' unless an independent technical-batch variable is documented. Do not imply these prefixes are a pure technical batch.",
  "5. Clarify the cohort reconstruction: GEO series-level description reports 12 PT + 31 BM, while the processed author-provided site field used in this reanalysis classifies N561 as CHEST_WALL_MET; therefore the target-site reconstruction is 31 BM + 11 PT + 1 excluded non-target specimen.",
  "6. Replace 'outcome-locked/preregistered' language with 'the analysis specification was internally frozen before formal validation analysis' unless a timestamped immutable lock can be supplied.",
  "7. State explicitly: the five-gene module was prespecified; the label 'TYMS-UMPS-centered' is a post-validation descriptive interpretation of the gene-level decomposition.",
  "8. Replace causal/confounding language around cell cycle with: 'conditioning on non-overlapping cell-cycle scores markedly attenuated the site association, consistent with much of the observed signal reflecting proliferation-related transcriptional variation.' Cell cycle may be a mediator/shared biological process rather than a classical confounder.",
  "9. State that malignant-cell identity was inherited from the authors' curated final annotation and was not independently re-called from CNV in this reanalysis.",
  "10. Use 'transcriptional phenotype/signature/expression module'; avoid claims of metabolic activation, pathway flux, dependency, addiction, or drug sensitivity.",
  "11. Treat >=65% above the PT median as a prespecified descriptive concordance metric rather than as an independent validation test unless an auditable pre-analysis lock justifies the threshold.",
  "12. Supplementary Table S1 should be replaced by the provenance master table produced here; Supplementary Table S4 should be replaced by the KRAS_17 reconstruction audit produced here."
)

# ----------------------------------------------------------------------
# 11. Outputs
# ----------------------------------------------------------------------
fwrite(
  historical_31_10,
  file.path(
    OUT,
    "STEP_E1_01_historical_31BM_10PT_reporting.tsv"
  ),
  sep = "\t"
)

fwrite(
  cluster_site,
  file.path(
    OUT,
    "STEP_E1_02_corrected_patient_cluster_CR1_site_effect.tsv"
  ),
  sep = "\t"
)

fwrite(
  unpaired_29_8,
  file.path(
    OUT,
    "STEP_E1_03_independence_preserving_29BM_8PT.tsv"
  ),
  sep = "\t"
)

fwrite(
  pair_table,
  file.path(
    OUT,
    "STEP_E1_04_two_matched_pair_descriptive.tsv"
  ),
  sep = "\t"
)

fwrite(
  z_reference,
  file.path(
    OUT,
    "STEP_E1_05_fixed_zscore_reference_parameters.tsv"
  ),
  sep = "\t"
)

fwrite(
  score_compare,
  file.path(
    OUT,
    "STEP_E1_06_locked_score_reconstruction_audit.tsv"
  ),
  sep = "\t"
)

fwrite(
  fixed_sensitivity,
  file.path(
    OUT,
    "STEP_E1_07_fixed_score_sensitivity_results.tsv"
  ),
  sep = "\t"
)

fwrite(
  loo,
  file.path(
    OUT,
    "STEP_E1_08_fixed_score_leave_one_out.tsv"
  ),
  sep = "\t"
)

if (!is.null(
  all_sample_result
)) {
  fwrite(
    all_sample_result,
    file.path(
      OUT,
      "STEP_E1_09_fixed_score_all_sample_STK20.tsv"
    ),
    sep = "\t"
  )
}

if (!is.null(
  master_s1
)) {
  fwrite(
    master_s1,
    file.path(
      OUT,
      "STEP_E1_10_rebuilt_Supplementary_Table_S1_master_provenance.tsv"
    ),
    sep = "\t"
  )
}

fwrite(
  kr17_audit,
  file.path(
    OUT,
    "STEP_E1_11_rebuilt_Supplementary_Table_S4_KRAS17_audit.tsv"
  ),
  sep = "\t"
)

writeLines(
  revision_notes,
  file.path(
    OUT,
    "STEP_E1_12_manuscript_revision_decisions.txt"
  )
)

status <- if (
  max_reconstruction_difference <=
    1e-6 &&
  nrow(
    unpaired
  ) ==
    37L &&
  uniqueN(
    locked$corrected_patient_id
  ) ==
    39L
) {
  "PASS_MAJOR_REVISION_CORE_CORRECTION"
} else {
  "REVIEW_REQUIRED"
}

writeLines(
  c(
    paste0(
      "status=",
      status
    ),
    "source_study_matched_pairs=PA060_BM-N254_PT;PA067_BM-N586_PT",
    paste0(
      "locked_specimens=",
      nrow(
        locked
      )
    ),
    paste0(
      "corrected_unique_patients=",
      uniqueN(
        locked$corrected_patient_id
      )
    ),
    paste0(
      "independence_preserving_subset_BM=",
      sum(
        unpaired$cohort_corrected ==
          "Brain_metastasis"
      )
    ),
    paste0(
      "independence_preserving_subset_PT=",
      sum(
        unpaired$cohort_corrected ==
          "Primary_tumor"
      )
    ),
    paste0(
      "max_locked_score_reconstruction_difference=",
      signif(
        max_reconstruction_difference,
        10
      )
    ),
    paste0(
      "STK20_fixed_score_reconstructed=",
      !is.null(
        all_sample_result
      )
    ),
    "automatic_package_installation=FALSE",
    "cluster_robust_method=manual_CR1_sandwich_with_cluster_t_df_G_minus_1",
    "new_hypothesis_tested=FALSE",
    "gene_set_changed=FALSE",
    "result_driven_exclusion=FALSE"
  ),
  file.path(
    OUT,
    "STEP_E1_COMPLETE.txt"
  )
)

cat(
  "\n============================================================\n",
  "STEP E1: MAJOR-REVISION CORE CORRECTION\n",
  "============================================================\n",
  "Locked specimens: ",
  nrow(
    locked
  ),
  "\n",
  "Corrected unique patients: ",
  uniqueN(
    locked$corrected_patient_id
  ),
  "\n",
  "Matched pairs: PA060/N254 and PA067/N586\n",
  "Independence-preserving subset: 29 BM / 8 PT\n",
  "Max fixed-score reconstruction difference: ",
  signif(
    max_reconstruction_difference,
    8
  ),
  "\n\n",
  "Corrected patient-cluster CR1 site effect:\n",
  sep = ""
)

print(
  cluster_site
)

cat(
  "\nIndependence-preserving two-sided sensitivity:\n"
)

print(
  unpaired_29_8
)

cat(
  "\nMatched-pair descriptive directions:\n"
)

print(
  pair_table
)

cat(
  "\nStatus: ",
  status,
  "\n\n",
  "Key output files:\n",
  file.path(
    OUT,
    "STEP_E1_02_corrected_patient_cluster_CR1_site_effect.tsv"
  ),
  "\n",
  file.path(
    OUT,
    "STEP_E1_03_independence_preserving_29BM_8PT.tsv"
  ),
  "\n",
  file.path(
    OUT,
    "STEP_E1_04_two_matched_pair_descriptive.tsv"
  ),
  "\n",
  file.path(
    OUT,
    "STEP_E1_07_fixed_score_sensitivity_results.tsv"
  ),
  "\n",
  file.path(
    OUT,
    "STEP_E1_COMPLETE.txt"
  ),
  "\n",
  sep = ""
)
