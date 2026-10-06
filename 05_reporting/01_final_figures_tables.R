# Final manuscript figures and supplementary figures
# Rebuilds the final figure set from completed analysis outputs. This script performs plotting,
# table assembly, and display-only transformations; it does not introduce a new inferential endpoint.
# Set ROOT to the local project directory before running.

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260818)

ROOT <- "D:/LUAD_LM_PYRIMIDINE"
MR <- file.path(ROOT, "results", "major_revision")
E8 <- file.path(MR, "E8_critical_reviewer_defense")
E9 <- file.path(MR, "E9A_target_panel_specificity")
V18 <- file.path(
  ROOT, "results", "final_figures_tables",
  "major_revision_figures_v18_complete_suite"
)
OUT <- file.path(
  ROOT, "results", "final_figures_tables",
  "major_revision_figures_v20_2_1_final_polish"
)
PANELS <- file.path(OUT, "vector_panels")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(PANELS, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "data.table", "ggplot2", "patchwork", "svglite", "scales"
)
miss <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(miss) > 0L) {
  stop(
    "Missing package(s): ", paste(miss, collapse = ", "),
    "\nInstall once with:\ninstall.packages(c(",
    paste(sprintf('"%s"', miss), collapse = ", "),
    "), dependencies = NA)"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(scales)
})

# ----------------------------------------------------------------------
# Visual system
# ----------------------------------------------------------------------
FONT <- "Arial"

COL_PRIMARY <- "#6F91B3"
COL_BM <- "#C96E70"
COL_TEAL <- "#4F9A91"
COL_GOLD <- "#B89436"
COL_GREEN <- "#6E9D7A"
COL_GREY <- "#B8BDC3"
COL_DARK <- "#242424"
COL_MID <- "#666666"
COL_LIGHT <- "#E4E7EA"

theme_v19 <- function(base_size = 7.0, margin_pt = 3) {
  theme_classic(
    base_size = base_size,
    base_family = FONT
  ) +
    theme(
      plot.title = element_text(
        size = 8.0, face = "plain",
        colour = "black", hjust = 0
      ),
      plot.subtitle = element_text(
        size = 6.1, colour = COL_MID,
        hjust = 0, lineheight = 0.96
      ),
      plot.caption = element_text(
        size = 5.5, colour = COL_MID,
        hjust = 0, lineheight = 0.95
      ),
      axis.title = element_text(
        size = 7.0, colour = "black"
      ),
      axis.text = element_text(
        size = 6.2, colour = "black"
      ),
      axis.line = element_line(
        linewidth = 0.32, colour = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.25, colour = "black"
      ),
      legend.title = element_text(
        size = 6.0, colour = "black"
      ),
      legend.text = element_text(
        size = 5.7, colour = COL_DARK
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 6.5, colour = "black"
      ),
      panel.background = element_rect(
        fill = "white", colour = NA
      ),
      plot.background = element_rect(
        fill = "white", colour = NA
      ),
      plot.margin = margin(
        margin_pt, margin_pt,
        margin_pt, margin_pt,
        unit = "pt"
      )
    )
}

tag_theme <- theme(
  plot.tag = element_text(
    family = FONT,
    face = "bold",
    size = 9,
    colour = "black"
  )
)

save_pub <- function(p, stub, width_mm, height_mm) {
  ggsave(
    file.path(OUT, paste0(stub, ".pdf")),
    p, width = width_mm, height = height_mm,
    units = "mm", device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggsave(
    file.path(OUT, paste0(stub, ".svg")),
    p, width = width_mm, height = height_mm,
    units = "mm", device = svglite::svglite,
    bg = "white"
  )
  ggsave(
    file.path(OUT, paste0(stub, ".tiff")),
    p, width = width_mm, height = height_mm,
    units = "mm", dpi = 600,
    compression = "lzw", bg = "white"
  )
  ggsave(
    file.path(OUT, paste0(stub, ".png")),
    p, width = width_mm, height = height_mm,
    units = "mm", dpi = 600, bg = "white"
  )
}

save_panel <- function(p, stub, width_mm, height_mm) {
  ggsave(
    file.path(PANELS, paste0(stub, ".pdf")),
    p, width = width_mm, height = height_mm,
    units = "mm", device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggsave(
    file.path(PANELS, paste0(stub, ".svg")),
    p, width = width_mm, height = height_mm,
    units = "mm", device = svglite::svglite,
    bg = "white"
  )
}

copy_accepted <- function(old_stub, new_stub) {
  for (ext in c("pdf", "svg", "tiff", "png")) {
    src <- file.path(V18, paste0(old_stub, ".", ext))
    dst <- file.path(OUT, paste0(new_stub, ".", ext))
    if (!file.exists(src)) {
      stop("Accepted V18 figure missing: ", src)
    }
    ok <- file.copy(src, dst, overwrite = TRUE)
    if (!ok) stop("Could not copy: ", src)
  }
}

fmt_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.0001) return(format(x, scientific = TRUE, digits = 2))
  if (x < 0.001) return(sprintf("%.3g", x))
  sprintf("%.3f", x)
}

# ----------------------------------------------------------------------
# Source files
# ----------------------------------------------------------------------
read_req <- function(path) {
  if (!file.exists(path)) stop("Required source missing: ", path)
  fread(path, showProgress = FALSE)
}

# Some legacy robustness artifacts predate the major_revision folder.
# Resolve those by exact basename across known project result roots rather than
# assuming that every source was copied into results/major_revision.
find_project_file <- function(basename_target) {
  roots <- c(
    file.path(ROOT, "results", "major_revision"),
    file.path(ROOT, "results", "reviewer_defense"),
    file.path(ROOT, "results")
  )

  direct_candidates <- file.path(roots, basename_target)
  direct_hits <- direct_candidates[file.exists(direct_candidates)]

  if (length(direct_hits) >= 1L) {
    return(direct_hits[1L])
  }

  recursive_hits <- character(0)

  for (rr in roots) {
    if (!dir.exists(rr)) next

    hh <- list.files(
      rr,
      recursive = TRUE,
      full.names = TRUE
    )

    hh <- hh[
      basename(hh) ==
        basename_target
    ]

    recursive_hits <- unique(
      c(
        recursive_hits,
        hh
      )
    )
  }

  if (length(recursive_hits) == 0L) {
    stop(
      "Required project artifact not found: ",
      basename_target,
      "\nSearched under:\n",
      paste(
        roots,
        collapse = "\n"
      )
    )
  }

  # Prefer major_revision if both copies exist, then reviewer_defense.
  ord <- order(
    !grepl(
      "[/\\\\]major_revision[/\\\\]",
      recursive_hits
    ),
    !grepl(
      "[/\\\\]reviewer_defense[/\\\\]",
      recursive_hits
    ),
    nchar(recursive_hits)
  )

  recursive_hits[ord][1L]
}

read_project_artifact <- function(basename_target) {
  f <- find_project_file(
    basename_target
  )

  cat(
    "Resolved project artifact: ",
    f,
    "\n",
    sep = ""
  )

  fread(
    f,
    showProgress = FALSE
  )
}

e8_ident <- read_req(
  file.path(E8, "STEP_E8_07_site_identification_CR2_models.tsv")
)
e8_loo_patient <- read_req(
  file.path(E8, "STEP_E8_09_leave_one_patient_out.tsv")
)
e8_loo_prefix <- read_req(
  file.path(E8, "STEP_E8_10_leave_one_prefix_out.tsv")
)
e8_prolif <- read_req(
  file.path(E8, "STEP_E8_12_module_score_conditioned_on_cycling_CR2.tsv")
)
e8_cycle <- read_req(
  file.path(E8, "STEP_E8_13_cycling_fraction_covariate_CR2.tsv")
)
e8_cycle_support <- read_req(
  file.path(E8, "STEP_E8_15_cycling_fraction_common_support_CR2.tsv")
)
e8_delta <- read_req(
  file.path(E8, "STEP_E8_17_delta_beta_cluster_bootstrap.tsv")
)
e8_spatial_bh <- read_req(
  file.path(E8, "STEP_E8_18_spatial_toroidal_shift_BH_FDR.tsv")
)

e9_gene <- read_req(
  file.path(E9, "STEP_E9A_02_gene_cycle_and_site_metrics.tsv")
)
e9_comp <- read_req(
  file.path(E9, "STEP_E9A_03_fixed_comparator_panels_CR2.tsv")
)
e9_nested <- read_req(
  file.path(E9, "STEP_E9A_04_nested_pyrimidine_panels_CR2.tsv")
)
e9_null_dist <- read_req(
  file.path(E9, "STEP_E9A_06_proliferation_matched_random_set_null_distribution.tsv")
)
e9_null <- read_req(
  file.path(E9, "STEP_E9A_07_proliferation_matched_random_set_null_summary.tsv")
)
e9_balance <- read_req(
  file.path(E9, "STEP_E9A_08_matching_feature_balance.tsv")
)

# Existing spatial data for supplementary figure.
# STEP_C5Bv2 predates the major_revision directory and may reside under
# results/reviewer_defense. Resolve by exact basename.
spatial_grid <- read_project_artifact(
  "STEP_C5Bv2_01_grid_level_scores.tsv"
)

e5_shift <- read_project_artifact(
  "STEP_E5_04_toroidal_shift_null_summary.tsv"
)

e5_rep <- read_project_artifact(
  "STEP_E5_07_representative_map_selection.tsv"
)

# Additional compact source-of-truth tables used only for visualization.
# No new statistical tests are run in this figure script.
e8_smd <- read_project_artifact(
  "STEP_E8_05_covariate_standardized_mean_differences.tsv"
)
e8_analysis41 <- read_project_artifact(
  "STEP_E8_01_analysis_dataset_41.tsv"
)

e2b_comp <- read_project_artifact(
  "STEP_E2B_04_specimen_cellcycle_composition.tsv"
)
e2b_stats <- read_project_artifact(
  "STEP_E2B_07_cellcycle_and_G0G1_module_statistics.tsv"
)
e2b_g1 <- read_project_artifact(
  "STEP_E2B_06_G0G1_like_module_score.tsv"
)
e2b_g1_cr1 <- read_project_artifact(
  "STEP_E2B_08_G0G1_patient_cluster_CR1.tsv"
)
e2b_g1_gene <- read_project_artifact(
  "STEP_E2B_09_G0G1_gene_level_statistics.tsv"
)

e9_neighbors <- read_project_artifact(
  "STEP_E9A_05_proliferation_matched_candidate_neighborhoods.tsv"
)

e5_dup <- read_project_artifact(
  "STEP_E5_01_duplicate_payload_sensitivity.tsv"
)
e5_epi_infer <- read_project_artifact(
  "STEP_E5_03_epithelial_enriched_proxy_inference.tsv"
)
e5_scale <- read_project_artifact(
  "STEP_E5_06_coordinate_scale_transparency.tsv"
)

# ----------------------------------------------------------------------
# FIGURE 1 — cohort, patient-aware endpoint and nonredundant robustness
# ----------------------------------------------------------------------

primary_locked <- read_project_artifact(
  "STEP_A3_locked_primary_31BM_10Primary_patient_table.tsv"
)
e1_cluster_v193 <- read_project_artifact(
  "STEP_E1_02_corrected_patient_cluster_CR1_site_effect.tsv"
)
e1_pairs_v193 <- read_project_artifact(
  "STEP_E1_04_two_matched_pair_descriptive.tsv"
)
e1_sens_v193 <- read_project_artifact(
  "STEP_E1_07_fixed_score_sensitivity_results.tsv"
)
e1_stk20_v193 <- read_project_artifact(
  "STEP_E1_09_fixed_score_all_sample_STK20.tsv"
)
e1_ind_v193 <- read_project_artifact(
  "STEP_E1_03_independence_preserving_29BM_8PT.tsv"
)
flow_v193 <- read_project_artifact(
  "STEP_E4_05_COHORT_FLOW.tsv"
)

canon_group_v193 <- function(x) {
  z <- tolower(as.character(x))
  fifelse(
    grepl("brain|bm", z),
    "Brain metastasis",
    "Primary tumor"
  )
}

if (nrow(primary_locked) != 41L) {
  stop("V20 expected 41 rows in the locked primary table.")
}

primary_locked[, group := factor(
  canon_group_v193(cohort),
  levels = c("Primary tumor", "Brain metastasis")
)]
primary_locked[, sample := as.character(samples)]

primary_map_v193 <- primary_locked[, .(
  sample,
  group,
  module_score = as.numeric(module_score)
)]

get_flow_n_v193 <- function(stage_label) {
  z <- flow_v193[stage == stage_label, n]
  if (length(z) != 1L) {
    stop("Could not resolve cohort-flow stage: ", stage_label)
  }
  as.integer(z)
}

n_all <- get_flow_n_v193("Retrieved processed specimens")
n_bm_proc <- get_flow_n_v193("Processed site = BRAIN_METS")
n_pt_proc <- get_flow_n_v193("Processed site = PRIMARY")
n_cw <- get_flow_n_v193("Processed site = CHEST_WALL_MET")
n_locked <- get_flow_n_v193("Locked primary specimens total")
n_patients <- get_flow_n_v193(
  "Locked primary unique patients after matched-pair correction"
)

flow_box <- data.table(
  y = c(4, 3, 2, 1),
  label = c(
    paste0("Processed author dataset\n", n_all, " specimens"),
    paste0(
      "Processed sites\n",
      n_bm_proc, " BM  |  ",
      n_pt_proc, " PT  |  ",
      n_cw, " chest-wall metastasis"
    ),
    "Target-site comparison\n31 BM + 11 PT specimens",
    paste0(
      "Locked primary analysis\n",
      n_locked, " specimens / ", n_patients,
      " patients\n31 BM + 10 PT; ≥50 malignant nuclei"
    )
  ),
  fill = c("#F3F4F5", "#ECEFF2", "#E8EEF4", "#E6F0ED")
)

callout_box <- data.table(
  ymin = c(2.66, 1.25),
  ymax = c(3.34, 1.93),
  label = c(
    "N561\nCHEST_WALL_MET\nOutside PT-vs-BM target-site comparison",
    "STK_20\n42 malignant nuclei\nIncluded only in no-threshold sensitivity"
  )
)

p1a <- ggplot() +
  geom_rect(
    data = flow_box,
    aes(
      xmin = 0, xmax = 6.35,
      ymin = y - 0.33, ymax = y + 0.33,
      fill = fill
    ),
    colour = "#6F6F6F",
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  geom_text(
    data = flow_box,
    aes(x = 3.175, y = y, label = label),
    family = FONT,
    size = 1.90,
    lineheight = 0.95,
    colour = COL_DARK
  ) +
  geom_segment(
    data = data.table(
      y1 = c(3.65, 2.65, 1.65),
      y2 = c(3.35, 2.35, 1.35)
    ),
    aes(x = 3.175, xend = 3.175, y = y1, yend = y2),
    arrow = grid::arrow(length = grid::unit(1.45, "mm")),
    linewidth = 0.35,
    colour = COL_DARK
  ) +
  geom_rect(
    data = callout_box,
    aes(
      xmin = 7.15, xmax = 11.85,
      ymin = ymin, ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "#FAFAFA",
    colour = "#B4B4B4",
    linewidth = 0.30
  ) +
  geom_text(
    data = callout_box,
    aes(
      x = 7.40,
      y = (ymin + ymax) / 2,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    family = FONT,
    size = 1.30,
    lineheight = 0.94,
    colour = COL_MID
  ) +
  annotate(
    "segment",
    x = 6.35, xend = 7.15,
    y = 3, yend = 3,
    linewidth = 0.30,
    colour = "#A8A8A8"
  ) +
  annotate(
    "segment",
    x = 6.35, xend = 7.15,
    y = 1.59, yend = 1.59,
    linewidth = 0.30,
    colour = "#A8A8A8"
  ) +
  coord_cartesian(
    xlim = c(-0.10, 12.05),
    ylim = c(0.55, 4.45),
    clip = "off"
  ) +
  theme_void(base_family = FONT) +
  labs(title = "Cohort reconstruction and locked analysis set") +
  theme(
    plot.title = element_text(
      family = FONT,
      size = 7.0,
      colour = "black",
      hjust = 0,
      margin = margin(b = 3, unit = "pt")
    ),
    plot.margin = margin(4, 4, 2, 2, unit = "pt")
  )

pairs <- copy(e1_pairs_v193)

pair_samples <- unique(c(pairs$PT_sample, pairs$BM_sample))
pair_draw <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  data.table(
    matched_pair_id = pairs$matched_pair_id[i],
    group = factor(
      c("Primary tumor", "Brain metastasis"),
      levels = c("Primary tumor", "Brain metastasis")
    ),
    x = c(1, 2),
    score = c(pairs$PT_score[i], pairs$BM_score[i])
  )
}))

score_plot <- copy(primary_map_v193)
score_plot[, x := fifelse(group == "Primary tumor", 1, 2)]

p1b <- ggplot(score_plot, aes(x = x, y = module_score)) +
  geom_violin(
    aes(fill = group, colour = group),
    width = 0.72,
    alpha = 0.07,
    linewidth = 0.35,
    trim = FALSE
  ) +
  geom_boxplot(
    aes(fill = group, colour = group),
    width = 0.30,
    alpha = 0.14,
    outlier.shape = NA,
    linewidth = 0.40
  ) +
  geom_point(
    data = score_plot[!sample %chin% pair_samples],
    aes(colour = group),
    position = position_jitter(
      width = 0.055, height = 0, seed = 15
    ),
    size = 0.95,
    alpha = 0.72
  ) +
  geom_line(
    data = pair_draw,
    aes(x = x, y = score, group = matched_pair_id),
    inherit.aes = FALSE,
    colour = "#8C8C8C",
    linewidth = 0.45
  ) +
  geom_point(
    data = pair_draw,
    aes(x = x, y = score, fill = group),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.75,
    stroke = 0.45,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "Primary tumor" = COL_PRIMARY,
      "Brain metastasis" = COL_BM
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      "Primary tumor" = COL_PRIMARY,
      "Brain metastasis" = COL_BM
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c(
      "Primary tumor\nn = 10 specimens",
      "Brain metastasis\nn = 31 specimens"
    ),
    limits = c(0.62, 2.38)
  ) +
  annotate(
    "text",
    x = 0.70,
    y = max(score_plot$module_score) +
      0.17 * diff(range(score_plot$module_score)),
    hjust = 0,
    label = paste0(
      "Patient-cluster CR1: β = ",
      sprintf("%.2f", e1_cluster_v193$estimate[1]),
      " (95% CI ",
      sprintf("%.2f", e1_cluster_v193$CI95_low_CR1[1]),
      "–",
      sprintf("%.2f", e1_cluster_v193$CI95_high_CR1[1]),
      ")\nTwo-sided P = ",
      fmt_p(e1_cluster_v193$two_sided_p_CR1[1])
    ),
    family = FONT,
    size = 1.62,
    lineheight = 0.95,
    colour = COL_DARK
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.06, 0.24))
  ) +
  labs(
    x = NULL,
    y = "Fixed five-gene score",
    title = "Primary patient-aware comparison"
  ) +
  theme_v19()

pair_long <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  data.table(
    pair = pairs$matched_pair_id[i],
    site = factor(c("PT", "BM"), levels = c("PT", "BM")),
    score = c(pairs$PT_score[i], pairs$BM_score[i]),
    delta = pairs$BM_minus_PT[i]
  )
}))
pair_delta <- unique(pair_long[, .(pair, delta)])

p1c <- ggplot(pair_long, aes(x = score, y = pair)) +
  geom_line(
    aes(group = pair),
    colour = "#8F8F8F",
    linewidth = 0.60
  ) +
  geom_point(
    aes(fill = site),
    shape = 21,
    size = 2.40,
    stroke = 0.48,
    colour = COL_DARK
  ) +
  geom_text(
    data = pair_delta,
    aes(
      x = Inf,
      y = pair,
      label = paste0("BM − PT = ", sprintf("%.3f", delta))
    ),
    inherit.aes = FALSE,
    hjust = 1.02,
    family = FONT,
    size = 1.55,
    colour = "#333333"
  ) +
  scale_fill_manual(
    values = c("PT" = COL_PRIMARY, "BM" = COL_BM),
    name = NULL
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.12, 0.45))
  ) +
  labs(
    x = "Fixed five-gene score",
    y = NULL,
    title = "Two matched PT–BM pairs",
    subtitle = "Descriptive only; n = 2 pairs"
  ) +
  theme_v19() +
  theme(legend.position = "top")

get_e1_sens_v193 <- function(pattern) {
  z <- e1_sens_v193[
    grepl(pattern, analysis, ignore.case = TRUE, perl = TRUE)
  ]
  if (nrow(z) != 1L) {
    stop("Could not uniquely resolve E1 sensitivity: ", pattern)
  }
  z
}

s_primary <- get_e1_sens_v193("primary.*50")
s_500 <- get_e1_sens_v193("500")
s_1000 <- get_e1_sens_v193("1000")
s_k17 <- get_e1_sens_v193("KRAS_17")

rob1 <- rbindlist(list(
  data.table(
    label = "Primary ≥50 nuclei",
    hedges_g = s_primary$hedges_g,
    low = s_primary$hedges_g_95CI_low,
    high = s_primary$hedges_g_95CI_high
  ),
  data.table(
    label = "No threshold + STK_20",
    hedges_g = e1_stk20_v193$hedges_g[1],
    low = e1_stk20_v193$hedges_g_95CI_low[1],
    high = e1_stk20_v193$hedges_g_95CI_high[1]
  ),
  data.table(
    label = "KRAS_17 excluded",
    hedges_g = s_k17$hedges_g,
    low = s_k17$hedges_g_95CI_low,
    high = s_k17$hedges_g_95CI_high
  ),
  data.table(
    label = "Matched pairs fully excluded",
    hedges_g = e1_ind_v193$hedges_g[1],
    low = e1_ind_v193$hedges_g_95CI_low[1],
    high = e1_ind_v193$hedges_g_95CI_high[1]
  ),
  data.table(
    label = "≥500 nuclei",
    hedges_g = s_500$hedges_g,
    low = s_500$hedges_g_95CI_low,
    high = s_500$hedges_g_95CI_high
  ),
  data.table(
    label = "≥1000 nuclei",
    hedges_g = s_1000$hedges_g,
    low = s_1000$hedges_g_95CI_low,
    high = s_1000$hedges_g_95CI_high
  )
))

rob1[, label := factor(
  label,
  levels = rev(c(
    "Primary ≥50 nuclei",
    "No threshold + STK_20",
    "KRAS_17 excluded",
    "Matched pairs fully excluded",
    "≥500 nuclei",
    "≥1000 nuclei"
  ))
)]

p1d <- ggplot(rob1, aes(y = label, x = hedges_g)) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.28,
    colour = "#A5A5A5"
  ) +
  geom_segment(
    aes(x = low, xend = high, y = label, yend = label),
    linewidth = 0.55,
    colour = COL_DARK
  ) +
  geom_point(
    shape = 21,
    size = 2.30,
    fill = COL_TEAL,
    colour = COL_DARK,
    stroke = 0.45
  ) +
  labs(
    x = "Hedges g (95% CI)",
    y = NULL,
    title = "Fixed-score robustness"
  ) +
  theme_v19()

top_row1 <- p1a | p1b +
  plot_layout(widths = c(1.24, 0.76))
bottom_row1 <- p1c | p1d +
  plot_layout(widths = c(1.0, 1.0))

fig1 <- top_row1 / bottom_row1 +
  plot_layout(heights = c(1.10, 0.90)) +
  plot_annotation(tag_levels = "A") &
  tag_theme

save_panel(p1a, "Figure1A_cohort_flow_V20_2", 89, 72)
save_panel(p1b, "Figure1B_primary_score_V20_2", 89, 72)
save_panel(p1c, "Figure1C_matched_pairs_V20_2", 89, 52)
save_panel(p1d, "Figure1D_nonredundant_robustness_V20_2", 89, 52)
save_pub(
  fig1,
  "Figure1_MAJOR_REVISION_V20_2",
  183, 126
)

# ----------------------------------------------------------------------
# FIGURE 2 — SPECIMEN/GENE ARCHITECTURE + TWO ORTHOGONAL EFFECT ESTIMATORS
# ----------------------------------------------------------------------

gene_file_v20 <- find_project_file(
  "GSE223499_STEP07C_patient_gene_expression_long.tsv"
)
e3_nb_v20 <- read_project_artifact(
  "STEP_E3_02_count_based_NB_gene_statistics.tsv"
)

gene_long <- fread(gene_file_v20, showProgress = FALSE)
fixed_genes <- c("DHFR", "DHODH", "SHMT1", "TYMS", "UMPS")

gene_long[, patient := as.character(patient)]
gene_long[, gene := toupper(as.character(gene))]
gene_long[, log2_CPM_plus1 := as.numeric(log2_CPM_plus1)]
gene_long <- gene_long[gene %chin% fixed_genes]

id_candidates <- list()
if ("patient" %in% names(primary_locked)) {
  id_candidates[["patient"]] <- as.character(primary_locked$patient)
}
id_candidates[["samples"]] <- as.character(primary_locked$samples)

expr_ids <- unique(gene_long$patient)
overlap <- vapply(
  id_candidates,
  function(z) sum(unique(z) %chin% expr_ids),
  integer(1)
)
best_name <- names(overlap)[which.max(overlap)]
if (max(overlap) != 41L) {
  stop("V20 Figure 2 could not resolve exactly 41 locked specimen identifiers.")
}

locked_ids <- id_candidates[[best_name]]
lock_map <- data.table(
  patient = locked_ids,
  group = as.character(primary_locked$group),
  module_score = as.numeric(primary_locked$module_score)
)

gene_long <- unique(
  gene_long[
    patient %chin% lock_map$patient,
    .(patient, gene, log2_CPM_plus1)
  ]
)
if (nrow(gene_long) != 41L * 5L) {
  stop("V20 Figure 2 failed the 41 x 5 audit.")
}

gene_long <- merge(
  gene_long,
  lock_map,
  by = "patient",
  all.x = TRUE,
  sort = FALSE
)
gene_long[, group := factor(
  group,
  levels = c("Primary tumor", "Brain metastasis")
)]
gene_long[, gene := factor(gene, levels = fixed_genes)]

hedges_g_v20 <- function(x, y) {
  x <- as.numeric(x[is.finite(x)])
  y <- as.numeric(y[is.finite(y)])
  nx <- length(x)
  ny <- length(y)
  df <- nx + ny - 2L
  sp <- sqrt(
    ((nx - 1) * var(x) + (ny - 1) * var(y)) / df
  )
  d <- (mean(x) - mean(y)) / sp
  J <- 1 - 3 / (4 * df - 1)
  J * d
}

hedges_ci_v20 <- function(g, n1, n0) {
  N <- n1 + n0
  se <- sqrt(
    N / (n1 * n0) +
      g^2 / (2 * (N - 2))
  )
  c(low = g - 1.96 * se, high = g + 1.96 * se)
}

gene_effect <- rbindlist(lapply(fixed_genes, function(gg) {
  d <- gene_long[gene == gg]
  bm <- d[group == "Brain metastasis", log2_CPM_plus1]
  pt <- d[group == "Primary tumor", log2_CPM_plus1]
  g <- hedges_g_v20(bm, pt)
  ci <- hedges_ci_v20(g, length(bm), length(pt))
  data.table(
    gene = gg,
    hedges_g = g,
    low = ci["low"],
    high = ci["high"]
  )
}))

e3main <- e3_nb_v20[
  grepl("40 technically concordant", analysis, fixed = TRUE)
]
if (nrow(e3main) != 5L) {
  stop("V20 expected five E3 main gene rows.")
}

# A. Specimen-by-gene heat map. Order is descriptive and based on the
# already-computed fixed score within site; no new inferential operation.
heat_order <- unique(
  lock_map[
    order(
      factor(
        group,
        levels = c("Primary tumor", "Brain metastasis")
      ),
      module_score
    ),
    patient
  ]
)
heat_index <- data.table(
  patient = heat_order,
  x = seq_along(heat_order)
)

heat_dt <- merge(
  copy(gene_long),
  heat_index,
  by = "patient",
  all.x = TRUE
)
heat_dt[
  ,
  z_expr := {
    ss <- sd(log2_CPM_plus1)
    if (!is.finite(ss) || ss == 0) {
      rep(0, .N)
    } else {
      (log2_CPM_plus1 - mean(log2_CPM_plus1)) / ss
    }
  },
  by = gene
]
heat_dt[, gene_label := as.character(gene)]
heat_dt[, gene_y := match(
  gene_label,
  rev(fixed_genes)
)]

p2a <- ggplot(
  heat_dt,
  aes(
    x = x,
    y = gene_y,
    fill = z_expr
  )
) +
  geom_tile(
    width = 0.94,
    height = 0.86,
    colour = "white",
    linewidth = 0.08
  ) +
  geom_vline(
    xintercept = 10.5,
    linewidth = 0.45,
    colour = "white"
  ) +
  annotate(
    "text",
    x = 5.5,
    y = 5.67,
    label = "Primary tumor · 10 specimens",
    family = FONT,
    size = 1.55,
    colour = COL_PRIMARY
  ) +
  annotate(
    "text",
    x = 26,
    y = 5.67,
    label = "Brain metastasis · 31 specimens",
    family = FONT,
    size = 1.55,
    colour = COL_BM
  ) +
  scale_fill_gradient2(
    low = COL_PRIMARY,
    mid = "white",
    high = COL_BM,
    midpoint = 0,
    limits = c(
      -max(abs(heat_dt$z_expr), na.rm = TRUE),
      max(abs(heat_dt$z_expr), na.rm = TRUE)
    ),
    oob = scales::squish,
    name = "Within-gene\nz score"
  ) +
  scale_x_continuous(
    breaks = NULL,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_continuous(
    breaks = 1:5,
    labels = rev(fixed_genes),
    limits = c(0.5, 5.95),
    expand = c(0, 0)
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Specimen-level architecture of the fixed five-gene module",
    subtitle = "Specimens are ordered within site by the fixed five-gene score"
  ) +
  theme_v19() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.justification = "left"
  )

# B. Distribution view.
p2b <- ggplot(
  gene_long,
  aes(
    x = group,
    y = log2_CPM_plus1,
    fill = group,
    colour = group
  )
) +
  geom_violin(
    width = 0.82,
    alpha = 0.07,
    linewidth = 0.32,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.30,
    alpha = 0.13,
    outlier.shape = NA,
    linewidth = 0.34
  ) +
  geom_point(
    position = position_jitter(
      width = 0.055,
      height = 0,
      seed = 31
    ),
    size = 0.62,
    alpha = 0.55
  ) +
  facet_wrap(
    ~ gene,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Primary tumor" = COL_PRIMARY,
      "Brain metastasis" = COL_BM
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      "Primary tumor" = COL_PRIMARY,
      "Brain metastasis" = COL_BM
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = c(
      "Primary tumor" = "PT",
      "Brain metastasis" = "BM"
    )
  ) +
  labs(
    x = NULL,
    y = "log2(CPM + 1)",
    title = "Expression distributions"
  ) +
  theme_v19(margin_pt = 1.5)

# C. Standardized expression effects.
gene_effect[, gene := factor(
  gene,
  levels = rev(fixed_genes)
)]
gene_effect[, focus := as.character(gene) %chin% c("TYMS", "UMPS")]

p2c <- ggplot(
  gene_effect,
  aes(y = gene, x = hedges_g)
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.28,
    colour = "#A5A5A5"
  ) +
  geom_segment(
    aes(
      x = low,
      xend = high,
      y = gene,
      yend = gene
    ),
    linewidth = 0.58,
    colour = COL_DARK
  ) +
  geom_point(
    aes(fill = focus),
    shape = 21,
    size = 2.25,
    stroke = 0.44,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = COL_GREY,
      "TRUE" = COL_BM
    ),
    guide = "none"
  ) +
  labs(
    x = "Hedges g (95% CI)",
    y = NULL,
    title = "Standardized expression effects"
  ) +
  theme_v19()

# D. Count-based NB effects.
count_effect <- e3main[, .(
  gene,
  ratio = CR1_ratio_BM_vs_PT,
  low = CR1_CI95_low_ratio,
  high = CR1_CI95_high_ratio,
  FDR = BH_FDR_CR1_across_5_genes
)]
count_effect[, gene := factor(
  gene,
  levels = rev(fixed_genes)
)]
count_effect[, focus := as.character(gene) %chin% c("TYMS", "UMPS")]
count_effect[, log2_ratio := log2(ratio)]
count_effect[, log2_low := log2(low)]
count_effect[, log2_high := log2(high)]

p2d <- ggplot(
  count_effect,
  aes(y = gene, x = log2_ratio)
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.28,
    colour = "#A5A5A5"
  ) +
  geom_segment(
    aes(
      x = log2_low,
      xend = log2_high,
      y = gene,
      yend = gene
    ),
    linewidth = 0.58,
    colour = COL_DARK
  ) +
  geom_point(
    aes(fill = focus),
    shape = 21,
    size = 2.25,
    stroke = 0.44,
    colour = COL_DARK
  ) +
  geom_text(
    aes(
      x = Inf,
      label = paste0(
        "FDR ",
        vapply(FDR, fmt_p, character(1))
      )
    ),
    hjust = 1.02,
    family = FONT,
    size = 1.28,
    colour = "#333333"
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = COL_GREY,
      "TRUE" = COL_GOLD
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.08, 0.18))
  ) +
  labs(
    x = expression(
      log[2]("BM/PT NB count-rate ratio")
    ),
    y = NULL,
    title = "Count-based pseudobulk"
  ) +
  theme_v19()

right2 <- p2c / p2d +
  plot_layout(heights = c(1, 1))

bottom2 <- p2b | right2 +
  plot_layout(widths = c(1.36, 0.64))

fig2 <- p2a / bottom2 +
  plot_layout(heights = c(0.62, 1.38)) +
  plot_annotation(tag_levels = "A") &
  tag_theme

save_panel(p2a, "Figure2A_specimen_gene_heatmap_V20_2", 183, 47)
save_panel(p2b, "Figure2B_gene_distributions_V20_2", 118, 90)
save_panel(p2c, "Figure2C_expression_effects_V20_2", 61, 43)
save_panel(p2d, "Figure2D_count_based_effects_V20_2", 61, 43)
save_pub(
  fig2,
  "Figure2_MAJOR_REVISION_V20_2",
  183, 150
)

# ======================================================================
# FIGURE 3 — IDENTIFICATION / COMMON SUPPORT / INFLUENCE / IMBALANCE
# ======================================================================

id_order <- c(
  "All 41: site only",
  "All 41: site + chemistry + age + sex",
  "Common chemistry only: 5' V2",
  "Common-prefix subset: KRAS/STK + prefix adjustment",
  grep(
    "^Empirical age common support",
    e8_ident$analysis,
    value = TRUE
  )[1]
)

id_lab <- c(
  "All 41 · site only",
  "All 41 · chemistry + age + sex",
  "5′ V2 only",
  "KRAS/STK common-prefix + prefix",
  "Age common-support + age"
)

id_plot <- e8_ident[
  match(id_order, analysis)
]
id_plot[, label := factor(
  id_lab,
  levels = rev(id_lab)
)]
id_plot[, class := c(
  "Primary",
  "Fully adjusted",
  "Common support",
  "Common support",
  "Common support"
)]

# A. CR2/Satterthwaite identification models.
p3a <- ggplot(
  id_plot,
  aes(
    y = label,
    x = estimate_site_BM
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.30,
    colour = "#A6A6A6"
  ) +
  geom_segment(
    aes(
      x = CI95_low_CR2,
      xend = CI95_high_CR2,
      y = label,
      yend = label
    ),
    linewidth = 0.68,
    colour = COL_DARK
  ) +
  geom_point(
    aes(fill = class),
    shape = 21,
    size = 2.70,
    stroke = 0.48,
    colour = COL_DARK
  ) +
  geom_text(
    aes(
      x = Inf,
      label = paste0(
        "P ",
        vapply(
          two_sided_p_CR2,
          fmt_p,
          character(1)
        )
      )
    ),
    hjust = 1.03,
    family = FONT,
    size = 1.45,
    colour = "#333333"
  ) +
  scale_fill_manual(
    values = c(
      "Primary" = COL_PRIMARY,
      "Full model" = COL_GREY,
      "Common support" = COL_TEAL
    ),
    name = NULL
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.07, 0.24))
  ) +
  labs(
    x = "Site coefficient β (BM vs PT; CR2 95% CI)",
    y = NULL,
    title = "Site association under small-sample and common-support inference",
    subtitle = "Common-support restrictions were defined from covariate/source overlap, not expression outcomes"
  ) +
  theme_v19() +
  theme(legend.position = "bottom")

# B. Leave-one-patient influence shown as a distribution, rather than an
# arbitrary-order dot sequence.
loo_pat <- copy(e8_loo_patient)
primary_beta_v20 <- e8_ident[
  analysis == "All 41: site only",
  estimate_site_BM
][1]

p3b <- ggplot(
  loo_pat,
  aes(x = estimate_site_BM)
) +
  geom_histogram(
    bins = 12,
    fill = "#DCE4EB",
    colour = "white",
    linewidth = 0.28
  ) +
  geom_rug(
    sides = "b",
    colour = COL_PRIMARY,
    alpha = 0.75,
    linewidth = 0.42
  ) +
  geom_vline(
    xintercept = primary_beta_v20,
    linetype = 2,
    linewidth = 0.48,
    colour = COL_BM
  ) +
  annotate(
    "text",
    x = primary_beta_v20,
    y = Inf,
    hjust = -0.08,
    vjust = 1.15,
    label = paste0(
      "Full β = ",
      sprintf("%.3f", primary_beta_v20),
      "\nLOO range\n",
      sprintf("%.3f", min(loo_pat$estimate_site_BM)),
      "–",
      sprintf("%.3f", max(loo_pat$estimate_site_BM))
    ),
    family = FONT,
    size = 1.45,
    colour = COL_MID,
    lineheight = 0.94
  ) +
  labs(
    x = "Leave-one-patient-out site coefficient β",
    y = "Omissions",
    title = "Patient-level influence"
  ) +
  theme_v19()

# C. Prefix influence as horizontal bars.
loo_pref <- copy(e8_loo_prefix)
loo_pref[, omitted_prefix := factor(
  omitted_prefix,
  levels = omitted_prefix[
    order(estimate_site_BM)
  ]
)]

p3c <- ggplot(
  loo_pref,
  aes(
    y = omitted_prefix,
    x = estimate_site_BM
  )
) +
  geom_vline(
    xintercept = primary_beta_v20,
    linetype = 2,
    linewidth = 0.36,
    colour = COL_MID
  ) +
  geom_col(
    width = 0.58,
    fill = "#D9E9E6"
  ) +
  geom_point(
    shape = 21,
    size = 2.25,
    fill = COL_TEAL,
    colour = COL_DARK,
    stroke = 0.42
  ) +
  geom_text(
    aes(
      x = estimate_site_BM,
      label = sprintf("%.3f", estimate_site_BM)
    ),
    hjust = -0.22,
    family = FONT,
    size = 1.40,
    colour = COL_MID
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.18))
  ) +
  labs(
    x = "Site coefficient β",
    y = NULL,
    title = "Processing-prefix influence",
    subtitle = "Prefix denotes a processing/source stratum"
  ) +
  theme_v19()

# D. Standardized mean differences reveal cohort composition.
smd <- copy(e8_smd)
smd[, label := gsub("^Prefix: ", "", covariate)]
smd[, class := fifelse(
  grepl("^Prefix:", covariate),
  "Processing/source",
  fifelse(
    grepl("Chemistry", covariate),
    "Technical",
    "Clinical"
  )
)]
smd[, label := factor(
  label,
  levels = rev(label[order(SMD_BM_minus_PT)])
)]

p3d <- ggplot(
  smd,
  aes(
    y = label,
    x = SMD_BM_minus_PT,
    fill = class
  )
) +
  annotate(
    "rect",
    xmin = -0.10,
    xmax = 0.10,
    ymin = -Inf,
    ymax = Inf,
    fill = "#F3F5F6",
    alpha = 0.7
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.30,
    colour = "#8E8E8E"
  ) +
  geom_col(
    width = 0.58,
    alpha = 0.88
  ) +
  geom_text(
    aes(
      x = SMD_BM_minus_PT,
      label = sprintf("%+.2f", SMD_BM_minus_PT),
      hjust = ifelse(
        SMD_BM_minus_PT >= 0,
        -0.15,
        1.15
      )
    ),
    family = FONT,
    size = 1.28,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "Clinical" = "#C6A7B8",
      "Technical" = "#83B7AC",
      "Processing/source" = "#9BAFC2"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(-2, -1, 0, 1, 2),
    limits = c(-2.2, 2.2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Standardized mean difference (BM − PT)",
    y = NULL,
    title = "Measured cohort imbalance"
  ) +
  theme_v19() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 4.7),
    legend.key.width = grid::unit(4.2, "mm"),
    legend.key.height = grid::unit(3.5, "mm")
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )

bottom3 <- (p3b | p3c | p3d) +
  plot_layout(
    widths = c(0.88, 0.76, 1.36)
  )

fig3 <- p3a /
  bottom3 +
  plot_layout(
    heights = c(0.88, 1.12)
  ) +
  plot_annotation(tag_levels = "A") &
  tag_theme

save_panel(p3a, "Figure3A_identification_models_V20_2", 183, 58)
save_panel(p3b, "Figure3B_patient_influence_distribution_V20_2", 55, 68)
save_panel(p3c, "Figure3C_prefix_influence_V20_2", 48, 68)
save_panel(p3d, "Figure3D_covariate_SMD_V20_2", 82, 70)
save_pub(
  fig3,
  "Figure3_MAJOR_REVISION_V20_2",
  183, 142
)

# ======================================================================
# FIGURE 4 — SIX-PANEL FORMAL PROLIFERATION-COUPLING ANALYSIS
# ======================================================================

canon_e2b_site <- function(x) {
  fifelse(
    grepl("Brain", x, ignore.case = TRUE),
    "BM",
    "PT"
  )
}

# A. Marker-defined cell-state composition across the 40 technically
# concordant specimens. Sequential x indexing removes artificial gaps.
comp40 <- copy(e2b_comp)
comp40[, site := canon_e2b_site(cohort)]
comp40[, site_order := fifelse(site == "PT", 0L, 1L)]
setorder(
  comp40,
  site_order,
  cycling_like_fraction,
  sample
)
comp40[, x := seq_len(.N)]

state_long <- melt(
  comp40,
  id.vars = c(
    "sample", "site", "x",
    "cycling_like_fraction"
  ),
  measure.vars = c(
    "G0G1_like_fraction",
    "S_like_fraction",
    "G2M_like_fraction"
  ),
  variable.name = "state",
  value.name = "fraction"
)

state_long[, state := factor(
  state,
  levels = c(
    "G2M_like_fraction",
    "S_like_fraction",
    "G0G1_like_fraction"
  ),
  labels = c(
    "G2/M-like",
    "S-like",
    "G0/G1-like"
  )
)]

n_pt40 <- comp40[site == "PT", .N]

p4a <- ggplot(
  state_long,
  aes(
    x = x,
    y = fraction,
    fill = state
  )
) +
  geom_col(
    width = 0.92,
    colour = "white",
    linewidth = 0.05
  ) +
  geom_vline(
    xintercept = n_pt40 + 0.5,
    linewidth = 0.55,
    colour = "white"
  ) +
  annotate(
    "text",
    x = (n_pt40 + 1) / 2,
    y = 1.055,
    label = "Primary tumor",
    family = FONT,
    size = 1.65,
    colour = COL_DARK
  ) +
  annotate(
    "text",
    x = n_pt40 + (40 - n_pt40 + 1) / 2,
    y = 1.055,
    label = "Brain metastasis",
    family = FONT,
    size = 1.65,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "G0/G1-like" = "#B3BAC1",
      "S-like" = "#E3A25F",
      "G2/M-like" = "#9675A8"
    ),
    name = "Marker-defined state"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 25),
    limits = c(0, 1.08),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  labs(
    x = "Specimens ordered by cycling-like fraction within site",
    y = "Malignant-nucleus composition",
    title = "Marker-defined malignant-cell states"
  ) +
  theme_v19() +
  theme(
    legend.position = "top",
    axis.line.x = element_blank()
  )

# B. Cycling-like fraction distribution.
cyc_stat <- e2b_stats[
  outcome == "cycling_like_fraction"
][1]

comp40[, site_factor := factor(
  site,
  levels = c("PT", "BM")
)]

p4b <- ggplot(
  comp40,
  aes(
    x = site_factor,
    y = cycling_like_fraction,
    fill = site_factor,
    colour = site_factor
  )
) +
  geom_violin(
    width = 0.84,
    alpha = 0.07,
    linewidth = 0.36,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.30,
    alpha = 0.14,
    outlier.shape = NA,
    linewidth = 0.38
  ) +
  geom_point(
    position = position_jitter(
      width = 0.055,
      height = 0,
      seed = 41
    ),
    size = 0.78,
    alpha = 0.66
  ) +
  scale_fill_manual(
    values = c(
      "PT" = COL_PRIMARY,
      "BM" = COL_BM
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      "PT" = COL_PRIMARY,
      "BM" = COL_BM
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 25),
    expand = expansion(mult = c(0.02, 0.16))
  ) +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    vjust = 1.25,
    label = paste0(
      "Δmedian ",
      sprintf("%.3f", cyc_stat$median_difference),
      "\nP ",
      fmt_p(cyc_stat$two_sided_wilcoxon_p),
      " · g ",
      sprintf("%.2f", cyc_stat$hedges_g)
    ),
    family = FONT,
    size = 1.45,
    colour = COL_MID,
    lineheight = 0.93
  ) +
  labs(
    x = NULL,
    y = "Cycling-like fraction",
    title = "Cycling-like fraction"
  ) +
  theme_v19()

# C. Marker-defined G0/G1-like module score.
g1_plot <- copy(e2b_g1)
g1_plot[, site := factor(
  canon_e2b_site(cohort),
  levels = c("PT", "BM")
)]

g1_stat <- e2b_stats[
  outcome == "module_score_G0G1_fixed" &
    grepl(">=50 nuclei", analysis, fixed = TRUE) &
    !grepl("independence", analysis, ignore.case = TRUE)
][1]

p4c <- ggplot(
  g1_plot,
  aes(
    x = site,
    y = module_score_G0G1_fixed,
    fill = site,
    colour = site
  )
) +
  geom_violin(
    width = 0.84,
    alpha = 0.07,
    linewidth = 0.36,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.30,
    alpha = 0.14,
    outlier.shape = NA,
    linewidth = 0.38
  ) +
  geom_point(
    position = position_jitter(
      width = 0.055,
      height = 0,
      seed = 43
    ),
    size = 0.78,
    alpha = 0.66
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.30,
    colour = "#A8A8A8"
  ) +
  scale_fill_manual(
    values = c(
      "PT" = COL_PRIMARY,
      "BM" = COL_BM
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      "PT" = COL_PRIMARY,
      "BM" = COL_BM
    ),
    guide = "none"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    vjust = 1.25,
    label = paste0(
      "CR1 β ",
      sprintf("%.3f", e2b_g1_cr1$estimate[1]),
      "\n95% CI ",
      sprintf("%.2f", e2b_g1_cr1$CI95_low_CR1[1]),
      " to ",
      sprintf("%.2f", e2b_g1_cr1$CI95_high_CR1[1]),
      "\nP ",
      fmt_p(e2b_g1_cr1$two_sided_p_CR1[1])
    ),
    family = FONT,
    size = 1.40,
    colour = COL_MID,
    lineheight = 0.93
  ) +
  labs(
    x = NULL,
    y = "Fixed five-gene score",
    title = "G0/G1-like compartment"
  ) +
  theme_v19()

# D. Formal attenuation summary. Horizontal bars avoid the cramped
# vertical labels seen in V20 while preserving the same estimates/CIs.
att <- rbindlist(list(
  e8_prolif[
    analysis == "40 concordant specimens: site only",
    .(
      state = "All malignant",
      beta = estimate_site_BM,
      low = CI95_low_CR2,
      high = CI95_high_CR2
    )
  ],
  e8_prolif[
    analysis ==
      "40 concordant: site + continuous cycling fraction",
    .(
      state = "+ cycling fraction",
      beta = estimate_site_BM,
      low = CI95_low_CR2,
      high = CI95_high_CR2
    )
  ],
  data.table(
    state = "G0/G1-like",
    beta = e2b_g1_cr1$estimate[1],
    low = e2b_g1_cr1$CI95_low_CR1[1],
    high = e2b_g1_cr1$CI95_high_CR1[1]
  )
))
att[, state := factor(
  state,
  levels = rev(c(
    "All malignant",
    "+ cycling fraction",
    "G0/G1-like"
  ))
)]
att[, class := factor(
  as.character(state),
  levels = c(
    "All malignant",
    "+ cycling fraction",
    "G0/G1-like"
  )
)]

p4d <- ggplot(
  att,
  aes(
    y = state,
    x = beta,
    fill = class
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.30,
    colour = "#8E8E8E"
  ) +
  geom_col(
    width = 0.56,
    alpha = 0.88
  ) +
  geom_segment(
    aes(
      x = low,
      xend = high,
      y = state,
      yend = state
    ),
    linewidth = 0.62,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "All malignant" = COL_BM,
      "+ cycling fraction" = "#D9C8A6",
      "G0/G1-like" = COL_TEAL
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(-0.5, 0, 0.5, 1.0),
    labels = c("-0.5", "0", "0.5", "1"),
    limits = c(-0.6, 1.1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    x = "Site β (95% CI)",
    y = NULL,
    title = "Site-effect attenuation",
    subtitle = paste0(
      "Cluster bootstrap: Δβ = ",
      sprintf("%.3f", e8_delta$delta_beta_all_minus_G0G1[1]),
      "\n95% CI ",
      sprintf("%.3f", e8_delta$bootstrap_CI95_low[1]),
      "–",
      sprintf("%.3f", e8_delta$bootstrap_CI95_high[1]),
      "; P ",
      fmt_p(e8_delta$bootstrap_two_sided_p_delta[1])
    )
  ) +
  theme_v19() +
  theme(
    axis.text.x = element_text(size = 6)
  )

# E. Cycling-fraction site coefficients across the robustness models,
# displayed as bars to diversify the visual grammar.
cycle_comb <- rbindlist(
  list(
    e8_cycle[, .(
      analysis,
      estimate_site_BM,
      CI95_low_CR2,
      CI95_high_CR2,
      two_sided_p_CR2
    )],
    e8_cycle_support[, .(
      analysis,
      estimate_site_BM,
      CI95_low_CR2,
      CI95_high_CR2,
      two_sided_p_CR2
    )]
  ),
  fill = TRUE
)

cycle_comb[, label := fcase(
  grepl("site only$", analysis), "Site only",
  grepl("site \\+ chemistry$", analysis), "+ chemistry",
  grepl("site \\+ age$", analysis), "+ age",
  grepl("chemistry \\+ age \\+ sex", analysis),
    "+ chemistry + age + sex",
  grepl("V2-only", analysis), "5′ V2 only",
  default = "KRAS/STK common-prefix + prefix"
)]
cycle_comb[, class := fifelse(
  grepl("V2-only|common-prefix", analysis),
  "Common support",
  fifelse(
    grepl("chemistry \\+ age \\+ sex", analysis),
    "Fully adjusted",
    "Covariate"
  )
)]
cycle_comb[, label := factor(
  label,
  levels = rev(c(
    "Site only",
    "+ chemistry",
    "+ age",
    "+ chemistry + age + sex",
    "5′ V2 only",
    "KRAS/STK common-prefix + prefix"
  ))
)]

p4e <- ggplot(
  cycle_comb,
  aes(
    y = label,
    x = estimate_site_BM,
    fill = class
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.30,
    colour = "#8E8E8E"
  ) +
  geom_col(
    width = 0.58,
    alpha = 0.88
  ) +
  geom_segment(
    aes(
      x = CI95_low_CR2,
      xend = CI95_high_CR2,
      y = label,
      yend = label
    ),
    linewidth = 0.50,
    colour = COL_DARK
  ) +
  geom_text(
    aes(
      x = Inf,
      label = paste0(
        "P ",
        vapply(
          two_sided_p_CR2,
          fmt_p,
          character(1)
        )
      )
    ),
    hjust = 1.02,
    family = FONT,
    size = 1.30,
    colour = "#333333"
  ) +
  scale_fill_manual(
    values = c(
      "Covariate" = COL_GOLD,
      "Fully adjusted" = COL_GREY,
      "Common support" = COL_TEAL
    ),
    name = NULL
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.10, 0.30))
  ) +
  labs(
    x = "Cycling-like fraction site coefficient β",
    y = NULL,
    title = "Cycling-fraction robustness"
  ) +
  theme_v19() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 5.0),
    legend.key.width = grid::unit(4.2, "mm")
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

# F. Gene-level G0/G1-like sensitivity.
g1g <- copy(e2b_g1_gene)
g1g[, gene := factor(
  gene,
  levels = rev(c(
    "DHFR", "DHODH", "SHMT1", "TYMS", "UMPS"
  ))
)]
g1g[, focus := gene == "UMPS"]

p4f <- ggplot(
  g1g,
  aes(
    y = gene,
    x = hedges_g
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.30,
    colour = "#A5A5A5"
  ) +
  geom_segment(
    aes(
      x = hedges_g_95CI_low,
      xend = hedges_g_95CI_high,
      y = gene,
      yend = gene
    ),
    linewidth = 0.62,
    colour = COL_DARK
  ) +
  geom_point(
    aes(fill = focus),
    shape = 21,
    size = 2.55,
    stroke = 0.46,
    colour = COL_DARK
  ) +
  geom_text(
    aes(
      x = Inf,
      label = paste0(
        "FDR ",
        vapply(
          BH_FDR_across_5_genes,
          fmt_p,
          character(1)
        )
      )
    ),
    hjust = 1.02,
    family = FONT,
    size = 1.28,
    colour = "#333333"
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = COL_GREY,
      "TRUE" = COL_GOLD
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.10, 0.27))
  ) +
  labs(
    x = "Hedges g (95% CI)",
    y = NULL,
    title = "Gene-level G0/G1 sensitivity"
  ) +
  theme_v19()

middle4 <- p4b | p4c | p4d +
  plot_layout(widths = c(0.86, 0.86, 1.28))
bottom4 <- p4e | p4f +
  plot_layout(widths = c(1.16, 0.84))

fig4 <- p4a /
  middle4 /
  bottom4 +
  plot_layout(
    heights = c(0.84, 1.00, 1.04)
  ) +
  plot_annotation(tag_levels = "A") &
  tag_theme

save_panel(p4a, "Figure4A_cell_state_composition_V20_2", 183, 54)
save_panel(p4b, "Figure4B_cycling_fraction_V20_2", 58, 62)
save_panel(p4c, "Figure4C_G0G1_score_V20_2", 58, 62)
save_panel(p4d, "Figure4D_formal_attenuation_V20_2", 58, 62)
save_panel(p4e, "Figure4E_cycling_robustness_models_V20_2", 104, 67)
save_panel(p4f, "Figure4F_G0G1_gene_sensitivity_V20_2", 75, 67)
save_pub(
  fig4,
  "Figure4_MAJOR_REVISION_V20_2",
  183, 190
)

# ======================================================================
# FIGURE 5 — SIX-PANEL TARGET-PANEL SPECIFICITY CALIBRATION
# ======================================================================

# A. Comparator programs are shown as a coefficient heat map rather than
# connected point-whisker lines. This avoids the visual implication that
# site-only and cycling-conditioned estimates are paired trajectories.
comp <- copy(e9_comp)
comp[, model := fifelse(
  grepl("site \\+ cycling", analysis),
  "+ cycling fraction",
  "Site only"
)]

panel_name_map <- c(
  "fixed_five_gene" = "Fixed five-gene",
  "nonoverlap_S_phase" = "Non-overlap S phase",
  "nonoverlap_G2M" = "Non-overlap G2/M",
  "DNA_replication_machinery" = "DNA-replication machinery"
)

panel_order_v20 <- c(
  "fixed_five_gene",
  "nonoverlap_S_phase",
  "nonoverlap_G2M",
  "DNA_replication_machinery"
)
panel_meta_v20 <- unique(
  comp[, .(panel, n_genes)]
)[match(panel_order_v20, panel)]
panel_meta_v20[, panel_label := paste0(
  panel_name_map[panel],
  "  (", n_genes, " genes)"
)]

comp <- merge(
  comp,
  panel_meta_v20[, .(panel, panel_label)],
  by = "panel",
  all.x = TRUE,
  sort = FALSE
)
comp[, panel_label := factor(
  panel_label,
  levels = rev(panel_meta_v20$panel_label)
)]
comp[, model := factor(
  model,
  levels = c(
    "Site only",
    "+ cycling fraction"
  )
)]

comp_lim <- max(abs(comp$estimate_site_BM), na.rm = TRUE)

p5a <- ggplot(
  comp,
  aes(
    x = model,
    y = panel_label,
    fill = estimate_site_BM
  )
) +
  geom_tile(
    width = 0.92,
    height = 0.76,
    colour = "white",
    linewidth = 0.65
  ) +
  geom_text(
    aes(
      label = paste0(
        "β ",
        sprintf("%.2f", estimate_site_BM),
        "\nP ",
        vapply(
          two_sided_p_CR2,
          fmt_p,
          character(1)
        )
      )
    ),
    family = FONT,
    size = 1.42,
    lineheight = 0.92,
    colour = COL_DARK
  ) +
  scale_fill_gradient2(
    low = COL_PRIMARY,
    mid = "white",
    high = COL_BM,
    midpoint = 0,
    limits = c(-comp_lim, comp_lim),
    oob = scales::squish,
    name = "Site β"
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      direction = "horizontal",
      barwidth = grid::unit(28, "mm"),
      barheight = grid::unit(2.6, "mm")
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Comparator programs",
    subtitle = "Cycling/replication comparators show larger site effects"
  ) +
  theme_v19() +
  theme(
    axis.text.x = element_text(
      size = 6.2,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.justification = "left"
  )

# B. Nested fixed-module decomposition.
nest <- copy(e9_nested)
nest <- nest[
  !grepl("site \\+ cycling", analysis)
]
nest[, panel_label := factor(
  panel,
  levels = rev(c(
    "five_gene_fixed",
    "TYMS_UMPS_only",
    "minus_TYMS",
    "minus_UMPS",
    "minus_TYMS_and_UMPS"
  )),
  labels = rev(c(
    "Five-gene fixed score",
    "TYMS + UMPS",
    "Five-gene minus TYMS",
    "Five-gene minus UMPS",
    "DHFR + DHODH + SHMT1"
  ))
)]

p5b <- ggplot(
  nest,
  aes(
    y = panel_label,
    x = estimate_site_BM
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.30,
    colour = "#A6A6A6"
  ) +
  geom_segment(
    aes(
      x = CI95_low_CR2,
      xend = CI95_high_CR2,
      y = panel_label,
      yend = panel_label
    ),
    linewidth = 0.65,
    colour = COL_DARK
  ) +
  geom_point(
    aes(
      fill = panel == "TYMS_UMPS_only"
    ),
    shape = 21,
    size = 2.65,
    stroke = 0.48,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = COL_GREY,
      "TRUE" = COL_GOLD
    ),
    guide = "none"
  ) +
  geom_text(
    aes(
      x = Inf,
      label = paste0(
        "P ",
        vapply(
          two_sided_p_CR2,
          fmt_p,
          character(1)
        )
      )
    ),
    hjust = 1.03,
    family = FONT,
    size = 1.18,
    colour = COL_MID
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.07, 0.34))
  ) +
  labs(
    x = "Site coefficient β (CR2 95% CI)",
    y = NULL,
    title = "TYMS–UMPS dominance",
    subtitle = "Secondary nested decomposition of the fixed module"
  ) +
  theme_v19()

# C. Proliferation-matched random-set null.
obs_beta <- e9_null$observed_site_beta[1]
null_med <- e9_null$null_beta_median[1]
null_q025 <- e9_null$null_beta_q025[1]
null_q975 <- e9_null$null_beta_q975[1]
emp_p <- e9_null$empirical_two_sided_p[1]
pct <- e9_null$observed_beta_percentile[1]

p5c <- ggplot(
  e9_null_dist,
  aes(x = site_beta_standardized_score)
) +
  annotate(
    "rect",
    xmin = null_q025,
    xmax = null_q975,
    ymin = -Inf,
    ymax = Inf,
    fill = "#D8E9E6",
    alpha = 0.48
  ) +
  geom_histogram(
    bins = 48,
    fill = "#D7DCE1",
    colour = "white",
    linewidth = 0.18
  ) +
  geom_vline(
    xintercept = null_med,
    colour = COL_DARK,
    linetype = 2,
    linewidth = 0.48
  ) +
  geom_vline(
    xintercept = obs_beta,
    colour = COL_BM,
    linewidth = 0.90
  ) +
  annotate(
    "text",
    x = obs_beta + 0.012,
    y = Inf,
    hjust = 0,
    vjust = 1.10,
    label = paste0(
      "Observed β ",
      sprintf("%.3f", obs_beta),
      " (",
      sprintf("%.2f", pct),
      "th percentile)",
      "\nMatched sets ≥ observed: ",
      sprintf("%.1f%%", 100 * emp_p)
    ),
    family = FONT,
    size = 1.25,
    lineheight = 0.93,
    colour = COL_BM
  ) +
  annotate(
    "text",
    x = null_med + 0.012,
    y = Inf,
    hjust = 0,
    vjust = 1.10,
    label = paste0(
      "Null median ",
      sprintf("%.3f", null_med)
    ),
    family = FONT,
    size = 1.22,
    colour = COL_DARK
  ) +
  labs(
    x = "Site coefficient β of matched random sets",
    y = "Random sets",
    title = "Proliferation-matched null",
    subtitle = "10,000 matched five-gene sets; site labels excluded"
  ) +
  theme_v19()

# D. Gene-level target-panel landscape.
gene_plot <- copy(e9_gene)
gene_plot[, is_module := gene %chin% c(
  "DHFR", "DHODH", "SHMT1", "TYMS", "UMPS"
)]

gene_plot[, label_x := spearman_rho_with_cycling_fraction]
gene_plot[, label_y := site_beta_standardized]
gene_plot[, label_hjust := 0]

gene_plot[gene == "TYMS", `:=`(
  label_x = spearman_rho_with_cycling_fraction - 0.025,
  label_y = site_beta_standardized + 0.085,
  label_hjust = 1
)]
gene_plot[gene == "UMPS", `:=`(
  label_x = spearman_rho_with_cycling_fraction - 0.020,
  label_y = site_beta_standardized - 0.080,
  label_hjust = 1
)]
gene_plot[gene == "DHFR", `:=`(
  label_x = spearman_rho_with_cycling_fraction + 0.025,
  label_y = site_beta_standardized + 0.045,
  label_hjust = 0
)]
gene_plot[gene == "DHODH", `:=`(
  label_x = spearman_rho_with_cycling_fraction + 0.025,
  label_y = site_beta_standardized - 0.045,
  label_hjust = 0
)]
gene_plot[gene == "SHMT1", `:=`(
  label_x = spearman_rho_with_cycling_fraction + 0.025,
  label_y = site_beta_standardized + 0.050,
  label_hjust = 0
)]

p5d <- ggplot(
  gene_plot,
  aes(
    x = spearman_rho_with_cycling_fraction,
    y = site_beta_standardized
  )
) +
  annotate(
    "rect",
    xmin = 0.55,
    xmax = Inf,
    ymin = 0.75,
    ymax = Inf,
    fill = "#F8EDEC",
    alpha = 0.55
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.28,
    colour = "#B8B8B8"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.28,
    colour = "#B8B8B8"
  ) +
  geom_point(
    data = gene_plot[is_module == FALSE],
    size = 1.15,
    alpha = 0.48,
    colour = "#9AA1A7"
  ) +
  geom_point(
    data = gene_plot[is_module == TRUE],
    aes(fill = gene),
    shape = 21,
    size = 2.35,
    stroke = 0.45,
    colour = COL_DARK
  ) +
  geom_text(
    data = gene_plot[
      gene %chin% c(
        "DHFR", "DHODH", "SHMT1", "TYMS", "UMPS"
      )
    ],
    aes(
      x = label_x,
      y = label_y,
      label = gene,
      hjust = label_hjust
    ),
    inherit.aes = FALSE,
    family = FONT,
    size = 1.35,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "DHFR" = "#AAB1B7",
      "DHODH" = "#AAB1B7",
      "SHMT1" = "#AAB1B7",
      "TYMS" = COL_BM,
      "UMPS" = COL_GOLD
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.08, 0.20))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  labs(
    x = "Spearman ρ with cycling-like fraction",
    y = "Standardized BM–PT site effect",
    title = "Gene-level landscape"
  ) +
  theme_v19()

# E. Matching-feature audit in the main figure.
bal <- melt(
  e9_balance,
  id.vars = "quantity",
  variable.name = "group",
  value.name = "value"
)
bal[, quantity_short := fcase(
  grepl("rho", quantity, ignore.case = TRUE),
    "|ρ(gene, cycling)|",
  grepl("mean", quantity, ignore.case = TRUE),
    "Mean log2CPM",
  default = "SD log2CPM"
)]
bal[, group := factor(
  group,
  levels = c(
    "observed_module_median",
    "candidate_pool_median"
  ),
  labels = c(
    "Fixed module",
    "Cycle-gene pool"
  )
)]

p5e <- ggplot(
  bal,
  aes(
    x = group,
    y = value,
    fill = group
  )
) +
  geom_col(
    width = 0.58,
    alpha = 0.90
  ) +
  facet_wrap(
    ~ quantity_short,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Fixed module" = COL_BM,
      "Cycle-gene pool" = COL_GREY
    ),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "Median",
    title = "Matching audit"
  ) +
  theme_v19() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

# F. Distribution of proliferation coupling in matched random sets.
p5f <- ggplot(
  e9_null_dist,
  aes(
    x = spearman_score_vs_cycling
  )
) +
  geom_histogram(
    bins = 45,
    fill = "#D8DDE1",
    colour = "white",
    linewidth = 0.18
  ) +
  geom_vline(
    xintercept =
      e9_null$observed_spearman_score_vs_cycling[1],
    colour = COL_BM,
    linewidth = 0.90
  ) +
  geom_vline(
    xintercept =
      e9_null$null_score_cycle_rho_median[1],
    colour = COL_DARK,
    linetype = 2,
    linewidth = 0.48
  ) +
  annotate(
    "text",
    x =
      e9_null$observed_spearman_score_vs_cycling[1],
    y = Inf,
    hjust = -0.05,
    vjust = 1.13,
    label = paste0(
      "Fixed score ρ ",
      sprintf(
        "%.3f",
        e9_null$observed_spearman_score_vs_cycling[1]
      )
    ),
    family = FONT,
    size = 1.42,
    colour = COL_BM
  ) +
  labs(
    x = "Spearman ρ of random-set score with cycling fraction",
    y = "Random sets",
    title = "Random-set cycling"
  ) +
  theme_v19()

top5 <- p5a | p5b +
  plot_layout(widths = c(1.12, 0.88))
mid5 <- p5c | p5d +
  plot_layout(widths = c(1.12, 0.88))
bot5 <- p5e | p5f +
  plot_layout(widths = c(1.10, 0.90))

fig5 <- top5 /
  mid5 /
  bot5 +
  plot_layout(
    heights = c(0.92, 1.08, 0.80)
  ) +
  plot_annotation(
    tag_levels = "A",
    caption = NULL
  ) &
  tag_theme

save_panel(p5a, "Figure5A_comparator_heatmap_V20_2", 89, 62)
save_panel(p5b, "Figure5B_nested_decomposition_V20_2", 89, 62)
save_panel(p5c, "Figure5C_matched_null_V20_2", 89, 72)
save_panel(p5d, "Figure5D_gene_landscape_V20_2", 89, 72)
save_panel(p5e, "Figure5E_matching_feature_audit_V20_2", 89, 55)
save_panel(p5f, "Figure5F_random_set_cycling_V20_2", 89, 55)
save_pub(
  fig5,
  "Figure5_MAJOR_REVISION_V20_2",
  183, 198
)

# ======================================================================
# SUPPLEMENTARY FIGURE 1 — NONREDUNDANT FIXED-SCORE ROBUSTNESS
# ======================================================================

s1_plot <- rbindlist(list(
  data.table(
    label = "Primary ≥50 nuclei",
    n_BM = 31L,
    n_PT = 10L,
    g = s_primary$hedges_g,
    low = s_primary$hedges_g_95CI_low,
    high = s_primary$hedges_g_95CI_high,
    p = s_primary$two_sided_wilcoxon_p
  ),
  data.table(
    label = "≥500 nuclei",
    n_BM = 28L,
    n_PT = 10L,
    g = s_500$hedges_g,
    low = s_500$hedges_g_95CI_low,
    high = s_500$hedges_g_95CI_high,
    p = s_500$two_sided_wilcoxon_p
  ),
  data.table(
    label = "≥1000 nuclei",
    n_BM = 25L,
    n_PT = 10L,
    g = s_1000$hedges_g,
    low = s_1000$hedges_g_95CI_low,
    high = s_1000$hedges_g_95CI_high,
    p = s_1000$two_sided_wilcoxon_p
  ),
  data.table(
    label = "KRAS_17 excluded",
    n_BM = 30L,
    n_PT = 10L,
    g = s_k17$hedges_g,
    low = s_k17$hedges_g_95CI_low,
    high = s_k17$hedges_g_95CI_high,
    p = s_k17$two_sided_wilcoxon_p
  ),
  data.table(
    label = "Matched pairs fully excluded",
    n_BM = 29L,
    n_PT = 8L,
    g = e1_ind_v193$hedges_g[1],
    low = e1_ind_v193$hedges_g_95CI_low[1],
    high = e1_ind_v193$hedges_g_95CI_high[1],
    p = e1_ind_v193$two_sided_wilcoxon_p[1]
  ),
  data.table(
    label = "No threshold + STK_20",
    n_BM = 31L,
    n_PT = 11L,
    g = e1_stk20_v193$hedges_g[1],
    low = e1_stk20_v193$hedges_g_95CI_low[1],
    high = e1_stk20_v193$hedges_g_95CI_high[1],
    p = e1_stk20_v193$two_sided_wilcoxon_p[1]
  )
))

s1_levels <- c(
  "Primary ≥50 nuclei",
  "≥500 nuclei",
  "≥1000 nuclei",
  "KRAS_17 excluded",
  "Matched pairs fully excluded",
  "No threshold + STK_20"
)

s1_plot[, label := factor(
  label,
  levels = rev(s1_levels)
)]
s1_plot[, minus_log10_p := -log10(p)]

ps1a <- ggplot(
  s1_plot,
  aes(y = label, x = g)
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "#A6A6A6",
    linewidth = 0.30
  ) +
  geom_segment(
    aes(
      x = low,
      xend = high,
      y = label,
      yend = label
    ),
    linewidth = 0.62,
    colour = COL_DARK
  ) +
  geom_point(
    shape = 21,
    size = 2.55,
    fill = COL_TEAL,
    colour = COL_DARK,
    stroke = 0.45
  ) +
  labs(
    x = "Hedges g (95% CI)",
    y = NULL,
    title = "Effect-size robustness",
    subtitle = "≥200 omitted because it is identical to the primary ≥50 analysis"
  ) +
  theme_v19()

ps1b <- ggplot(
  s1_plot,
  aes(y = label, x = minus_log10_p)
) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = 2,
    colour = COL_GOLD,
    linewidth = 0.36
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = minus_log10_p,
      y = label,
      yend = label
    ),
    colour = COL_LIGHT,
    linewidth = 0.70
  ) +
  geom_point(
    shape = 21,
    size = 2.45,
    fill = COL_PRIMARY,
    colour = COL_DARK,
    stroke = 0.45
  ) +
  labs(
    x = expression(-log[10]("two-sided P")),
    y = NULL,
    title = "Two-sided sensitivity P values"
  ) +
  theme_v19() +
  theme(axis.text.y = element_blank())

s1_n <- melt(
  s1_plot[, .(
    label,
    `Brain metastasis` = n_BM,
    `Primary tumor` = n_PT
  )],
  id.vars = "label",
  variable.name = "site",
  value.name = "n"
)

ps1c <- ggplot(
  s1_n,
  aes(
    y = label,
    x = n,
    fill = site
  )
) +
  geom_col(
    width = 0.60
  ) +
  scale_fill_manual(
    values = c(
      "Primary tumor" = COL_PRIMARY,
      "Brain metastasis" = COL_BM
    ),
    name = NULL
  ) +
  labs(
    x = "Specimens retained",
    y = NULL,
    title = "Sample-set size"
  ) +
  theme_v19() +
  theme(
    axis.text.y = element_blank(),
    legend.position = "top"
  )

supp1 <- ps1a | ps1b | ps1c +
  plot_layout(
    widths = c(1.25, 0.78, 0.72)
  ) +
  plot_annotation(tag_levels = "A") &
  tag_theme

save_panel(ps1a, "SupplementaryFigure1A_effect_size_V20_2", 82, 78)
save_panel(ps1b, "SupplementaryFigure1B_pvalues_V20_2", 51, 78)
save_panel(ps1c, "SupplementaryFigure1C_sample_counts_V20_2", 46, 78)
save_pub(
  supp1,
  "Supplementary_Figure1_MAJOR_REVISION_V20_2",
  183, 88
)

# ======================================================================
# SUPPLEMENTARY FIGURE 2 — COHORT COMPOSITION + IMBALANCE AUDIT
# ======================================================================

base41 <- copy(e8_analysis41)
base41[, site := factor(
  fifelse(site_BM == 1L, "Brain metastasis", "Primary tumor"),
  levels = c("Primary tumor", "Brain metastasis")
)]

chem <- base41[
  ,
  .N,
  by = .(site, snrna_seq_chemistry)
]
chem[, fraction := N / sum(N), by = site]

chem_levels <- sort(unique(chem$snrna_seq_chemistry))
chem_palette <- setNames(
  if (length(chem_levels) == 2L) {
    c("#2F7F73", "#87B8AE")
  } else {
    scales::hue_pal()(length(chem_levels))
  },
  chem_levels
)

ps2a <- ggplot(
  chem,
  aes(
    x = site,
    y = fraction,
    fill = snrna_seq_chemistry
  )
) +
  geom_col(
    width = 0.58,
    colour = "white",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = chem_palette,
    name = "snRNA-seq chemistry"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 25),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    labels = c(
      "Primary tumor" = "PT",
      "Brain metastasis" = "BM"
    )
  ) +
  labs(
    x = NULL,
    y = "Within-site fraction",
    title = "Documented snRNA-seq chemistry"
  ) +
  theme_v19() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom"
  )

ps2b <- ggplot(
  base41,
  aes(
    x = site,
    y = age,
    fill = site,
    colour = site
  )
) +
  geom_violin(
    width = 0.82,
    alpha = 0.07,
    linewidth = 0.35,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.30,
    alpha = 0.14,
    outlier.shape = NA,
    linewidth = 0.40
  ) +
  geom_point(
    position = position_jitter(
      width = 0.055,
      height = 0,
      seed = 77
    ),
    size = 0.78,
    alpha = 0.65
  ) +
  scale_fill_manual(
    values = c(
      "Primary tumor" = COL_PRIMARY,
      "Brain metastasis" = COL_BM
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      "Primary tumor" = COL_PRIMARY,
      "Brain metastasis" = COL_BM
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = c(
      "Primary tumor" = "PT",
      "Brain metastasis" = "BM"
    )
  ) +
  labs(
    x = NULL,
    y = "Age at resection (years)",
    title = "Age distribution"
  ) +
  theme_v19() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

sex_dt <- base41[
  ,
  .N,
  by = .(site, sex)
]
sex_dt[, fraction := N / sum(N), by = site]

sex_levels <- sort(unique(as.character(sex_dt$sex)))
sex_palette <- setNames(
  c("#C4A5B7", "#7E9DB8")[seq_along(sex_levels)],
  sex_levels
)

ps2c <- ggplot(
  sex_dt,
  aes(
    x = site,
    y = fraction,
    fill = sex
  )
) +
  geom_col(
    width = 0.58,
    colour = "white",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = sex_palette,
    name = "Sex"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 25),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    labels = c(
      "Primary tumor" = "PT",
      "Brain metastasis" = "BM"
    )
  ) +
  labs(
    x = NULL,
    y = "Within-site fraction",
    title = "Sex composition"
  ) +
  theme_v19() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom"
  )

smd_supp <- copy(e8_smd)
smd_supp[, abs_smd := abs(SMD_BM_minus_PT)]
smd_supp[, label := gsub("^Prefix: ", "", covariate)]
smd_supp[, class := fcase(
  grepl("^Prefix:", covariate), "Processing/source",
  grepl("Chemistry", covariate), "Technical",
  default = "Clinical"
)]
smd_supp[, label := factor(
  label,
  levels = label[order(abs_smd)]
)]

ps2d <- ggplot(
  smd_supp,
  aes(
    y = label,
    x = abs_smd,
    fill = class
  )
) +
  annotate(
    "rect",
    xmin = 0,
    xmax = 0.10,
    ymin = -Inf,
    ymax = Inf,
    fill = "#F1F3F5",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = 0.10,
    linetype = 2,
    linewidth = 0.32,
    colour = "#9C9C9C"
  ) +
  geom_col(
    width = 0.58,
    alpha = 0.90
  ) +
  geom_text(
    aes(
      x = abs_smd,
      label = sprintf("%.2f", abs_smd)
    ),
    hjust = -0.15,
    family = FONT,
    size = 1.30,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "Clinical" = "#C6A7B8",
      "Technical" = "#83B7AC",
      "Processing/source" = "#9BAFC2"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.18))
  ) +
  labs(
    x = "|standardized mean difference|",
    y = NULL,
    title = "Magnitude of measured imbalance"
  ) +
  theme_v19() +
  theme(legend.position = "bottom")

top_s2 <- ps2a | ps2b | ps2c +
  plot_layout(
    widths = c(1.08, 1.12, 0.80)
  )

supp2 <- top_s2 /
  ps2d +
  plot_layout(
    heights = c(1.0, 0.90)
  ) +
  plot_annotation(tag_levels = "A") &
  tag_theme

save_panel(ps2a, "SupplementaryFigure2A_chemistry_V20_2", 58, 68)
save_panel(ps2b, "SupplementaryFigure2B_age_V20_2", 58, 68)
save_panel(ps2c, "SupplementaryFigure2C_sex_V20_2", 58, 68)
save_panel(ps2d, "SupplementaryFigure2D_SMD_V20_2", 183, 58)
save_pub(
  supp2,
  "Supplementary_Figure2_MAJOR_REVISION_V20_2",
  183, 132
)

# ======================================================================
# SUPPLEMENTARY FIGURE 3 — SAME-STUDY SPATIAL CONTEXT + AUTOCORRELATION NULL
# ======================================================================

rep_bm <- e5_rep[
  grepl("Brain", cohort, ignore.case = TRUE),
  sample
][1]
rep_pt <- e5_rep[
  grepl("Primary", cohort, ignore.case = TRUE),
  sample
][1]

map_data <- spatial_grid[
  sample %chin% c(rep_bm, rep_pt),
  .(
    sample, cohort, gx, gy,
    pyrimidine_score,
    proliferation_score
  )
]

map_long <- melt(
  map_data,
  id.vars = c(
    "sample", "cohort", "gx", "gy"
  ),
  measure.vars = c(
    "pyrimidine_score",
    "proliferation_score"
  ),
  variable.name = "signal",
  value.name = "score"
)

map_long[, site_lab := ifelse(
  grepl("Brain", cohort, ignore.case = TRUE),
  paste0(sample, " · BM"),
  paste0(sample, " · PT")
)]
map_long[, sig_lab := ifelse(
  signal == "pyrimidine_score",
  "TYMS–UMPS",
  "Proliferation"
)]

map_lim <- quantile(
  abs(map_long$score),
  0.98,
  na.rm = TRUE
)

ps3a <- ggplot(
  map_long,
  aes(
    x = gx,
    y = gy,
    fill = score
  )
) +
  geom_tile(colour = NA) +
  coord_fixed(expand = FALSE) +
  facet_grid(site_lab ~ sig_lab) +
  scale_fill_gradient2(
    low = COL_PRIMARY,
    mid = "white",
    high = COL_BM,
    midpoint = 0,
    limits = c(-map_lim, map_lim),
    oob = scales::squish,
    name = "Within-specimen\nstandardized score"
  ) +
  labs(
    title = "Same-study Slide-seq maps",
    subtitle = "Representative specimens selected by eligible-grid count; grid width = 400 deposited coordinate units"
  ) +
  theme_void(
    base_family = FONT,
    base_size = 7
  ) +
  theme(
    plot.title = element_text(
      size = 8,
      colour = "black",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 5.8,
      colour = COL_MID,
      hjust = 0
    ),
    strip.text = element_text(
      size = 6.3,
      colour = "black"
    ),
    legend.position = "bottom"
  )

sp_adj <- e5_shift[
  association == "partial_rank_residual_surface"
]
sp_adj[, site := fifelse(
  grepl("Brain", cohort, ignore.case = TRUE),
  "BM",
  "PT"
)]
sp_adj[, sample_order := factor(
  sample,
  levels = sample[
    order(observed_rho)
  ]
)]

# B. Observed adjusted rho with the specimen-specific toroidal null interval.
ps3b <- ggplot(
  sp_adj,
  aes(
    y = sample_order,
    x = observed_rho
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.30,
    colour = "#A8A8A8"
  ) +
  geom_segment(
    aes(
      x = toroidal_null_q025,
      xend = toroidal_null_q975,
      y = sample_order,
      yend = sample_order
    ),
    linewidth = 1.10,
    colour = "#D6D9DC"
  ) +
  geom_point(
    aes(fill = site),
    shape = 21,
    size = 2.25,
    stroke = 0.44,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "PT" = COL_PRIMARY,
      "BM" = COL_BM
    ),
    name = NULL
  ) +
  labs(
    x = "Observed adjusted spatial ρ",
    y = NULL,
    title = "Observed local co-variation versus toroidal null",
    subtitle = "Grey segments: specimen-specific 95% null intervals"
  ) +
  theme_v19() +
  theme(legend.position = "top")

# C. Specimen-specific toroidal-shift P values after multiplicity correction.
spbh <- copy(e8_spatial_bh)
spbh[, sample_order := factor(
  sample,
  levels = sample[
    order(toroidal_two_sided_p)
  ]
)]

ps3c <- ggplot(
  spbh,
  aes(
    x = -log10(toroidal_two_sided_p),
    y = sample_order
  )
) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = 2,
    colour = COL_GOLD,
    linewidth = 0.35
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = -log10(toroidal_two_sided_p),
      y = sample_order,
      yend = sample_order
    ),
    colour = COL_LIGHT,
    linewidth = 0.70
  ) +
  geom_point(
    aes(
      fill = BH_FDR_12 < 0.05
    ),
    shape = 21,
    size = 2.15,
    stroke = 0.42,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = COL_GREY,
      "TRUE" = COL_BM
    ),
    guide = "none"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.02,
    vjust = 1.10,
    label = paste0(
      "Nominal P<0.05: ",
      sum(spbh$toroidal_two_sided_p < 0.05),
      "/12\nBH-FDR<0.05: ",
      sum(spbh$BH_FDR_12 < 0.05),
      "/12\nminimum FDR ",
      sprintf("%.3f", min(spbh$BH_FDR_12))
    ),
    family = FONT,
    size = 1.38,
    colour = COL_MID,
    lineheight = 0.93
  ) +
  labs(
    x = expression(
      -log[10]("toroidal-shift P")
    ),
    y = NULL,
    title = "Autocorrelation-aware local-alignment tests"
  ) +
  theme_v19()

bottom3s <- ps3b | ps3c +
  plot_layout(widths = c(1, 1))

supp3 <- ps3a /
  bottom3s +
  plot_layout(heights = c(1.15, 0.85)) +
  plot_annotation(tag_levels = "A") &
  tag_theme

save_panel(ps3a, "SupplementaryFigure3A_spatial_maps_V20_2", 183, 78)
save_panel(ps3b, "SupplementaryFigure3B_observed_rho_null_V20_2", 89, 68)
save_panel(ps3c, "SupplementaryFigure3C_toroidal_pvalues_V20_2", 89, 68)
save_pub(
  supp3,
  "Supplementary_Figure3_MAJOR_REVISION_V20_2",
  183, 150
)

# ======================================================================
# SUPPLEMENTARY FIGURE 4 — SPATIAL DUPLICATE / EPITHELIAL / SUPPORT AUDIT
# ======================================================================

dup <- copy(e5_dup)
dup[, association_label := fifelse(
  association == "raw_spearman",
  "Raw",
  "Adjusted"
)]
dup[, scenario_label := fcase(
  scenario == "exclude_both_primary_n12",
    "Exclude both duplicates",
  scenario == "retain_PA056_only_n13",
    "Retain PA056 only",
  scenario == "retain_KRAS11_only_n13",
    "Retain KRAS_11 only",
  default = "Retain both (descriptive)"
)]
dup[, scenario_label := factor(
  scenario_label,
  levels = rev(c(
    "Exclude both duplicates",
    "Retain PA056 only",
    "Retain KRAS_11 only",
    "Retain both (descriptive)"
  ))
)]

ps4a <- ggplot(
  dup,
  aes(
    y = scenario_label,
    x = median_rho,
    fill = association_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.30,
    colour = "#A8A8A8"
  ) +
  geom_col(
    position = position_dodge(width = 0.68),
    width = 0.58,
    alpha = 0.90
  ) +
  scale_fill_manual(
    values = c(
      "Raw" = COL_PRIMARY,
      "Adjusted" = COL_TEAL
    ),
    name = NULL
  ) +
  labs(
    x = "Median specimen-level spatial ρ",
    y = NULL,
    title = "Duplicate-payload sensitivity"
  ) +
  theme_v19() +
  theme(legend.position = "top")

epi <- copy(e5_epi_infer)
epi[, association_label := fifelse(
  association == "raw_spearman",
  "Raw",
  "Adjusted"
)]

ps4b <- ggplot(
  epi,
  aes(
    x = association_label,
    y = median_rho,
    fill = association_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.30,
    colour = "#A8A8A8"
  ) +
  geom_col(
    width = 0.54,
    alpha = 0.90
  ) +
  geom_text(
    aes(
      label = paste0(
        n_positive, "/", n_specimens,
        " positive\nP ",
        vapply(
          exact_two_sided_p,
          fmt_p,
          character(1)
        )
      )
    ),
    vjust = -0.45,
    family = FONT,
    size = 1.45,
    colour = COL_DARK,
    lineheight = 0.93
  ) +
  scale_fill_manual(
    values = c(
      "Raw" = COL_PRIMARY,
      "Adjusted" = COL_GOLD
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.28))
  ) +
  labs(
    x = NULL,
    y = "Median spatial ρ",
    title = "Epithelial-enriched proxy"
  ) +
  theme_v19()

scale_dt <- copy(e5_scale)
scale_dt[, site := fifelse(
  grepl("Brain", cohort, ignore.case = TRUE),
  "BM",
  "PT"
)]
scale_dt[, sample_order := factor(
  sample,
  levels = sample[
    order(n_eligible_grids)
  ]
)]

ps4c <- ggplot(
  scale_dt,
  aes(
    y = sample_order,
    x = n_eligible_grids,
    fill = site
  )
) +
  geom_col(
    width = 0.60,
    alpha = 0.90
  ) +
  scale_fill_manual(
    values = c(
      "PT" = COL_PRIMARY,
      "BM" = COL_BM
    ),
    name = NULL
  ) +
  labs(
    x = "Eligible grids",
    y = NULL,
    title = "Eligible spatial grids",
    subtitle = "Grid width = 400 deposited coordinate units"
  ) +
  theme_v19() +
  theme(legend.position = "top")

supp4 <- ps4a | ps4b | ps4c +
  plot_layout(
    widths = c(1.10, 0.72, 1.02)
  ) +
  plot_annotation(tag_levels = "A") &
  tag_theme

save_panel(ps4a, "SupplementaryFigure4A_duplicate_sensitivity_V20_2", 75, 76)
save_panel(ps4b, "SupplementaryFigure4B_epithelial_proxy_V20_2", 47, 76)
save_panel(ps4c, "SupplementaryFigure4C_spatial_support_V20_2", 57, 76)
save_pub(
  supp4,
  "Supplementary_Figure4_MAJOR_REVISION_V20_2",
  183, 86
)

# ======================================================================
# SUPPLEMENTARY FIGURE 5 — DETAILED E9A MATCHING DIAGNOSTICS
# ======================================================================

bal_s5 <- melt(
  e9_balance,
  id.vars = "quantity",
  variable.name = "group",
  value.name = "value"
)
bal_s5[, quantity_short := fcase(
  grepl("rho", quantity, ignore.case = TRUE),
    "|ρ(gene, cycling)|",
  grepl("mean", quantity, ignore.case = TRUE),
    "Mean log2CPM",
  default = "SD log2CPM"
)]
bal_s5[, group := factor(
  group,
  levels = c(
    "observed_module_median",
    "candidate_pool_median"
  ),
  labels = c(
    "Fixed module",
    "Cycle-gene candidate pool"
  )
)]

ps5a <- ggplot(
  bal_s5,
  aes(
    x = group,
    y = value,
    fill = group
  )
) +
  geom_col(
    width = 0.58,
    alpha = 0.90
  ) +
  facet_wrap(
    ~ quantity_short,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Fixed module" = COL_BM,
      "Cycle-gene candidate pool" = COL_GREY
    ),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "Median",
    title = "Matching-feature audit",
    subtitle = "Site labels were not used to construct candidate neighborhoods"
  ) +
  theme_v19() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

neighbors <- copy(e9_neighbors)
neighbors[, module_gene := factor(
  module_gene,
  levels = c(
    "DHFR", "DHODH", "SHMT1", "TYMS", "UMPS"
  )
)]

ps5b <- ggplot(
  neighbors,
  aes(
    x = module_gene,
    y = squared_feature_distance,
    fill = module_gene
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.28,
    linewidth = 0.35
  ) +
  geom_point(
    position = position_jitter(
      width = 0.08,
      height = 0,
      seed = 91
    ),
    size = 0.65,
    alpha = 0.55,
    colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "DHFR" = "#AAB1B7",
      "DHODH" = "#AAB1B7",
      "SHMT1" = "#AAB1B7",
      "TYMS" = COL_BM,
      "UMPS" = COL_GOLD
    ),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = "Squared matching-feature distance",
    title = "Nearest candidate neighborhoods",
    subtitle = "20 candidates per fixed module gene"
  ) +
  theme_v19()

ps5c <- ggplot(
  e9_null_dist,
  aes(
    x = spearman_score_vs_cycling
  )
) +
  geom_histogram(
    bins = 45,
    fill = "#D8DDE1",
    colour = "white",
    linewidth = 0.18
  ) +
  geom_vline(
    xintercept =
      e9_null$observed_spearman_score_vs_cycling[1],
    colour = COL_BM,
    linewidth = 0.90
  ) +
  geom_vline(
    xintercept =
      e9_null$null_score_cycle_rho_median[1],
    colour = COL_DARK,
    linetype = 2,
    linewidth = 0.48
  ) +
  annotate(
    "text",
    x =
      e9_null$observed_spearman_score_vs_cycling[1],
    y = Inf,
    hjust = -0.05,
    vjust = 1.12,
    label = paste0(
      "Fixed five-gene ρ ",
      sprintf(
        "%.3f",
        e9_null$observed_spearman_score_vs_cycling[1]
      ),
      "\nNull median ",
      sprintf(
        "%.3f",
        e9_null$null_score_cycle_rho_median[1]
      )
    ),
    family = FONT,
    size = 1.42,
    colour = COL_BM,
    lineheight = 0.93
  ) +
  labs(
    x = "Spearman ρ of random-set score with cycling fraction",
    y = "Random sets",
    title = "Random-set proliferation coupling"
  ) +
  theme_v19()

top_s5 <- ps5a | ps5b +
  plot_layout(
    widths = c(1.18, 0.82)
  )

supp5 <- top_s5 /
  ps5c +
  plot_layout(
    heights = c(0.95, 1.05)
  ) +
  plot_annotation(tag_levels = "A") &
  tag_theme

save_panel(ps5a, "SupplementaryFigure5A_matching_features_V20_2", 104, 62)
save_panel(ps5b, "SupplementaryFigure5B_candidate_distances_V20_2", 75, 62)
save_panel(ps5c, "SupplementaryFigure5C_random_set_cycling_V20_2", 183, 65)
save_pub(
  supp5,
  "Supplementary_Figure5_MAJOR_REVISION_V20_2",
  183, 132
)

# ======================================================================
# COMPLETION MANIFEST
# ======================================================================

writeLines(
  c(
    paste0(
      "completed_at=",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ),
    "status=PASS_MAJOR_REVISION_FIGURES_V20_2_FINAL_FREEZE",
    "rendering_engine=R_ggplot2_patchwork_only",
    "generative_image_content_used=FALSE",
    "synthetic_data_used=FALSE",
    "new_hypothesis_tested_in_figure_script=FALSE",
    "primary_endpoint_changed=FALSE",
    "gene_set_changed=FALSE",
    "result_driven_exclusion=FALSE",
    "Figure1=cohort_patient_aware_primary_endpoint_plus_nonredundant_robustness",
    "Figure2=specimen_gene_heatmap_plus_expression_and_count_effects",
    "Figure3=CR2_common_support_plus_patient_prefix_and_covariate_imbalance_audit",
    "Figure4=six_panel_formal_proliferation_coupling",
    "Figure5=six_panel_target_panel_specificity_calibration",
    "SupplementaryFigure1=three_panel_nonredundant_fixed_score_robustness",
    "SupplementaryFigure2=four_panel_cohort_composition_and_SMD_audit",
    "SupplementaryFigure3=spatial_maps_observed_rho_and_toroidal_multiplicity",
    "SupplementaryFigure4=duplicate_epithelial_proxy_and_spatial_support_audit",
    "SupplementaryFigure5=detailed_matching_feature_candidate_distance_and_null_diagnostics",
    "spatial_role=supplementary_same_study_contextual_support",
    "specificity_decision=TARGET_PANEL_SPECIFICITY_NOT_ESTABLISHED",
    "visual_design=high_density_multi_geometry_journal_style",
    "V20_1_visual_QA=shortened_panel_titles_fixed_clipping_rebalanced_layouts_no_statistical_changes",
    "V20_2_final_freeze=resolved_remaining_Figure3_legend_Figure4D_E_and_Figure5_annotation_clipping",
    "figure5A_design=coefficient_heatmap_no_sloped_connected_CI_lines",
    "figure4D_design=horizontal_coefficient_bars_to_avoid_cramped_labels",
    "figure5_layout=short_titles_horizontal_heatmap_legend_manual_gene_labels",
    "supplement_layout=short_PT_BM_labels_and_shortened_spatial_support_text",
    "statistical_language_guard=site_associated_not_independent_site_specific_program",
    paste0("output_directory=", OUT),
    paste0("panel_directory=", PANELS)
  ),
  file.path(
    OUT,
    "FIGURE_V20_2_COMPLETE.txt"
  )
)

cat("\n============================================================\n")
cat("V20.2 FINAL FREEZE JOURNAL FIGURE SUITE COMPLETE\n")
cat("============================================================\n")
cat("Output: ", OUT, "\n", sep = "")
cat(
  "Status: PASS_MAJOR_REVISION_FIGURES_V20_2_FINAL_FREEZE\n\n"
)

cat("Main figures:\n")
for (i in 1:5) {
  cat(
    "Figure",
    i,
    "_MAJOR_REVISION_V20_2.png\n",
    sep = ""
  )
}

cat("\nSupplementary figures:\n")
for (i in 1:5) {
  cat(
    "Supplementary_Figure",
    i,
    "_MAJOR_REVISION_V20_2.png\n",
    sep = ""
  )
}

cat("\nManifest:\nFIGURE_V20_2_COMPLETE.txt\n")
