
rm(list = ls())
options(stringsAsFactors = FALSE)

# ----------------------------
# 0. Reproducibility settings
# ----------------------------
SEED <- 20260918
N_PERM <- 9999

# Main prevalence threshold:
# A pathway must be detected (>0) in at least this proportion of samples.
PREVALENCE_MAIN <- 0.50

# Sensitivity thresholds. These are evaluated independently.
PREVALENCE_SENS <- c(0.10, 0.20, 0.50)

# Optional adjusted sensitivity model.
# These should be baseline covariates, not downstream symptoms/mediators.
ADJUSTED_COVARIATES <- c("Age", "BMI", "Education")
RUN_ADJUSTED_MODEL <- TRUE

# Plot settings
ELLIPSE_LEVEL <- 0.95
TOP_N_PATHWAYS <- 8

# Input files
DEMOGRAPHY_FILE <- "Data/Metobolimics_data2.xlsx"
KO_FILE <- "Data/Unigenes_sample_relative_level3.xlsx"

# Output directory
OUTDIR <- "Metagenomics/KO_beta_diversity"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1. Load packages
# ----------------------------
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "stringr",
  "ggplot2", "vegan", "RColorBrewer", "ggrepel"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(vegan)
  library(RColorBrewer)
  library(ggrepel)
})

# Save session information for reproducibility
capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))

# ----------------------------
# 2. Helper functions
# ----------------------------

check_unique_ids <- function(x, label) {
  dup <- unique(x[duplicated(x)])
  if (length(dup) > 0) {
    stop(
      label, " contains duplicated Subid values: ",
      paste(head(dup, 20), collapse = ", "),
      if (length(dup) > 20) " ..." else ""
    )
  }
}

extract_adonis_row <- function(adonis_obj, term = "group") {
  tab <- as.data.frame(adonis_obj)
  if (!term %in% rownames(tab)) {
    stop("Term '", term, "' not found in adonis2 output.")
  }
  out <- tab[term, , drop = FALSE]
  out$Term <- term
  out
}

run_permanova_permdisp <- function(abundance_mat, meta, n_perm = 9999, seed = 1) {
  if (!all(rownames(abundance_mat) == meta$Subid)) {
    stop("Sample order mismatch between abundance matrix and metadata.")
  }

  # Bray-Curtis dissimilarity on transformed abundance data
  bray <- vegan::vegdist(abundance_mat, method = "bray")

  # PERMANOVA: global group difference
  set.seed(seed)
  permanova <- vegan::adonis2(
    bray ~ group,
    data = meta,
    permutations = n_perm
  )

  # PERMDISP: test homogeneity of multivariate dispersion
  # Median is the default/recommended robust center in betadisper.
  bd <- vegan::betadisper(bray, group = meta$group, type = "median")

  set.seed(seed + 1)
  bd_perm <- vegan::permutest(bd, permutations = n_perm)

  permanova_row <- extract_adonis_row(permanova, term = "group")

  dispersion_table <- as.data.frame(bd_perm$tab)
  dispersion_table$Term <- rownames(dispersion_table)

  list(
    bray = bray,
    permanova = permanova,
    permanova_row = permanova_row,
    betadisper = bd,
    betadisper_test = bd_perm,
    dispersion_table = dispersion_table
  )
}

make_site_score_table <- function(ord, meta) {
  # With one binary constraint, the first constrained axis is dbRDA1.
  # The next displayed axis is typically the first unconstrained axis (MDS1).
  sc <- as.data.frame(
    vegan::scores(ord, display = "sites", choices = 1:2, scaling = 1)
  )

  if (ncol(sc) < 2) {
    stop("Fewer than two ordination axes were returned.")
  }

  original_axis_names <- colnames(sc)[1:2]

  sc <- sc[, 1:2, drop = FALSE]
  colnames(sc) <- c("Axis1", "Axis2")
  sc$Subid <- rownames(sc)

  # Explicit keyed join instead of assuming row order
  sc <- sc %>%
    left_join(meta %>% select(Subid, group), by = "Subid")

  if (any(is.na(sc$group))) {
    stop("Some ordination site scores could not be matched to metadata.")
  }

  attr(sc, "axis_names") <- original_axis_names
  sc
}

get_axis_labels <- function(ord, site_scores) {
  # Use the public eigenvals() API rather than directly accessing internal slots.
  eig_all <- vegan::eigenvals(ord, model = "all")

  # dbrda may contain negative eigenvalues for non-Euclidean distances.
  # For a transparent plotting label, percentages are calculated relative
  # to the sum of positive eigenvalues only.
  eig_positive <- eig_all[eig_all > 0]
  denom <- sum(eig_positive)

  axis_names <- attr(site_scores, "axis_names")

  pct1 <- NA_real_
  pct2 <- NA_real_

  if (length(eig_all) >= 1 && denom > 0) {
    pct1 <- 100 * eig_all[1] / denom
  }
  if (length(eig_all) >= 2 && denom > 0) {
    pct2 <- 100 * eig_all[2] / denom
  }

  label1 <- if (is.finite(pct1)) {
    sprintf("%s (%.2f%% of positive inertia)", axis_names[1], pct1)
  } else {
    axis_names[1]
  }

  label2 <- if (is.finite(pct2)) {
    sprintf("%s (%.2f%% of positive inertia)", axis_names[2], pct2)
  } else {
    axis_names[2]
  }

  list(x = label1, y = label2)
}

# Descriptive pathway vectors:
# correlations between each transformed pathway and the two displayed site axes.
# These are for annotation only and are NOT used as differential-abundance tests.
get_top_pathway_vectors <- function(abundance_mat, site_scores, n_top = 8) {
  if (!all(rownames(abundance_mat) == site_scores$Subid)) {
    # Reorder explicitly if needed.
    abundance_mat <- abundance_mat[match(site_scores$Subid, rownames(abundance_mat)), , drop = FALSE]
  }

  x <- site_scores$Axis1
  y <- site_scores$Axis2

  vec_x <- apply(abundance_mat, 2, function(z) {
    suppressWarnings(cor(z, x, use = "pairwise.complete.obs", method = "pearson"))
  })

  vec_y <- apply(abundance_mat, 2, function(z) {
    suppressWarnings(cor(z, y, use = "pairwise.complete.obs", method = "pearson"))
  })

  vec <- data.frame(
    Pathway = colnames(abundance_mat),
    Axis1 = as.numeric(vec_x),
    Axis2 = as.numeric(vec_y),
    stringsAsFactors = FALSE
  )

  vec <- vec %>%
    filter(is.finite(Axis1), is.finite(Axis2)) %>%
    mutate(VectorLength = sqrt(Axis1^2 + Axis2^2)) %>%
    arrange(desc(VectorLength)) %>%
    slice_head(n = n_top)

  vec
}

# Scale annotation arrows to the spread of site scores.
scale_vectors_for_plot <- function(vectors, site_scores, fraction = 0.55) {
  if (nrow(vectors) == 0) return(vectors)

  xr <- diff(range(site_scores$Axis1, na.rm = TRUE))
  yr <- diff(range(site_scores$Axis2, na.rm = TRUE))

  max_vx <- max(abs(vectors$Axis1), na.rm = TRUE)
  max_vy <- max(abs(vectors$Axis2), na.rm = TRUE)

  sx <- if (max_vx > 0) fraction * xr / (2 * max_vx) else 1
  sy <- if (max_vy > 0) fraction * yr / (2 * max_vy) else 1
  s <- min(sx, sy)

  vectors %>%
    mutate(
      Axis1_plot = Axis1 * s,
      Axis2_plot = Axis2 * s
    )
}

run_prevalence_pipeline <- function(raw_abundance, meta, prevalence_threshold,
                                    n_perm = 9999, seed = 1) {
  # Prevalence is calculated on the untransformed nonnegative abundance matrix.
  prevalence <- colMeans(raw_abundance > 0, na.rm = TRUE)
  keep <- names(prevalence)[prevalence >= prevalence_threshold]

  if (length(keep) < 2) {
    stop(
      "Too few pathways retained at prevalence threshold ",
      prevalence_threshold
    )
  }

  filtered <- raw_abundance[, keep, drop = FALSE]

  # Remove invariant zero columns defensively.
  filtered <- filtered[, colSums(filtered, na.rm = TRUE) > 0, drop = FALSE]

  # Square-root transformation reduces dominance of highly abundant pathways.
  transformed <- sqrt(filtered)

  res <- run_permanova_permdisp(
    abundance_mat = transformed,
    meta = meta,
    n_perm = n_perm,
    seed = seed
  )

  p_row <- res$permanova_row

  # Extract PERMDISP p robustly from the first model row.
  disp_tab <- res$dispersion_table
  disp_p_col <- grep("^Pr\\(", colnames(disp_tab), value = TRUE)[1]
  disp_p <- if (!is.na(disp_p_col)) disp_tab[[disp_p_col]][1] else NA_real_

  data.frame(
    PrevalenceThreshold = prevalence_threshold,
    N_Pathways = ncol(transformed),
    PERMANOVA_F = p_row$F,
    PERMANOVA_R2 = p_row$R2,
    PERMANOVA_P = p_row$`Pr(>F)`,
    PERMDISP_P = disp_p,
    stringsAsFactors = FALSE
  )
}

# ----------------------------
# 3. Read input data
# ----------------------------
data_demography <- read_excel(DEMOGRAPHY_FILE, sheet = 1)
metagenomics <- read_excel(KO_FILE, sheet = 1)

required_demo_cols <- c(
  "Subid", "group", "Stage", "Age", "Education", "BMI",
  "PSQI", "MoCA", "SAS", "KMI", "MRS", "SOL",
  "TST", "TIB", "SE", "WASO"
)

missing_demo <- setdiff(required_demo_cols, colnames(data_demography))
if (length(missing_demo) > 0) {
  stop(
    "Demography file is missing columns: ",
    paste(missing_demo, collapse = ", ")
  )
}

required_ko_cols <- c("KO_Pathway_Level3", "Description")
missing_ko <- setdiff(required_ko_cols, colnames(metagenomics))
if (length(missing_ko) > 0) {
  stop(
    "KO file is missing columns: ",
    paste(missing_ko, collapse = ", ")
  )
}

# ----------------------------
# 4. Build annotation table
# ----------------------------
taxonomy_table <- stringr::str_split_fixed(
  as.character(metagenomics$Description),
  ";",
  3
)

colnames(taxonomy_table) <- c(
  "Pathway_first_class",
  "Pathway_second_class",
  "Pathway_name"
)

pathway_annotation <- data.frame(
  KO_Pathway_Level3 = as.character(metagenomics$KO_Pathway_Level3),
  taxonomy_table,
  stringsAsFactors = FALSE
)

# KO names must be unique for use as matrix column names.
# Preserve original KO identifier in pathway_annotation.
ko_original <- as.character(metagenomics$KO_Pathway_Level3)

if (any(is.na(ko_original) | ko_original == "")) {
  stop("KO_Pathway_Level3 contains missing/empty pathway identifiers.")
}

ko_unique <- make.unique(ko_original, sep = "__dup")
pathway_annotation$KO_Unique <- ko_unique

# ----------------------------
# 5. Convert KO table to samples x pathways
# ----------------------------
sample_cols <- setdiff(
  colnames(metagenomics),
  c("KO_Pathway_Level3", "Description")
)

if (length(sample_cols) < 2) {
  stop("No sample abundance columns found in KO file.")
}

ko_numeric <- metagenomics[, sample_cols, drop = FALSE]

# Convert all abundance values to numeric with explicit failure check.
ko_numeric[] <- lapply(ko_numeric, function(x) suppressWarnings(as.numeric(x)))

if (anyNA(ko_numeric)) {
  bad_n <- sum(is.na(as.matrix(ko_numeric)))
  stop(
    "KO abundance table contains ", bad_n,
    " missing/non-numeric values after conversion. ",
    "Resolve these values explicitly rather than replacing them with zero."
  )
}

if (any(as.matrix(ko_numeric) < 0)) {
  stop("Negative abundance values detected; Bray-Curtis abundance analysis requires nonnegative values.")
}

raw_abundance_all <- t(as.matrix(ko_numeric))
colnames(raw_abundance_all) <- ko_unique
rownames(raw_abundance_all) <- sample_cols
storage.mode(raw_abundance_all) <- "numeric"

# ----------------------------
# 6. Metadata QC and sample matching
# ----------------------------
data_demography$Subid <- as.character(data_demography$Subid)
rownames(raw_abundance_all) <- as.character(rownames(raw_abundance_all))

check_unique_ids(data_demography$Subid, "Demography table")
check_unique_ids(rownames(raw_abundance_all), "KO abundance table")

# Keep only the two intended study groups.
meta <- data_demography %>%
  filter(group %in% c("HCs", "PMI")) %>%
  mutate(
    group = factor(group, levels = c("HCs", "PMI"))
  )

# Match by ID and preserve metadata order.
common_ids <- meta$Subid[meta$Subid %in% rownames(raw_abundance_all)]

if (length(common_ids) < 3) {
  stop("Too few matched samples between metadata and KO abundance data.")
}

meta <- meta %>%
  filter(Subid %in% common_ids) %>%
  arrange(match(Subid, common_ids))

raw_abundance <- raw_abundance_all[meta$Subid, , drop = FALSE]

if (!identical(rownames(raw_abundance), meta$Subid)) {
  stop("Internal sample-order alignment failed.")
}

# Exclude samples with all-zero pathway profiles.
sample_sum <- rowSums(raw_abundance, na.rm = TRUE)

sample_qc <- data.frame(
  Subid = meta$Subid,
  group = as.character(meta$group),
  TotalAbundance = sample_sum,
  Include = sample_sum > 0,
  ExclusionReason = ifelse(sample_sum > 0, "", "All-zero KO pathway profile"),
  stringsAsFactors = FALSE
)

write.csv(
  sample_qc,
  file.path(OUTDIR, "sample_QC.csv"),
  row.names = FALSE
)

keep_samples <- sample_qc$Include
meta <- meta[keep_samples, , drop = FALSE]
raw_abundance <- raw_abundance[keep_samples, , drop = FALSE]

if (nlevels(droplevels(meta$group)) != 2) {
  stop("Both HCs and PMI groups must remain after QC.")
}
meta$group <- droplevels(meta$group)

# Save final subject order.
write.csv(
  meta,
  file.path(OUTDIR, "analysis_metadata_and_subject_order.csv"),
  row.names = FALSE
)

# ----------------------------
# 7. Main prevalence filtering
# ----------------------------
pathway_prevalence <- colMeans(raw_abundance > 0)

prevalence_table <- data.frame(
  KO_Unique = names(pathway_prevalence),
  Prevalence = as.numeric(pathway_prevalence),
  Keep_Main = pathway_prevalence >= PREVALENCE_MAIN,
  stringsAsFactors = FALSE
) %>%
  left_join(pathway_annotation, by = "KO_Unique")

write.csv(
  prevalence_table,
  file.path(OUTDIR, "pathway_prevalence.csv"),
  row.names = FALSE
)

keep_pathways <- names(pathway_prevalence)[
  pathway_prevalence >= PREVALENCE_MAIN
]

if (length(keep_pathways) < 2) {
  stop("Too few pathways retained at the main prevalence threshold.")
}

abundance_main_raw <- raw_abundance[, keep_pathways, drop = FALSE]
abundance_main_raw <- abundance_main_raw[
  ,
  colSums(abundance_main_raw) > 0,
  drop = FALSE
]

# Primary transformation
abundance_main <- sqrt(abundance_main_raw)

# ----------------------------
# 8. Main PERMANOVA + PERMDISP
# ----------------------------
main_tests <- run_permanova_permdisp(
  abundance_mat = abundance_main,
  meta = meta,
  n_perm = N_PERM,
  seed = SEED
)

write.csv(
  as.data.frame(main_tests$permanova),
  file.path(OUTDIR, "PERMANOVA_main.csv"),
  row.names = TRUE
)

write.csv(
  main_tests$dispersion_table,
  file.path(OUTDIR, "PERMDISP_main.csv"),
  row.names = FALSE
)

# ----------------------------
# 9. Main dbRDA ordination
# ----------------------------
# vegan currently recommends dbrda() over the older capscale() implementation.
ord_main <- vegan::dbrda(
  abundance_main ~ group,
  data = meta,
  distance = "bray"
)

set.seed(SEED + 2)
dbrda_overall_test <- anova(
  ord_main,
  permutations = N_PERM
)

set.seed(SEED + 3)
dbrda_axis_test <- anova(
  ord_main,
  by = "axis",
  permutations = N_PERM
)

write.csv(
  as.data.frame(dbrda_overall_test),
  file.path(OUTDIR, "dbRDA_overall_permutation_test.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(dbrda_axis_test),
  file.path(OUTDIR, "dbRDA_axis_permutation_test.csv"),
  row.names = TRUE
)

# Site coordinates: untouched observed ordination scores.
site_scores <- make_site_score_table(ord_main, meta)

write.csv(
  site_scores,
  file.path(OUTDIR, "dbRDA_site_scores.csv"),
  row.names = FALSE
)

axis_labels <- get_axis_labels(ord_main, site_scores)

# ----------------------------
# 10. Descriptive top pathway vectors
# ----------------------------
# These vectors describe correlation with displayed ordination coordinates.
# They are NOT differential abundance tests and should be described as
# "pathways most strongly associated with the displayed ordination axes".
top_vectors <- get_top_pathway_vectors(
  abundance_mat = abundance_main,
  site_scores = site_scores,
  n_top = TOP_N_PATHWAYS
)

top_vectors <- top_vectors %>%
  left_join(pathway_annotation, by = c("Pathway" = "KO_Unique"))

write.csv(
  top_vectors,
  file.path(OUTDIR, "top_pathway_ordination_vectors.csv"),
  row.names = FALSE
)

top_vectors_plot <- scale_vectors_for_plot(
  vectors = top_vectors,
  site_scores = site_scores
)

# ----------------------------
# 11. Plot dbRDA
# ----------------------------
colors <- c(
  "HCs" = "#1f77b4",
  "PMI" = "#d62728"
)

permanova_main <- main_tests$permanova_row
permanova_p <- permanova_main$`Pr(>F)`
permanova_r2 <- permanova_main$R2

# PERMDISP p-value
disp_tab <- main_tests$dispersion_table
disp_p_col <- grep("^Pr\\(", colnames(disp_tab), value = TRUE)[1]
permdisp_p <- if (!is.na(disp_p_col)) disp_tab[[disp_p_col]][1] else NA_real_

# Group centers are computed from observed coordinates only.
group_centers <- site_scores %>%
  group_by(group) %>%
  summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2),
    .groups = "drop"
  )

p <- ggplot() +
  # Lines from observed sample locations to observed group centers.
  geom_segment(
    data = site_scores %>%
      left_join(
        group_centers,
        by = "group",
        suffix = c("", "_center")
      ),
    aes(
      x = Axis1_center, y = Axis2_center,
      xend = Axis1, yend = Axis2,
      color = group
    ),
    alpha = 0.25,
    linewidth = 0.25
  ) +

  # 95% data ellipse based on observed coordinates.
  stat_ellipse(
    data = site_scores,
    aes(x = Axis1, y = Axis2, fill = group),
    geom = "polygon",
    level = ELLIPSE_LEVEL,
    alpha = 0.15,
    color = NA,
    type = "t"
  ) +

  geom_point(
    data = site_scores,
    aes(x = Axis1, y = Axis2, color = group),
    size = 1.5,
    alpha = 0.85
  ) +

  geom_label(
    data = group_centers,
    aes(x = Axis1, y = Axis2, label = group),
    fill = "white",
    color = "black",
    size = 3.5,
    fontface = "bold",
    label.size = 0.3
  ) +

  geom_vline(
    xintercept = 0,
    linewidth = 0.3,
    color = "grey75"
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    color = "grey75"
  ) +

  scale_color_manual(values = colors, name = NULL) +
  scale_fill_manual(values = colors, name = NULL) +

  labs(
    x = axis_labels$x,
    y = axis_labels$y,
    title = "KEGG orthologue (KO) pathway profiles",
    subtitle = sprintf(
      "PERMANOVA: R² = %.3f, P = %.4g; PERMDISP: P = %.4g",
      permanova_r2,
      permanova_p,
      permdisp_p
    )
  ) +

  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Optional descriptive pathway arrows.
# Comment out this block if a cleaner manuscript figure is preferred.
if (nrow(top_vectors_plot) > 0) {
  label_var <- ifelse(
    !is.na(top_vectors_plot$Pathway_name) &
      top_vectors_plot$Pathway_name != "",
    top_vectors_plot$Pathway_name,
    top_vectors_plot$Pathway
  )

  top_vectors_plot$PlotLabel <- label_var

  p <- p +
    geom_segment(
      data = top_vectors_plot,
      aes(
        x = 0, y = 0,
        xend = Axis1_plot,
        yend = Axis2_plot
      ),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.15, "cm")),
      linewidth = 0.35,
      color = "grey35",
      alpha = 0.65
    ) +
    ggrepel::geom_text_repel(
      data = top_vectors_plot,
      aes(
        x = Axis1_plot,
        y = Axis2_plot,
        label = PlotLabel
      ),
      inherit.aes = FALSE,
      size = 2.8,
      color = "grey20",
      max.overlaps = Inf,
      min.segment.length = 0
    )
}

ggsave(
  filename = file.path(OUTDIR, "dbRDA_KO_main.pdf"),
  plot = p,
  width = 8,
  height = 6.5,
  device = cairo_pdf
)

ggsave(
  filename = file.path(OUTDIR, "dbRDA_KO_main.png"),
  plot = p,
  width = 8,
  height = 6.5,
  dpi = 600
)

# ----------------------------
# 12. Covariate-adjusted sensitivity analysis
# ----------------------------
if (RUN_ADJUSTED_MODEL) {

  missing_cov <- setdiff(ADJUSTED_COVARIATES, colnames(meta))

  if (length(missing_cov) > 0) {
    warning(
      "Adjusted model skipped. Missing covariates: ",
      paste(missing_cov, collapse = ", ")
    )
  } else {

    adjusted_vars <- c("Subid", "group", ADJUSTED_COVARIATES)

    adjusted_complete <- complete.cases(meta[, adjusted_vars, drop = FALSE])

    meta_adj <- meta[adjusted_complete, , drop = FALSE]
    abundance_adj <- abundance_main[meta_adj$Subid, , drop = FALSE]

    if (nrow(meta_adj) >= 10 &&
        nlevels(droplevels(meta_adj$group)) == 2) {

      meta_adj$group <- droplevels(meta_adj$group)

      # Construct formula programmatically:
      # Bray distance ~ group + Age + BMI + Education
      bray_adj <- vegan::vegdist(abundance_adj, method = "bray")

      rhs <- paste(
        c("group", ADJUSTED_COVARIATES),
        collapse = " + "
      )

      perm_formula <- as.formula(
        paste("bray_adj ~", rhs)
      )

      set.seed(SEED + 10)
      permanova_adj <- vegan::adonis2(
        perm_formula,
        data = meta_adj,
        permutations = N_PERM,
        by = "margin"
      )

      write.csv(
        as.data.frame(permanova_adj),
        file.path(OUTDIR, "PERMANOVA_covariate_adjusted_sensitivity.csv"),
        row.names = TRUE
      )

      # Partial dbRDA visualization:
      # isolate group-associated ordination after conditioning on covariates.
      cond <- paste(
        ADJUSTED_COVARIATES,
        collapse = " + "
      )

      dbrda_formula <- as.formula(
        paste(
          "abundance_adj ~ group + Condition(",
          cond,
          ")"
        )
      )

      ord_adj <- vegan::dbrda(
        dbrda_formula,
        data = meta_adj,
        distance = "bray"
      )

      set.seed(SEED + 11)
      dbrda_adj_test <- anova(
        ord_adj,
        permutations = N_PERM
      )

      write.csv(
        as.data.frame(dbrda_adj_test),
        file.path(OUTDIR, "dbRDA_covariate_adjusted_sensitivity.csv"),
        row.names = TRUE
      )

      write.csv(
        meta_adj,
        file.path(OUTDIR, "adjusted_model_subjects.csv"),
        row.names = FALSE
      )

    } else {
      warning(
        "Adjusted model skipped because too few complete cases ",
        "or one group was lost after complete-case filtering."
      )
    }
  }
}

# ----------------------------
# 13. Prevalence-threshold sensitivity analysis
# ----------------------------
sens_results <- lapply(
  seq_along(PREVALENCE_SENS),
  function(i) {
    thr <- PREVALENCE_SENS[i]

    run_prevalence_pipeline(
      raw_abundance = raw_abundance,
      meta = meta,
      prevalence_threshold = thr,
      n_perm = N_PERM,
      seed = SEED + 100 + i
    )
  }
)

sens_table <- bind_rows(sens_results)

write.csv(
  sens_table,
  file.path(OUTDIR, "prevalence_threshold_sensitivity.csv"),
  row.names = FALSE
)

# ----------------------------
# 14. Compact main-results summary
# ----------------------------
main_summary <- data.frame(
  Analysis = c(
    "Sample size",
    "HCs n",
    "PMI n",
    "Main prevalence threshold",
    "Retained pathways",
    "PERMANOVA R2",
    "PERMANOVA P",
    "PERMDISP P",
    "Permutations"
  ),
  Value = c(
    nrow(meta),
    sum(meta$group == "HCs"),
    sum(meta$group == "PMI"),
    PREVALENCE_MAIN,
    ncol(abundance_main),
    permanova_r2,
    permanova_p,
    permdisp_p,
    N_PERM
  ),
  stringsAsFactors = FALSE
)

write.csv(
  main_summary,
  file.path(OUTDIR, "main_analysis_summary.csv"),
  row.names = FALSE
)

# ----------------------------
# 15. Console output
# ----------------------------
cat("\n========================================\n")
cat("KO pathway analysis completed\n")
cat("========================================\n")
cat("N =", nrow(meta), "\n")
cat("HCs =", sum(meta$group == "HCs"), "\n")
cat("PMI =", sum(meta$group == "PMI"), "\n")
cat("Main prevalence threshold =", PREVALENCE_MAIN, "\n")
cat("Retained pathways =", ncol(abundance_main), "\n")
cat(
  sprintf(
    "PERMANOVA: R2 = %.4f, P = %.6g\n",
    permanova_r2,
    permanova_p
  )
)
cat(
  sprintf(
    "PERMDISP: P = %.6g\n",
    permdisp_p
  )
)
cat("Outputs saved to:", OUTDIR, "\n")
cat("No artificial jitter/noise was added to ordination scores.\n")
cat("========================================\n")
'''

path = Path("/mnt/data/KO_pathway_beta_diversity_review_ready.R")
path.write_text(code, encoding="utf-8")
print(f"Created: {path}")
print(f"Lines: {len(code.splitlines())}")
