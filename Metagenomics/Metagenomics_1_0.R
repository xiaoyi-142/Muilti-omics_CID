# ============================================================
# Metagenomics_1_0_review_ready.R
# Species-level metagenomic biomarker analysis: PMI vs HCs
#
# Purpose
#   1) Import and audit species-level abundance/taxonomy data
#   2) Match samples by Subid (no positional matching)
#   3) Apply a pre-specified prevalence filter
#   4) Run reproducible LEfSe
#   5) Independently compute Wilcoxon tests with BH-FDR across
#      ALL prevalence-filtered features
#   6) Export marker tables, abundances, QC/audit files, and
#      a reproducible family-level summary figure
#


rm(list = ls())
options(stringsAsFactors = FALSE)

# ----------------------------
# 0. Analysis parameters
# ----------------------------
SEED <- 20260918L
PREVALENCE_CUTOFF <- 0.10
LEFSE_LDA_CUTOFF <- 2.0
LEFSE_KW_CUTOFF <- 0.05
LEFSE_WILCOX_CUTOFF <- 0.05
LEFSE_BOOTSTRAP_N <- 100L
LEFSE_BOOTSTRAP_FRACTION <- 2 / 3

GROUP_LEVELS <- c("HCs", "PMI")

DEMOGRAPHY_FILE <- "Data/Metobolimics_data2.xlsx"
METAGENOMICS_FILE <- "Data/Unigenes_relative_s.xlsx"

OUTDIR <- "Results/Metagenomics_species"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1. Packages
# ----------------------------
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "stringr", "tibble",
  "phyloseq", "microbiomeMarker", "ggplot2", "scales"
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
  library(tibble)
  library(phyloseq)
  library(microbiomeMarker)
  library(ggplot2)
  library(scales)
})

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))

# ----------------------------
# 2. Helper functions
# ----------------------------
check_unique <- function(x, label) {
  dup <- unique(x[duplicated(x)])
  if (length(dup) > 0) {
    stop(
      label, " contains duplicated IDs: ",
      paste(head(dup, 20), collapse = ", "),
      if (length(dup) > 20) " ..." else ""
    )
  }
}

clean_taxon_rank <- function(x) {
  x <- as.character(x)
  x <- sub("^[A-Za-z]__", "", x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

safe_first <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & trimws(x) != ""]
  if (length(x) == 0) "Unclassified" else x[1]
}

safe_wilcox <- function(x, group) {
  x_hc <- x[group == "HCs"]
  x_pmi <- x[group == "PMI"]

  if (length(x_hc) == 0 || length(x_pmi) == 0) {
    return(NA_real_)
  }

  if (length(unique(c(x_hc, x_pmi))) <= 1) {
    return(1)
  }

  suppressWarnings(
    stats::wilcox.test(
      x_pmi,
      x_hc,
      alternative = "two.sided",
      exact = FALSE,
      correct = FALSE
    )$p.value
  )
}

# ----------------------------
# 3. Read data
# ----------------------------
data_demography <- read_excel(DEMOGRAPHY_FILE, sheet = 1)
metagenomics_raw <- read_excel(METAGENOMICS_FILE, sheet = 1)

required_demo <- c("Subid", "group")
missing_demo <- setdiff(required_demo, names(data_demography))
if (length(missing_demo) > 0) {
  stop(
    "Demography file is missing columns: ",
    paste(missing_demo, collapse = ", ")
  )
}

required_meta <- c("Taxonomy", "species")
missing_meta <- setdiff(required_meta, names(metagenomics_raw))
if (length(missing_meta) > 0) {
  stop(
    "Metagenomics file is missing columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

data_demography <- data_demography %>%
  mutate(
    Subid = trimws(as.character(Subid)),
    group = trimws(as.character(group))
  ) %>%
  filter(group %in% GROUP_LEVELS) %>%
  mutate(group = factor(group, levels = GROUP_LEVELS))

check_unique(data_demography$Subid, "Demography table")

# Use ONLY columns whose names are actual study Subids.
# This avoids accidentally treating annotation columns as abundance columns.
sample_cols <- names(metagenomics_raw)[
  names(metagenomics_raw) %in% data_demography$Subid
]

if (length(sample_cols) < 2) {
  stop(
    "Fewer than two abundance columns matched metadata Subids. ",
    "Check sample names in the metagenomics Excel file."
  )
}

# Audit unmatched subjects.
sample_match_audit <- data.frame(
  Subid = data_demography$Subid,
  group = as.character(data_demography$group),
  HasMetagenomics = data_demography$Subid %in% sample_cols,
  stringsAsFactors = FALSE
)

write.csv(
  sample_match_audit,
  file.path(OUTDIR, "sample_match_audit.csv"),
  row.names = FALSE
)

# ----------------------------
# 4. Parse taxonomy for EVERY input row
# ----------------------------
taxonomy_split <- stringr::str_split_fixed(
  as.character(metagenomics_raw$Taxonomy),
  ";",
  7
)

taxonomy_df <- as.data.frame(
  taxonomy_split,
  stringsAsFactors = FALSE
)

names(taxonomy_df) <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species_from_taxonomy"
)

taxonomy_df[] <- lapply(taxonomy_df, clean_taxon_rank)

# Prefer explicit s__ label from the `species` column.
species_from_species_col <- stringr::str_extract(
  as.character(metagenomics_raw$species),
  "(?<=s__)[^;]*$"
)
species_from_species_col <- trimws(species_from_species_col)

# Fall back to the species rank parsed from Taxonomy.
species_name <- ifelse(
  !is.na(species_from_species_col) &
    species_from_species_col != "",
  species_from_species_col,
  taxonomy_df$Species_from_taxonomy
)

species_name <- trimws(as.character(species_name))

# Explicit feature exclusion rule.
# "uncultured X sp." is retained; only fully unassigned generic labels are removed.
generic_unassigned <- grepl(
  "^(other|others|unknown|unclassified|unassigned|not assigned)$",
  species_name,
  ignore.case = TRUE
)

valid_feature <- !is.na(species_name) &
  species_name != "" &
  !generic_unassigned

feature_audit <- data.frame(
  InputRow = seq_len(nrow(metagenomics_raw)),
  OriginalSpeciesField = as.character(metagenomics_raw$species),
  Species = species_name,
  Include = valid_feature,
  ExclusionReason = ifelse(
    valid_feature,
    "",
    "Missing or generic unassigned species label"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  feature_audit,
  file.path(OUTDIR, "feature_exclusion_audit.csv"),
  row.names = FALSE
)

if (sum(valid_feature) < 2) {
  stop("Too few valid species-level features remain after explicit filtering.")
}

# Keep original rows as distinct features.
# If display species names are duplicated, assign stable unique FeatureIDs
# rather than silently summing biologically ambiguous entries.
feature_id <- make.unique(species_name[valid_feature], sep = "__dup")

taxonomy_valid <- taxonomy_df[valid_feature, , drop = FALSE] %>%
  mutate(
    FeatureID = feature_id,
    Species = species_name[valid_feature]
  ) %>%
  select(
    FeatureID, Kingdom, Phylum, Class, Order,
    Family, Genus, Species
  )

# Fill missing taxonomy ranks transparently for phyloseq compatibility.
taxonomy_for_ps <- taxonomy_valid
for (nm in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  taxonomy_for_ps[[nm]][
    is.na(taxonomy_for_ps[[nm]]) |
      trimws(taxonomy_for_ps[[nm]]) == ""
  ] <- "Unclassified"
}

# Duplicate display-name audit.
duplicate_species_audit <- taxonomy_valid %>%
  count(Species, name = "N_Features") %>%
  filter(N_Features > 1)

write.csv(
  duplicate_species_audit,
  file.path(OUTDIR, "duplicate_species_name_audit.csv"),
  row.names = FALSE
)

# ----------------------------
# 5. Build feature x sample abundance matrix
# ----------------------------
abundance_feature_sample <- metagenomics_raw[
  valid_feature,
  sample_cols,
  drop = FALSE
]

abundance_feature_sample[] <- lapply(
  abundance_feature_sample,
  function(x) suppressWarnings(as.numeric(x))
)

abundance_feature_sample <- as.matrix(abundance_feature_sample)
rownames(abundance_feature_sample) <- feature_id
storage.mode(abundance_feature_sample) <- "numeric"

if (anyNA(abundance_feature_sample)) {
  stop(
    "Missing/non-numeric abundance values were detected. ",
    "Resolve them explicitly; this script will not replace them with zero."
  )
}

if (any(!is.finite(abundance_feature_sample))) {
  stop("Non-finite abundance values (Inf/-Inf) were detected.")
}

if (any(abundance_feature_sample < 0)) {
  stop("Negative abundance values were detected.")
}

# Convert to samples x features.
abundance_sample_feature_all <- t(abundance_feature_sample)

# ----------------------------
# 6. Match samples by Subid and preserve explicit order
# ----------------------------
meta <- data_demography %>%
  filter(Subid %in% rownames(abundance_sample_feature_all))

if (nrow(meta) < 3) {
  stop("Too few matched samples for analysis.")
}

abundance_sample_feature_all <- abundance_sample_feature_all[
  meta$Subid,
  ,
  drop = FALSE
]

if (!identical(rownames(abundance_sample_feature_all), meta$Subid)) {
  stop("Sample order mismatch after Subid matching.")
}

# Exclude all-zero sample profiles explicitly.
sample_total <- rowSums(abundance_sample_feature_all)

sample_qc <- data.frame(
  Subid = meta$Subid,
  group = as.character(meta$group),
  TotalInputAbundance = sample_total,
  Include = sample_total > 0,
  ExclusionReason = ifelse(
    sample_total > 0,
    "",
    "All-zero species abundance profile"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sample_qc,
  file.path(OUTDIR, "sample_QC.csv"),
  row.names = FALSE
)

keep_samples <- sample_qc$Include
meta <- meta[keep_samples, , drop = FALSE]
abundance_sample_feature_all <- abundance_sample_feature_all[
  keep_samples,
  ,
  drop = FALSE
]
meta$group <- droplevels(meta$group)

if (!all(GROUP_LEVELS %in% levels(meta$group))) {
  stop("Both HCs and PMI must remain after sample QC.")
}

write.csv(
  meta,
  file.path(OUTDIR, "analysis_subject_order.csv"),
  row.names = FALSE
)

# ----------------------------
# 7. Prevalence filter
# ----------------------------
prevalence <- colMeans(abundance_sample_feature_all > 0)

prevalence_table <- data.frame(
  FeatureID = names(prevalence),
  Prevalence = as.numeric(prevalence),
  Keep = prevalence >= PREVALENCE_CUTOFF,
  stringsAsFactors = FALSE
) %>%
  left_join(taxonomy_valid, by = "FeatureID")

write.csv(
  prevalence_table,
  file.path(OUTDIR, "species_prevalence.csv"),
  row.names = FALSE
)

keep_features <- names(prevalence)[
  prevalence >= PREVALENCE_CUTOFF
]

if (length(keep_features) < 2) {
  stop(
    "Too few features retained at prevalence cutoff ",
    PREVALENCE_CUTOFF
  )
}

abundance_filtered_raw <- abundance_sample_feature_all[
  ,
  keep_features,
  drop = FALSE
]

# Remove zero-total features defensively.
abundance_filtered_raw <- abundance_filtered_raw[
  ,
  colSums(abundance_filtered_raw) > 0,
  drop = FALSE
]

# TSS relative abundance is used for descriptive plots/statistical summaries.
# LEfSe below performs its own explicit CPM normalization.
sample_sums <- rowSums(abundance_filtered_raw)

if (any(sample_sums <= 0)) {
  stop("Zero-total samples remain after prevalence filtering.")
}

relative_abundance <- sweep(
  abundance_filtered_raw,
  1,
  sample_sums,
  "/"
)

# ----------------------------
# 8. Construct phyloseq object
# ----------------------------
otu_mat <- t(abundance_filtered_raw)

sample_df <- data.frame(
  group = meta$group,
  row.names = meta$Subid,
  check.names = FALSE
)

taxonomy_kept <- taxonomy_for_ps %>%
  filter(FeatureID %in% colnames(abundance_filtered_raw)) %>%
  arrange(match(FeatureID, colnames(abundance_filtered_raw)))

if (!identical(taxonomy_kept$FeatureID, colnames(abundance_filtered_raw))) {
  stop("Taxonomy and abundance feature order mismatch.")
}

tax_mat <- as.matrix(
  taxonomy_kept[
    ,
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    drop = FALSE
  ]
)
rownames(tax_mat) <- taxonomy_kept$FeatureID

ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sample_df)
)

# ----------------------------
# 9. Reproducible LEfSe
# ----------------------------
# taxa_rank = "none": features are already species-level rows and are not
# re-summarized by microbiomeMarker.
#
# CPM is stated explicitly to match the standard LEfSe-style normalization.
# A fixed seed makes the bootstrap LDA step reproducible.
set.seed(SEED)

lefse_res <- microbiomeMarker::run_lefse(
  ps,
  group = "group",
  taxa_rank = "none",
  transform = "identity",
  norm = "CPM",
  kw_cutoff = LEFSE_KW_CUTOFF,
  wilcoxon_cutoff = LEFSE_WILCOX_CUTOFF,
  lda_cutoff = LEFSE_LDA_CUTOFF,
  bootstrap_n = LEFSE_BOOTSTRAP_N,
  bootstrap_fraction = LEFSE_BOOTSTRAP_FRACTION,
  multigrp_strat = FALSE,
  strict = "0"
)

marker_obj <- microbiomeMarker::marker_table(lefse_res)

# IMPORTANT:
# In current microbiomeMarker LEfSe implementation, package `padj` is not
# a BH-adjusted q-value. Do not use it as FDR. We calculate BH-FDR
# independently across all prevalence-filtered features in Step 10.
if (is.null(marker_obj) || nrow(as.data.frame(marker_obj)) == 0) {
  warning("LEfSe identified no markers at the specified thresholds.")

  marker_df <- data.frame(
    FeatureID = character(0),
    LEfSe_EnrichedGroup = character(0),
    LEfSe_LDA = numeric(0),
    LEfSe_KW_P = numeric(0),
    stringsAsFactors = FALSE
  ) %>%
    left_join(taxonomy_valid, by = "FeatureID")
} else {
  marker_raw <- as.data.frame(marker_obj)

  if ("padj" %in% names(marker_raw)) {
    marker_raw$padj <- NULL
  }

  marker_df <- marker_raw %>%
    rename(
      FeatureID = feature,
      LEfSe_EnrichedGroup = enrich_group,
      LEfSe_LDA = ef_lda,
      LEfSe_KW_P = pvalue
    ) %>%
    left_join(taxonomy_valid, by = "FeatureID")
}

# ----------------------------
# 10. Independent all-feature Wilcoxon + BH-FDR
# ----------------------------
# This is intentionally performed on ALL prevalence-filtered features,
# not only on LEfSe-selected markers, to avoid post-selection p-value correction.
feature_stats <- lapply(
  colnames(relative_abundance),
  function(fid) {
    x <- relative_abundance[, fid]

    data.frame(
      FeatureID = fid,
      Median_HCs = median(x[meta$group == "HCs"], na.rm = TRUE),
      Median_PMI = median(x[meta$group == "PMI"], na.rm = TRUE),
      Mean_HCs = mean(x[meta$group == "HCs"], na.rm = TRUE),
      Mean_PMI = mean(x[meta$group == "PMI"], na.rm = TRUE),
      Prevalence_HCs = mean(x[meta$group == "HCs"] > 0, na.rm = TRUE),
      Prevalence_PMI = mean(x[meta$group == "PMI"] > 0, na.rm = TRUE),
      Wilcoxon_P = safe_wilcox(x, meta$group),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

feature_stats$Wilcoxon_FDR_BH <- p.adjust(
  feature_stats$Wilcoxon_P,
  method = "BH"
)

feature_stats <- feature_stats %>%
  left_join(taxonomy_valid, by = "FeatureID")

write.csv(
  feature_stats,
  file.path(OUTDIR, "all_filtered_species_wilcoxon_BH.csv"),
  row.names = FALSE
)

# Join valid BH q-values to LEfSe marker output.
marker_df <- marker_df %>%
  left_join(
    feature_stats %>%
      select(
        FeatureID,
        Median_HCs, Median_PMI,
        Mean_HCs, Mean_PMI,
        Prevalence_HCs, Prevalence_PMI,
        Wilcoxon_P, Wilcoxon_FDR_BH
      ),
    by = "FeatureID"
  ) %>%
  mutate(
    FDR_supported = !is.na(Wilcoxon_FDR_BH) &
      Wilcoxon_FDR_BH < 0.05
  ) %>%
  arrange(
    LEfSe_EnrichedGroup,
    desc(LEfSe_LDA),
    Wilcoxon_FDR_BH
  )

write.csv(
  marker_df,
  file.path(OUTDIR, "LEfSe_marker_results_with_taxonomy_and_FDR.csv"),
  row.names = FALSE
)

write.csv(
  marker_df %>% filter(FDR_supported),
  file.path(OUTDIR, "LEfSe_markers_FDR_supported.csv"),
  row.names = FALSE
)

# ----------------------------
# 11. Export abundance matrices for downstream analysis
# ----------------------------
# LEfSe marker set (standard LEfSe definition).
lefse_features <- marker_df$FeatureID

lefse_abundance <- data.frame(
  Subid = meta$Subid,
  group = as.character(meta$group),
  relative_abundance[
    ,
    lefse_features,
    drop = FALSE
  ],
  check.names = FALSE
)

write.csv(
  lefse_abundance,
  file.path(OUTDIR, "metagenomics_relative_abundance_LEfSe_markers.csv"),
  row.names = FALSE
)

# More conservative FDR-supported LEfSe subset.
fdr_features <- marker_df %>%
  filter(FDR_supported) %>%
  pull(FeatureID)

fdr_abundance <- data.frame(
  Subid = meta$Subid,
  group = as.character(meta$group),
  relative_abundance[
    ,
    fdr_features,
    drop = FALSE
  ],
  check.names = FALSE
)

write.csv(
  fdr_abundance,
  file.path(OUTDIR, "metagenomics_relative_abundance_LEfSe_FDR_supported.csv"),
  row.names = FALSE
)

# ----------------------------
# 12. Programmatic family-level marker summary
# ----------------------------
# No hand-edited CSV is used.
family_count <- marker_df %>%
  mutate(
    Family = ifelse(
      is.na(Family) | Family == "",
      "Unclassified",
      Family
    ),
    Phylum = ifelse(
      is.na(Phylum) | Phylum == "",
      "Unclassified",
      Phylum
    )
  ) %>%
  count(
    Phylum,
    Family,
    LEfSe_EnrichedGroup,
    name = "N_LEfSe_markers"
  ) %>%
  complete(
    Phylum,
    Family,
    LEfSe_EnrichedGroup = GROUP_LEVELS,
    fill = list(N_LEfSe_markers = 0)
  )

write.csv(
  family_count,
  file.path(OUTDIR, "LEfSe_marker_count_by_family.csv"),
  row.names = FALSE
)

if (nrow(marker_df) > 0) {
  family_plot_df <- family_count %>%
    mutate(
      SignedCount = ifelse(
        LEfSe_EnrichedGroup == "PMI",
        -N_LEfSe_markers,
        N_LEfSe_markers
      ),
      GroupLabel = factor(
        LEfSe_EnrichedGroup,
        levels = GROUP_LEVELS
      )
    )

  group_colors <- c(
    "HCs" = "#457B9D",
    "PMI" = "#E76F51"
  )

  p_family <- ggplot(
    family_plot_df,
    aes(
      x = SignedCount,
      y = Family,
      fill = GroupLabel
    )
  ) +
    geom_vline(
      xintercept = 0,
      color = "grey45",
      linewidth = 0.4
    ) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = group_colors,
      name = "LEfSe enriched group"
    ) +
    scale_x_continuous(
      labels = abs,
      expand = expansion(mult = c(0.10, 0.10))
    ) +
    facet_grid(
      Phylum ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    labs(
      x = "Number of LEfSe marker features",
      y = NULL,
      title = "Species-level LEfSe markers by taxonomic family"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey95"),
      axis.text.y = element_text(face = "italic"),
      legend.position = "top"
    )

  ggsave(
    file.path(OUTDIR, "LEfSe_marker_family_summary.pdf"),
    p_family,
    width = 9,
    height = max(6, 0.22 * length(unique(family_plot_df$Family)) + 2),
    limitsize = FALSE
  )
}

# ----------------------------
# 13. Save a single analysis object for Script 1_1
# ----------------------------
analysis_object <- list(
  parameters = list(
    seed = SEED,
    prevalence_cutoff = PREVALENCE_CUTOFF,
    lefse_lda_cutoff = LEFSE_LDA_CUTOFF,
    lefse_kw_cutoff = LEFSE_KW_CUTOFF,
    lefse_wilcox_cutoff = LEFSE_WILCOX_CUTOFF,
    lefse_bootstrap_n = LEFSE_BOOTSTRAP_N,
    lefse_bootstrap_fraction = LEFSE_BOOTSTRAP_FRACTION,
    group_levels = GROUP_LEVELS
  ),
  meta = meta,
  taxonomy = taxonomy_valid,
  relative_abundance = relative_abundance,
  raw_filtered_abundance = abundance_filtered_raw,
  prevalence = prevalence_table,
  all_feature_stats = feature_stats,
  lefse_markers = marker_df
)

saveRDS(
  analysis_object,
  file.path(OUTDIR, "metagenomics_species_analysis_object.rds")
)

# ----------------------------
# 14. Analysis parameter manifest
# ----------------------------
parameter_manifest <- data.frame(
  Parameter = c(
    "SEED",
    "PREVALENCE_CUTOFF",
    "LEFSE_LDA_CUTOFF",
    "LEFSE_KW_CUTOFF",
    "LEFSE_WILCOX_CUTOFF",
    "LEFSE_BOOTSTRAP_N",
    "LEFSE_BOOTSTRAP_FRACTION",
    "LEFSE_NORMALIZATION",
    "LEFSE_TAXA_RANK",
    "MULTIPLE_TEST_CORRECTION"
  ),
  Value = c(
    SEED,
    PREVALENCE_CUTOFF,
    LEFSE_LDA_CUTOFF,
    LEFSE_KW_CUTOFF,
    LEFSE_WILCOX_CUTOFF,
    LEFSE_BOOTSTRAP_N,
    LEFSE_BOOTSTRAP_FRACTION,
    "CPM",
    "none (input features already species-level)",
    "BH-FDR across all prevalence-filtered species for independent Wilcoxon tests"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  parameter_manifest,
  file.path(OUTDIR, "analysis_parameters.csv"),
  row.names = FALSE
)

# ----------------------------
# 15. Console summary
# ----------------------------
cat("\n========================================\n")
cat("Species-level metagenomics analysis complete\n")
cat("========================================\n")
cat("Subjects:", nrow(meta), "\n")
cat("HCs:", sum(meta$group == "HCs"), "\n")
cat("PMI:", sum(meta$group == "PMI"), "\n")
cat("Input valid features:", sum(valid_feature), "\n")
cat("Retained after prevalence filter:", ncol(relative_abundance), "\n")
cat("LEfSe markers:", nrow(marker_df), "\n")
cat(
  "LEfSe markers also supported at Wilcoxon BH-FDR < 0.05:",
  sum(marker_df$FDR_supported, na.rm = TRUE),
  "\n"
)
cat("Output directory:", OUTDIR, "\n")
cat("No manual intermediate edits or artificial observations were used.\n")
cat("========================================\n")
