# ============================================================
# metabolomics_3_review_ready.R
# Differential analysis of KEGG-annotated metabolites: HCs vs PMI
#
# Primary rule (prespecified):
#   BH-FDR < 0.05 AND absolute log2 fold change >= 1
#   (equivalent to FC >= 2 or <= 0.5; PMI / HCs)
#
# Purpose:
#   - Explicit Subid matching
#   - Auditable feature mapping
#   - BH correction across ALL eligible KEGG-annotated metabolites
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ----------------------------
# 0. Parameters
# ----------------------------
INPUT_FILE <- "Data/Metobolimics_data2.xlsx"
META_SHEET <- 1
METAB_SHEET <- 6
MAP_SHEET <- 5

GROUP_LEVELS <- c("HCs", "PMI")
FDR_CUTOFF <- 0.05
ABS_LOG2FC_CUTOFF <- 1.0
MIN_NONMISSING_PER_GROUP <- 3L

OUTDIR <- "Results/Metabolomics"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1. Packages
# ----------------------------
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "purrr",
  "ggplot2", "scales"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(scales)
})

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo_metabolomics_3.txt"))

# ----------------------------
# 2. Helpers
# ----------------------------
check_unique_ids <- function(x, label) {
  dup <- unique(x[duplicated(x)])
  if (length(dup) > 0) {
    stop(label, " contains duplicated Subid values: ",
         paste(head(dup, 20), collapse = ", "))
  }
}

strict_numeric_df <- function(df) {
  out <- df
  for (nm in names(out)) {
    original <- out[[nm]]
    numeric_x <- suppressWarnings(as.numeric(as.character(original)))
    bad <- !is.na(original) & is.na(numeric_x)
    if (any(bad)) {
      stop("Non-numeric values detected in feature column: ", nm)
    }
    out[[nm]] <- numeric_x
  }
  out
}

build_feature_map <- function(raw_feature_names, map_df) {
  map_df <- as.data.frame(map_df)

  required <- c("Name")
  if (!all(required %in% names(map_df))) {
    stop("Mapping sheet must contain Name.")
  }

  if ("KEGG_ID" %in% names(map_df) &&
      all(raw_feature_names %in% as.character(map_df$KEGG_ID))) {
    matched_map <- map_df[as.character(map_df$KEGG_ID) %in% raw_feature_names, , drop = FALSE]
    if (anyDuplicated(as.character(matched_map$KEGG_ID))) {
      stop("KEGG_ID is not unique in mapping sheet for one or more data columns.")
    }
    idx <- match(raw_feature_names, as.character(map_df$KEGG_ID))
    out <- map_df[idx, , drop = FALSE]
    out$DataColumn <- raw_feature_names
    out$MappingMethod <- "Matched by KEGG_ID"
  } else if (all(raw_feature_names %in% as.character(map_df$Name))) {
    matched_map <- map_df[as.character(map_df$Name) %in% raw_feature_names, , drop = FALSE]
    if (anyDuplicated(as.character(matched_map$Name))) {
      stop("Name is not unique in mapping sheet for one or more data columns.")
    }
    idx <- match(raw_feature_names, as.character(map_df$Name))
    out <- map_df[idx, , drop = FALSE]
    out$DataColumn <- raw_feature_names
    out$MappingMethod <- "Matched by Name"
  } else if (nrow(map_df) == length(raw_feature_names)) {
    warning(
      "Using documented positional mapping because explicit IDs did not match. ",
      "The mapping is exported for audit."
    )
    out <- map_df
    out$DataColumn <- raw_feature_names
    out$MappingMethod <- "Positional fallback"
  } else {
    stop(
      "Unable to map sheet-6 features safely to sheet-5 annotation. ",
      "Provide an explicit DataColumn-to-metabolite mapping."
    )
  }

  out$Name <- as.character(out$Name)
  out$Name[is.na(out$Name) | out$Name == ""] <-
    out$DataColumn[is.na(out$Name) | out$Name == ""]
  out$AnalysisName <- make.unique(out$Name, sep = "__dup")
  out
}

safe_wilcox <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- droplevels(g[keep])

  if (nlevels(g) != 2) return(c(W = NA_real_, P = NA_real_))
  if (sum(g == "HCs") < MIN_NONMISSING_PER_GROUP ||
      sum(g == "PMI") < MIN_NONMISSING_PER_GROUP) {
    return(c(W = NA_real_, P = NA_real_))
  }
  if (length(unique(x)) <= 1) return(c(W = NA_real_, P = 1))

  tst <- suppressWarnings(
    wilcox.test(x ~ g, exact = FALSE, correct = FALSE)
  )
  c(W = unname(tst$statistic), P = tst$p.value)
}

safe_fold_change <- function(mean_hc, mean_pmi) {
  # Input data are required to be non-negative.
  if (is.na(mean_hc) || is.na(mean_pmi)) return(NA_real_)
  if (mean_hc > 0) return(mean_pmi / mean_hc)
  if (mean_hc == 0 && mean_pmi > 0) return(Inf)
  if (mean_hc == 0 && mean_pmi == 0) return(NA_real_)
  NA_real_
}

# ----------------------------
# 3. Read data and map features
# ----------------------------
meta <- read_excel(INPUT_FILE, sheet = META_SHEET)
metab_raw <- read_excel(INPUT_FILE, sheet = METAB_SHEET)
map_df <- read_excel(INPUT_FILE, sheet = MAP_SHEET)

if (!all(c("Subid", "group") %in% names(meta))) {
  stop("Metadata must contain Subid and group.")
}
if (!"Subid" %in% names(metab_raw)) {
  stop("Metabolomics sheet must contain Subid.")
}

meta <- meta %>%
  mutate(
    Subid = trimws(as.character(Subid)),
    group = trimws(as.character(group))
  ) %>%
  filter(group %in% GROUP_LEVELS) %>%
  mutate(group = factor(group, levels = GROUP_LEVELS))

metab_raw$Subid <- trimws(as.character(metab_raw$Subid))

check_unique_ids(meta$Subid, "Metadata")
check_unique_ids(metab_raw$Subid, "Metabolomics")

raw_feature_names <- setdiff(names(metab_raw), "Subid")
feature_map <- build_feature_map(raw_feature_names, map_df)

# KEGG annotation eligibility.
if (!"KEGG_MapID" %in% names(feature_map)) {
  stop("Mapping sheet must contain KEGG_MapID for this analysis.")
}

feature_map <- feature_map %>%
  mutate(
    KEGG_MapID = as.character(KEGG_MapID),
    KEGG_eligible = !is.na(KEGG_MapID) &
      trimws(KEGG_MapID) != "" &
      trimws(KEGG_MapID) != "-"
  )

write.csv(
  feature_map,
  file.path(OUTDIR, "metabolomics_3_feature_mapping_audit.csv"),
  row.names = FALSE
)

eligible_map <- feature_map %>% filter(KEGG_eligible)
if (nrow(eligible_map) == 0) stop("No KEGG-annotated metabolites available.")

# ----------------------------
# 4. Match samples by Subid
# ----------------------------
meta_analysis <- meta %>%
  filter(Subid %in% metab_raw$Subid)

metab_analysis <- metab_raw %>%
  filter(Subid %in% meta_analysis$Subid) %>%
  arrange(match(Subid, meta_analysis$Subid))

meta_analysis <- meta_analysis %>%
  arrange(match(Subid, metab_analysis$Subid))

if (!identical(meta_analysis$Subid, metab_analysis$Subid)) {
  stop("Sample order mismatch.")
}

write.csv(
  meta_analysis %>% select(Subid, group),
  file.path(OUTDIR, "metabolomics_3_analysis_subject_order.csv"),
  row.names = FALSE
)

# ----------------------------
# 5. Build unique analysis matrix
# ----------------------------
abundance_df <- metab_analysis %>%
  select(all_of(eligible_map$DataColumn))

abundance_df <- strict_numeric_df(abundance_df)
colnames(abundance_df) <- eligible_map$AnalysisName

if (any(as.matrix(abundance_df) < 0, na.rm = TRUE)) {
  stop(
    "Negative values detected. Mean-ratio fold change used here requires ",
    "non-negative abundance/intensity values."
  )
}

# ----------------------------
# 6. Differential testing
# ----------------------------
long_df <- abundance_df %>%
  mutate(
    Subid = metab_analysis$Subid,
    group = meta_analysis$group
  ) %>%
  pivot_longer(
    cols = all_of(eligible_map$AnalysisName),
    names_to = "AnalysisName",
    values_to = "Value"
  ) %>%
  left_join(
    eligible_map %>%
      select(
        AnalysisName, DataColumn, Name,
        any_of(c("KEGG_ID", "KEGG_MapID"))
      ),
    by = "AnalysisName"
  )

group_cols <- c(
  "AnalysisName", "DataColumn", "Name",
  intersect(c("KEGG_ID", "KEGG_MapID"), names(long_df))
)

results_all <- long_df %>%
  group_by(across(all_of(group_cols))) %>%
  group_modify(~{
    dat <- .x
    w <- safe_wilcox(dat$Value, dat$group)

    hc <- dat$Value[dat$group == "HCs"]
    pmi <- dat$Value[dat$group == "PMI"]

    mean_hc <- mean(hc, na.rm = TRUE)
    mean_pmi <- mean(pmi, na.rm = TRUE)
    fc <- safe_fold_change(mean_hc, mean_pmi)

    tibble(
      N_HCs = sum(!is.na(hc)),
      N_PMI = sum(!is.na(pmi)),
      Median_HCs = median(hc, na.rm = TRUE),
      IQR_HCs = IQR(hc, na.rm = TRUE),
      Median_PMI = median(pmi, na.rm = TRUE),
      IQR_PMI = IQR(pmi, na.rm = TRUE),
      Mean_HCs = mean_hc,
      Mean_PMI = mean_pmi,
      FoldChange_PMI_vs_HCs = fc,
      log2FC_PMI_vs_HCs = ifelse(
        is.na(fc),
        NA_real_,
        ifelse(is.infinite(fc), Inf, log2(fc))
      ),
      W = w["W"],
      P_value = w["P"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    FDR_BH = p.adjust(P_value, method = "BH"),
    FDR_significant = !is.na(FDR_BH) & FDR_BH < FDR_CUTOFF,
    Effect_size_pass = !is.na(log2FC_PMI_vs_HCs) &
      abs(log2FC_PMI_vs_HCs) >= ABS_LOG2FC_CUTOFF,
    Primary_marker = FDR_significant & Effect_size_pass,
    Direction = case_when(
      Primary_marker & log2FC_PMI_vs_HCs > 0 ~ "PMI",
      Primary_marker & log2FC_PMI_vs_HCs < 0 ~ "HCs",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(FDR_BH, desc(abs(log2FC_PMI_vs_HCs)))

results_primary <- results_all %>%
  filter(Primary_marker)

write.csv(
  results_all,
  file.path(OUTDIR, "metabolomics_3_all_KEGG_metabolites_Wilcoxon_FDR_FC.csv"),
  row.names = FALSE
)

write.csv(
  results_primary,
  file.path(OUTDIR, "Filtered_Metabolites_Wilcoxon_FDR_FC.csv"),
  row.names = FALSE
)

# ----------------------------
# 7. Export abundance for downstream analyses
# ----------------------------
# All KEGG-annotated metabolites.
abundance_export_all <- abundance_df
colnames(abundance_export_all) <- eligible_map$AnalysisName

abundance_export_all <- data.frame(
  Subid = metab_analysis$Subid,
  abundance_export_all,
  check.names = FALSE
)

write.csv(
  abundance_export_all,
  file.path(OUTDIR, "metabolomics_KEGG_abundance_all.csv"),
  row.names = FALSE
)

# Primary marker subset.
primary_names <- results_primary$AnalysisName
abundance_export_sig <- data.frame(
  Subid = metab_analysis$Subid,
  abundance_df[, primary_names, drop = FALSE],
  check.names = FALSE
)

write.csv(
  abundance_export_sig,
  file.path(OUTDIR, "Significant_Metabolites_Abundance.csv"),
  row.names = FALSE
)

# If KEGG IDs are available, export a KEGG-ID-keyed matrix without losing
# duplicate IDs: make.unique() preserves every original feature.
if ("KEGG_ID" %in% names(results_primary) && nrow(results_primary) > 0) {
  kegg_ids <- as.character(results_primary$KEGG_ID)
  kegg_ids[is.na(kegg_ids) | kegg_ids == ""] <-
    results_primary$AnalysisName[is.na(kegg_ids) | kegg_ids == ""]
  kegg_ids_unique <- make.unique(kegg_ids, sep = "__dup")

  kegg_mat <- abundance_df[, results_primary$AnalysisName, drop = FALSE]
  colnames(kegg_mat) <- kegg_ids_unique

  write.csv(
    data.frame(
      Subid = metab_analysis$Subid,
      kegg_mat,
      check.names = FALSE
    ),
    file.path(OUTDIR, "sig_metabolite_abundance_by_KEGGID.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      AnalysisName = results_primary$AnalysisName,
      Metabolite = results_primary$Name,
      KEGG_ID = results_primary$KEGG_ID,
      ExportColumn = kegg_ids_unique,
      stringsAsFactors = FALSE
    ),
    file.path(OUTDIR, "sig_KEGGID_export_mapping.csv"),
    row.names = FALSE
  )
}

# ----------------------------
# 8. Scientific visualization
# ----------------------------
# Use observed raw values and a pseudo-log display transform so zeros are
# not silently dropped. Statistical tests above use the original values.
if (nrow(results_primary) > 0) {

  plot_df <- long_df %>%
    filter(AnalysisName %in% results_primary$AnalysisName) %>%
    left_join(
      results_primary %>%
        select(
          AnalysisName,
          FDR_BH,
          log2FC_PMI_vs_HCs,
          Direction
        ),
      by = "AnalysisName"
    ) %>%
    mutate(
      PlotLabel = paste0(
        Name,
        "\nq=", formatC(FDR_BH, format = "g", digits = 2),
        "; log2FC=", formatC(log2FC_PMI_vs_HCs, format = "f", digits = 2)
      )
    )

  positive_values <- plot_df$Value[
    is.finite(plot_df$Value) & plot_df$Value > 0
  ]

  pseudo_sigma <- if (length(positive_values) > 0) {
    min(positive_values) / 2
  } else {
    1
  }

  colors <- c("HCs" = "#4575B4", "PMI" = "#D73027")

  make_direction_plot <- function(direction_value, filename) {
    dd <- plot_df %>% filter(Direction == direction_value)
    if (nrow(dd) == 0) return(invisible(NULL))

    p <- ggplot(dd, aes(x = group, y = Value, fill = group)) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.55,
        linewidth = 0.45
      ) +
      geom_jitter(
        aes(color = group),
        width = 0.12,
        size = 0.8,
        alpha = 0.55,
        show.legend = FALSE
      ) +
      facet_wrap(~PlotLabel, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_y_continuous(
        trans = scales::pseudo_log_trans(
          base = 10,
          sigma = pseudo_sigma
        )
      ) +
      labs(
        x = NULL,
        y = "Metabolite abundance/intensity",
        fill = "Group",
        title = ifelse(
          direction_value == "PMI",
          "Metabolites higher in PMI",
          "Metabolites higher in HCs"
        ),
        subtitle = sprintf(
          "Primary rule: BH-FDR < %.2f and |log2FC| >= %.1f",
          FDR_CUTOFF, ABS_LOG2FC_CUTOFF
        )
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        legend.position = "top",
        strip.text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )

    ggsave(
      file.path(OUTDIR, filename),
      p,
      width = 10,
      height = max(5, ceiling(length(unique(dd$AnalysisName)) / 4) * 2.5)
    )
  }

  make_direction_plot("PMI", "metabolomics_3_higher_in_PMI.pdf")
  make_direction_plot("HCs", "metabolomics_3_higher_in_HCs.pdf")
}

# ----------------------------
# 9. Summary
# ----------------------------
summary_df <- data.frame(
  Metric = c(
    "N_HCs",
    "N_PMI",
    "N_KEGG_eligible_metabolites",
    "N_tested",
    "N_FDR_significant",
    "N_primary_FDR_plus_FC",
    "FDR_cutoff",
    "abs_log2FC_cutoff"
  ),
  Value = c(
    sum(meta_analysis$group == "HCs"),
    sum(meta_analysis$group == "PMI"),
    nrow(eligible_map),
    sum(!is.na(results_all$P_value)),
    sum(results_all$FDR_significant),
    sum(results_all$Primary_marker),
    FDR_CUTOFF,
    ABS_LOG2FC_CUTOFF
  )
)

write.csv(
  summary_df,
  file.path(OUTDIR, "metabolomics_3_analysis_summary.csv"),
  row.names = FALSE
)

cat("\nKEGG-annotated metabolite differential analysis complete.\n")
cat("Primary markers:", nrow(results_primary), "\n")
