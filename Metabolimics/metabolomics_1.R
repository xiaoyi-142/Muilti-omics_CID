# ============================================================
# metabolomics_1_review_ready.R
# Cross-sectional univariate metabolite comparison: HCs vs PMI
#
# Purpose:
#   - Match samples by Subid
#   - Test every testable metabolite once using Wilcoxon rank-sum
#   - Apply BH-FDR across the full tested metabolite set
#   - Export complete QC and result tables
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ----------------------------
# 0. Parameters
# ----------------------------
INPUT_FILE <- "Data/Metobolimics_data2.xlsx"
META_SHEET <- 1
METAB_SHEET <- 7

GROUP_LEVELS <- c("HCs", "PMI")
FDR_CUTOFF <- 0.05
MIN_NONMISSING_PER_GROUP <- 3L

OUTDIR <- "Results/Metabolomics"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1. Packages
# ----------------------------
required_pkgs <- c("readxl", "dplyr", "tidyr", "purrr")
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
})

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo_metabolomics_1.txt"))

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

strict_numeric <- function(x, feature_name) {
  original <- x
  numeric_x <- suppressWarnings(as.numeric(as.character(x)))
  bad <- !is.na(original) & is.na(numeric_x)
  if (any(bad)) {
    stop(
      "Non-numeric values detected in metabolite '", feature_name,
      "'. Resolve them explicitly before analysis."
    )
  }
  numeric_x
}

safe_wilcox <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- droplevels(g[keep])

  if (nlevels(g) != 2) return(c(W = NA_real_, P = NA_real_))

  n1 <- sum(g == GROUP_LEVELS[1])
  n2 <- sum(g == GROUP_LEVELS[2])
  if (n1 < MIN_NONMISSING_PER_GROUP || n2 < MIN_NONMISSING_PER_GROUP) {
    return(c(W = NA_real_, P = NA_real_))
  }

  if (length(unique(x)) <= 1) {
    return(c(W = NA_real_, P = 1))
  }

  tst <- suppressWarnings(
    wilcox.test(
      x ~ g,
      exact = FALSE,
      correct = FALSE,
      conf.int = FALSE
    )
  )
  c(W = unname(tst$statistic), P = tst$p.value)
}

# ----------------------------
# 3. Read and validate metadata
# ----------------------------
meta <- read_excel(INPUT_FILE, sheet = META_SHEET)
metab <- read_excel(INPUT_FILE, sheet = METAB_SHEET)

if (!all(c("Subid", "group") %in% names(meta))) {
  stop("Metadata sheet must contain Subid and group.")
}
if (!"Subid" %in% names(metab)) {
  stop("Metabolomics sheet must contain Subid.")
}

meta <- meta %>%
  mutate(
    Subid = trimws(as.character(Subid)),
    group = trimws(as.character(group))
  ) %>%
  filter(group %in% GROUP_LEVELS) %>%
  mutate(group = factor(group, levels = GROUP_LEVELS))

metab$Subid <- trimws(as.character(metab$Subid))

check_unique_ids(meta$Subid, "Metadata")
check_unique_ids(metab$Subid, "Metabolomics data")

# ----------------------------
# 4. Explicit sample matching
# ----------------------------
sample_audit <- meta %>%
  transmute(
    Subid,
    group = as.character(group),
    HasMetabolomics = Subid %in% metab$Subid
  )

write.csv(
  sample_audit,
  file.path(OUTDIR, "metabolomics_1_sample_match_audit.csv"),
  row.names = FALSE
)

meta_analysis <- meta %>%
  filter(Subid %in% metab$Subid)

metab_analysis <- metab %>%
  filter(Subid %in% meta_analysis$Subid) %>%
  arrange(match(Subid, meta_analysis$Subid))

meta_analysis <- meta_analysis %>%
  arrange(match(Subid, metab_analysis$Subid))

if (!identical(metab_analysis$Subid, meta_analysis$Subid)) {
  stop("Sample order mismatch after Subid matching.")
}

write.csv(
  meta_analysis %>% select(Subid, group),
  file.path(OUTDIR, "metabolomics_1_analysis_subject_order.csv"),
  row.names = FALSE
)

# ----------------------------
# 5. Strict numeric conversion
# ----------------------------
feature_cols <- setdiff(names(metab_analysis), "Subid")
if (length(feature_cols) == 0) stop("No metabolite columns found.")

for (nm in feature_cols) {
  metab_analysis[[nm]] <- strict_numeric(metab_analysis[[nm]], nm)
}

# ----------------------------
# 6. Long table and one test per metabolite
# ----------------------------
long_df <- metab_analysis %>%
  left_join(
    meta_analysis %>% select(Subid, group),
    by = "Subid",
    relationship = "one-to-one"
  ) %>%
  pivot_longer(
    cols = all_of(feature_cols),
    names_to = "Metabolite",
    values_to = "Value"
  )

results <- long_df %>%
  group_by(Metabolite) %>%
  group_modify(~{
    dat <- .x
    w <- safe_wilcox(dat$Value, dat$group)

    x_hc <- dat$Value[dat$group == "HCs"]
    x_pmi <- dat$Value[dat$group == "PMI"]

    tibble(
      N_HCs = sum(!is.na(x_hc)),
      N_PMI = sum(!is.na(x_pmi)),
      Median_HCs = median(x_hc, na.rm = TRUE),
      IQR_HCs = IQR(x_hc, na.rm = TRUE),
      Median_PMI = median(x_pmi, na.rm = TRUE),
      IQR_PMI = IQR(x_pmi, na.rm = TRUE),
      Mean_HCs = mean(x_hc, na.rm = TRUE),
      Mean_PMI = mean(x_pmi, na.rm = TRUE),
      W = w["W"],
      P_value = w["P"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    FDR_BH = p.adjust(P_value, method = "BH"),
    FDR_significant = !is.na(FDR_BH) & FDR_BH < FDR_CUTOFF
  ) %>%
  arrange(FDR_BH, P_value)

# Feature audit
feature_audit <- results %>%
  transmute(
    Metabolite,
    N_HCs,
    N_PMI,
    Testable = !is.na(P_value),
    Reason = case_when(
      !Testable ~ paste0(
        "Insufficient non-missing data (<", MIN_NONMISSING_PER_GROUP,
        " per group) or invalid two-group comparison"
      ),
      TRUE ~ ""
    )
  )

write.csv(
  feature_audit,
  file.path(OUTDIR, "metabolomics_1_feature_testability_audit.csv"),
  row.names = FALSE
)

write.csv(
  results,
  file.path(OUTDIR, "metabolomics_1_all_metabolites_Wilcoxon_BH.csv"),
  row.names = FALSE
)

write.csv(
  results %>% filter(FDR_significant),
  file.path(OUTDIR, "metabolomics_1_FDR_significant_metabolites.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  Metric = c(
    "N_total_samples",
    "N_HCs",
    "N_PMI",
    "N_metabolites_input",
    "N_metabolites_tested",
    "N_FDR_significant",
    "FDR_cutoff"
  ),
  Value = c(
    nrow(meta_analysis),
    sum(meta_analysis$group == "HCs"),
    sum(meta_analysis$group == "PMI"),
    length(feature_cols),
    sum(!is.na(results$P_value)),
    sum(results$FDR_significant),
    FDR_CUTOFF
  )
)

write.csv(
  summary_table,
  file.path(OUTDIR, "metabolomics_1_analysis_summary.csv"),
  row.names = FALSE
)

cat("\nMetabolomics univariate analysis complete.\n")
cat("HCs:", sum(meta_analysis$group == "HCs"),
    " PMI:", sum(meta_analysis$group == "PMI"), "\n")
cat("Tested metabolites:", sum(!is.na(results$P_value)), "\n")
cat("BH-FDR significant:", sum(results$FDR_significant), "\n")
