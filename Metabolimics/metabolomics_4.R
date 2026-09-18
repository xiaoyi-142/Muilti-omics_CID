# ============================================================
# metabolomics_4_review_ready.R
# Direction-specific metabolite-to-pathway mapping
#
# Purpose:
#   - Map the significant metabolites from script 3 to pathway tables
#   - Preserve many-to-many metabolite-pathway relationships
#   - Keep PMI- and HCs-direction pathway sources explicit
#   - Export unmatched cases for audit
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ----------------------------
# 0. Paths
# ----------------------------
INPUT_FILE <- "Data/Metobolimics_data2.xlsx"
PMI_PATHWAY_SHEET <- "PMI_Pathway"
HCS_PATHWAY_SHEET <- "HCs_Pathway"

RESULT_FILE <- "Results/Metabolomics/Filtered_Metabolites_Wilcoxon_FDR_FC.csv"
OUTDIR <- "Results/Metabolomics"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1. Packages
# ----------------------------
required_pkgs <- c("readxl", "dplyr", "tidyr", "stringr")
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
  library(stringr)
})

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo_metabolomics_4.txt"))

# ----------------------------
# 2. Helpers
# ----------------------------
normalize_key <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_squish() %>%
    stringr::str_to_lower()
}

pathway_to_long <- function(df, direction, source_sheet) {
  required <- c("Pathway", "Metabolimics")
  if (!all(required %in% names(df))) {
    stop(
      "Pathway sheet '", source_sheet,
      "' must contain columns Pathway and Metabolimics."
    )
  }

  df %>%
    transmute(
      Pathway = stringr::str_squish(as.character(Pathway)),
      Metabolimics = as.character(Metabolimics)
    ) %>%
    filter(
      !is.na(Pathway), Pathway != "",
      !is.na(Metabolimics), Metabolimics != ""
    ) %>%
    separate_rows(Metabolimics, sep = ";") %>%
    mutate(
      Metabolite_in_pathway_table = stringr::str_squish(Metabolimics),
      MetaboliteKey = normalize_key(Metabolite_in_pathway_table),
      Direction = direction,
      SourceSheet = source_sheet
    ) %>%
    select(
      Direction,
      SourceSheet,
      Pathway,
      Metabolite_in_pathway_table,
      MetaboliteKey
    ) %>%
    distinct()
}

# ----------------------------
# 3. Read inputs
# ----------------------------
if (!file.exists(RESULT_FILE)) {
  stop(
    "Primary differential-metabolite result not found: ", RESULT_FILE,
    "\nRun metabolomics_3_review_ready.R first."
  )
}

pmi_pathway <- read_excel(INPUT_FILE, sheet = PMI_PATHWAY_SHEET)
hcs_pathway <- read_excel(INPUT_FILE, sheet = HCS_PATHWAY_SHEET)
sig_metabolites <- read.csv(
  RESULT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_result_cols <- c("Name", "Direction")
if (!all(required_result_cols %in% names(sig_metabolites))) {
  stop(
    "Result file must contain columns: ",
    paste(required_result_cols, collapse = ", ")
  )
}

# ----------------------------
# 4. Build auditable long pathway map
# ----------------------------
pathway_map <- bind_rows(
  pathway_to_long(
    pmi_pathway,
    direction = "PMI",
    source_sheet = PMI_PATHWAY_SHEET
  ),
  pathway_to_long(
    hcs_pathway,
    direction = "HCs",
    source_sheet = HCS_PATHWAY_SHEET
  )
)

write.csv(
  pathway_map,
  file.path(OUTDIR, "metabolomics_4_direction_specific_pathway_map_long.csv"),
  row.names = FALSE
)

# ----------------------------
# 5. Join by metabolite + direction
# ----------------------------
sig_for_join <- sig_metabolites %>%
  mutate(
    ResultRowID = row_number(),
    Metabolite = stringr::str_squish(as.character(Name)),
    MetaboliteKey = normalize_key(Metabolite),
    Direction = as.character(Direction)
  )

invalid_direction <- setdiff(
  unique(na.omit(sig_for_join$Direction)),
  c("PMI", "HCs")
)
if (length(invalid_direction) > 0) {
  stop(
    "Unexpected Direction values: ",
    paste(invalid_direction, collapse = ", ")
  )
}

matches_long <- sig_for_join %>%
  left_join(
    pathway_map,
    by = c("Direction", "MetaboliteKey")
  )

write.csv(
  matches_long,
  file.path(OUTDIR, "Filtered_Metabolites_Pathway_matches_long.csv"),
  row.names = FALSE
)

# ----------------------------
# 6. Collapsed one-row-per-metabolite result
# ----------------------------
collapsed_pathways <- matches_long %>%
  group_by(ResultRowID) %>%
  summarise(
    Pathway = if (all(is.na(Pathway))) {
      NA_character_
    } else {
      paste(sort(unique(na.omit(Pathway))), collapse = "; ")
    },
    N_Pathways = n_distinct(Pathway, na.rm = TRUE),
    .groups = "drop"
  )

collapsed <- sig_for_join %>%
  select(ResultRowID, all_of(names(sig_metabolites))) %>%
  left_join(collapsed_pathways, by = "ResultRowID") %>%
  select(-ResultRowID)

write.csv(
  collapsed,
  file.path(OUTDIR, "Filtered_Metabolites_with_Pathway.csv"),
  row.names = FALSE
)

# ----------------------------
# 7. Unmatched audit
# ----------------------------
unmatched <- collapsed %>%
  filter(is.na(Pathway)) %>%
  select(
    any_of(c(
      "AnalysisName", "Name", "KEGG_ID",
      "KEGG_MapID", "Direction"
    ))
  )

# Distinguish "not found in direction-specific sheet" from "not found anywhere".
union_keys <- pathway_map %>%
  distinct(MetaboliteKey)

unmatched_audit <- sig_for_join %>%
  filter(!MetaboliteKey %in%
           (matches_long %>%
              filter(!is.na(Pathway)) %>%
              pull(MetaboliteKey))) %>%
  mutate(
    FoundInEitherPathwaySheet = MetaboliteKey %in% union_keys$MetaboliteKey,
    AuditInterpretation = ifelse(
      FoundInEitherPathwaySheet,
      "Metabolite exists in a pathway sheet, but not in the sheet matching its differential direction",
      "Metabolite not found in either pathway sheet"
    )
  ) %>%
  select(
    any_of(c(
      "AnalysisName", "Name", "KEGG_ID",
      "Direction",
      "FoundInEitherPathwaySheet",
      "AuditInterpretation"
    ))
  )

write.csv(
  unmatched_audit,
  file.path(OUTDIR, "metabolomics_4_unmatched_metabolites_audit.csv"),
  row.names = FALSE
)

# ----------------------------
# 8. Unique pathway summary
# ----------------------------
unique_pathways <- matches_long %>%
  filter(!is.na(Pathway)) %>%
  distinct(Direction, Pathway) %>%
  arrange(Direction, Pathway)

write.csv(
  unique_pathways,
  file.path(OUTDIR, "metabolomics_4_unique_pathways.csv"),
  row.names = FALSE
)

summary_df <- data.frame(
  Metric = c(
    "N_significant_metabolites",
    "N_with_at_least_one_direction_matched_pathway",
    "N_unmatched",
    "N_unique_PMI_pathways",
    "N_unique_HCs_pathways"
  ),
  Value = c(
    nrow(sig_metabolites),
    sum(!is.na(collapsed$Pathway)),
    sum(is.na(collapsed$Pathway)),
    sum(unique_pathways$Direction == "PMI"),
    sum(unique_pathways$Direction == "HCs")
  )
)

write.csv(
  summary_df,
  file.path(OUTDIR, "metabolomics_4_pathway_mapping_summary.csv"),
  row.names = FALSE
)

cat("\nPathway mapping complete.\n")
cat("Significant metabolites:", nrow(sig_metabolites), "\n")
cat("Mapped:", sum(!is.na(collapsed$Pathway)), "\n")
cat("Unmatched:", sum(is.na(collapsed$Pathway)), "\n")
