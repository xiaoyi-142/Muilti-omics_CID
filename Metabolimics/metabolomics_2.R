# ============================================================
# metabolomics_2_review_ready.R
# Global metabolomic profile analysis: PERMANOVA + PERMDISP + dbRDA
#
# Purpose:
#   - Match samples explicitly by Subid
#   - Transform non-negative metabolite abundance/intensity profiles
#   - Test overall HCs vs PMI profile difference with PERMANOVA
#   - Check dispersion with PERMDISP
#   - Visualize observed ordination coordinates with dbRDA
#   - Run exploratory univariable associations between profile structure
#     and prespecified clinical/demographic variables
#
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ----------------------------
# 0. Parameters
# ----------------------------
SEED <- 20260918L
N_PERM <- 9999L

INPUT_FILE <- "Data/Metobolimics_data2.xlsx"
META_SHEET <- 1
METAB_SHEET <- 6
MAP_SHEET <- 5

GROUP_LEVELS <- c("HCs", "PMI")
ELLIPSE_LEVEL <- 0.95
TOP_N_VECTORS <- 8L

EXPLORATORY_VARIABLES <- c(
  "Stage", "Age", "BMI", "Education",
  "TST", "SE", "TIB", "SOL", "WASO",
  "PSQI", "MoCA", "KMI", "MRS", "SAS"
)

OUTDIR <- "Results/Metabolomics"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1. Packages
# ----------------------------
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "purrr",
  "ggplot2", "ggrepel", "vegan"
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
  library(ggrepel)
  library(vegan)
})

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo_metabolomics_2.txt"))

# ----------------------------
# 2. Helpers
# ----------------------------
check_unique_ids <- function(x, label) {
  dup <- unique(x[duplicated(x)])
  if (length(dup) > 0) {
    stop(label, " contains duplicated IDs: ",
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
      stop("Non-numeric values detected in metabolite column: ", nm)
    }
    out[[nm]] <- numeric_x
  }
  out
}

build_feature_map <- function(raw_feature_names, map_df) {
  map_df <- as.data.frame(map_df)

  if (!"Name" %in% names(map_df)) {
    stop("Mapping sheet must contain a Name column.")
  }

  # Preferred: explicit matching by current column name to KEGG_ID.
  if ("KEGG_ID" %in% names(map_df) &&
      all(raw_feature_names %in% as.character(map_df$KEGG_ID))) {
    matched_map <- map_df[as.character(map_df$KEGG_ID) %in% raw_feature_names, , drop = FALSE]
    if (anyDuplicated(as.character(matched_map$KEGG_ID))) {
      stop("KEGG_ID is not unique in mapping sheet for one or more data columns.")
    }
    idx <- match(raw_feature_names, as.character(map_df$KEGG_ID))
    out <- data.frame(
      DataColumn = raw_feature_names,
      Metabolite = as.character(map_df$Name[idx]),
      KEGG_ID = as.character(map_df$KEGG_ID[idx]),
      MappingMethod = "Matched by KEGG_ID",
      stringsAsFactors = FALSE
    )
  } else if (all(raw_feature_names %in% as.character(map_df$Name))) {
    matched_map <- map_df[as.character(map_df$Name) %in% raw_feature_names, , drop = FALSE]
    if (anyDuplicated(as.character(matched_map$Name))) {
      stop("Name is not unique in mapping sheet for one or more data columns.")
    }
    idx <- match(raw_feature_names, as.character(map_df$Name))
    out <- data.frame(
      DataColumn = raw_feature_names,
      Metabolite = as.character(map_df$Name[idx]),
      KEGG_ID = if ("KEGG_ID" %in% names(map_df))
        as.character(map_df$KEGG_ID[idx]) else NA_character_,
      MappingMethod = "Matched by Name",
      stringsAsFactors = FALSE
    )
  } else if (nrow(map_df) == length(raw_feature_names)) {
    warning(
      "Feature names could not be explicitly matched to mapping IDs. ",
      "Using the workbook's documented positional mapping. ",
      "This fallback is recorded in the audit table."
    )
    out <- data.frame(
      DataColumn = raw_feature_names,
      Metabolite = as.character(map_df$Name),
      KEGG_ID = if ("KEGG_ID" %in% names(map_df))
        as.character(map_df$KEGG_ID) else NA_character_,
      MappingMethod = "Positional fallback",
      stringsAsFactors = FALSE
    )
  } else {
    stop(
      "Could not map metabolite columns safely. Add an explicit column-ID ",
      "mapping to sheet 5."
    )
  }

  out$Metabolite[is.na(out$Metabolite) | out$Metabolite == ""] <-
    out$DataColumn[is.na(out$Metabolite) | out$Metabolite == ""]
  out$AnalysisName <- make.unique(out$Metabolite, sep = "__dup")
  out
}

# ----------------------------
# 3. Read data
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
check_unique_ids(metab_raw$Subid, "Metabolomics data")

# ----------------------------
# 4. Match samples and metadata
# ----------------------------
meta_analysis <- meta %>%
  filter(Subid %in% metab_raw$Subid)

metab_analysis <- metab_raw %>%
  filter(Subid %in% meta_analysis$Subid) %>%
  arrange(match(Subid, meta_analysis$Subid))

meta_analysis <- meta_analysis %>%
  arrange(match(Subid, metab_analysis$Subid))

if (!identical(meta_analysis$Subid, metab_analysis$Subid)) {
  stop("Sample order mismatch after Subid matching.")
}

write.csv(
  meta_analysis,
  file.path(OUTDIR, "metabolomics_2_analysis_metadata.csv"),
  row.names = FALSE
)

# ----------------------------
# 5. Feature mapping and numeric QC
# ----------------------------
raw_feature_names <- setdiff(names(metab_analysis), "Subid")
feature_map <- build_feature_map(raw_feature_names, map_df)

write.csv(
  feature_map,
  file.path(OUTDIR, "metabolomics_2_feature_mapping_audit.csv"),
  row.names = FALSE
)

abundance_df <- metab_analysis %>%
  select(all_of(raw_feature_names))
abundance_df <- strict_numeric_df(abundance_df)

if (anyNA(abundance_df)) {
  stop(
    "Missing metabolite values detected. ",
    "Handle missingness with a prespecified strategy before global profile analysis."
  )
}

abundance_mat <- as.matrix(abundance_df)
colnames(abundance_mat) <- feature_map$AnalysisName
rownames(abundance_mat) <- metab_analysis$Subid
storage.mode(abundance_mat) <- "numeric"

if (any(!is.finite(abundance_mat))) stop("Inf/-Inf values detected.")
if (any(abundance_mat < 0)) {
  stop(
    "Negative metabolite values detected. Bray-Curtis requires non-negative input. ",
    "Use an appropriate preprocessing strategy before this analysis."
  )
}

# Remove zero-total features and samples explicitly.
feature_keep <- colSums(abundance_mat) > 0
sample_keep <- rowSums(abundance_mat) > 0

feature_qc <- data.frame(
  AnalysisName = colnames(abundance_mat),
  TotalAbundance = colSums(abundance_mat),
  Include = feature_keep
)
write.csv(
  feature_qc,
  file.path(OUTDIR, "metabolomics_2_feature_QC.csv"),
  row.names = FALSE
)

sample_qc <- data.frame(
  Subid = rownames(abundance_mat),
  group = as.character(meta_analysis$group),
  TotalAbundance = rowSums(abundance_mat),
  Include = sample_keep
)
write.csv(
  sample_qc,
  file.path(OUTDIR, "metabolomics_2_sample_QC.csv"),
  row.names = FALSE
)

abundance_mat <- abundance_mat[sample_keep, feature_keep, drop = FALSE]
meta_analysis <- meta_analysis[sample_keep, , drop = FALSE]

if (!identical(rownames(abundance_mat), meta_analysis$Subid)) {
  stop("Order mismatch after QC.")
}

# Prespecified square-root transform before Bray-Curtis.
abundance_sqrt <- sqrt(abundance_mat)

# ----------------------------
# 6. PERMANOVA and PERMDISP
# ----------------------------
bray <- vegan::vegdist(abundance_sqrt, method = "bray")

set.seed(SEED)
permanova <- vegan::adonis2(
  bray ~ group,
  data = meta_analysis,
  permutations = N_PERM
)

set.seed(SEED + 1L)
disp_obj <- vegan::betadisper(
  bray,
  group = meta_analysis$group,
  type = "median"
)
permdisp <- vegan::permutest(
  disp_obj,
  permutations = N_PERM
)

write.csv(
  as.data.frame(permanova),
  file.path(OUTDIR, "metabolomics_2_PERMANOVA.csv"),
  row.names = TRUE
)
write.csv(
  as.data.frame(permdisp$tab),
  file.path(OUTDIR, "metabolomics_2_PERMDISP.csv"),
  row.names = TRUE
)

# ----------------------------
# 7. dbRDA for visualization
# ----------------------------
ord <- vegan::dbrda(
  abundance_sqrt ~ group,
  data = meta_analysis,
  distance = "bray"
)

set.seed(SEED + 2L)
ord_overall_test <- anova(ord, permutations = N_PERM)

set.seed(SEED + 3L)
ord_axis_test <- anova(ord, by = "axis", permutations = N_PERM)

write.csv(
  as.data.frame(ord_overall_test),
  file.path(OUTDIR, "metabolomics_2_dbRDA_overall_test.csv"),
  row.names = TRUE
)
write.csv(
  as.data.frame(ord_axis_test),
  file.path(OUTDIR, "metabolomics_2_dbRDA_axis_test.csv"),
  row.names = TRUE
)

site_scores <- as.data.frame(
  vegan::scores(ord, display = "sites", choices = 1:2, scaling = 1)
)
if (ncol(site_scores) < 2) stop("dbRDA returned fewer than two display axes.")

axis_names <- colnames(site_scores)[1:2]
colnames(site_scores)[1:2] <- c("Axis1", "Axis2")
site_scores$Subid <- rownames(site_scores)

site_scores <- site_scores %>%
  left_join(
    meta_analysis %>% select(Subid, group),
    by = "Subid",
    relationship = "one-to-one"
  )

# Deterministically orient Axis1 so PMI is on the positive side.
if (mean(site_scores$Axis1[site_scores$group == "PMI"]) <
    mean(site_scores$Axis1[site_scores$group == "HCs"])) {
  site_scores$Axis1 <- -site_scores$Axis1
}

write.csv(
  site_scores,
  file.path(OUTDIR, "metabolomics_2_dbRDA_site_scores.csv"),
  row.names = FALSE
)

# Axis percentage labels: denominator is positive inertia only.
eig <- vegan::eigenvals(ord, model = "all")
positive_eig <- eig[eig > 0]
denom <- sum(positive_eig)

pct1 <- if (length(eig) >= 1 && denom > 0) 100 * eig[1] / denom else NA_real_
pct2 <- if (length(eig) >= 2 && denom > 0) 100 * eig[2] / denom else NA_real_

x_lab <- if (is.finite(pct1))
  sprintf("%s (%.2f%% of positive inertia)", axis_names[1], pct1)
else axis_names[1]

y_lab <- if (is.finite(pct2))
  sprintf("%s (%.2f%% of positive inertia)", axis_names[2], pct2)
else axis_names[2]

# Descriptive metabolite vectors: correlation with observed ordination scores.
abund_for_scores <- abundance_sqrt[
  match(site_scores$Subid, rownames(abundance_sqrt)),
  ,
  drop = FALSE
]

vec <- data.frame(
  AnalysisName = colnames(abund_for_scores),
  Axis1 = apply(
    abund_for_scores, 2,
    function(z) suppressWarnings(cor(z, site_scores$Axis1, method = "pearson"))
  ),
  Axis2 = apply(
    abund_for_scores, 2,
    function(z) suppressWarnings(cor(z, site_scores$Axis2, method = "pearson"))
  ),
  stringsAsFactors = FALSE
) %>%
  filter(is.finite(Axis1), is.finite(Axis2)) %>%
  mutate(VectorLength = sqrt(Axis1^2 + Axis2^2)) %>%
  arrange(desc(VectorLength)) %>%
  slice_head(n = TOP_N_VECTORS) %>%
  left_join(feature_map, by = "AnalysisName")

write.csv(
  vec,
  file.path(OUTDIR, "metabolomics_2_top_ordination_vectors.csv"),
  row.names = FALSE
)

# Scale vectors for annotation only.
if (nrow(vec) > 0) {
  xr <- diff(range(site_scores$Axis1))
  yr <- diff(range(site_scores$Axis2))
  max_vx <- max(abs(vec$Axis1))
  max_vy <- max(abs(vec$Axis2))
  sx <- if (max_vx > 0) 0.55 * xr / (2 * max_vx) else 1
  sy <- if (max_vy > 0) 0.55 * yr / (2 * max_vy) else 1
  s <- min(sx, sy)
  vec <- vec %>%
    mutate(
      Axis1_plot = Axis1 * s,
      Axis2_plot = Axis2 * s
    )
}

# Extract primary PERMANOVA values.
perm_tab <- as.data.frame(permanova)
perm_r2 <- perm_tab["group", "R2"]
perm_p <- perm_tab["group", "Pr(>F)"]

disp_tab <- as.data.frame(permdisp$tab)
disp_p_col <- grep("^Pr\\(", colnames(disp_tab), value = TRUE)[1]
disp_p <- if (length(disp_p_col) == 1) disp_tab[1, disp_p_col] else NA_real_

group_centers <- site_scores %>%
  group_by(group) %>%
  summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2),
    .groups = "drop"
  )

colors <- c("HCs" = "#4575B4", "PMI" = "#D73027")

p_ord <- ggplot(site_scores, aes(Axis1, Axis2, color = group)) +
  stat_ellipse(
    aes(fill = group),
    geom = "polygon",
    level = ELLIPSE_LEVEL,
    type = "t",
    alpha = 0.15,
    color = NA
  ) +
  geom_point(size = 1.8, alpha = 0.85) +
  geom_label(
    data = group_centers,
    aes(Axis1, Axis2, label = group),
    inherit.aes = FALSE,
    fill = "white",
    color = "black",
    fontface = "bold",
    size = 3.5
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey75") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75") +
  scale_color_manual(values = colors, name = NULL) +
  scale_fill_manual(values = colors, name = NULL) +
  labs(
    x = x_lab,
    y = y_lab,
    title = "dbRDA of KEGG-annotated metabolomic profiles",
    subtitle = sprintf(
      "PERMANOVA: R² = %.3f, P = %.4g; PERMDISP: P = %.4g",
      perm_r2, perm_p, disp_p
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

if (nrow(vec) > 0) {
  p_ord <- p_ord +
    geom_segment(
      data = vec,
      aes(
        x = 0, y = 0,
        xend = Axis1_plot, yend = Axis2_plot
      ),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.14, "cm")),
      linewidth = 0.35,
      color = "grey35",
      alpha = 0.7
    ) +
    ggrepel::geom_text_repel(
      data = vec,
      aes(Axis1_plot, Axis2_plot, label = Metabolite),
      inherit.aes = FALSE,
      size = 2.8,
      min.segment.length = 0,
      max.overlaps = Inf,
      color = "grey20"
    )
}

ggsave(
  file.path(OUTDIR, "metabolomics_2_dbRDA.pdf"),
  p_ord,
  width = 8,
  height = 6.5
)

# ----------------------------
# 8. Exploratory univariable profile associations
# ----------------------------
# Each variable is tested in a separate model. These results describe
# univariable association with the global metabolomic profile and should
# NOT be interpreted as mutually adjusted/independent contributions.
available_vars <- intersect(EXPLORATORY_VARIABLES, names(meta_analysis))

assoc_results <- map_dfr(seq_along(available_vars), function(i) {
  var_name <- available_vars[i]

  keep <- !is.na(meta_analysis[[var_name]])
  meta_sub <- meta_analysis[keep, , drop = FALSE]
  A_sub <- abundance_sqrt[keep, , drop = FALSE]

  if (nrow(meta_sub) < 5 || length(unique(meta_sub[[var_name]])) < 2) {
    return(tibble(
      Variable = var_name,
      N = nrow(meta_sub),
      R2 = NA_real_,
      R2_adj = NA_real_,
      P_value = NA_real_
    ))
  }

  pred <- meta_sub[[var_name]]
  if (is.character(pred)) pred <- factor(pred)
  meta_model <- data.frame(predictor = pred)

  ord_var <- vegan::dbrda(
    A_sub ~ predictor,
    data = meta_model,
    distance = "bray"
  )

  set.seed(SEED + 100L + i)
  ptab <- as.data.frame(anova(ord_var, permutations = N_PERM))

  r2_obj <- vegan::RsquareAdj(ord_var)

  tibble(
    Variable = var_name,
    N = nrow(meta_sub),
    R2 = unname(r2_obj$r.squared),
    R2_adj = unname(r2_obj$adj.r.squared),
    P_value = ptab[1, "Pr(>F)"]
  )
}) %>%
  mutate(
    FDR_BH = p.adjust(P_value, method = "BH")
  ) %>%
  arrange(FDR_BH, desc(R2_adj))

write.csv(
  assoc_results,
  file.path(OUTDIR, "metabolomics_2_exploratory_profile_associations.csv"),
  row.names = FALSE
)

if (nrow(assoc_results) > 0) {
  category_map <- function(x) {
    case_when(
      x %in% c("Stage", "Age", "BMI", "Education") ~ "Demographics",
      x %in% c("TST", "SE", "TIB", "SOL", "WASO") ~ "Sleep diary",
      x %in% c("PSQI", "MoCA", "KMI", "MRS", "SAS") ~ "Rating scale",
      TRUE ~ "Other"
    )
  }

  assoc_plot <- assoc_results %>%
    mutate(
      Category = category_map(Variable),
      Variable = factor(Variable, levels = rev(Variable))
    )

  p_assoc <- ggplot(
    assoc_plot,
    aes(x = R2_adj, y = Variable, fill = Category)
  ) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_col() +
    geom_text(
      aes(label = ifelse(is.na(R2_adj), "", sprintf("%.3f", R2_adj))),
      hjust = ifelse(assoc_plot$R2_adj >= 0, -0.1, 1.1),
      size = 3
    ) +
    labs(
      x = "Adjusted R² (separate univariable dbRDA models)",
      y = NULL,
      fill = NULL,
      title = "Exploratory associations with global metabolomic profile"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank())

  ggsave(
    file.path(OUTDIR, "metabolomics_2_exploratory_R2adj.pdf"),
    p_assoc,
    width = 7,
    height = 6
  )
}

cat("\nGlobal metabolomic profile analysis complete.\n")
cat(sprintf("PERMANOVA R2 = %.4f; P = %.6g\n", perm_r2, perm_p))
cat(sprintf("PERMDISP P = %.6g\n", disp_p))
