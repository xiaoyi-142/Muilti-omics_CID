# # 进行Kegg ID分析
rm(list = ls()) #2.清理运行环境*`
library(tableone) #1.加载R包*`
library(xlsx)
library(readxl)
library(tidyverse)
library(rstatix)
library(dplyr)
library(tidyr)
library(broom)
library(readxl)
# 加载所需包
library(RLdbRDA)
library(vegan)
library(RColorBrewer)
library(ggrepel)

# 读取宏基因数据
metagenomics_abundance_data <- read.csv("Results/metagenomics_abundance_data.csv")

# 读取代谢组学数据
metabolimics_abundance_data <- read.csv("Results/metabolimics_abundance_data.csv")


# 求共有的 Subid
common_subid <- intersect(metagenomics_abundance_data$Subid, metabolimics_abundance_data$Subid)

# 筛选出共有Subid的子集
metagenomics_common <- metagenomics_abundance_data %>% filter(Subid %in% common_subid)
metabolimics_common <- metabolimics_abundance_data %>% filter(Subid %in% common_subid)

# 可选：保证Subid顺序一致
metagenomics_common <- metagenomics_common %>% arrange(Subid)
metabolimics_common <- metabolimics_common %>% arrange(Subid)

# 检查一下
dim(metagenomics_common)
dim(metabolimics_common)
head(metagenomics_common[, 1:5])
head(metabolimics_common[, 1:5])


library(vegan)
library(ggplot2)
library(dplyr)

# 假设 metagenomics_common, metabolimics_common 已经配对并排序一致
otu_data <- metagenomics_common %>% select(-Subid, -group)
metab_data <- metabolimics_common %>% select(-Subid)
group <- metagenomics_common$group  # or metabolimics_common$group

# PCA
otu_pca <- prcomp(otu_data, scale. = TRUE)
metab_pca <- prcomp(metab_data, scale. = TRUE)

otu_scores <- scores(otu_pca, choices = 1:900)
metab_scores <- scores(metab_pca, choices = 1:900)

# Procrustes分析和Permutation检验
set.seed(123)
proc <- procrustes(otu_scores, metab_scores, symmetric = TRUE)
prot <- protest(otu_scores, metab_scores, permutations = 999)

M2 <- round(proc$ss, 4)
r_approx <- round(sqrt(1 - M2), 3)
pval <- signif(prot$signif, 3)
title_txt <- sprintf("Procrustes: M² = %.3f, r ≈ %.3f, p = %s", M2, r_approx, pval)


n <- nrow(otu_scores)
df <- data.frame(
  Subid = metagenomics_common$Subid,
  group = group,
  x_micro = proc$Yrot[,1],   # 微生物PCA对齐后的x
  y_micro = proc$Yrot[,2],   # 微生物PCA对齐后的y
  x_metab = proc$X[,1],      # 代谢物PCA对齐后的x
  y_metab = proc$X[,2]       # 代谢物PCA对齐后的y
)


# 分组配色
group_colors <- c("NA" = "#e76f51", "NNA" = "#457b9d")
# 自动适配实际分组名
names(group_colors) <- unique(df$group)

# 画图
ggplot(df) +
  # 连线
  geom_segment(aes(x = x_micro, y = y_micro, xend = x_metab, yend = y_metab), color = "grey70") +
  # 微生物点
  geom_point(aes(x = x_micro, y = y_micro, color = group), shape = 16, size = 2) +
  # 代谢物点
  geom_point(aes(x = x_metab, y = y_metab, color = group), shape = 15, size = 2) +
  scale_color_manual(values = group_colors, name = "") +
  labs(
    x = NULL, y = NULL,
    title = title_txt
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0, vjust = 1, size = 16),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  # 中心虚线
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  # 限制纵坐标范围
  ylim(c(-0.05, 0.04))


ggsave("Metagenomics/Procrustes_plot.pdf", width = 10, height = 8)


