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


# 如未安装请先运行：install.packages("ade4")
library(ade4)
library(ggplot2)
library(dplyr)

# 数据已按Subid配对，顺序一致
otu_data <- metagenomics_common %>% select(-Subid, -group)
metab_data <- metabolimics_common %>% select(-Subid)
group <- metagenomics_common$group


otu_pca <- dudi.pca(otu_data, scannf = FALSE, nf = 2)
metab_pca <- dudi.pca(metab_data, scannf = FALSE, nf = 2)

# 进行Co-inertia分析
coin_res <- coinertia(otu_pca, metab_pca, scannf = FALSE, nf = 2)

# RV系数（全局相关系数）
RV <- round(coin_res$RV, 4)
cat("RV coefficient: ", RV, "\n")

# # Permutation检验（默认999次）
set.seed(123)
coin_perm <- randtest(coin_res, nrepet = 50)
pval <- signif(coin_perm$pvalue, 3)
cat("Permutation p-value: ", pval, "\n")

df_cia <- data.frame(
  Subid = metagenomics_common$Subid,
  group = group,
  x_micro = coin_res$lX[,1],   # 微生物组（样本，coinertia1）
  y_micro = coin_res$lX[,2],
  x_metab = coin_res$lY[,1],   # 代谢组（样本，coinertia1）
  y_metab = coin_res$lY[,2]
)

# 分组配色（如有NA/NNA请调整）
group_colors <- c("NA" = "#e76f51", "NNA" = "#457b9d")
names(group_colors) <- unique(df_cia$group)

# 绘图
# 1. 解释度
otu_var_exp <- otu_pca$eig / sum(otu_pca$eig)
metab_var_exp <- metab_pca$eig / sum(metab_pca$eig)
otu_var_exp_1 <- round(otu_var_exp[1] * 100, 1)
otu_var_exp_2 <- round(otu_var_exp[2] * 100, 1)
metab_var_exp_1 <- round(metab_var_exp[1] * 100, 1)
metab_var_exp_2 <- round(metab_var_exp[2] * 100, 1)

# 横纵轴标签（任选其一，建议第二种）
# xlab <- sprintf("CIA1 (Microbiome PC1: %.1f%%, Metabolome PC1: %.1f%%)", otu_var_exp_1, metab_var_exp_1)
# ylab <- sprintf("CIA2 (Microbiome PC2: %.1f%%, Metabolome PC2: %.1f%%)", otu_var_exp_2, metab_var_exp_2)

xlab <- sprintf("CIA1 (mean explained variance: %.1f%%)", (otu_var_exp_1 + metab_var_exp_1) / 2)
ylab <- sprintf("CIA2 (mean explained variance: %.1f%%)", (otu_var_exp_2 + metab_var_exp_2) / 2)



# 2. 作图
# 计算四个变量的2.5%和97.5%分位点
filter_by_quantile <- function(vec) {
  q <- quantile(vec, probs = c(0.025, 0.975), na.rm = TRUE)
  (vec >= q[1]) & (vec <= q[2])
}

# 同时在四个坐标都在95%区间内的样本保留
keep_idx <- 
  filter_by_quantile(df_cia$x_micro) &
  filter_by_quantile(df_cia$y_micro) &
  filter_by_quantile(df_cia$x_metab) &
  filter_by_quantile(df_cia$y_metab)

df_cia_filt <- df_cia[keep_idx, ]
cat("保留样本数：", nrow(df_cia_filt), "\n")

# 作图
ggplot(df_cia_filt) +
  geom_segment(aes(x = x_micro, y = y_micro, xend = x_metab, yend = y_metab), color = "grey70") +
  geom_point(aes(x = x_micro, y = y_micro, color = group), shape = 16, size = 2) +
  geom_point(aes(x = x_metab, y = y_metab, color = group), shape = 15, size = 2) +
  scale_color_manual(values = group_colors, name = "") +
  labs(
    x = xlab,
    y = ylab,
    title = paste0("Co-inertia analysis: RV = ", RV, ", p-value = ", pval)
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0, vjust = 1, size = 16),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")

ggsave("Metagenomics/Coinertia_plot_filtered.pdf", width = 8, height = 6)




