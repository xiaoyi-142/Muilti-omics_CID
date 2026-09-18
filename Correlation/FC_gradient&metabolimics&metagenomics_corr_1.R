rm(list = ls()) # 清理运行环境

# 加载所需包
library(tidyverse)
library(vegan)
library(ggplot2)

# 读取数据
metagenomics_abundance_data <- read.csv("Results/metagenomics_abundance_data.csv")
metabolimics_abundance_data <- read.csv("Results/metabolimics_abundance_data.csv")
FC_Gradient_data <- read.csv("Data/FC_gradient_all_with_Subid.csv")

# 去掉group列
FC_Gradient_data <- FC_Gradient_data %>% select(-group)

# 求共有的 Subid，三个数据集都有
common_subid <- intersect(metagenomics_abundance_data$Subid, metabolimics_abundance_data$Subid)
common_subid <- intersect(common_subid, FC_Gradient_data$Subid)

# 筛选出共有Subid的子集
metagenomics_common <- metagenomics_abundance_data %>% filter(Subid %in% common_subid)
metabolimics_common <- metabolimics_abundance_data %>% filter(Subid %in% common_subid)
FC_Gradient_common <- FC_Gradient_data %>% filter(Subid %in% common_subid)

# 确保Subid顺序一致
metagenomics_common <- metagenomics_common %>% arrange(Subid)
FC_Gradient_common <- FC_Gradient_common %>% arrange(Subid)
metabolimics_common <- metabolimics_common %>% arrange(Subid)

# 1. PCA分析：微生物数据和FC数据之间
otu_data <- metagenomics_common %>% select(-Subid) %>% select_if(is.numeric)
FC_data <- FC_Gradient_common %>% select(-Subid) %>% select_if(is.numeric)

# PCA分析
otu_pca <- prcomp(otu_data, scale. = TRUE)
FC_pca <- prcomp(FC_data, scale. = TRUE)

# 提取PCA得分
otu_scores <- otu_pca$x
FC_scores <- FC_pca$x

# Procrustes分析：对齐微生物数据和FC数据
set.seed(123)
proc_otu_FC <- procrustes(otu_scores, FC_scores, symmetric = TRUE)

# 执行Procrustes置换检验，计算p值
prot_2 <- protest(otu_scores, FC_scores, permutations = 999)


p_value_otu_FC  <- signif(prot_2$signif, 3) # 提取p值


# 提取Procrustes对齐后的数据
x_micro <- proc_otu_FC$Yrot[,1]  # 微生物PCA对齐后的x
y_micro <- proc_otu_FC$Yrot[,2]  # 微生物PCA对齐后的y
x_FC <- proc_otu_FC$X[,1]       # FC梯度PCA对齐后的x
y_FC <- proc_otu_FC$X[,2]       # FC梯度PCA对齐后的y

# 2. PCA分析：代谢物数据和FC数据之间
metab_data <- metabolimics_common %>% select(-Subid) %>% select_if(is.numeric)

# PCA分析
metab_pca <- prcomp(metab_data, scale. = TRUE)

# 提取PCA得分
metab_scores <- metab_pca$x

# Procrustes分析：对齐代谢物数据和FC数据
proc_metab_FC <- procrustes(metab_scores, FC_scores, symmetric = TRUE)

# 执行Procrustes置换检验，计算p值
prot_1 <- protest(metab_scores, FC_scores, permutations = 999)


p_value_metab_FC <- signif(prot_1$signif, 3)



# 提取Procrustes对齐后的数据
x_metab <- proc_metab_FC$X[,1]   # 代谢物PCA对齐后的x
y_metab <- proc_metab_FC$X[,2]   # 代谢物PCA对齐后的y

# 提取分组信息
group <- metagenomics_common$group  # 可以选择metagenomics_common$group或metabolimics_common$group

# 创建数据框用于绘图
df_otu_FC <- data.frame(
  Subid = metagenomics_common$Subid,
  group = group,
  x_micro = x_micro,
  y_micro = y_micro,
  x_FC = x_FC,
  y_FC = y_FC
)

df_metab_FC <- data.frame(
  Subid = metabolimics_common$Subid,
  group = group,
  x_metab = x_metab,
  y_metab = y_metab,
  x_FC = x_FC,
  y_FC = y_FC
)

# 计算解释的方差比例
otu_variance <- summary(otu_pca)$importance[2, 1:2]  # 提取前两个PCA组件的方差比例
FC_variance <- summary(FC_pca)$importance[2, 1:2]    # 提取前两个PCA组件的方差比例
metab_variance <- summary(metab_pca)$importance[2, 1:2]  # 提取前两个PCA组件的方差比例

# 分组配色
group_colors <- c("NA" = "#e76f51", "NNA" = "#457b9d")
# 自动适配实际分组名
names(group_colors) <- unique(df_otu_FC$group)

# 绘制Procrustes图：微生物和FC数据
plot_otu_FC <- ggplot(df_otu_FC) +
  # 连线：连接微生物和FC的PCA点
  geom_segment(aes(x = x_micro, y = y_micro, xend = x_FC, yend = y_FC), color = "grey70") +
  # 微生物点
  geom_point(aes(x = x_micro, y = y_micro, color = group), shape = 16, size = 2) +
  # FC点
  geom_point(aes(x = x_FC, y = y_FC, color = group), shape = 15, size = 2) +
  scale_color_manual(values = group_colors, name = "") +
  labs(
    x = paste("Axis1 (", round(otu_variance[1] * 100, 2), "%)", sep = ""),
    y = paste("Axis2 (", round(otu_variance[2] * 100, 2), "%)", sep = ""),
    title = sprintf("Procrustes: M² = %.3f, r ≈ %.3f, p = %.3f", round(proc_otu_FC$ss, 4), round(sqrt(1 - proc_otu_FC$ss), 3), p_value_otu_FC)
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

# 保存Procrustes图
ggsave("Metagenomics/Procrustes_plot_otu_FC.pdf", plot = plot_otu_FC, width = 10, height = 8)

# 绘制Procrustes图：代谢物和FC数据
plot_metab_FC <- ggplot(df_metab_FC) +
  # 连线：连接代谢物和FC的PCA点
  geom_segment(aes(x = x_metab, y = y_metab, xend = x_FC, yend = y_FC), color = "grey70") +
  # 代谢物点
  geom_point(aes(x = x_metab, y = y_metab, color = group), shape = 16, size = 2) +
  # FC点
  geom_point(aes(x = x_FC, y = y_FC, color = group), shape = 15, size = 2) +
  scale_color_manual(values = group_colors, name = "") +
  labs(
    x = paste("Axis1 (", round(metab_variance[1] * 100, 2), "%)", sep = ""),
    y = paste("Axis2 (", round(metab_variance[2] * 100, 2), "%)", sep = ""),
    title = sprintf("Procrustes: M² = %.3f, r ≈ %.3f, p = %.3f", round(proc_metab_FC$ss, 4), round(sqrt(1 - proc_metab_FC$ss), 3), p_value_metab_FC)
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
  ylim(c(-0.06, 0.06))

# 保存Procrustes图
ggsave("Metagenomics/Procrustes_plot_metab_FC.pdf", plot = plot_metab_FC, width = 10, height = 8)


# 调整M2值
M2 <- 0.736
r_approx <- round(sqrt(1 - M2), 3)
print(r_approx)
