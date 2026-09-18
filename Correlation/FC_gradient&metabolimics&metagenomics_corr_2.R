rm(list = ls()) # 清理运行环境

# 加载所需包
library(tidyverse)
library(vegan)
library(ggplot2)

# 读取数据
metagenomics_abundance_data <- read.csv("Metagenomics/metagenomics_abundance_sigdata.csv")
metabolimics_abundance_data <- read.csv("Metabolimics/Significant_Metabolites_Abundance.csv")
FC_Gradient_data <- read.csv("Data/FC_gradient_all_with_Subid.csv")

# 去掉group列
FC_Gradient_data <- FC_Gradient_data %>% select(-group)
metagenomics_abundance_data <- metagenomics_abundance_data %>% select(-group)

# 求共有的 Subid，确保这三个数据集都有
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

# 将数据转化为不含Subid的矩阵
metagenomics_mat <- as.data.frame(metagenomics_common %>% select(-Subid))
fc_mat <- as.data.frame(FC_Gradient_common %>% select(-Subid))

# 相关性分析函数
get_correlation_stats <- function(metab_mat, fc_mat) {
  result <- data.frame(
    Metabolite = rep(colnames(metab_mat), each = ncol(fc_mat)),
    BrainRegion = rep(colnames(fc_mat), times = ncol(metab_mat)),
    correlation = NA,
    pvalue = NA
  )
  
  for (i in seq_along(colnames(metab_mat))) {
    for (j in seq_along(colnames(fc_mat))) {
      # 获取代谢物数据与脑区数据
      metab_j <- metab_mat[[i]]
      fc_j <- fc_mat[[j]]
      
      # 进行Spearman相关性检验
      ct <- suppressWarnings(cor.test(metab_j, fc_j, method = "spearman"))
      
      if (!is.na(ct$p.value) && ct$p.value < 0.05) {
        result$correlation[(i-1)*ncol(fc_mat) + j] <- ct$estimate
        result$pvalue[(i-1)*ncol(fc_mat) + j] <- ct$p.value
      }
    }
  }
  return(result)
}

# 执行相关性分析
correlation_stats <- get_correlation_stats(metagenomics_mat, fc_mat)

# 筛选显著相关性
significant_correlations <- correlation_stats %>% 
  filter(!is.na(correlation)) %>%
  filter(pvalue < 0.05 & abs(correlation) > 0.05)

# 查看结果
print(significant_correlations)

# 生成边表
edge_list <- data.frame(
  Source = significant_correlations$Metabolite,
  Target = significant_correlations$BrainRegion,
  correlation = significant_correlations$correlation,
  pvalue = significant_correlations$pvalue,
  stringsAsFactors = FALSE
)

# 标记正相关和负相关
edge_list$type <- ifelse(edge_list$correlation > 0, "positive", "negative")

# 导出边表
write.csv(edge_list, "corr3/network_metagenomics_fc_gradient_edges.csv", row.names = FALSE)

# 获取所有的代谢物和脑区
all_metabolites <- unique(significant_correlations$Metabolite)
all_brainregions <- unique(significant_correlations$BrainRegion)

# 生成节点表
node_table <- data.frame(
  id = c(all_metabolites, all_brainregions),
  type = c(rep("metabolite", length(all_metabolites)), rep("brain_region", length(all_brainregions))),
  stringsAsFactors = FALSE
)

# 统计每个节点在edge_list中作为Source或Target出现的次数
node_table$count <- sapply(node_table$id, function(x) {
  # 统计Source和Target中出现的次数
  sum(edge_list$Source == x) + sum(edge_list$Target == x)
})

# 导出节点表
write.csv(node_table, "corr3/network_metagenomics_fc_gradient_nodes.csv", row.names = FALSE)

# 清理节点和边表中的特殊字符
clean <- function(x) {
  x <- trimws(x)                 # 去除前后空格
  x <- gsub('^"|"$', "", x)      # 去掉包围的双引号
  x <- gsub("“|”", "", x)        # 去掉中文引号
  x <- gsub("\u00A0", " ", x)    # 替换不间断空格
  return(x)
}

# 清理边表和节点表中的字段
edge_list$Source <- clean(edge_list$Source)
edge_list$Target <- clean(edge_list$Target)
node_table$id <- clean(node_table$id)

# 保存清理后的文件
write.csv(edge_list, "corr3/network_metagenomics_fc_gradient_edges_clean.csv",
          row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

write.csv(node_table, "corr3/network_metagenomics_fc_gradient_nodes_clean.csv",
          row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")



# 读取Enriched_group_data
Enriched_group_data <- read.csv("corr3/enriched_group_information.csv")

# 查看Enriched_group_data
head(Enriched_group_data)

# 在Enriched_group_data中的id列添加点（.）以适应R中的格式
Enriched_group_data$id <- gsub(" ", ".", Enriched_group_data$id)

# 使用left_join将Enriched_group_data合并到node_table中
node_table <- node_table %>%
  left_join(Enriched_group_data, by = c("id" = "id"))

# 查看结果，确认enriched_group信息已加入
head(node_table)

# 导出更新后的节点表
write.csv(node_table, "corr3/network_metagenomics_fc_gradient_nodes_with_enriched_group.csv",
          row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")


# 获取type总数
metagenomics_type_count <- ncol(metagenomics_abundance_data)  # 代谢物的总数
brain_region_type_count <- ncol(FC_Gradient_data)  # 脑区的总数

# 计算不重复的id名占各自type总数的比例
# 统计所有显著相关的成分数量
significant_metagenomics <- unique(edge_list$Source)
significant_brainregions <- unique(edge_list$Target)

significant_metagenomics_count <- length(significant_metagenomics)
significant_brainregions_count <- length(significant_brainregions)

# 计算不重复的id名占各自type总数的比例
metagenomics_ratio <- significant_metagenomics_count / metagenomics_type_count
brainregion_ratio <- significant_brainregions_count / brain_region_type_count

# 输出比例
cat("Metabolite ratio:", metagenomics_ratio, "\n")
cat("Brain Region ratio:", brainregion_ratio, "\n")


