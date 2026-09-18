
rm(list = ls()) #2.清理运行环境*`
library(tableone) #1.加载R包*`
library(xlsx)
library(tidyverse)
library(rstatix)
library(dplyr)
library(tidyr)
library(broom)
library(readxl)

# data_baseline <- read.xlsx("Data/N2PC_Demography.xlsx", 2)
data_baseline <- read_excel("Data/Metobolimics_data2.xlsx", 1)


head(data_baseline)
str(data_baseline)  


x <- colnames(data_baseline)[3:17]
x
x_cloname<- as.character(unlist(strsplit(x, split = ",")))
data_baseline$group <- as.factor(data_baseline$group)
x_cloname

# # 检验数据正态性
# 合起来不符合正态分布，单独则符合

x_cloname
gather_data <- data_baseline %>%
  gather(key = "type", value = "score", x_cloname[3:15]) %>% 
  convert_as_factor(group,type)

write.xlsx(gather_data,"data/Metobolimics_Baseline_results.xlsx",sheetName = "正态检验准备",append = T)
# 直接到excel中对数据进行更改，而后再导入此数据
gather_data <- read.xlsx("data/Metobolimics_Baseline_results.xlsx",sheetName = "正态检验准备")
# 需要去除NA值后进行正态性检验
data_shapiro <- drop_na(gather_data)

shapiro_results <- data_shapiro%>%
  group_by(group,type)%>%
  shapiro_test(score)

shapiro_results
# 保存正态性检验结果
write.xlsx(shapiro_results,"data/Metobolimics_Baseline_results.xlsx",sheetName = "正态检验结果",append = T)


myVars <- x_cloname
myVars 
# 分类变量名称
# catVars <- c( "Gender", "group", "Stage" )
catVars <- c( "Stage", "group")

table_bas <- CreateTableOne(vars = myVars, 
                            factorVars = catVars,#条件2*
                            strata = "group", #条件4*
                            data = data_baseline, #原始数据*
)
table_bas

table0 <- print(table_bas, #构建的table函数（带条件1.2.3）
                catDigits = 2,contDigits = 3,pDigits = 3, #附加条件
                showAllLevels=TRUE, #显示所有变量
                quote = FALSE, # 不显示引号
                # addOverall = TRUE,
                noSpaces = TRUE, # #删除用于对齐的空格
                printToggle = TRUE) #展示输出结果*`
table0

write.xlsx(table0,"data/Metobolimics_Baseline_results.xlsx",sheetName = "基线情况",append = T)
#不是正态分布的变量
nonvar <- c("Education", "MoCA")
# 条件5新加入 假如有T<5或n<40变量应使用Fisher精确检验，本文数量大，无需;例如再次加入性别变量
# exactvars <- c("Gender")
exactvars <- c()
table1<- print(table_bas, #构建的table函数（带条件1.2.3）
               # nonnormal = nonvar,#条件4
               exact = exactvars, #条件5
               catDigits = 2,contDigits = 3,pDigits = 3, #附加条件
               showAllLevels=TRUE, #显示所有变量
               quote = FALSE, # 不显示引号
               # addOverall = TRUE,
               noSpaces = TRUE, # #删除用于对齐的空格
               printToggle = TRUE) #展示输出结果*`
write.xlsx(table1,"data/Metobolimics_Baseline_results.xlsx",sheetName = "基线情况",append = T)



# 初始化一个空的数据框来存储统计检验结果
test_results <- data.frame(Variable = character(),
                           Statistic = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

x_cloname
myVars <-x_cloname[3:15] 

# 对连续变量进行t检验
for (var in myVars) {
  if (var %in% catVars) {
    next
  }
  formula_str <- paste(var, "~ group")
  t_test_result <- t.test(as.formula(formula_str), data = data_baseline)
  
  new_row <- data.frame(Variable = var,
                        Statistic = t_test_result$statistic,
                        p_value = t_test_result$p.value)
  
  test_results <- rbind(test_results, new_row)
}

# 对分类变量（如性别）进行卡方检验
for (var in catVars) {
  chisq_test_result <- chisq.test(table(data_baseline[, var], data_baseline$group))
  
  new_row <- data.frame(Variable = var,
                        Statistic = chisq_test_result$statistic,
                        p_value = chisq_test_result$p.value)
  
  test_results <- rbind(test_results, new_row)
}



# # 首先，确保我们的Stage列是字符类型，以便我们可以进行替换
# data_baseline$Stage <- as.character(data_baseline$Stage)

# 将 Stage 设为有序因子
data_baseline$Stage<- factor(data_baseline$Stage, levels = c("-2", "-1", "+1a", "+1b", "+1c"), ordered = TRUE)
data_baseline$Stage_num <- as.numeric(data_baseline$Stage)  # 转成有序数值
wilcox.test(Stage_num ~ group, data = data_baseline)

# 打印测试结果
print(test_result)
# 将table1转换为数据框
table1_df <- as.data.frame(table1)

# # 将统计检验结果添加到table1数据框中
# table1_df <- merge(table1_df, test_results, by = "Variable", all.x = TRUE)
# 
# # 保存新的数据框
# write.xlsx(table1_df, "data/N2PC_Demography-updated.xlsx", sheetName = "基线情况（更新）")
# 

# 计算Cohen's d
calculate_cohens_d <- function(group1, group2) {
  mean1 <- mean(group1)
  mean2 <- mean(group2)
  sd1 <- sd(group1)
  sd2 <- sd(group2)
  n1 <- length(group1)
  n2 <- length(group2)
  
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  cohens_d <- (mean1 - mean2) / pooled_sd
  return(cohens_d)
}

# 初始化一个空的数据框来存储Cohen's d
cohens_d_results <- data.frame(Variable = character(),
                               Cohens_d = numeric(),
                               stringsAsFactors = FALSE)

# 对连续变量计算Cohen's d
for (var in myVars) {
  if (var %in% catVars) {
    next
  }
  group1_data <- data_baseline[data_baseline$group == "HCs", var]
  group2_data <- data_baseline[data_baseline$group == "PMI", var]
  
  cohens_d <- calculate_cohens_d(unlist(group1_data),unlist(group2_data))
  
  new_row <- data.frame(Variable = var, Cohens_d = cohens_d)
  cohens_d_results <- rbind(cohens_d_results, new_row)
}

# 将Cohen's d结果合并到table1_df
table1_df <- merge(table1_df, cohens_d_results, by = "Variable", all.x = TRUE)


# 进行卡方检验
chisq_result <- chisq.test(table(data_baseline$Gender, data_baseline$group))

# 计算Cramér的V
n <- sum(chisq_result$observed)  # 总样本数
k <- min(nrow(chisq_result$observed), ncol(chisq_result$observed))  # 行或列的最小数
cramers_v <- sqrt(chisq_result$statistic / (n * (k - 1)))

cramers_v

