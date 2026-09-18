
rm(list = ls()) #2.清理运行环境*`
library(tableone) #1.加载R包*`
library(xlsx)
library(tidyverse)
library(rstatix)
library(dplyr)
library(tidyr)
library(broom)
library(readxl)


data_baseline <- read_excel("Data/Acupuncture_Demo_data.xlsx", 2)

head(data_baseline)
str(data_baseline)  



x <- colnames(data_baseline)[4:33]
x
x_cloname<- as.character(unlist(strsplit(x, split = ",")))
data_baseline$group <- as.factor(data_baseline$group)
x_cloname


# # 检验数据正态性
# 合起来不符合正态分布，单独则符合
x_cloname
gather_data <- data_baseline %>%
  gather(key = "type", value = "score", x_cloname[1:30]) %>% 
  convert_as_factor(group,type)

write.xlsx(gather_data,"RCT/Metobolimics_Baseline_results.xlsx",sheetName = "正态检验准备",append = T)
# 直接到excel中对数据进行更改，而后再导入此数据
gather_data <- read.xlsx("RCT/Metobolimics_Baseline_results.xlsx",sheetName = "正态检验准备")
# 需要去除NA值后进行正态性检验
data_shapiro <- drop_na(gather_data)

shapiro_results <- data_shapiro%>%
  group_by(group,type)%>%
  shapiro_test(score)

shapiro_results
# 保存正态性检验结果
write.xlsx(shapiro_results,"RCT/Metobolimics_Baseline_results.xlsx",sheetName = "正态检验结果",append = T)


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

write.xlsx(table0,"RCT/Metobolimics_Baseline_results.xlsx",sheetName = "基线情况",append = T)
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
               printToggle = FALSE) #展示输出结果*`

write.xlsx(table1,"RCT/Acupuncture_Demo_data_results.xlsx",sheetName = "基线情况")
# 保存为 CSV
write.csv(table1,
          file = "RCT/Acupuncture_Demo_data_results.csv",
          row.names = TRUE,    # 保留变量名在第一列
          fileEncoding = "UTF-8")

# 初始化一个空的数据框来存储统计检验结果
test_results <- data.frame(Variable = character(),
                           Statistic = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

x_cloname
myVars <-x_cloname[1:30] 

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

# 打印测试结果
print(test_results)



# 清理环境
# ================================
# 0) 环境与依赖
# ================================
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readxl)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(openxlsx)
})

# 小助手：若对象不存在/为空/全NA，则给默认值
`%OR%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (is.atomic(x) && all(is.na(x))) return(y)
  x
}

# ================================
# 1) 读取数据并筛选
# ================================
data_baseline <- read_excel("Data/Acupuncture_Demo_data.xlsx", sheet = 2)
data_RCT      <- read_excel("Data/Acupuncture_Demo_data.xlsx", sheet = 1)

# 统一关键列名（兼容大小写）
names(data_baseline) <- sub("^subid$", "Subid", names(data_baseline), ignore.case = TRUE)
names(data_RCT)      <- sub("^rawid$", "Rawid", names(data_RCT), ignore.case = TRUE)
names(data_RCT)      <- sub("^group$", "Group", names(data_RCT), ignore.case = TRUE)
names(data_RCT)      <- sub("^condition$", "Condition", names(data_RCT), ignore.case = TRUE)
names(data_RCT)      <- sub("^age$", "Age", names(data_RCT), ignore.case = TRUE)

dat <- data_RCT %>% filter(Rawid %in% data_baseline$Subid)
# 保存结果为新的 Excel
library(openxlsx)
write.xlsx(dat, "RCT/Acupuncture_Demo_data.xlsx", sheetName = "data_RCT_filtered",append = T)

# ================================
# 2) 因子与协变量设置（按你数据实际值）
# ================================
# 确保 Group 含有 "Sham" "Acupuncture"；Condition 含 "Pre" "Post"
if (!all(c("Sham","Acupuncture") %in% unique(dat$Group))) {
  stop("Group 中未包含期望水平：'Sham','Acupuncture'")
}
if (!all(c("Pre","Post") %in% unique(dat$Condition))) {
  stop("Condition 中未包含期望水平：'Pre','Post'")
}

dat <- dat %>%
  mutate(
    Rawid     = factor(Rawid),
    Group     = factor(Group, levels = c("Sham","Acupuncture")),  # 参照 Sham
    Condition = factor(Condition, levels = c("Pre","Post")),      # 参照 Pre
    age_c     = as.numeric(scale(Age, center = TRUE, scale = FALSE))
  )

# 结局列表（按你 head(dat) 的列）
clinical_outcomes <- c("PSQI","SAS","KMI","MRS","SOL","TST","TIB","SE","WASO")
clinical_outcomes <- clinical_outcomes[clinical_outcomes %in% names(dat)]

# 输出目录（解决 Permission denied）
dir.create("RCT", showWarnings = FALSE, recursive = TRUE)

# ================================
# 3) 定义单结局的 LMM + EMM + DiD 流程
# ================================
fit_one_outcome <- function(outcome) {
  message(sprintf("==> Fitting LMM for %s", outcome))
  # 仅保留该结局有观测的数据行，避免奇异拟合
  dsub <- dat %>% filter(!is.na(.data[[outcome]]))
  
  # 至少每组每时点要有数据
  chk <- dsub %>% count(Group, Condition)
  if (nrow(chk) < 4 || any(chk$n == 0)) {
    return(list(
      outcome = outcome,
      error = sprintf("数据不足：%s（某些 Group×Condition 组合缺失）", outcome)
    ))
  }
  
  # 拟合 LMM
  f <- as.formula(paste0(outcome, " ~ Group * Condition + age_c + (1|Rawid)"))
  m <- lmer(f, data = dsub, REML = TRUE)
  
  # 四格 EMM
  emm <- emmeans(m, ~ Group * Condition)
  
  # 两时点组间差 (Acupuncture − Sham)，分层按 Condition
  gd <- emmeans(m, ~ Group | Condition)
  gd_con <- contrast(gd, method = "revpairwise")  # (Acupuncture − Sham)
  gd_sum <- summary(gd_con, infer = c(TRUE, TRUE)) %>% as.data.frame()
  
  # DiD：先在每组内部 Post-Pre，再(Acu − Sham)
  wdiff_by_grp <- contrast(emmeans(m, ~ Condition | Group), method = "revpairwise")  # (Post − Pre)
  did <- contrast(wdiff_by_grp, method = list("DiD (Acupuncture vs Sham)" = c(-1, 1)), by = NULL)
  did_sum <- summary(did, infer = c(TRUE, TRUE)) %>% as.data.frame()
  
  # Type III 主效应与交互
  a3 <- anova(m, type = 3) %>% as.data.frame()
  p_group <- a3["Group","Pr(>F)"]            %OR% NA_real_
  p_time  <- a3["Condition","Pr(>F)"]        %OR% NA_real_
  p_int   <- a3["Group:Condition","Pr(>F)"]  %OR% NA_real_
  
  list(
    outcome   = outcome,
    model     = m,
    emm       = emm,
    gd_sum    = gd_sum,
    did_sum   = did_sum,
    p_group   = p_group,
    p_time    = p_time,
    p_int     = p_int
  )
}

# 批量拟合（容错）
fits <- map(clinical_outcomes, ~tryCatch(fit_one_outcome(.x), error = function(e) list(outcome=.x, error=e$message)))
names(fits) <- clinical_outcomes

# ================================
# 4) 生成主结果表
# ================================
fmt_mean_sd <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (is.finite(m) && is.finite(s)) sprintf("%.1f (%.1f)", m, s) else NA_character_
}
cond_map <- c(Pre = "Baseline", Post = "4 wk")

make_row <- function(res) {
  oc <- res$outcome
  if (!is.null(res$error)) {
    return(tibble::tibble(
      Outcome = oc,
      Note = res$error
    ))
  }
  
  # 观测均值(SD)
  obs_ms <- dat %>%
    filter(!is.na(.data[[oc]])) %>%
    complete(Group, Condition) %>%
    group_by(Group, Condition) %>%
    summarise(ms = if (all(is.na(.data[[oc]]))) NA_character_ else fmt_mean_sd(.data[[oc]]),
              .groups = "drop") %>%
    pivot_wider(names_from = c(Group, Condition), names_sep = "_", values_from = ms)
  
  exp_pre  <- obs_ms[["Acupuncture_Pre"]] %OR% NA_character_
  exp_post <- obs_ms[["Acupuncture_Post"]] %OR% NA_character_
  sha_pre  <- obs_ms[["Sham_Pre"]]        %OR% NA_character_
  sha_post <- obs_ms[["Sham_Post"]]       %OR% NA_character_
  
  # 参与者数（按受试者）
  n_tab <- dat %>%
    filter(!is.na(.data[[oc]])) %>%
    distinct(Rawid, Group) %>%
    count(Group) %>%
    pivot_wider(names_from = Group, values_from = n, values_fill = 0)
  n_exp  <- n_tab[["Acupuncture"]] %OR% 0
  n_sham <- n_tab[["Sham"]]        %OR% 0
  
  # 组间差（Acu − Sham），分时点
  mdiff <- res$gd_sum %>%
    mutate(Time = dplyr::recode(Condition, !!!cond_map)) %>%
    transmute(
      Time,
      MeanDiff_CI = sprintf("%.2f (%.2f to %.2f)", estimate, lower.CL, upper.CL)
    )
  mdiff_baseline <- mdiff$MeanDiff_CI[mdiff$Time == "Baseline"] %OR% NA_character_
  mdiff_4wk      <- mdiff$MeanDiff_CI[mdiff$Time == "4 wk"]     %OR% NA_character_
  
  # DiD
  did_row <- res$did_sum[1, , drop = FALSE]
  did_est <- if (nrow(did_row)) sprintf("%.3f", did_row$estimate) else NA_character_
  did_ci  <- if (nrow(did_row)) sprintf("%.3f to %.3f", did_row$lower.CL, did_row$upper.CL) else NA_character_
  did_p   <- if (nrow(did_row)) signif(did_row$p.value, 3) else NA
  
  tibble::tibble(
    Outcome = oc,
    `No. of participants (Experimental)` = n_exp,
    `No. of participants (Sham)`         = n_sham,
    
    `Experimental mean (SD) at Baseline` = exp_pre,
    `Sham mean (SD) at Baseline`         = sha_pre,
    `Mean difference at Baseline (95% CI)` = mdiff_baseline,
    
    `Experimental mean (SD) at 4 wk`     = exp_post,
    `Sham mean (SD) at 4 wk`             = sha_post,
    `Mean difference at 4 wk (95% CI)`   = mdiff_4wk,
    
    `P (Group effect)`       = signif(res$p_group, 3),
    `P (Time effect)`        = signif(res$p_time, 3),
    `P (Interaction effect)` = signif(res$p_int, 3),
    
    `DiD estimate`           = did_est,
    `DiD 95% CI`             = did_ci,
    `DiD p-value`            = did_p
  )
}

final_table <- map_dfr(fits, make_row) %>% relocate(Outcome)
print(final_table, n = Inf)

# ================================
# 5) 导出所有结果
# ================================
# 5.1 主表
write.xlsx(final_table,
           file = "RCT/primary_outcome_table_LMM.xlsx",
           sheetName = "Primary",
           overwrite = TRUE)

head(final_table)

# 5.2 每个结局：四格 EMM、两时点组间差、DiD
wb <- createWorkbook()
addWorksheet(wb, "Index")
writeData(wb, "Index", data.frame(Outcome = names(fits)), withFilter = TRUE)

for (nm in names(fits)) {
  res <- fits[[nm]]
  if (!is.null(res$error)) next
  
  # EMM 四格
  emm_df <- as.data.frame(res$emm)
  addWorksheet(wb, paste0(nm, "_EMM"))
  writeData(wb, paste0(nm, "_EMM"), emm_df)
  
  # 两时点组间差
  addWorksheet(wb, paste0(nm, "_Diff_byTime"))
  writeData(wb, paste0(nm, "_Diff_byTime"), res$gd_sum)
  
  # DiD
  addWorksheet(wb, paste0(nm, "_DiD"))
  writeData(wb, paste0(nm, "_DiD"), res$did_sum)
}
saveWorkbook(wb, "RCT/emm_and_DiD_details.xlsx", overwrite = TRUE)



