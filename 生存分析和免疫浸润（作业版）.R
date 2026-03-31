rm(list=ls())
# 加载必要的包
library(TCGAbiolinks)
library(survival)
library(survminer)
library(glmnet)
library(caret)
library(dplyr)
#下载TCGA-LUAD临床信息（肺腺癌）
TCGAbiolinks:::getGDCprojects()$project_id
cancer_type="TCGA-LUAD"
clinical_data=GDCquery_clinic(project=cancer_type,type="clinical")
dim(clinical)
clinical[1:4,1:4]
head(colnames(clinical))
# 查询和下载表达数据（RNA-seq）
query_exp  <- GDCquery(project = cancer_type, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification")

GDCdownload(query_exp)
exp_data <- GDCprepare(query_exp)
# 保存原始数据
save(clinical_data, exp_data, file = "TCGA-LUAD_rawdata.RData")

load("E:\\生信案例分析与实践\\专题四 生存分析\\TCGA-LUAD_rawdata.RData")
#数据预处理
# 安装并加载SummarizedExperiment包
if(!require(SummarizedExperiment)) {
  if(!require(BiocManager)) install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
  library(SummarizedExperiment)
}
# 提取表达矩阵（count数据）
count_matrix <- assay(exp_data)
#步骤1：提取表达数据的样本类型并筛选肿瘤样本
# 提取样本类型代码（第4部分的前两位）
sample_names <- colnames(count_matrix)
sample_type <- substr(sample_names, 14, 15)

# 筛选肿瘤样本（01表示原发肿瘤）
#tumor_samples <- sample_type == "01"
#count_matrix_tumor <- count_matrix[, tumor_samples]
# 筛选肿瘤样本 (01-09: 原发肿瘤, 10-19: 正常组织)
tumor_samples <- colnames(count_matrix)[sample_type %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09")]
normal_samples<-colnames(count_matrix)[sample_type%in% 10:19]
count_matrix_tumor <- count_matrix[, tumor_samples]
#步骤2：提取肿瘤样本的患者ID
# 提取患者ID（前三个部分）
patient_ids_tumor <- sapply(strsplit(colnames(count_matrix_tumor), "-"), 
                            function(x) paste(x[1:3], collapse = "-"))
#步骤3：匹配临床数据中的患者ID
#确保表达数据中的患者ID存在于临床数据中：
# 过滤表达数据，仅保留有临床信息的样本
valid_patients <- patient_ids_tumor %in% clinical_data$submitter_id
count_matrix_filtered <- count_matrix_tumor[, valid_patients]
patient_ids_filtered <- patient_ids_tumor[valid_patients]

# 检查匹配数量
ncol(count_matrix_filtered)  # 应接近临床数据样本数（585）
#步骤4：处理同一患者的多个样本（若有）
# 去除重复患者ID，保留每个患者的第一个样本
keep <- !duplicated(patient_ids_filtered)
count_matrix_final <- count_matrix_filtered[, keep]
patient_ids_final <- patient_ids_filtered[keep]

# 设置列名为患者ID
colnames(count_matrix_final) <- patient_ids_final
#步骤5：验证与临床数据的匹配
all(colnames(count_matrix_final) %in% clinical_data$submitter_id)  # 应为TRUE
#步骤6：合并数据进行分析:将临床数据与表达数据按患者ID对齐
# 按submitter_id排序临床数据
clinical_data <- clinical_data[match(colnames(count_matrix_final), clinical_data$submitter_id), ]

# 确认行顺序一致
all(colnames(count_matrix_final) == clinical_data$submitter_id)  # 应为TRUE

#tips:若临床数据中存在未在表达数据中出现的患者，可进一步筛选临床数据
rownames(count_matrix_final)
# 基因注释
gene_info <- rowData(exp_data)
dim(gene_info)
###将表达矩阵的行名（ENSEMBL ID）转换为基因名（gene symbol）
#步骤1：构建基因ID与基因名的映射关系
# 从gene_info中提取基因ID和基因名
gene_id_to_name <- data.frame(
  ensembl_id = gene_info$gene_id,
  gene_name = gene_info$gene_name
)

# 查看前5个映射关系
head(gene_id_to_name)
#步骤2：替换行名为基因名
# 1. 匹配行名与gene_id
matched_indices <- match(rownames(count_matrix_final), gene_id_to_name$ensembl_id)

# 2. 提取对应的gene_name
new_gene_names <- gene_id_to_name$gene_name[matched_indices]

# 3. 处理NA值（未匹配的ENSEMBL ID保留原始ID）
new_gene_names <- ifelse(is.na(new_gene_names), 
                         rownames(count_matrix), 
                         new_gene_names)

# 4. 处理重复基因名（添加后缀保证唯一性）
rownames(count_matrix_final) <- make.unique(new_gene_names)

# === 新增：加载编码基因列表并过滤 ===
load("E:\\大作业\\gene.Rdata")  # 加载编码基因列表
coding_genes <- uniqe_gene  # 假设对象名为uniqe_gene
# 过滤非编码基因（只保留编码蛋白质基因）
count_matrix_final <- count_matrix_final[rownames(count_matrix_final) %in% coding_genes, ]

# 检查过滤结果
cat("过滤后基因数量:", nrow(count_matrix_final), "\n")
cat("编码基因占比:", mean(rownames(count_matrix_final) %in% coding_genes), "\n")

#表达数据标准化
# 基因表达标准化（DESeq2方差稳定变换）
# 加载DESeq2包
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_matrix_final,
                              colData = clinical_data,# # 临床数据作为样本注释
                              design = ~1)# 若无分组，使用~1
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=FALSE)
dataExp <- assay(vsd)
# 保存数据
save(clinical_data, exp_data, file = "TCGA-LUAD_rawdata.RData")
# 准备生存数据
# 创建生存分析数据集 ------------------------------------------------------------
sur_data<-NULL
sur_data <- data.frame(
  # 基础信息
  patient.id = clinical_data$bcr_patient_barcode,
  
  # ============== 生存数据 ==============
  # 总生存时间（天）：死亡患者取死亡时间，存活患者取最后随访时间
  OS.time = ifelse(clinical_data$vital_status == "Dead",
                   clinical_data$days_to_death,
                   clinical_data$days_to_last_follow_up),
  stage = clinical_data$ajcc_pathologic_stage,
  t_stage = clinical_data$ajcc_pathologic_t,
  n_stage = clinical_data$ajcc_pathologic_n,
  m_stage = clinical_data$ajcc_pathologic_m,
  # 生存状态：死亡=1，存活=0
  OS = as.numeric(clinical_data$vital_status == "Dead"),
  # ============ 人口统计学 ==============
  # 年龄处理：将天数转换为年（保留1位小数）
  age = round(clinical_data$age_at_diagnosis / 365, 1),
  # 性别处理：转换为因子
  gender = factor(clinical_data$gender,
                  levels = c("male", "female"),
                  labels = c("Male", "Female")),
  
  # 种族处理：简化分类（合并小样本）
  race = factor(dplyr::case_when(
    clinical_data$race %in% c("white", "White") ~ "White",
    clinical_data$race %in% c("black or african american", "Black") ~ "Black",
    clinical_data$race %in% c("asian", "Asian") ~ "Asian",
    TRUE ~ "Other/Unknown"
  )),
  # ============ 生活方式 ================
  smoking_status = factor(
    case_when(
      # 精确匹配终生不吸烟者
      grepl("^Lifelong\\s*Non[-\\s]?Smoker$", 
            clinical_data$tobacco_smoking_status, 
            ignore.case = TRUE) ~ "Never Smoker",
      
      # 捕获所有Reformed变体（包括不同戒烟时长）
      grepl("Current Reformed Smoker", 
            clinical_data$tobacco_smoking_status) ~ "Former Smoker",
      
      # 精确匹配当前吸烟者
      clinical_data$tobacco_smoking_status == "Current Smoker" ~ "Current Smoker",
      
      # 处理特殊类别
      clinical_data$tobacco_smoking_status %in% 
        c("Not Reported", "Unknown") ~ "Missing/Unknown",
      
      # 安全处理未预见类别
      TRUE ~ "Missing/Unknown"
    ),
    levels = c("Never Smoker", "Former Smoker", "Current Smoker", "Missing/Unknown")
  ),
  # 吸烟包年数处理：保留原始数值（NA表示缺失）
  pack_years = clinical_data$pack_years_smoked)
save(sur_data, dataExp, file = "TCGA-LUAD_norm_data.RData")

#########生存分析与免疫浸润#####################################################
load("E:\\生信案例分析与实践\\专题四 生存分析\\TCGA-LUAD_norm_data.RData")
Index1<-which(is.na(sur_data$OS))
Index2<-which(is.na(sur_data$OS.time)|sur_data$OS.time<=0)
Index<-unique(c(Index1,Index2))

data_Exp<-dataExp[,-Index]#数据框[1] 19938   503
dim(data_Exp)
Sur_data<-sur_data[-Index,]#数据框,[1] 503  12
dim(Sur_data)
###生存分析#########
#（1）利用 LASSO 挑选标记构建风险模型
# 设置随机种子确保结果可重复
set.seed(123)

# （1）随机划分训练集和验证集 (70%训练集，30%验证集)
train_indices <- sample(1:ncol(data_Exp), 0.7 * ncol(data_Exp))
dataExp.Train <- data_Exp[, train_indices]
dim(dataExp.Train)#[1] 19938   352

dataExp.Test <- data_Exp[, -train_indices]
dim(dataExp.Test)#[1] 19938   151

sur_data.Train <- Sur_data[train_indices, ]#[1] 352  12
dim(sur_data.Train)

sur_data.Test <- Sur_data[-train_indices, ]#[1] 151  12
dim(sur_data.Test)

# 准备训练集的生存标签
Label1 <- data.frame(
  OS.time = sur_data.Train$OS.time,
  OS = sur_data.Train$OS
)

# 准备验证集的生存标签
dataLabel.Test <- data.frame(
  OS.time = sur_data.Test$OS.time,
  OS = sur_data.Test$OS
)

# 识别与预后相关的基因(是否是生存相关因素）
library("survival")
m <- dim(dataExp.Train)[1]
sur_ass_gene <- matrix(nrow = m, ncol = 2)
rownames(sur_ass_gene) <- rownames(dataExp.Train)

for (i in 1:m) {
  Sur_Obj <- summary(coxph(Surv(Label1[,'OS.time'], Label1[,'OS']) ~ dataExp.Train[i,]))
  sur_ass_gene[i,] <- c(rownames(dataExp.Train)[i], Sur_Obj$logtest[3])
}

# 添加FDR校正
sur_ass_gene <- cbind(sur_ass_gene, p.adjust(as.numeric(sur_ass_gene[,2]), method = "BH"))
colnames(sur_ass_gene) <- c("gene", "pvalue", "padj")
sur_ass_gene<-data.frame(sur_ass_gene)

# 选择前500个最显著的基因
a <- order(as.numeric(sur_ass_gene[,3]), decreasing = FALSE)
Top500surgene <- a[1:500]
dataExp_Top500 <- dataExp.Train[Top500surgene,]
sur_ass_gene1 <- sur_ass_gene[Top500surgene,]

# 构建LASSO风险模型
library(glmnet)
Y <- Surv(Label1[,'OS.time'], Label1[,'OS'])
cvfit <- cv.glmnet(t(dataExp_Top500), Y, family = "cox", nlambda = 1000, alpha = 1)
print(cvfit)#列是样本、行是特征
plot(cvfit)#  min 46 ise 0
title(main = "LASSO Cox模型\n交叉验证曲线", line = 2.5, cex.main = 0.95)

# LASSO选择基因和系数
coefficients <- coef(cvfit, s = cvfit$lambda.min)
coefficients1 <- coef(cvfit, s = cvfit$lambda.1se)
#查看基因和系数
coef_min<-coef(cvfit, s = cvfit$lambda.min)
risk_model1111 <- data.frame(
  gene = coef_min@Dimnames[[1]][which(coef_min != 0)],  # 基因名
  coef = coef_min[which(coef_min != 0)]                 # 对应系数
)

# 选择非零系数的基因 ise
Active.Index <- which(as.numeric(coefficients) != 0)
Active.coef <- as.numeric(coefficients)[Active.Index]
Sig_gene_cox <- rownames(coefficients)[Active.Index]

# 计算训练集风险评分
riskscore1 <- predict(cvfit, t(dataExp_Top500), s = c(cvfit$lambda.min, cvfit$lambda.1se))
thrhold <- median(riskscore1[,1])
High_lowIndex <- ifelse(riskscore1[,1] >= thrhold, "High_risk", "Low_risk")

# 训练集生存分析
fit31_1 <- survfit(Surv(Label1[,'OS.time'], Label1[,'OS']) ~ High_lowIndex)
# 绘制训练集KM曲线
library(survminer)
ggsurvplot(fit31_1, 
           data = data.frame(High_lowIndex), 
           conf.int = TRUE, 
           pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, 
           title = "LASSO Cox模型训练集生存分析",  # 主标题
           ggtheme = theme_survminer() +          # 自定义主题
             theme(plot.title = element_text(
               size = 14,                        # 字体大小
               face = "bold",                    # 加粗
               hjust = 0.5,                      # 水平居中
               margin = margin(b = 10))),         # 下边距10pt
           xlab = "Time (days)", 
           legend = c(0.8, 0.9), 
           legend.title = "", 
           palette = c("#E7B800", "#2E9FDF"))

# 验证集风险评分预测
dataExp_Top500.test <- dataExp.Test[Top500surgene,]
riskscore2 <- predict(cvfit, t(dataExp_Top500.test), s = c(cvfit$lambda.min, cvfit$lambda.1se))
High_lowIndex.test <- ifelse(riskscore2[,1] >= thrhold, "High_risk", "Low_risk")

# 验证集生存分析
fit31_1.test <- survfit(Surv(dataLabel.Test[,'OS.time'], dataLabel.Test[,'OS']) ~ High_lowIndex.test)

# 绘制验证集KM曲线
ggsurvplot(fit31_1.test, 
           data = data.frame(High_lowIndex.test), 
           conf.int = TRUE, 
           pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, 
           title = "LASSO Cox模型测试集生存分析",  # 主标题
           ggtheme = theme_survminer() +          # 自定义主题
             theme(plot.title = element_text(
               size = 14,                        # 字体大小
               face = "bold",                    # 加粗
               hjust = 0.5,                      # 水平居中
               margin = margin(b = 10))),         # 下边距10pt
           xlab = "Time (days)", 
           legend = c(0.8, 0.8), 
           legend.title = "", 
           palette = c("#E7B800", "#2E9FDF"))

##逐步回归进一步挑选#####
f11<-coxph(Surv(Label1[,'OS.time'],Label1[,'OS'])~.,data=as.data.frame(t(dataExp_Top500[Active.Index,])))
library(MASS)
f22<-stepAIC(f11)
###--------------------------------------------------
riskscore33 <- predict(f22, as.data.frame(t(dataExp_Top500[Active.Index,])))
thrhold33<-median(riskscore33)
High_lowIndex<-ifelse(riskscore33>=thrhold33,"High_risk","Low_risk")
fit331_2 <- survfit(Surv(Label1[,'OS.time'],Label1[,'OS']) ~ High_lowIndex)
#逐步回归训练集KM曲线
ggsurvplot(fit331_2, # 创建的拟合对象
           data = as.factor(High_lowIndex),  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           xlab = "tme(d)", # 指定x轴标签
           legend = c(0.8,0.9), # 指定图例位置
           title = "逐步 Lasso 回归模型训练集 KM 曲线",  # 主标题
           ggtheme = theme_survminer() +          # 自定义主题
             theme(plot.title = element_text(
               size = 14,                        # 字体大小
               face = "bold",                    # 加粗
               hjust = 0.5,                      # 水平居中
               margin = margin(b = 10))),         # 下边距10pt
           #break.x.by = 100)  # 设置x轴刻度间距
           palette=c("#E7B800", "#2E9FDF"))


riskscore44 <- predict(f22, as.data.frame(t(dataExp_Top500.test[Active.Index,])))
High_lowIndex<-ifelse(riskscore44>=thrhold33,"High_risk","Low_risk")
fit331_2.test <- survfit(Surv(dataLabel.Test[,'OS.time'],dataLabel.Test[,'OS']) ~ High_lowIndex)
#逐步回归验证集KM曲线
ggsurvplot(fit331_2.test, # 创建的拟合对象
           data = as.factor(High_lowIndex),  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           xlab = "time(d)", # 指定x轴标签
           legend = c(0.8,0.9), # 指定图例位置
           title = "逐步 Lasso 回归模型测试集 KM 曲线",  # 主标题
           ggtheme = theme_survminer() +          # 自定义主题
             theme(plot.title = element_text(
               size = 14,                        # 字体大小
               face = "bold",                    # 加粗
               hjust = 0.5,                      # 水平居中
               margin = margin(b = 10))),         # 下边距10pt
           #break.x.by = 100)  # 设置x轴刻度间距
           palette=c("#E7B800", "#2E9FDF"))
#森林图
ggforest(model = f22,#cox逐步回归方程
         main = "senlintu",#标题
         cpositions = c(0.05, 0.15, 0.35),#前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 )#保留HR值以及95%CI的小数位数


############（2）加入临床信息，重新构建模型#####
#将基因表达数据与临床特征结合，构建多变量模型，并评估临床特征的影响
# 准备临床特征数据
# 处理分类变量
sur_data.Train$stage <- factor(sur_data.Train$stage)
sur_data.Train$gender <- factor(sur_data.Train$gender)
sur_data.Train$race <- factor(sur_data.Train$race)
sur_data.Train$smoking_status <- factor(sur_data.Train$smoking_status)

# 创建包含基因表达和临床特征的数据集
# 提取LASSO选择的基因表达
gene_expression_train <- t(dataExp_Top500[Active.Index, ])
gene_expression_test <- t(dataExp_Top500.test[Active.Index, ])

# 合并基因表达和临床特征
clinical_features <- c("age", "gender", "stage", "t_stage", "n_stage", "m_stage", "smoking_status")
combined_data_train <- data.frame(gene_expression_train, sur_data.Train[, clinical_features])
combined_data_test <- data.frame(gene_expression_test, sur_data.Test[, clinical_features])

# 缺失值：使用多重插补（更复杂但更可靠）
# 检查各变量的缺失值情况
missing_summary <- sapply(cbind(sur_data.Train, gene_expression_train), function(x) sum(is.na(x)))
print(missing_summary[missing_summary > 0])  # 显示有缺失值的变量
#BiocManager::install("mice")
library(mice)
imputed_data <- mice(cbind(sur_data.Train, gene_expression_train), m=5, maxit=50, method='pmm')
data_imputed <- complete(imputed_data, 1)  # 使用第一次插补结果

# 使用删除缺失值后的数据
combined_formula <- as.formula(paste("Surv(OS.time, OS) ~", 
                                     paste(c(colnames(gene_expression_train), clinical_features), collapse = " + ")))
full_model <- coxph(combined_formula, data = data_imputed)  # 使用处理后的完整数据

#使用逐步回归选择重要特征
library(MASS)
# 再次尝试逐步回归
step_model <- stepAIC(full_model, direction = "both")
# 查看最终模型
summary(step_model)

# 计算多变量模型风险评分（使用插补后的数据）
riskscore_combined <- predict(step_model, newdata = data_imputed)
thrhold_combined <- median(riskscore_combined)
High_lowIndex_combined <- ifelse(riskscore_combined >= thrhold_combined, "High_risk", "Low_risk")

# 合并风险分组与生存数据
surv_data_with_risk <- data.frame(
  OS.time = data_imputed$OS.time,
  OS = data_imputed$OS,
  Risk_Group = High_lowIndex_combined
)

# 多变量模型生存分析
fit_combined <- survfit(Surv(OS.time, OS) ~ Risk_Group, data = surv_data_with_risk)

# 绘制多变量模型KM曲线（修正数据参数）
ggsurvplot(fit_combined, 
           data = surv_data_with_risk, 
           conf.int = TRUE, 
           pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, 
           xlab = "Time (days)", 
           legend = c(0.8, 0.9), 
           legend.title = "", 
           palette = c("#E7B800", "#2E9FDF"),
           title = "基于基因和临床特征的多变量模型生存分析(训练集)",  # 主标题
           ggtheme = theme_survminer() +          # 自定义主题
             theme(plot.title = element_text(
               size = 14,                        # 字体大小
               face = "bold",                    # 加粗
               hjust = 0.5,                      # 水平居中
               margin = margin(b = 10)))  )       # 下边距10pt)

# 在验证集上评估多变量模型（使用插补逻辑处理验证集缺失值）
# 假设验证集已进行相同的缺失值处理
# 先对验证集进行多重插补
imputed_data_test <- mice(cbind(sur_data.Test, gene_expression_test), m=5, maxit=50, method='pmm')
data_imputed_test <- complete(imputed_data_test, 1)

# 计算验证集风险评分
riskscore_combined_test <- predict(step_model, newdata = data_imputed_test)
High_lowIndex_combined_test <- ifelse(riskscore_combined_test >= thrhold_combined, "High_risk", "Low_risk")

# 合并验证集风险分组与生存数据
surv_data_test_with_risk <- data.frame(
  OS.time = data_imputed_test$OS.time,
  OS = data_imputed_test$OS,
  Risk_Group = High_lowIndex_combined_test
)

# 绘制验证集多变量模型KM曲线
fit_combined_test <- survfit(Surv(OS.time, OS) ~ Risk_Group, data = surv_data_test_with_risk)
ggsurvplot(fit_combined_test, 
           data = surv_data_test_with_risk, 
           conf.int = TRUE, 
           pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, 
           xlab = "Time (days)", 
           legend = c(0.8, 0.9), 
           legend.title = "", 
           palette = c("#E7B800", "#2E9FDF"),
           title = "基于基因和临床特征的多变量模型生存分析(验证集)",  # 主标题
           ggtheme = theme_survminer() +          # 自定义主题
             theme(plot.title = element_text(
               size = 14,                        # 字体大小
               face = "bold",                    # 加粗
               hjust = 0.5,                      # 水平居中
               margin = margin(b = 10)))  )     # 下边距10pt))

#模型性能比较（时间依赖 ROC）修改
# 比较基因模型和多变量模型的性能
library(timeROC)
time_points <- c(1, 3, 5) * 365  # 1年、3年和5年生存率
riskscore2
# 假设riskscore2是基因模型的风险评分
# 基因模型的时间依赖ROC（使用验证集插补后的数据）
roc_gene <- timeROC(T = data_imputed_test$OS.time, 
                    delta = data_imputed_test$OS, 
                    marker = riskscore2[,2],  # 假设riskscore2是二维矩阵
                    cause = 1, 
                    times = time_points,
                    ROC = TRUE)

# 多变量模型的时间依赖ROC
roc_combined <- timeROC(T = data_imputed_test$OS.time, 
                        delta = data_imputed_test$OS, 
                        marker = riskscore_combined_test,
                        cause = 1, 
                        times = time_points,
                        ROC = TRUE)

# 打印ROC曲线下面积
cat("基因模型AUC:\n")
print(roc_gene$AUC)
#    t=365    t=1095    t=1825 
#0.6476714 0.5688504 0.4384183 
cat("多变量模型AUC:\n")
print(roc_combined$AUC)
#t=365    t=1095    t=1825 
#0.7131725 0.6100808 0.5097586 

# 加载必要的包
library(ggplot2)
library(gridExtra)

# 整理基因模型数据
roc_gene_df <- data.frame(
  time = rep(c(365, 1095, 1825), each = length(roc_gene$FP[,1])),
  type = "基因模型",
  specificity = 1 - roc_gene$FP[,1],
  sensitivity = roc_gene$TP[,1],
  auc = rep(roc_gene$AUC, each = length(roc_gene$FP[,1]))
)

# 整理多变量模型数据
roc_combined_df <- data.frame(
  time = rep(c(365, 1095, 1825), each = length(roc_combined$FP[,1])),
  type = "多变量模型",
  specificity = 1 - roc_combined$FP[,1],
  sensitivity = roc_combined$TP[,1],
  auc = rep(roc_combined$AUC, each = length(roc_combined$FP[,1]))
)

# 合并数据
all_roc_df <- rbind(roc_gene_df, roc_combined_df)

# 创建时间标签映射
time_labels <- c("365" = "1-Year ROC", "1095" = "3-Year ROC", "1825" = "5-Year ROC")

# 创建三个时间点的ROC图
p1 <- ggplot(subset(all_roc_df, time == 365), aes(x = specificity, y = sensitivity, color = type)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "gray") +
  labs(title = time_labels["365"], x = "1-Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  annotate("text", x = 0.3, y = 0.2, 
           label = paste0("AUC:\n基因模型 = ", round(roc_gene$AUC[1], 3),
                          "\n多变量模型 = ", round(roc_combined$AUC[1], 3)),
           size = 3)

p2 <- ggplot(subset(all_roc_df, time == 1095), aes(x = specificity, y = sensitivity, color = type)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "gray") +
  labs(title = time_labels["1095"], x = "1-Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  annotate("text", x = 0.3, y = 0.2, 
           label = paste0("AUC:\n基因模型 = ", round(roc_gene$AUC[2], 3),
                          "\n多变量模型 = ", round(roc_combined$AUC[2], 3)),
           size = 3)

p3 <- ggplot(subset(all_roc_df, time == 1825), aes(x = specificity, y = sensitivity, color = type)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "gray") +
  labs(title = time_labels["1825"], x = "1-Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  annotate("text", x = 0.3, y = 0.2, 
           label = paste0("AUC:\n基因模型 = ", round(roc_gene$AUC[3], 3),
                          "\n多变量模型 = ", round(roc_combined$AUC[3], 3)),
           size = 3)

# 使用gridExtra包组合三个图
grid.arrange(p1, p2, p3, ncol = 3, 
             top = "Time-dependent ROC Comparison for Gene and Multivariable Models")


#####生存分析可视化（森林图和列线图）####################
ggforest(model = f22,#cox回归方程
         main = "senlintu",#标题
         cpositions = c(0.05, 0.15, 0.35),#前三列距离
         fontsize = 1, #字体大小
         refLabel = 'reference', #相对变量的数值标签，也可改为1
         noDigits = 3 )#保留HR值以及95%CI的小数位数

#列线图
##多因素cox回归
surage_chose<-sur_data.Train[,c(2,3,7,8,9,11)]
NewLabel<-cbind(surage_chose,riskS=riskscore33)

fit3 <- coxph(Surv(OS.time,OS) ~riskS+age+stage+gender,data=NewLabel)
# 提取stage列并创建新变量
stage_new <- NewLabel$stage

# 逐次替换为简化格式
stage_new <- gsub("Stage I", "I", stage_new)
stage_new <- gsub("Stage IA", "IA", stage_new)
stage_new <- gsub("Stage IB", "IB", stage_new)
stage_new <- gsub("Stage IIA", "IIA", stage_new)
stage_new <- gsub("Stage IIB", "IIB", stage_new)
stage_new <- gsub("Stage IIIA", "IIIA", stage_new)
stage_new <- gsub("Stage IIIB", "IIIB", stage_new)
stage_new <- gsub("Stage IV", "IV", stage_new)

# 将处理后的值重新赋值给原数据框
NewLabel$stage <- stage_new
#BiocManager::install("rms")
library(rms)
library(survival)
dd=datadist(NewLabel)
options(datadist='dd') #非常重要
fit55 <- psm(Surv(OS.time,OS) ~riskS+age+stage+gender,data =NewLabel, dist='lognormal') 
med <- Quantile(fit55) # 计算中位生存时间
surv <- Survival(fit55) # 构建生存概率函数
## 绘制COX回归中位生存时间的Nomogram图
nom <- nomogram(fit55, fun=function(x) med(lp=x),funlabel="Median Survival Time")
plot(nom)








################################################################################
# 免疫浸润分析与可视化
# 加载必要的包

library(ESTIMATE)
library(xCell)
library(limma)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(gridExtra)

###########ssGSEA分析################
library("GSVA")
geneSet <- read.csv("E:\\生信案例分析与实践\\专题五 免疫浸润分析\\CellReports.txt",header = F,sep = "\t")
# CellReports.txt 就是TILs.txt
geneSet=t(geneSet)
colnames(geneSet)=geneSet[1,]
geneSet<-geneSet[-1,]
f<-function(X)
{ temp<-which(X=="")
if (length(temp)>0)
{aa<-X[-temp]}else
{aa<-X}
aa}
#注意geneSet中有空值（“”）
cellmarker<-list()
for (i in 1:(dim(geneSet)[2]))
{
  temp1<-geneSet[,i]
  cellmarker[[i]]<-f(temp1)
}
names(cellmarker)<-colnames(geneSet)

# 创建 ssGSEA 参数对象
param <- ssgseaParam(
  exprData = as.matrix(dataExp.Train),
  geneSets = cellmarker,
)
# 确保基因集里的基因都能在表达矩阵中找到
cellmarker <- lapply(cellmarker, function(genes) {
  intersect(genes, rownames(dataExp.Train))
})

# 移除空的基因集
cellmarker <- cellmarker[lengths(cellmarker) > 0]
########运行 GSVA
data_immue_inflation <- gsva(param)
library(pheatmap)
# 检查样本数是否一致
length(High_lowIndex) == ncol(data_immue_inflation)

# 如果不一致，可能是样本顺序不同或部分样本被过滤
# 确保 High_lowIndex 的顺序与 data_immue_inflation 的列名一致
High_lowIndex <- High_lowIndex[colnames(data_immue_inflation)]
# 准备注释信息（以之前的高低风险组为例）
annotation_col <- data.frame(
  RiskGroup = factor(High_lowIndex, levels = c("Low_risk", "High_risk")),
  row.names = colnames(data_immue_inflation)
)

# 绘制热图
pheatmap(
  mat = data_immue_inflation,
  scale = "row",          # 按行标准化
  annotation_col = annotation_col,
  show_colnames = FALSE,  # 样本多时不显示列名
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Immune Infiltration Patterns (ssGSEA)"
)
# 5. 准备箱线图数据 --------------------------------------------------
plot_data <- as.data.frame(t(data_immue_inflation))
plot_data$Sample <- rownames(plot_data)
plot_data$RiskGroup <- annotation_col$RiskGroup

data_melt <- melt(plot_data, 
                  id.vars = c("Sample", "RiskGroup"),
                  variable.name = "CellType",
                  value.name = "Score")
# ------------------------- 箱线图绘制 -------------------------
# 1. 转换数据格式
plot_data <- as.data.frame(t(data_immue_inflation))
plot_data$Sample <- rownames(plot_data)
plot_data$RiskGroup <- annotation_col$RiskGroup

# 2. 使用 melt 转换长格式
data_melt <- melt(
  plot_data,
  id.vars = c("Sample", "RiskGroup"),
  variable.name = "CellType",
  value.name = "Score"
)

# 3. 绘制箱线图并添加统计检验
ggplot(data_melt, aes(x = CellType, y = Score, fill = RiskGroup)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, position = position_dodge(0.8)) +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),
    method = "wilcox.test",
    label.x = 1.5,
    size = 3
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # 调整x轴标签角度
  labs(x = "Cell Type", y = "ssGSEA Score", fill = "Risk Group")

###########xCell分析######################
#xCell安装
BiocManager::install("devtools")
devtools::install_github('dviraran/xCell')#最好用R4.2版本
library(xCell)
tscores <- xCellAnalysis(dataExp.Train,rnaseq = T,scale = TRUE, alpha = 0.5, 
                          parallel.sz = 4, parallel.type = "SOCK")
#二代测序数据rnaseq=T，芯片 rnaseq=F
#xCell biomarker gene:64种细胞
makers<-xCell.data$signature
library(pheatmap)

# 准备分组信息
annotation_col <- data.frame(
  RiskGroup = factor(High_lowIndex[colnames(tscores)], 
                     levels = c("Low_risk", "High_risk")),
  row.names = colnames(tscores)
)

# 绘制热图（显示前30个变异最大的细胞类型）
cell_var <- apply(tscores, 1, var)
top_cells <- names(sort(cell_var, decreasing = TRUE))[1:30]

pheatmap(
  mat = tscores[top_cells, ],
  scale = "row",
  annotation_col = annotation_col,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  show_colnames = FALSE,
  main = "xCell Immune Infiltration (Top 30 Variable Cell Types)",
  fontsize_row = 8,
  clustering_method = "complete"
)
library(ggplot2)
library(ggpubr)
library(reshape2)

# 准备绘图数据
plot_data <- as.data.frame(t(tscores[top_cells, ]))
# 1. 检查并匹配样本 -------------------------------------------------
# 确保 xCell 结果和分组信息的样本完全一致
common_samples <- intersect(colnames(tscores), names(High_lowIndex))

# 过滤数据
tscores <- tscores[, common_samples]
High_lowIndex <- High_lowIndex[common_samples]

# 2. 重新创建注释数据框 ---------------------------------------------
annotation_col <- data.frame(
  RiskGroup = factor(High_lowIndex, levels = c("Low_risk", "High_risk")),
  row.names = common_samples
)

# 3. 准备绘图数据 --------------------------------------------------
plot_data <- as.data.frame(t(tscores[top_cells, ]))
stopifnot(nrow(plot_data) == length(annotation_col$RiskGroup)) # 验证行数

plot_data$RiskGroup <- annotation_col$RiskGroup  # 现在长度一致

# 4. 继续后续绘图流程 ----------------------------------------------
data_melt <- melt(plot_data, id.vars = "RiskGroup")
plot_data$RiskGroup <- annotation_col$RiskGroup
data_melt <- melt(plot_data, id.vars = "RiskGroup")

# 绘制箱线图
ggplot(data_melt, aes(x = variable, y = value, fill = RiskGroup)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), 
              size = 0.8, alpha = 0.4) +
  stat_compare_means(
    aes(group = RiskGroup),
    method = "wilcox.test",
    label = "p.signif",
    label.y = max(data_melt$value) * 1.1,
    size = 4
  ) +
  scale_fill_manual(values = c("Low_risk" = "blue", "High_risk" = "red")) +
  labs(
    title = "xCell Immune Cell Infiltration Scores",
    x = "Cell Type",
    y = "xCell Score",
    fill = "Risk Group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(min(data_melt$value), max(data_melt$value) * 1.15))







#########Cibersort分析###################

# BiocManager::install("devtools")
# devtools::install_github("Moonerss/CIBERSORT")## github装包
##改造数据集

Exp2<-dataExp.Train
temp<-row.names(dataExp.Train)
Exp2_1<-cbind(GeneSymbol=temp,Exp2)
write.table(Exp2_1,file = "tcga_mianyi_expr_data.txt", col.names = TRUE,sep = "\t ",row.names = F)
source("E:\\生信案例分析与实践\\专题四 生存分析\\cibersort.R")
# 尝试大写函数名
results <- CIBERSORT(
  sig_matrix = system.file("extdata", "LM22.txt", package = "CIBERSORT"),
  mixture_file = "tcga_mianyi_expr_data.txt",
  perm = 100,
  QN = F
)


library(pheatmap)

# 假设 results 是 CIBERSORT 的输出（22 种免疫细胞比例）
# 提取免疫细胞比例（去掉 P-value 和 Correlation 列）
immune_matrix <- results[, 1:22]

# 按样本分组（假设有分组信息 risk_groups）
annotation_col <- data.frame(
  RiskGroup = factor(High_lowIndex, levels = c("Low_risk", "High_risk")),
  row.names = colnames(data_immue_inflation)
)

# 绘制热图
pheatmap(
  t(immune_matrix),  # 转置矩阵（行为细胞，列为样本）
  scale = "row",     # 按行标准化（Z-score）
  annotation_col = annotation_col,  # 添加分组注释
  color = colorRampPalette(c("blue", "white", "red"))(100),  # 颜色映射
  show_colnames = FALSE,  # 不显示样本名（避免重叠）
  main = "Immune Cell Infiltration (CIBERSORT)"
)
# 假设 results 是 CIBERSORT 的输出
head(results)

# 检查样本是否匹配
stopifnot(rownames(results) == colnames(dataExp.Train))

# 合并 CIBERSORT 结果与风险分组
cibersort_df <- as.data.frame(results[, 1:22])  # 提取22种免疫细胞比例
cibersort_df$RiskGroup <- High_lowIndex  # 添加风险分组列
library(tidyr)
library(ggplot2)

# 转换为长格式（适合 ggplot2）
cibersort_long <- pivot_longer(
  cibersort_df,
  cols = -RiskGroup,
  names_to = "CellType",
  values_to = "Proportion"
)
ggplot(cibersort_long, aes(x = CellType, y = Proportion, fill = RiskGroup)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  scale_fill_manual(values = c("Low_risk" = "blue", "High_risk" = "red")) +
  labs(
    title = "cibersort Immune Cell Infiltration by Risk Group",
    x = "Cell Type",
    y = "Proportion",
    fill = "Risk Group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
