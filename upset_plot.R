# 加载必要包
if (!require("UpSetR")) install.packages("UpSetR")
if (!require("dplyr")) install.packages("dplyr")
library(UpSetR)
library(dplyr)
#install.packages("RColorBrewer")
# ===================== 1. 读取数据（替换为你的数据读取方式） =====================
# 方式1：直接读取合并后的Excel文件（之前生成的merged_gene_id_matrix.xlsx）
# 注意：若已在环境中，可跳过此步
if (!require("readxl")) install.packages("readxl")
library(readxl)
merged_df <- read_excel("combined_gene_ids.xlsx")  # 替换为你的文件路径


# ===================== 2. 转换为UpSetR所需的集合列表格式 =====================
# 将每一列（比较组）转换为“集合”，元素为该组的gene_id
gene_sets <- lapply(merged_df, function(col) {
  # 过滤NA值（仅保留该组存在的gene_id）
  na.omit(col) %>% as.character()
})

# 命名集合（列名即比较组名称）
names(gene_sets) <- colnames(merged_df)

n <-15
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n)


# ===================== 3. 绘制UpSet图（匹配例图样式） =====================
p<- upset(
  fromList(gene_sets),  # 将列表转换为UpSetR格式
  sets = colnames(merged_df),  # 手动指定全部15个组
  # 图形布局与样式
  mb.ratio = c(0.6, 0.4),  # 上下图比例（上方柱状图:下方矩阵+条形图）
  order.by = "freq",        # 按交集频率排序（与例图一致）
  decreasing = TRUE,        # 降序排列
  keep.order = F,        # 保留集合的原始顺序
  # 关键添加：显示所有数值
  show.numbers = "yes",     # 在上方条形图显示交集大小
  set_size.show = TRUE,     # 在左侧显示集合大小
  set_size.scale_max = max(sapply(gene_sets, length)) * 1.1,  # 自动缩放
  set_size.numbers_size = 8,  # 集合大小数字的字体大小
  set_size.angles = 0,      # 集合大小数字的角度（水平）
  # 元素大小调整
  point.size = 2,           # 矩阵中点的大小
  line.size = 1,            # 矩阵中线的粗细
  #
  sets.bar.color = my_colors,
  # 文本样式
  text.scale = c(1.2, 1.2, 1, 1, 1.5, 0.8),  # 文本缩放（对应：轴标签、刻度、集合名、数值等）
  number.angles = 0,        # 柱状图数值的角度（0=水平）
  # 轴标签（可选）
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Total gene number",
  nintersects=NA
  
)
p
