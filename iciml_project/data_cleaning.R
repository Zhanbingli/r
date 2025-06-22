# 加载核心包
library(GEOquery)

# 下载 GSE78220 数据集。GSEMatrix = TRUE 确保下载的是处理好的矩阵形式
# 这会返回一个 ExpressionSet 对象的列表，通常我们只关心第一个元素
gse <- getGEO("GSE78220", GSEMatrix = TRUE)[[1]]

# 从 ExpressionSet 对象中提取表型（样本）信息
# pData 是一个专门用于提取这部分信息的函数
pheno_raw <- pData(gse)

# 从 ExpressionSet 对象中提取表达矩阵
# exprs 是一个专门用于提取表达矩阵的函数
expr_raw <- exprs(gse)


# 加载相关包
library(GEOquery)
library(readxl) # 用于读取 .xlsx 文件
library(R.utils) # 用于解压 .gz 文件

# 只下载补充文件到指定目录 data/GSE78220/
getGEOSuppFiles("GSE78220", baseDir = "data", makeDirectory = TRUE)

# 在下载的文件夹中找到 FPKM 文件
# full.names = TRUE 返回完整路径
fpkm_file_gz <- list.files("data/GSE78220", pattern = "FPKM.*gz$", full.names = TRUE)

# 解压文件 (fread 也可以直接读取 .gz 文件，但解压出来更方便)
gunzip(fpkm_file_gz, overwrite = FALSE)
fpkm_file_xlsx <- sub("\\.gz$", "", fpkm_file_gz) # 获取解压后的文件名

# 使用 readxl 包读取 Excel 文件
# 这个 Excel 文件包含了我们需要的 FPKM 表达矩阵
expr_raw_from_file <- read_xlsx(fpkm_file_xlsx)

# 加载 dplyr 包用于数据处理
library(dplyr)

# 假设我们使用方法一获取的 pheno_raw
# 查看所有列名，找到包含响应信息的列
# names(pheno_raw) 或 colnames(pheno_raw)
# 经过检查，信息在 'characteristics_ch1.3' 或 'anti-pd-1 response:ch1' 列中
# (不同版本的 GEOquery 或 GEO 数据可能导致列名稍有不同)

pheno_clean <- pheno_raw %>%
  mutate(
    # 提取响应信息文本，并转为小写以统一格式
    # `sub()` 函数用于删除 "response: " 前缀
    response_raw = tolower(sub("response: ", "", characteristics_ch1.3)),
    
    # 根据文本创建二元响应变量 (0 或 1)
    # 这是机器学习或统计建模的标准做法
    # "responder" 被编码为 1, "nonresponder" 被编码为 0
    response_bin = ifelse(response_raw == "responder", 1, 0)
  ) %>%
  # 为了方便后续匹配，将样本的 GEO ID (例如 GSM2079940) 设置为行名
  # 这通常存储在 'geo_accession' 列
  tibble::column_to_rownames(var = "geo_accession")

# 检查清洗结果
# table() 函数可以统计每个分类的数量
table(pheno_clean$response_raw, pheno_clean$response_bin)

讲解:
  
  mutate(): dplyr 的核心函数，用于创建或修改列。

tolower(sub(...)): 这是一个组合操作。sub("response: ", "", ...) 先把多余的文字 "response: " 去掉，然后 tolower() 把结果（如 "Responder"）全部转为小写（"responder"），这样可以避免因大小写问题导致的匹配错误。

ifelse(condition, value_if_true, value_if_false): 这是一个条件判断。如果 response_raw 是 "responder"，那么 response_bin 列的值就是 1，否则就是 0。这种 0/1 编码对于后续的统计分析至关重要。

tibble::column_to_rownames(): 一个方便的函数，它将某一列（这里是 geo_accession）的内容设置为数据框的行名，并删除原始列。将样本ID设为行名是匹配表达矩阵和表型数据的关键步骤。

第三步：清洗和整理表达矩阵 (Expression Matrix)
从文件（方法二）或 ExpressionSet（方法一）中得到的原始表达矩阵也需要处理，主要是确保其行名是基因，列名是与表型数据一致的样本ID。

# 假设我们使用从补充文件读取的 expr_raw_from_file
# 它是一个数据框，第一列是基因名，其余列是样本

# 检查数据维度和前几行
# dim(expr_raw_from_file)
# head(expr_raw_from_file[, 1:5])

# 将第一列 'Gene' 设置为行名，并移除该列，得到一个纯数值矩阵
expr_matrix <- expr_raw_from_file %>%
  tibble::column_to_rownames(var = "Gene")

# 清理列名：样本列名可能是 "Ptid1.FPKM" 的形式，我们需要把它变成 "Ptid1"
# 但在这个数据集中，列名已经是样本ID了，所以这步可能不需要
# 如果需要，可以使用： colnames(expr_matrix) <- sub("\\..*", "", colnames(expr_matrix))

# ！！！关键对齐步骤！！！
# 检查表达矩阵的列名 (样本ID) 和表型数据的行名 (样本ID) 是否完全一致
# 这是确保数据可以正确分析的前提
all(colnames(expr_matrix) %in% rownames(pheno_clean))
# > TRUE  (如果结果是 TRUE，说明匹配成功！)

# 如果有不匹配的，可以用 setdiff 查看是哪些样本出了问题
# setdiff(colnames(expr_matrix), rownames(pheno_clean))

# 为了安全，只保留那些在两边都存在的样本
common_samples <- intersect(colnames(expr_matrix), rownames(pheno_clean))
expr_final <- expr_matrix[, common_samples]
pheno_final <- pheno_clean[common_samples, ]

# 再次确认
all(colnames(expr_final) == rownames(pheno_final))
# > TRUE


# 创建一个目录来存放处理好的数据
if (!dir.exists("data_clean")) {
  dir.create("data_clean")
}

# 保存清洗后的表型数据
write.csv(pheno_final, "data_clean/GSE78220_pheno_clean.csv", row.names = TRUE)

# 保存清洗后的表达矩阵
# quote = FALSE 可以避免在列名和行名上添加不必要的引号，让文件更干净
# row.names = TRUE 会把行名（基因名）写在文件的第一列
write.csv(expr_final, "data_clean/GSE78220_expr_clean.csv", quote = FALSE, row.names = TRUE)

# (可选) 清理内存
# 在处理大数据集时，及时清理不再需要的变量是个好习惯
rm(gse, expr_raw, pheno_raw, expr_matrix, pheno_clean)
gc() # 触发垃圾回收