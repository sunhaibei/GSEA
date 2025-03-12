# 安装和加载所需的包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")  # 小鼠的基因注释数据库
BiocManager::install("enrichplot")    # 用于可视化

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# 安装并加载 readxl 包
install.packages("readxl")
library(readxl)


# 打开 XLSX 文件
degs_all <- read_excel("degs_all.xlsx")

# 检查数据
head(degs_all)

# 使用 bitr 函数进行 ID 转换
id_mapping <- bitr(degs_all$id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# 检查转换结果
head(id_mapping)

# 合并 ID 转换结果到原始数据
degs_all <- merge(degs_all, id_mapping, by.x = "id", by.y = "ENSEMBL")

# 检查合并后的数据
head(degs_all)

# 提取 log2fc 和 Entrez ID
geneList <- degs_all$log2fc
names(geneList) <- degs_all$ENTREZID

# 按 log2fc 降序排序
geneList <- sort(geneList, decreasing = TRUE)

# 检查基因列表
head(geneList)

# 使用 GO 数据库进行 GSEA 分析
gsea_result <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",  # 使用转换后的 Entrez ID
  ont          = "BP",        # 选择生物学过程（BP）
  minGSSize    = 10,          # 最小基因集大小
  maxGSSize    = 500,         # 最大基因集大小
  pvalueCutoff = 0.05,        # p值阈值
  verbose      = FALSE
)

# 查看结果
head(gsea_result)

# 绘制前 10 个显著富集的基因集
dotplot(gsea_result, showCategory = 10)

# 绘制第一个显著富集的基因集的 GSEA 图
gseaplot(gsea_result, geneSetID = 1)

# 使用 KEGG 数据库进行 GSEA 分析
gsea_result_kegg <- gseKEGG(
  geneList     = geneList,
  organism     = "mmu",       # 小鼠的 KEGG 代码
  pvalueCutoff = 0.05
)

# 查看结果
head(gsea_result_kegg)

# 可视化 KEGG 结果
dotplot(gsea_result_kegg, showCategory = 10)
gseaplot(gsea_result_kegg, geneSetID = 1)
