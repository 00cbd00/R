setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_07",recursive = T,showWarnings = F)



library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)



# reads_count.csv
#################################################################################################
setwd("F:/code/biotrainee2017/result_04")

temp = read.csv("reads_count.csv",header = T,encoding = "UTF-8")
reads_count = temp
setwd("F:/code/biotrainee2017/data")
#################################################################################################
# kegg_extract.csv
#################################################################################################
setwd("F:/code/biotrainee2017/result_06")

temp = read.csv("kegg_extract.csv",header = T,encoding = "UTF-8")
kegg_extract = temp
setwd("F:/code/biotrainee2017/data")
#################################################################################################


# DEseq
#################################################################################################
# 前期数据准备
rownames(reads_count) = reads_count[,1]
reads_count = reads_count[,-1]
# 奇数为before,偶数为after
group_list = factor(rep(c("before","after"),c(3,3)))
colData = data.frame(row.names=colnames(reads_count), group_list=group_list)

# 开始执行 DEseq
dds = DESeqDataSetFromMatrix(countData = reads_count,
                             colData = colData,
                             design = ~ group_list)
dds2 = DESeq(dds)
resultsNames(dds2)
result_DEseq = results(dds2, contrast=c("group_list","after","before"))
result_DEseq = result_DEseq[order(result_DEseq$padj),]
result_DEseq = as.data.frame(result_DEseq)
result_DEseq = result_DEseq[!is.na(result_DEseq[,5]),]
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_07")

write.csv(result_DEseq,"result_DEseq.csv",row.names = F)
temp = read.csv("result_DEseq.csv",header = T,encoding = "UTF-8")
# FALSE   TRUE 
# 185232  11070 小数点位偏差，影响不大
table(temp==result_DEseq)
# result_DEseq = temp
setwd("F:/code/biotrainee2017/data")
#################################################################################################


#################################################################################################
# 提取差异表达基因
ID = which(result_DEseq[,5]<=0.05 & abs(result_DEseq[,2])>=1.5)
dif_expression = result_DEseq[ID,]
# 将 ensg id 转换成 gene id
keytypes(org.Hs.eg.db)
gene_id = bitr(rownames(dif_expression),
               fromType="ENSEMBL",
               toType=keytypes(org.Hs.eg.db)[6],
               OrgDb="org.Hs.eg.db")

# 用包做富集分析 提取结果
enrich_kegg = enrichKEGG(gen = gene_id[,2],organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1)
head(enrich_kegg)
result_kegg = as.data.frame(enrich_kegg)
# N 总数:8108   n 及格数:(136)   M 总抽取数:32   m 总抽取及格数:(6)
result_kegg[1,]
# 总共得到 100 个通路
nrow(result_kegg)


# N 总数:8104 与kegg提取的结果偏差不大
N = length(unique(kegg_extract[,1]))

# M 总抽取数:32 与kegg提取的结果一致
M = length(which(gene_id[,2]%in%kegg_extract[,1]))

# m 总抽取及格数
m = kegg_extract[which(kegg_extract[,1]%in%gene_id[,2]),]
m = as.data.frame(table(m[,2]))
# 总共得到 100 个通路,与kegg提取的结果一致
nrow(m)

# n 总及格数
n = as.data.frame(table(kegg_extract[,2]))
n = n[which(n[,1]%in%m[,1]),]
table(as.matrix(m[,1])==as.matrix(n[,1]))

# hsa05200
m[-which(m[,1]%in%result_kegg[,1]),]
# hsa01230 有偏差
result_kegg[-which(result_kegg[,1]%in%m[,1]),]

# 结果与kegg富集结果接近
phyper(m[,2]-1,n[,2],N-n[,2],M,lower.tail = F)
# 结果与kegg富集结果一致
phyper(m[,2]-1,n[,2],8108-n[,2],M,lower.tail = F)


#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_07")

write.csv(result_kegg,"result_kegg.csv",row.names = F)
temp = read.csv("result_kegg.csv",header = T,encoding = "UTF-8")
# FALSE  TRUE 
# 292    608  小数点位偏差，影响不大
table(temp==result_kegg)
# result_kegg = temp
setwd("F:/code/biotrainee2017/data")
#################################################################################################
