setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_08",recursive = T,showWarnings = F)



library(org.Hs.eg.db)
library(clusterProfiler)




# gtf_grch38_extract.csv
#################################################################################################
setwd("F:/code/biotrainee2017/result_05")

temp = read.csv("gtf_grch38_extract.csv",header = T,encoding = "UTF-8")
gtf_grch38_extract = temp
temp = gtf_grch38_extract[1:1000,]
setwd("F:/code/biotrainee2017/data")
#################################################################################################


#################################################################################################
temp = gtf_grch38_extract[,2:3]
temp = temp[!duplicated(temp),]

keytypes(org.Hs.eg.db)
gene_id = bitr(temp[,1],
               fromType="ENSEMBL",
               toType=keytypes(org.Hs.eg.db)[6],
               OrgDb="org.Hs.eg.db")
gene_name = bitr(gtf_grch38_extract[,2],
                 fromType="ENSEMBL",
                 toType=keytypes(org.Hs.eg.db)[10],
                 OrgDb="org.Hs.eg.db")
# 26585
nrow(gene_id)
# 40609
length(which(temp[,2]!=""))
# 手动转换和用包转换还是存在偏差的，可能是数据来源与版本更新存在差别
#################################################################################################

