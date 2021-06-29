setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_03",recursive = T,showWarnings = F)


library(stringr)



# Homo_sapiens.GRCh38.104.gtf.gz
# ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
#################################################################################################
# 读取索引文件
data_1 = read.table(gzfile("Homo_sapiens.GRCh38.104.gtf.gz"),
                    sep = "\t",
                    stringsAsFactors=F,
                    encoding = "UTF-8",
                    comment.char = "#")
head(data_1)
temp = data_1[1:100,]
#################################################################################################


# 对Homo_sapiens.GRCh38.104.gtf.gz文件进行统计
#################################################################################################
# 3146132 全是 gene_id
table(substr(data_1[,9],1,7))

ID = which(data_1[,3]=="gene")
# 60664 个gene
length(ID)
# 60664 个ENSG号
length(unique(substr(data_1[,9],9,23)))
# 60664 说明一个gene对应一个ENSG号
length(unique(substr(data_1[ID,9],9,23)))

# 提取 染色体编号、ENSG号
temp = cbind(data_1[,1],substr(data_1[,9],9,23))
# 去重
temp = temp[!duplicated(temp),]
# 全为TRUE,说明同样的ENSG号只在一条染色体上没有问题
table(!duplicated(temp[,2]))


ID = which(data_1[,3]=="gene")
# 按照 "; " 分隔
temp = str_split(data_1[ID,9],"; ",simplify = T)
# 提取 gene_name ,没有 gene_named 的显示空白
temp[,3] = ifelse(substr(temp[,3],1,9)=="gene_name",substr(temp[,3],11,99999),"")

# ensg_summary
ensg_summary = data.frame(chr=data_1[ID,1],
                          gene_id=substr(data_1[ID,9],9,23),
                          gene_name=temp[,3],
                          start=data_1[ID,4],
                          end=data_1[ID,5],
                          transcript=0,
                          exon=0,
                          cds=0
                          )

# 提取ENSG id 和id对应的类型做统计
temp = cbind(substr(data_1[,9],9,23),data_1[,3])
temp = as.data.frame(table(paste(temp[,1],temp[,2],sep = ";")))
temp = cbind(substr(temp[,1],1,15),substr(temp[,1],17,99999),temp[,2])
# transcript
temp1 = temp[temp[,2]=="transcript",]
ID = match(temp1[,1],ensg_summary[,2])
ensg_summary[ID,6] = as.numeric(temp1[,3])
# exon
temp1 = temp[temp[,2]=="exon",]
ID = match(temp1[,1],ensg_summary[,2])
ensg_summary[ID,7] = as.numeric(temp1[,3])
# CDS
temp1 = temp[temp[,2]=="CDS",]
ID = match(temp1[,1],ensg_summary[,2])
ensg_summary[ID,8] = as.numeric(temp1[,3])
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_03")

write.csv(ensg_summary,"ensg_summary.csv",row.names = F)
temp = read.csv("ensg_summary.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==ensg_summary)
setwd("F:/code/biotrainee2017/data")
#################################################################################################


# 对ensg_summary的每条染色体进行统计
#################################################################################################
chr_summary = matrix(NA,0,6)
for (i in unique(ensg_summary[,1])) {
  ID = which(ensg_summary[,1]==i)
  # 提取指定一条染色体的数据
  temp = ensg_summary[ID,]
  temp = c(i,
           nrow(temp),
           length(which(temp[,3]!="")),
           sum(temp[,6]),
           sum(temp[,7]),
           sum(temp[,8])
           )
  chr_summary = rbind(chr_summary,temp)
}

# 提取指定染色体编号
ID = match(c(1:22,"X","Y","MT"),chr_summary[,1])
chr_summary = data.frame(chr=chr_summary[ID,1],
                         ensg_number=as.numeric(chr_summary[ID,2]),
                         gene_number=as.numeric(chr_summary[ID,3]),
                         transcript_number=as.numeric(chr_summary[ID,4]),
                         exon_number=as.numeric(chr_summary[ID,5]),
                         cds_number=as.numeric(chr_summary[ID,6])
                         )
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_03")

write.csv(chr_summary,"chr_summary.csv",row.names = F)
temp = read.csv("chr_summary.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==chr_summary)
setwd("F:/code/biotrainee2017/data")
#################################################################################################



