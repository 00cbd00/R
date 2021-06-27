http://www.biotrainee.com/thread-3-1-1.html


setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_01",recursive = T,showWarnings = F)


# overlap_clear:清除每行区域的重叠位点
# 要求输入的数据框为 染色体编号、起始、终止位点 3列数据
#################################################################################################
overlap_clear = function(data){
  # 记录运行起始时间
  time_1 = Sys.time()
  
  # 数字化数据的起始、终止位点
  data = data.frame(data[,1],as.numeric(data[,2]),as.numeric(data[,3]))

  result = matrix(NA,0,3)
  for (i in unique(data[,1])) {
    # 对每一个染色体编号分别处理
    ID = which(data[,1] == i)
    
    # 该染色体编号 存在多行
    if(length(ID) > 1){
      # 提取起始、终止位点
      temp = data[(data[,1] == i),2:3]
      # 按照起始位置、终止位置的升序排序
      temp = temp[order(temp[,1],temp[,2]),]
      # 开始清除重叠位点
      condition = T
      while (condition) {
        
        # 存在重叠的行(上一行的终止位置 大于 下一行的起始位置,即 存在重叠)
        ID = which((temp[1:(nrow(temp)-1),2] - temp[2:nrow(temp),1]) >= 0)
        
        if(length(ID) != 0){
          
          # 上一行终止位置 覆盖 下一行终止位置 的行
          ID_1 = which((temp[ID,2]-temp[ID+1,2]) >= 0)
          
          if(length(ID_1) != 0){
            # 下一行起始、终止位置均被上一行覆盖，可直接舍去
            temp = temp[-(ID[ID_1]+1),]
          }
          else{
            # 将上一行终止位置小于下一行终止位置的行 的终止位置 延长至 与下一行终止位置一致
            temp[ID,2] = temp[ID+1,2]
          }
        }
        else{
          # 已经不存在重叠的行，循环结束
          condition = F
        }
      }
      temp = cbind(i,temp)
      result = rbind(result,temp)
    }
    # 该染色体编号 仅存在1行
    else{
      result = rbind(result,c(data[ID,1],data[ID,2],data[ID,3]))
    }
  }
  result = data.frame(chr = result[,1],
                      start = as.numeric(result[,2]),
                      end = as.numeric(result[,3])
                      )
  # 显示运行总时间
  print(difftime(Sys.time(), time_1, units = 'secs'))
  return(result)
}
#################################################################################################


# count_length:统计 每条染色体对应序列总体长度
# 要求输入的数据框为 染色体编号、起始、终止位点 3列数据
#################################################################################################
count_length = function(data){
  data = data.frame(data[,1],as.numeric(data[,2]),as.numeric(data[,3]))
  result = matrix(NA,0,4)
  for (i in unique(data[,1])) {
    # 对每一个染色编号分别处理
    ID = which(data[,1] == i)
    if(length(ID) > 1){
      temp = data[(data[,1] == i),2:3]
      # 统计总体长度
      temp_length = sum(temp[,2]-temp[,1]) + nrow(temp)
      # 增加一行总结
      result = rbind(result,c(i,min(temp),max(temp),temp_length))
    }
    else{
      result = rbind(result,c(i,data[ID,2],data[ID,3],(data[ID,3]-data[ID,2]+1)))
    }
  }
  result = data.frame(chr = result[,1],
                      start = as.numeric(result[,2]),
                      end = as.numeric(result[,3]),
                      length = as.numeric(result[,4]))
  print(sum(result[,4]))
  result = rbind(result,c(nrow(result),min(result[,2]),max(result[,3]),sum(result[,4])))
  rownames(result)[nrow(result)] = "summary"
  return(result)
}
#################################################################################################


# Homo_sapiens.GRCh38.104.gtf.gz
# ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
#################################################################################################
# 预读100行
temp = read.table(gzfile("Homo_sapiens.GRCh38.104.gtf.gz"),
                  sep = "\t",
                  stringsAsFactors=F,
                  nrows=100,
                  encoding = "UTF-8",
                  comment.char = "#")
# 读取索引文件
data_1 = read.table(gzfile("Homo_sapiens.GRCh38.104.gtf.gz"),
                  sep = "\t",
                  stringsAsFactors=F,
                  encoding = "UTF-8",
                  comment.char = "#")
head(data_1)



# 提取 exon 起始、终止位置
ID = which(data_1[,3]=="exon")
data_exon_1 = data_1[ID,c(1,4,5)]
# 先剔除起始位置、终止位置完全一致的行
data_exon_2 = data_exon_1[!duplicated(data_exon_1),]
temp = data_exon_2[1:100,]

# Time difference of 8.632494 secs
overlap_clear_exon_1 = overlap_clear(data_exon_2)
# 150586755
count_length_exon_1 = count_length(overlap_clear_exon_1)



# 提取 cds 起始、终止位置
ID = which(data_1[,3]=="CDS")
data_cds_1 = data_1[ID,c(1,4,5)]
# 先剔除起始位置、终止位置完全一致的行
data_cds_2 = data_cds_1[!duplicated(data_cds_1),]
temp = data_cds_2[1:100,]

# Time difference of 3.614206 secs
overlap_clear_cds_1 = overlap_clear(data_cds_2)
# 35758443
count_length_cds_1 = count_length(overlap_clear_cds_1)
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_01")

write.csv(overlap_clear_exon_1,"overlap_clear_exon_1_Homo_sapiens.GRCh38.104.gtf.gz.csv",
          row.names = F)
temp = read.csv("overlap_clear_exon_1_Homo_sapiens.GRCh38.104.gtf.gz.csv",
                header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp == overlap_clear_exon_1)

write.csv(count_length_exon_1,"count_length_exon_1_Homo_sapiens.GRCh38.104.gtf.gz.csv",
          row.names = F)
temp = read.csv("count_length_exon_1_Homo_sapiens.GRCh38.104.gtf.gz.csv",
                header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp == count_length_exon_1)

write.csv(overlap_clear_cds_1,"overlap_clear_cds_1_Homo_sapiens.GRCh38.104.gtf.gz.csv",
          row.names = F)
temp = read.csv("overlap_clear_cds_1_Homo_sapiens.GRCh38.104.gtf.gz.csv",
                header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp == overlap_clear_cds_1)

write.csv(count_length_cds_1,"count_length_cds_1_Homo_sapiens.GRCh38.104.gtf.gz.csv",
          row.names = F)
temp = read.csv("count_length_cds_1_Homo_sapiens.GRCh38.104.gtf.gz.csv",
                header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp == count_length_cds_1)

setwd("F:/code/biotrainee2017/data")
#################################################################################################



# TxDb.Hsapiens.UCSC.hg19.knownGene
#################################################################################################
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
data_2 = TxDb.Hsapiens.UCSC.hg19.knownGene
data_2



data_exon_1 = exons(data_2)
data_exon_1 = as.data.frame(data_exon_1)
head(data_exon_1)
data_exon_1 = data_exon_1[,1:3]
# 先剔除起始位置、终止位置完全一致的行
data_exon_2 = data_exon_1[!duplicated(data_exon_1),]

# Time difference of 5.555318 secs
overlap_clear_exon_2 = overlap_clear(data_exon_2)
# 88341206
count_length_exon_2 = count_length(overlap_clear_exon_2)



data_cds_1 = cds(data_2)
data_cds_1 = as.data.frame(data_cds_1)
head(data_cds_1)
data_cds_1 = data_cds_1[,1:3]
# 先剔除起始位置、终止位置完全一致的行
data_cds_2 = data_cds_1[!duplicated(data_cds_1),]

# Time difference of 2.349134 secs
overlap_clear_cds_2 = overlap_clear(data_cds_2)
# 36747178
count_length_cds_2 = count_length(overlap_clear_cds_2)
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_01")

write.csv(overlap_clear_exon_2,"overlap_clear_exon_2_TxDb.Hsapiens.UCSC.hg19.knownGene.csv",
          row.names = F)
temp = read.csv("overlap_clear_exon_2_TxDb.Hsapiens.UCSC.hg19.knownGene.csv",
                header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp == overlap_clear_exon_2)

write.csv(count_length_exon_2,"count_length_exon_2_TxDb.Hsapiens.UCSC.hg19.knownGene.csv",
          row.names = F)
temp = read.csv("count_length_exon_2_TxDb.Hsapiens.UCSC.hg19.knownGene.csv",
                header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp == count_length_exon_2)

write.csv(overlap_clear_cds_2,"overlap_clear_cds_2_TxDb.Hsapiens.UCSC.hg19.knownGene.csv",
          row.names = F)
temp = read.csv("overlap_clear_cds_2_TxDb.Hsapiens.UCSC.hg19.knownGene.csv",
                header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp == overlap_clear_cds_2)

write.csv(count_length_cds_2,"count_length_cds_2_TxDb.Hsapiens.UCSC.hg19.knownGene.csv",
          row.names = F)
temp = read.csv("count_length_cds_2_TxDb.Hsapiens.UCSC.hg19.knownGene.csv",
                header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp == count_length_cds_2)

setwd("F:/code/biotrainee2017/data")
#################################################################################################



# CCDS.current.txt
# ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/
#################################################################################################
data_3 = read.table("CCDS.current.txt",
                  sep = "\t",
                  stringsAsFactors=F,
                  encoding = "UTF-8",
                  comment.char = "",
                  header = T)

# 提取染色体编号 和 编码的起始、终止位点
data_cds_1 = data_3[,c(1,10)]
# 剔除没有编码位点的行
data_cds_1 = data_cds_1[(data_cds_1[,2] != "-"),]
# 剔除编码位点的大括号
data_cds_1[,2] = gsub("[][]","",data_cds_1[,2])

data_cds_2 = matrix(NA,0,3)
for (i in unique(data_cds_1[,1])) {
  # 对每一个染色编号分别处理
  temp = data_cds_1[(data_cds_1[,1] == i),2]
  # 拆分位点
  temp = unlist(strsplit(temp,split = ', '))
  temp = t(matrix(unlist(strsplit(temp,split = '-')),2,length(temp)))
  temp = data.frame(chr = i,start = as.numeric(temp[,1]),end = as.numeric(temp[,2]))
  temp = temp[order(temp[,2],temp[,3]),]
  data_cds_2 = rbind(data_cds_2,temp)
}

# Time difference of 2.164124 secs
overlap_clear_cds_3 = overlap_clear(data_cds_2)
# 35185772
count_length_cds_3 = count_length(overlap_clear_cds_3)
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_01")

write.csv(overlap_clear_cds_3,"overlap_clear_cds_3_CCDS.current.txt.csv",row.names = F)
temp = read.csv("overlap_clear_cds_3_CCDS.current.txt.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==overlap_clear_cds_3)

write.csv(count_length_cds_3,"count_length_cds_3_CCDS.current.txt.csv",row.names = F)
temp = read.csv("count_length_cds_3_CCDS.current.txt.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==count_length_cds_3)

setwd("F:/code/biotrainee2017/data")
#################################################################################################




