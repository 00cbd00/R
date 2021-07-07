setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_11",recursive = T,showWarnings = F)




# count_length:统计 每条染色体对应序列总体长度
# 要求输入的数据框为 染色体编号、起始、终止位点 3列数据
#################################################################################################
chr_split_output = function(data){
  setwd("F:/code/biotrainee2017/result_11")
  for (i in unique(data_1[,1])) {
    # 对每一个染色编号分别处理
    ID = which(data_1[,1] == i)
    chr_temp =data_1[ID,]
    # 定义文件名
    chr_file = paste("chr_",i,".csv",sep = "")
    # 输出
    write.csv(chr_temp,chr_file,row.names = F)
    # temp = read.csv(chr_file,header = T,encoding = "UTF-8")
    # ##全为TRUE没有错
    # table(temp==chr_temp)
  }
  setwd("F:/code/biotrainee2017/data")
}
#################################################################################################




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
temp = data_1[1:1000,]
#################################################################################################

#################################################################################################
chr_split_output(data = data_1)
#################################################################################################
