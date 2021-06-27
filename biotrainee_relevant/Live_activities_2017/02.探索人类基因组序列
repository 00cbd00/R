http://www.biotrainee.com/thread-625-1-1.html


setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_02",recursive = T,showWarnings = F)


library(stringr)



# https://www.ncbi.nlm.nih.gov/datasets/genomes/?acc=GCF_000001405.39
# https://www.ncbi.nlm.nih.gov/genome/?term=GRCh38
#################################################################################################
# 预读100000行 观察数据
temp = read.table(gzfile("GCF_000001405.39_GRCh38.p13_genomic.fna.gz"),
                  sep = "\t",
                  stringsAsFactors=F,
                  nrows=100000,
                  encoding = "UTF-8",
                  comment.char = "")

# 读取基因组文件
# data_1 = read.table(gzfile("GCF_000001405.39_GRCh38.p13_genomic.fna.gz"),
#                     sep = "\t",
#                     stringsAsFactors=F,
#                     # nrows=100000,
#                     encoding = "UTF-8",
#                     comment.char = "")
# 读取后占内存7个多G,电脑不行，换分段读取
#################################################################################################


# seq_info:查看每条序列信息
# 要求输入的数据为fasta格式的染色体序列信息
#################################################################################################
seq_info = function(data){
  time_1 = Sys.time()
  con = file(data,"r")
  # 每10w行进行分段读取
  line=readLines(con,n=100000)
  i = 0
  result = matrix(NA,0,2)
  while( length(line) != 0 ) {
    # 读取以 > 为开头的行
    ID = grep("^>",line)
    if(length(ID) != 0){
      # 提取以 > 为开头的行所在行号以及整行信息
      result = rbind(result,matrix(c(ID+i*100000,line[ID]),length(ID),2))
    }
    if(length(line) != 100000){
      # 提取最后一行的行号并命名为： end_line
      result = rbind(result,matrix(c(length(line)+i*100000,"end_line"),length(ID),2))
    }
    line=readLines(con,n=100000)
    i = i+1
  }
  print(difftime(Sys.time(), time_1, units = 'secs'))
  close(con)
  return(result)
}
#################################################################################################


# base_count:统计各种碱基含量
# 要求输入的数据为fasta格式的染色体序列信息
# ,"GC_ratio","N_ratio"


# base_count:统计各种碱基含量
# 要求输入的数据:
# data:     fasta格式的染色体序列信息
# seq_info: seq_info运行的结果
# cds_info: cds相关的染色体编号、起始、终止位点 3列数据
#################################################################################################
base_count = function(data,seq_info,cds_info){
  time_1 = Sys.time()
  result = matrix(NA,0,7)
  colnames(result) = c("chr","Length","A","T","G","C","N")
  con = file(data,"r")
  for (i in 1:(nrow(seq_info)-1)) {
    time_2 = Sys.time()
    # 以 每段染色体的序列信息行 开始,终止位点 行结束,分段读取
    line=readLines(con,n=as.numeric(seq_info[i+1,1])-as.numeric(seq_info[i,1]))
    # 以 >NC 开头的完整的染色体序列信息
    if(substr(line[1],2,3) == "NC"){
      # 提取染色体编号
      chr = strsplit(line[1],split = "[ ,]")[[1]][5]
      # 如果提取不到染色体编码，则为线粒体
      chr = ifelse(chr == "","MT",chr)
      # 重复序列为小写，全部替换成大写
      line = paste(toupper(line[-1]),collapse = "")
      # 统计 A、T、G、C、N 含量
      temp = str_count(line,c("A","T","G","C","N"))
      # 合并
      result = rbind(result,c(chr,nchar(line),temp))
      
      # 匹配cds文件染色体编号
      ID = which(cds_info[,1] == chr)
      # 合并染色体上的所有cds序列
      line = paste(substring(line,cds_info[ID,2],cds_info[ID,3]),collapse = "")
      # 统计 A、T、G、C、N 含量
      temp = str_count(line,c("A","T","G","C","N"))
      # 合并
      result = rbind(result,c(paste("cds_",chr,sep = ""),nchar(line),temp))
      # 每读取一段运行后显示时间
      print(c(chr,difftime(Sys.time(), time_2, units = 'secs')))
    }
    else{
      # 每读取一段运行后显示时间
      print(c(substr(line[1],2,3),difftime(Sys.time(), time_2, units = 'secs')))
    }
  }
  # 转换数据类型
  result = data.frame(chr=result[,1],
                      length = as.numeric(result[,2]),
                      A=as.numeric(result[,3]),
                      T=as.numeric(result[,4]),
                      G=as.numeric(result[,5]),
                      C=as.numeric(result[,6]),
                      N=as.numeric(result[,7]))
  # 统计非 A、T、G、C、N 的碱基数量
  result$other = result[,2]-(result[,3]+result[,4]+result[,5]+result[,6]+result[,7])
  # 统计 N 占总长的比例
  result$N_ratio = round(result[,7]/result[,2],7)
  # 统计 GC 含量
  result$GC_ratio = round((result[,5]+result[,6])/(result[,2]-result[,7]),7)
  print(difftime(Sys.time(), time_1, units = 'secs'))
  close(con)
  return(result)
}
#################################################################################################


# F:/code/biotrainee2017/result_01
#################################################################################################
setwd("F:/code/biotrainee2017/result_01")
# 提取去重叠后的cds文件
cds_position = read.csv("overlap_clear_cds_1_Homo_sapiens.GRCh38.104.gtf.gz.csv",
                        header = T,encoding = "UTF-8")
setwd("F:/code/biotrainee2017/data")
#################################################################################################


# https://www.ncbi.nlm.nih.gov/datasets/genomes/?acc=GCF_000001405.39
#################################################################################################
# seq_info:查看每条序列信息
# Time difference of 148.1825 secs
seq_info_grch38 = seq_info("GCF_000001405.39_GRCh38.p13_genomic.fna.gz")
# >NC >NT >NW end 
# 25 358 256   1 
table(substr(seq_info_grch38[,2],1,3))

# Time difference of 571.9037  secs
# base_count:统计各种碱基含量
base_count_grch38 = base_count(data = "GCF_000001405.39_GRCh38.p13_genomic.fna.gz",
                               seq_info = seq_info_grch38,
                               cds_info = cds_position)
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_02")

write.csv(seq_info_grch38,"seq_info_grch38.csv",row.names = F)
temp = read.csv("seq_info_grch38.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==seq_info_grch38)

write.csv(base_count_grch38,"base_count_grch38.csv",row.names = F)
temp = read.csv("base_count_grch38.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==base_count_grch38)
setwd("F:/code/biotrainee2017/data")
#################################################################################################




