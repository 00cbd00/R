setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_04",recursive = T,showWarnings = F)


# 整合count结果
#################################################################################################
count_data = list.files(pattern=".count",recursive = T)
# 先读一个数据观察
temp = read.table(count_data[1],sep = "\t",header = F,
                  stringsAsFactors=F,encoding = "UTF-8",comment.char = "")


# 读取每个数据的第二列--count列
reads_count = temp
for (i in 2:length(count_data)) {
  temp = read.table(count_data[i],sep = "\t",header = F,
                    stringsAsFactors=F,encoding = "UTF-8",comment.char = "")
  # reads_count的每一个ensg id在temp中的位置
  ID = match(reads_count[,1],temp[,1])
  reads_count = cbind(reads_count,temp[ID,2])
}

length(reads_count[1,])==length(count_data)+1
##TRUE,没问题，数据全读进来了
temp = unlist(strsplit(count_data,"[.]count"))
colnames(reads_count) = c("ENSG_id",temp)

##治疗后的数据列
ID=grep("[02468]",substring(colnames(reads_count),10,10))
##按照治疗前治疗后排序
reads_count = cbind(reads_count[,-ID],reads_count[,ID])
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_04")

write.csv(reads_count,"reads_count.csv",row.names = F)
temp = read.csv("reads_count.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==reads_count)
setwd("F:/code/biotrainee2017/data")
#################################################################################################





