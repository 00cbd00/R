setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_10",recursive = T,showWarnings = F)







# position_info:查看每条序列的位置信息
# 要求输入的数据为 染色体编号、起始、终止位点
#################################################################################################
position_info = function(data,info=ensg_summary){
  info = data.frame(chr=info$chr,
                    start=info$start,
                    end=info$end,
                    gene_id=info$gene_id,
                    gene_name=info$gene_name)
  result = matrix(NA,0,7)
  for (i in 1:nrow(data)) {
    # 提取指定的染色体编号
    ID = which(info$chr==data$chr[i])
    temp_info = info[ID,]
    
    # 判断是否在某个基因中
    ID = which(temp_info$start<=data$start[i] & temp_info$end>=data$end[i])
    if(length(ID)!=0){
        temp = cbind(data[i,],type="in",temp_info[ID,c(4,5,2,3)],row.names = NULL)
        result = rbind(result,temp)
    }
    
    # 判断是否包含某个基因
    ID = which(temp_info$start>=data$start[i] & temp_info$end<=data$end[i])
    if(length(ID)!=0){
      temp = cbind(data[i,],type="include",temp_info[ID,c(4,5,2,3)],row.names = NULL)
      result = rbind(result,temp)
    }
    
    # 判断起始位点是否与某个基因重叠
    ID = which(temp_info$start<data$start[i] 
               & data$start[i]<=temp_info$end
               & temp_info$end<data$end[i])
    if(length(ID)!=0){
      temp = cbind(data[i,],type="start_overlap",temp_info[ID,c(4,5,2,3)],row.names = NULL)
      result = rbind(result,temp)
    }
    
    # 判断终止位点是否与某个基因重叠
    ID = which(data$end[i]<temp_info$end
               & temp_info$start<=data$end[i] 
               & data$start[i]<temp_info$start)
    if(length(ID)!=0){
      temp = cbind(data[i,],type="end_overlap",temp_info[ID,c(4,5,2,3)],row.names = NULL)
      result = rbind(result,temp)
    }
    
    # 提取与起始位点最接近的基因
    ID = which(temp_info$end<data$start[i])
    ID = which(temp_info$end==max(temp_info$end[ID]))
    temp = cbind(data[i,],type="start_close",temp_info[ID,c(4,5,2,3)],row.names = NULL)
    result = rbind(result,temp)
    
    # 提取与终止位点最接近的基因
    ID = which(temp_info$start>data$end[i])
    ID = which(temp_info$start==min(temp_info$start[ID]))
    temp = cbind(data[i,],type="end_close",temp_info[ID,c(4,5,2,3)],row.names = NULL)
    result = rbind(result,temp)
  }
  return(result)
}
#################################################################################################


# ensg_summary.csv
#################################################################################################
setwd("F:/code/biotrainee2017/result_03")

temp = read.csv("ensg_summary.csv",header = T,encoding = "UTF-8")
ensg_summary = temp
temp = ensg_summary[1:1000,]
setwd("F:/code/biotrainee2017/data")
#################################################################################################
# 10.data.test.txt
#################################################################################################
data_test = read.table("10.data.test.txt",sep = " ",stringsAsFactors=F,encoding = "UTF-8")
data_test = data.frame(chr=substr(data_test[,1],4,99999),
                       start=as.numeric(data_test[,9]),
                       end=as.numeric(data_test[,17]))
#################################################################################################


position_result = position_info(data = data_test)
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_10")

write.csv(position_result,"position_result.csv",row.names = F)
temp = read.csv("position_result.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==position_result)
# position_result = temp
setwd("F:/code/biotrainee2017/data")
#################################################################################################


