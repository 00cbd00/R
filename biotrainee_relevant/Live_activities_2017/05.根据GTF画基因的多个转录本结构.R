setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_05",recursive = T,showWarnings = F)


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


# 对gtf进行提取
#################################################################################################
data_2 = data_1
# 第3题已经对gtf文件做过了解,此处不再检验检查
ID = which(data_2[,3]=="gene")
# 按照 "; " 分隔
temp = str_split(data_2[ID,9],"; ",simplify = T)
# 提取ensg id
temp[,1] = str_split(temp[,1]," ",simplify = T)[,2]
# 提取 gene_name ,没有 gene_named 的显示空白
temp[,3] = ifelse(substr(temp[,3],1,9)=="gene_name",substr(temp[,3],11,99999),"")

# 提取ensg id
data_2[,9] = substr(data_2[,9],9,23)
# data_1 的ensg id 在temp中的位置
ID = match(data_2[,9],temp[,1])
# 增加 gene_name
data_2$gene_name = temp[ID,3]

gtf_grch38_extract = data.frame(chr=data_2[,1],
                                gene_id=data_2[,9],
                                gene_name=data_2[,10],
                                type=data_1[,3],
                                start=data_2[,4],
                                end=data_2[,5],
                                source=data_2[,2]
                                )
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_05")

write.csv(gtf_grch38_extract,"gtf_grch38_extract.csv",row.names = F)
temp = read.csv("gtf_grch38_extract.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==gtf_grch38_extract)
# gtf_grch38_extract = temp
setwd("F:/code/biotrainee2017/data")
#################################################################################################


# 大致画一下图，理解原理就行，不美化图片了
#################################################################################################
temp = gtf_grch38_extract[gtf_grch38_extract[,3]=="TP53",]
temp = temp[temp[,7]=="havana",]

num_transcript = length(which(temp[,4]=="transcript"))
min_start = min(temp[,5])
max_end = max(temp[,6])
max_length = max_end-min_start

plot(x = NULL,
     y = NULL,
     type = "n",
     xlab = "",
     ylab = "",
     xlim = c(min_start-max_length/20,max_end+max_length/20) ,
     ylim = c(0,24))
title(main = temp[1,3],
      sub = paste("chr",temp[1,1],":",min_start,"-",max_end,sep=""))
j=0
for (i in 1:nrow(temp)) {
  if(temp[i,4]=="transcript"){
    j=j+1
    lines(c(temp[i,5]-max_length/20,temp[i,6]+max_length/20),c(j,j),lwd=2)
  }
  if(temp[i,4]=="exon"){
    barplot(temp[i,6]-temp[i,5],
            width=0.8,
            horiz = T,
            add = T,
            space=j+0.25*j-0.45,
            offset=temp[i,5],
            col = "pink")
  }
  if(temp[i,4]=="CDS"){
    barplot(temp[i,6]-temp[i,5],
            width=0.8,
            horiz = T,
            add = T,
            space=j+0.25*j-0.45,
            offset=temp[i,5],
            col = "red")
  }
  if(temp[i,4]=="three_prime_utr" | temp[i,4]=="five_prime_utr"){
    barplot(temp[i,6]-temp[i,5],
            width=0.8,
            horiz = T,
            add = T,
            space=j+0.25*j-0.45,
            offset=temp[i,5],
            col = "yellow",
            legend.text = T)
  }
}
#################################################################################################

