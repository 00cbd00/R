setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_06",recursive = T,showWarnings = F)


library(stringr)



# hsa00001.keg
# https://www.genome.jp/kegg-bin/get_htext?hsa00001+3101
#################################################################################################
# 数据格式不是很规律，不好读取，直接用 readLines 函数
con = file("hsa00001.keg","r")
line = readLines(con,n=100000)
close(con)
#################################################################################################


#################################################################################################
# 观察数据结构
temp = as.data.frame(line[1:1000])
#   !     #     %     +     A     B     C     D 
#   2    11     1     1     8   112   534 50197 
table(substr(line,1,1))
#  5  3176  4799  9682 12326 20533 33647 50250
which(substr(line,1,1)=="A")
# [1] "A09100 Metabolism"                           "A09120 Genetic Information Processing"      
# [3] "A09130 Environmental Information Processing" "A09140 Cellular Processes"                  
# [5] "A09150 Organismal Systems"                   "A09160 Human Diseases"                      
# [7] "A09180 Brite Hierarchies"                    "A09190 Not Included in Pathway or Brite" 
line[which(substr(line,1,1)=="A")]

# 提取C通路和D基因
kegg_gene = line[which(substr(line,1,1)%in%c("C","D"))]
head(kegg_gene)

# 提取C通路
ID = which(substr(kegg_gene,1,1)=="C")
temp = kegg_gene[ID]
head(temp)
# 将每个C通路扩大至与对应的D基因数相同
kegg_gene_2 = rep(temp,c(ID[2:length(ID)]-1,length(kegg_gene))-ID)
# 提取C编号和对应的通路名
kegg_gene_2 = cbind(substr(kegg_gene_2,6,10),substr(kegg_gene_2,12,99999))
# 清除通路名中多余的通路编号
kegg_gene_2[,2] = str_split(kegg_gene_2[,2]," \\[",simplify = T)[,1]

# 提取D基因
ID = which(substr(kegg_gene,1,1)=="D")
temp = kegg_gene[ID]
head(temp)
# 提取基因编号
kegg_gene_2 = cbind(kegg_gene_2,str_split(temp,"[ ]",simplify = T)[,7])
# 提取基因名
temp = substr(temp,8,99999)
temp = str_split(temp,";",simplify = T)[,1]
temp = str_split(temp,"\t",simplify = T)[,1]
temp = gsub("^[0-9]* ","",temp)
kegg_gene_2 = cbind(kegg_gene_2,temp)

# 整合
kegg_extract = data.frame(gene_id = as.numeric(kegg_gene_2[,3]),
                          pathway_id = as.numeric(kegg_gene_2[,1]),
                          gene_name = kegg_gene_2[,4],
                          pathway_name = kegg_gene_2[,2])
#################################################################################################
# 输出结果
#################################################################################################
setwd("F:/code/biotrainee2017/result_06")

write.csv(kegg_extract,"kegg_extract.csv",row.names = F)
temp = read.csv("kegg_extract.csv",header = T,encoding = "UTF-8")
##全为TRUE没有错
table(temp==kegg_extract)
# kegg_extract = temp
setwd("F:/code/biotrainee2017/data")
#################################################################################################

