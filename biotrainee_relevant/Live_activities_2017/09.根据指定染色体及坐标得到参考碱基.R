setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_09",recursive = T,showWarnings = F)




# 根据指定染色体及坐标得到参考碱基
#################################################################################################
# 其实在 02.计算人类基因序列GC含量.R 中的思路已经能 根据指定染色体及坐标得到参考碱基

# 生成一条随机序列
temp = sample(c("A","T","G","C"),50,replace = T)
temp_seq = paste(temp,collapse = "")
# CAATAGGTACGCGGGACATATCCAGCTCTTGTTGATTGTCTCAGTCTGGT
temp_seq

#      G
substring(temp_seq,6,6)
# CAATAGGTAC
substring(temp_seq,1,10)
# CAATAGGTAC
#     AGGTACGCGGGACATATCCAGCTCTTGTTGATTGTCTCAGTCTGGT
substring(temp_seq,c(1,5),c(10,50))
#################################################################################################

