setwd("F:/code/biotrainee2017/data")
dir.create("F:/code/biotrainee2017/result_12",recursive = T,showWarnings = F)


library("rjson")

json_data = fromJSON(file = "metadata.cart.2021-07-07.json")

# JSON格式有专门的包读取，返回list
# 练习内容没啥意义，具体项目具体分析
