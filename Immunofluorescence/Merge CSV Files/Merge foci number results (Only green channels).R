rm(list=ls())
# 注意：引号里面的文件位置是Batch export TIFF images from LIF files脚本得到的TIFF文件所在的文件夹的路径；路径中的子文件夹之间需要使用/隔开，而不是\
# Note: The file path refers to the folder containing the TIFF files generated by the "Batch export TIFF images from LIF files" script; subfolders in the path should be separated by / instead of \.
folder_path <- "C:/Users/amiba-xqq/Desktop/confocol/IF_test_merged tif/"

subfolders <- list.dirs(folder_path, recursive = FALSE)
k=0
repeat{
  k=k+1
  if(k>length(subfolders)){
    break 
  }
  setwd(subfolders[k])
  getwd()
  
  # 查找当前工作目录下所有以"_green_foci.csv"结尾的CSV文件
  # Search for all CSV files ending with "_green_foci.csv" in the current working directory
  file_list <- list.files(pattern = "_green_foci.csv$")

  merged_data <- read.csv(file_list[1], header = TRUE)
  i=1
  repeat{
    i=i+1
    if(i>(length(file_list))){
      break
    }
    data2 <- read.csv(file_list[i], header = TRUE)
    merged_data <- rbind(merged_data, data2)
  }
  write.csv(merged_data, file = "green_foci_merge.csv", row.names = FALSE)
  
}

# 最后得到的green_foci_merge.csv文件展示了合并后各组green通道foci数目。
# The final green_foci_merge.csv file shows the merged number of green channel foci for each group.


