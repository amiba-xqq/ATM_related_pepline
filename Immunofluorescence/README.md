Immunofluorescence assay
=============================

这些ImageJ macro和R脚本是用于定量分析DNA修复蛋白形成的每个细胞的foci数目，核内相对荧光强度。同样也适用于分析PLA实验的分析，用于统计PLA foci。

免疫荧光图片通过Leica共聚焦荧光显微镜拍摄，文件是.lif格式的，每个lif文件包含多组图片。

分析步骤为：
* 使用`Batch export TIFF images from LIF files`里的ImageJ macro脚本导出TIFF图片；
* 使用`Counting the foci number`或`Relative fluorescence intensity within nuclei`里的ImageJ macro脚本分析对应的指标，得到对应的csv结果文件；
* 使用`Merge CSV Files`里的R脚本合并所有统计的结果。

所有的ImageJ macro脚本只适用于windows系统。


ImageJ macro and R scripts were used to quantify the number of foci, and relative fluorescence intensity of DNA repair proteins within nuclei for each cell. Additionally, these scripts can also be used to analyze PLA (Proximity Ligation Assay) assay to count PLA foci.

Immunofluorescence images were captured by a confocal fluorescence microscope (Leica). Files are in .LIF format, with each LIF file containing multiple sets of images.

The analysis workflow is as follows:
(1) Use the ImageJ macro script in `Batch export TIFF images from LIF files` to export TIFF images;
(2) Use the ImageJ macro scripts in `Counting the foci number` or `Relative fluorescence intensity within nuclei` to analyze the 
corresponding metrics and obtain the corresponding CSV result files;
(3) Use the R script in `Merge CSV Files` to merge all statistical results.

All ImageJ macro scripts are only compatible with Windows systems.
