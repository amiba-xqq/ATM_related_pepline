# Counting the foci number

ImageJ macro脚本用于统计每个细胞核内荧光蛋白的foci数目，最终得到的csv格式的结果文件，第一列是细胞核编号，第二列是对应核内foci的数目。

`Batch export merged images from LIF files (including DAPI, red, and green channels).ijm`得到的TIFF图片，使用`Counting the foci number (Both red and green channels).ijm`分析foci数目，其它类型的TIFF图片同理。

具体使用流程可以参考`Protocol of Counting the foci number.pptx`。

These ImageJ macro scripts count the number of fluorescent protein foci within each nucleus. The results are saved in a CSV file where the first column is the nucleus ID and the second column lists the number of foci per nucleus.

Use the `Batch export merged images from LIF files (including DAPI, red, and green channels).ijm` script to generate TIFF images. Then, use the `Counting the foci number (Both red and green channels).ijm` script to analyze foci number in these images. Similar methods apply to other TIFF image types.

For specifics on how to use these scripts, refer to the `Protocol of Counting the foci number.pptx`.
