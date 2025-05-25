//第2行和第32行是唯一需要改的地方，本代码其它地方不可改动
input_path = "C:\\Users\\amiba-xqq\\Desktop\\confocol\\niechen\\sample figure\\";
// 创建新的文件夹
input_path2 = replace(input_path, "(.*\\\\).*\\\\.*$", "$1");
program_name = replace(input_path, ".*\\\\([^\\\\]+)\\\\$", "$1");
output_dir = input_path2 + program_name + "_merged tif";
File.makeDirectory(output_dir);

mlist = getFileList(input_path);
for (j=0; j<mlist.length; j++) {
	mfolder = mlist[j];
	mfolder_nolif = replace(mfolder, "^(.*)\\.lif$", "$1");
	base_dir_below = input_path + "\\" + mfolder;
	mfolder_file = output_dir + "\\" + mfolder_nolif;
	File.makeDirectory(mfolder_file);
// 打开LIF文件并获取所有图像元数据
run("Bio-Formats Importer", "open=[" + base_dir_below + "] autoscale color_mode=Colorized open_all_series rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT contains=[] name=[" + base_dir_below + "]");
n_images = nImages/2;
width = getWidth();
height = getHeight();
// 循环处理每个图像并将其保存为TIFF文件
for (image=0; image < n_images; image++) {
  figure_m = image*2;
  selectImage(figure_m+1);
  merged_image = getTitle();
  merged_image0 = getTitle();
  selectImage(figure_m+2);
  merged_image1 = getTitle();
  //下面中c1是red通道，c3是blue通道
  //merged_image0指C=0结尾的图片，merged_image1指C=1结尾的图片
  //更改merged_image0/1的位置，使得它们和c1/c3的通道颜色相同

  run("Merge Channels...", "c1=[" + merged_image1 + "] c3=[" + merged_image0 + "] keep");
  figure_tile = replace(merged_image, ".*lif - (.*) - C=0", "$1");
  figure_tile = figure_tile + "_merge.tif";
  saveAs("Tiff", mfolder_file + "\\" +figure_tile);
  close("RGB.tif");
  close(figure_tile); 
}
close("*");
}

