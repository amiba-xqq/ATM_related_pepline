// 第5行和第38行是唯一需要改的地方
// Lines 5 and 38 are the only ones that need to be modified
// 第5行是多个lif文件所在文件夹的路径，注意子文件夹名称之间间隔用\\而不是/
// On line 5, specify the path to the folder containing multiple LIF files. Note that subfolder names should be separated by backslashes (\\) instead of forward slashes (/)
input_path = "C:\\Users\\Yourname\\Desktop\\confocol\\PLA assay\\";
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
// Open the LIF file to retrieve all its image metadata
run("Bio-Formats Importer", "open=[" + base_dir_below + "] autoscale color_mode=Colorized open_all_series rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT contains=[] name=[" + base_dir_below + "]");
n_images = nImages/2;
width = getWidth();
height = getHeight();
// 循环处理每个图像并将其保存为TIFF文件
// Loop through each image and save it as a TIFF file
for (image=0; image < n_images; image++) {
  figure_m = image*2;
  selectImage(figure_m+1);
  merged_image = getTitle();
  merged_image0 = getTitle();
  selectImage(figure_m+2);
  merged_image1 = getTitle();
  // 下面中c2是green通道，c3是blue通道; merged_image0指C=0结尾的图片，merged_image1指C=1结尾的图片
  // Below: c2 is the green channel, c3 is the blue channel; merged_image0 refers to images ending with C=0, merged_image1 refers to those ending with C=1
  // 更改merged_image0/1的位置，使得它们和c2/c3的通道颜色相同
  // Modify the sequence of merged_image0, and 1 so that their colors align with those of c2, and c3. For instance, merged_image0 corresponds to images ending with C=0 (DAPI channel), merged_image1 to C=1 (green channel)
  
  run("Merge Channels...", "c2=[" + merged_image1 + "] c3=[" + merged_image0 + "] keep");
  figure_tile = replace(merged_image, ".*lif - (.*) - C=0", "$1");
  figure_tile = figure_tile + "_merge.tif";
  saveAs("Tiff", mfolder_file + "\\" +figure_tile);
  close("RGB.tif");
  close(figure_tile); 
}
close("*");
}

