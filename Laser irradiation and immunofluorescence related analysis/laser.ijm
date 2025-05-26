// 设置参数:base_dir是laser图片的所在路径；roi_dir是图片中每个细胞对应的laser的ROI区域文件所在路径，需要提前在Fiji(ImageJ)中手动选择区域，并保存；Threshold_green是设置的阈值
// Set parameters: base_dir is the path where laser images are stored; roi_dir is the path to ROI files corresponding to each cell in the laser images, which need to be manually selected in Fiji (ImageJ) and saved beforehand; Threshold_green is the defined threshold
base_dir = "C:\\Users\\Yourname\\Laser\\NBS1-GFP IRAK1i\\";
roi_dir = "C:\\Users\\Yourname\\Laser\\roi\\cell1.zip";
Threshold_green = 60;

list = getFileList(base_dir);
for (i=0; i<list.length; i++) {
    folder = list[i];
    folder_path = base_dir + folder;
    green_path = folder_path + "_green_level.csv";    
    open(folder_path);
    selectImage(folder);
    setAutoThreshold("Default dark");
    run("Threshold...");
    setThreshold(Threshold_green, 255, "raw");
    roiManager("Open", roi_dir);
    roiManager("Measure");

    close("Threshold");
    close("ROI Manager");   
    close(folder);
}
saveAs("Results", green_path);
close("*");
close("Results");
 
       
