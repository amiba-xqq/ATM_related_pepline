// 设置参数，只需要更改3-8行代码的参数; 第3行是Batch export TIFF images from LIF files脚本得到的TIFF文件所在的文件夹
// Set parameters; only modify the parameters in lines 3–8 of the code;  Line 3 specifies the folder containing the TIFF files exported from 'Batch export TIFF images from LIF files' script
base_dir = "C:\\Users\\13103\\Desktop\\confocol\\niechen_20230609\\RNase rescue pATM_merged tif\\";
blue_min = 0;
blue_max = 255;
Threshold_green = 60;
min_size = 1000;
max_size = 20000;

// 遍历每个文件夹
// Loop through every folder
mlist = getFileList(base_dir);
for (j=0; j<mlist.length; j++) {
	mfolder = mlist[j];
	base_dir_below = base_dir + "\\" + mfolder + "\\";
	list = getFileList(base_dir_below);
for (i=0; i<list.length; i++) {
    folder = list[i];
    folder_path = base_dir_below + folder;
    folder_blue = "\\" + folder + " (blue)";
    folder_green = "\\" + folder + " (green)";
    green_path = folder_path + "_green_level.csv";
    figure_window = "Drawing of " + folder_blue;
    figure_path = folder_path + "_figure.tif" ;
    
    open(folder_path);
    run("Split Channels");
    // 打开DAPI图像并进行处理
    // Open the DAPI channel image and process it
    selectWindow(folder_blue);
    run("RGB Color");
    run("Brightness/Contrast...");
    setMinAndMax(blue_min, blue_max);
    run("8-bit");
    selectWindow(folder_blue);
    run("Convert to Mask");
    run("Convert to Mask");
    run("Fill Holes");
    run("Smooth");
    run("Smooth");
    run("Smooth");
    run("Smooth");
    run("Smooth");
    run("Smooth");
    run("Smooth");
    run("Smooth");
    run("Smooth");
    run("Smooth");
    run("Convert to Mask");
    run("Convert to Mask");
    run("Watershed");
   
   selectWindow(folder_green);
   setAutoThreshold("Default dark");
   run("Threshold...");
   setThreshold(Threshold_green, 255, "raw");
   
   selectWindow(folder_blue);
   run("Analyze Particles...", "size=" + min_size + "-" + max_size + " pixel show=Outlines exclude add");
 
   selectWindow(folder_green);
   roiManager("Measure");
   saveAs("Results", green_path);
   close("Results"); 
   selectWindow(figure_window);
   saveAs("Tiff", figure_path);
    // 关闭所有窗口
    // close all windows
    close("Threshold");
    close("ROI Manager");   
    close("B&C");  
    close("*");
}
}
    
 
       
