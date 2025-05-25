// 设置参数，只需要更改3-11行代码的参数; 第3行是Batch export TIFF images from LIF files脚本得到的TIFF文件所在的文件夹
// Set parameters; only modify the parameters in lines 3–11 of the code;  Line 3 specifies the folder containing the TIFF files exported from 'Batch export TIFF images from LIF files' script
base_dir = "C:\\Users\\Yourname\\Desktop\\confocol\\HeLa_merged tif\\";
blue_min = 0;
blue_max = 255;
red_min = 127;
red_max = 128;
min_size = 1000;
max_size = 20000;
min_speckle = 2;
max_speckle = 15;

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
    folder_red = "\\" +folder + " (red)";
    red_path = folder_path + "_red_foci.csv";
    red_path2 = folder + "_red_foci.csv";
    red_figure_path = folder_path + "_red_foci.tif";
    list_blue= "Speckle List " + folder_blue;
    
    // 打开并分离RGB通道
    // Open the image and split it into RGB channels
    open(folder_path);
    run("Split Channels");
    // 打开DAPI图像并进行处理
    // Open the DAPI channel image and process it
    selectWindow(folder_blue);
    run("RGB Color");
    run("Brightness/Contrast...");
    setMinAndMax(blue_min, blue_max);
    run("8-bit");
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
    
    // 打开red图像并进行处理
    // Open the green channel image and process it
    selectWindow(folder_red);
    run("RGB Color");
    run("Brightness/Contrast...");
    setMinAndMax(red_min, red_max);
    run("8-bit");
    run("Make Binary");
    run("Watershed");
    
    // 分析red foci数量
    // Analyze the number of red foci
    selectWindow(folder_red);
    run("Speckle Inspector", "primary=[" + folder_blue + "] secondary=[" + folder_red + "] redirect=[" + folder_red + "] min_primary_size=" + min_size + " max_primary_size=" + max_size + " min_secondary_size=" + min_speckle + " max_secondary_size=" + max_speckle + " show=primary exclude speckle");
    selectWindow(list_blue); 
    saveAs("Results", red_path);
    selectWindow("Inspector of "+folder_blue+"-1");
    saveAs("Tiff", red_figure_path);
    
    // 关闭所有窗口
    // close all windows
    close("B&C");
    close("ROI Manager"); 
    close(red_path2);
    close("*");
}
}
    
