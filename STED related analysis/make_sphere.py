#该脚本用于
import numpy as np
import pandas as pd
from skimage import io, measure, morphology
import tifffile as tiff
from scipy.stats import spearmanr

def calculate_integrated_density(image, roi_mask):
    """计算ROI区域的integrated density，即该区域像素值的总和"""
    return np.sum(image[roi_mask])

def create_sphere_mask(center, radius, shape):
    z, y, x = np.indices(shape)
    distance = np.sqrt((z - center[0])**2 + (y - center[1])**2 + (x - center[2])**2)
    sphere_mask = distance <= radius
    return sphere_mask
    
def process_image_sequence(image_a_path, image_b_path, threshold_a, radius):
    img_a = tiff.imread(image_a_path)  
    img_b = tiff.imread(image_b_path)  

    # 对图像A进行二值化处理
    binary_a = img_a > threshold_a

    # 标签化图像A中的连通区域
    labeled_a = measure.label(binary_a)
    tiff.imsave('binary_pATM.tif', labeled_a)
    # 提取每个foci的几何属性
    regions = measure.regionprops(labeled_a)
    print("pATM foci number:",len(np.array([region.label for region in regions])))
    
    results = []
    i = 0
    sphere_mask_total = np.zeros(img_a.shape, dtype=bool)

    for region in regions:
        # 获取foci的几何中心和体积
        center = region.centroid
        volume = region.area  
        if volume < 6:
            continue
        if radius < center[0] < ((labeled_a.shape[0])-radius) and radius < center[1] < ((labeled_a.shape[1])-radius)  and radius < center[2] < ((labeled_a.shape[2])-radius) :
        # 创建球体ROI mask  
            sphere_mask = create_sphere_mask(center, radius, img_a.shape)
            sphere_mask_total = sphere_mask_total | sphere_mask

        # 计算tif B图像在该ROI区域的integrated density
            integrated_density = calculate_integrated_density(img_b, sphere_mask)
        
        # 保存结果
            results.append([volume, integrated_density])
            i = i+1
            
    tiff.imsave('pATM_sphere.tif', sphere_mask_total)
    # 将结果保存为CSV
    df = pd.DataFrame(results, columns=['Foci Volume', 'Integrated Density'])
    df.to_csv('foci_statistics.csv', index=False)

# 设置参数
image_a_path = 'pATM.tif'
image_b_path = 'EU.tif'
threshold_a = 100 
radius = 10

# 调用函数处理图像序列
process_image_sequence(image_a_path, image_b_path, threshold_a, radius)
