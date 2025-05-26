# 该脚本用于基于原始的pATM/EU荧光TIFF图片，得到pATM foci的sphere模型，并计算EU在每个球体区域的相对荧光强度，即EU enrichment
# 最后导出的foci_statistics.csv，分别包含每个foci的体积与EU enrichment，接着分析两者的Spearman相关系数
# This script processes original pATM/EU fluorescence TIFF images to generate sphere models of pATM foci and calculates relative fluorescence intensity of EU within each spherical region, i.e., EU enrichment
# The final exported foci_statistics.csv includes the volume and EU enrichment for each focus. Subsequently, their Spearman correlation coefficient is analyzed
import numpy as np
import pandas as pd
from skimage import io, measure, morphology
import tifffile as tiff
from scipy.stats import spearmanr

def calculate_integrated_density(image, roi_mask):
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
    # Apply thresholding to Image A
    binary_a = img_a > threshold_a
    
    # 标签化图像A中的连通区域
    # Label the connected regions in Image A
    labeled_a = measure.label(binary_a)
    tiff.imsave('binary_pATM.tif', labeled_a)
    # 提取每个foci的几何属性
    # Extract the geometric properties of each focus
    regions = measure.regionprops(labeled_a)
    print("pATM foci number:",len(np.array([region.label for region in regions])))
    
    results = []
    i = 0

    for region in regions:
        # 获取foci的几何中心和体积
        # Get the center and volume of each focus
        center = region.centroid
        volume = region.area
        if volume < 6:
            continue
        if radius < center[0] < ((labeled_a.shape[0])-radius) and radius < center[1] < ((labeled_a.shape[1])-radius)  and radius < center[2] < ((labeled_a.shape[2])-radius) :
        # 创建球体ROI mask
        # Create a sphere ROI mask
            sphere_mask = create_sphere_mask(center, radius, img_a.shape)

        # 计算tif B图像在该ROI区域的integrated density
        # Compute integrated density of TIF B in ROI
            integrated_density = calculate_integrated_density(img_b, sphere_mask)
            results.append([volume, integrated_density])
            i = i+1
    
    df = pd.DataFrame(results, columns=['Foci Volume', 'Integrated Density'])
    df.to_csv('foci_statistics.csv', index=False)

# 设置参数:threshold_a是pATM荧光图片进行二值化的阈值;radius是每个pATM foci对应的球体的半径，单位为px，根据比例尺进行，更改需要确保这个长度是0.5 μm
# Set parameters: threshold_a is the binary threshold for pATM fluorescence images; radius specifies the sphere radius (in pixels) corresponding to each pATM focus, adjustments should ensure this length corresponds to 0.5 µm
image_a_path = 'pATM.tif'
image_b_path = 'EU.tif'
threshold_a = 122  
radius = 6.25

process_image_sequence(image_a_path, image_b_path, threshold_a, radius)

# foci的体积与EU enrichment的Spearman相关系数分析
# Spearman correlation coefficient analysis between the volume and EU enrichment of each focus
df = pd.read_csv('foci_statistics.csv')
df = df.iloc[1:].reset_index(drop=True)
column1 = df.iloc[:, 0]
column2 = df.iloc[:, 1]
correlation, p_value = spearmanr(column1, column2)
print(f"Spearman Correlation: {correlation}")
print(f"P-value: {p_value}")
