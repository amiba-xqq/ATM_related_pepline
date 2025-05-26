# 该脚本用于基于原始pATM/EU的荧光图片，得到它们的3D模型
# This script is used for processing the original pATM/EU fluorescence images to obtain their 3D models
import numpy as np
import tifffile as tiff
from skimage import morphology
from skimage.filters import threshold_otsu
import pyvista as pv
from skimage.measure import label, regionprops
import imageio

def process_and_save_3d_model(tiff_file_path):
    image = tiff.imread(tiff_file_path)
    
    # 对图像进行二值化
    # Apply thresholding to image
    binary_image = image > 100
          
    # 去除噪点
    # Remove noise
    binary_clean = morphology.remove_small_objects(binary_image, min_size=6)    
    labeled_image = label(binary_clean)

    # 获取区域属性并转换为3D模型
    # Get region properties and convert to 3D model
    meshes = []
    regions = regionprops(labeled_image)
    for region in regions:        
        coords = region.coords
        if len(coords) > 0:
            points = np.array(coords, dtype=np.float32)
            mesh = pv.PolyData(points)
            mesh.point_data["colors"] = np.linspace(0, 1, len(points))
            meshes.append(mesh)
    full_mesh = meshes[0] if len(meshes) == 1 else meshes[0].merge(meshes[1:])
    
    # 如果只需要可视化3D模型，而不保存3D模型，则只需要运行第38、67行代码，不运行39、68行代码；如果只需要保存3D模型，不需要可视化，则只需要运行第39、68行代码，不运行38、67行代码
    # If you only want to visualize the 3D model without saving it, run lines 38 and 67; do not run lines 39 and 68. If you only want to save the 3D model without visualization, run lines 39 and 68; do not run lines 38 and 67
    #plotter = pv.Plotter()
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(
        full_mesh,
        scalars=full_mesh.points[:, 0],
        cmap='plasma',
        opacity=0.8,
        show_edges=False,
        lighting=True,
        show_scalar_bar=False,
        metallic=0.2,                   
        roughness=0.6,                 
        ambient=0.3,                    
        diffuse=0.7,                   
        specular=0.5                    
    )
    outline = full_mesh.outline()
    plotter.add_mesh(outline, color='white', line_width=4)
    plotter.add_light(pv.Light(
        position=(10, 10, 10),         
        color='white',                
        intensity=0.8                  
    ))
    
    plotter.enable_anti_aliasing()
    plotter.set_background('black')
    plotter.camera_position = 'xy'
    plotter.camera.azimuth = 60
    plotter.camera.elevation = 20
    #plotter.show()
    plotter.screenshot("pATM_model.tif")
    
    plotter.open_gif("rotation.gif")
    for i in range(0, 360, 5):
        plotter.camera.azimuth = i
        plotter.write_frame()
    plotter.close()

# 这里指原始pATM荧光TIFF图片的路径
# This refers to the path of the original pATM fluorescence TIFF images
tiff_file_path = "pATM.tif"  
process_and_save_3d_model(tiff_file_path)
