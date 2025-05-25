# 读取make_sphere.py得到的pATM_sphere.tif文件，并构建3D模型
# Read the 'pATM_sphere.tif' file obtained from 'make_sphere.py' and build the 3D model
import numpy as np
import tifffile as tiff
from skimage import morphology
from skimage.filters import threshold_otsu
import pyvista as pv
from skimage.measure import label, regionprops
import imageio

def process_and_save_3d_model(tiff_file_path):
    image = tiff.imread(tiff_file_path)
    
    labeled_image = label(image)
    #tiff.imsave('binary_pATM.tif', labeled_image)

    # 获取区域属性并转换为3D模型
    # Get region properties and convert to 3D model
    meshes = []
    regions = regionprops(labeled_image)
    for region in regions:        
        coords = region.coords
        if len(coords) > 0:
            points = np.array(coords, dtype=np.float32)
            mesh = pv.PolyData(points)
            # 设置为颜色渐变
            # Set to color gradient
            mesh.point_data["colors"] = np.linspace(0, 1, len(points))
            meshes.append(mesh)    
    full_mesh = meshes[0] if len(meshes) == 1 else meshes[0].merge(meshes[1:])
    
    # 如果只需要可视化3D模型，而不保存3D模型，则只需要运行第34、63行代码，不运行35、64行代码；如果只需要保存3D模型，不需要可视化，则只需要运行第35、64行代码，不运行34、63行代码
    # If you only need to visualize the 3D model without saving it, run lines 34 and 63; do not run lines 35 and 64. If you only need to save the 3D model without visualization, run lines 35 and 64; do not run lines 34 and 63.
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
    plotter.screenshot("pATM_sphere_model.tif")
    
    plotter.open_gif("pATM_sphere_rotation.gif")
    for i in range(0, 360, 5):
        plotter.camera.azimuth = i
        plotter.write_frame()
    plotter.close()

# 这里是make_sphere.py得到的pATM_sphere.tif文件
# This is the 'pATM_sphere.tif' file obtained from 'make_sphere.py'
tiff_file_path = "pATM_sphere.tif"  
process_and_save_3d_model(tiff_file_path)
