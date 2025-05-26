# STED related analysis

EU_enrichment_analysis.py基于原始的pATM/EU荧光TIFF图片，得到pATM foci的sphere模型，并计算EU在每个球体区域的相对荧光强度（EU enrichment），最后计算每个foci的体积与EU enrichment的Spearman相关系数。

make_sphere.py基于原始pATM荧光图片，对于每个pATM foci，以它的几何中心为球心，做直径为1 μm的球体，并转化为可视化的TIFF图片。

sphere_model.py用于将make_sphere.py得到的pATM foci的球体模型可视化，构建3D模型。
