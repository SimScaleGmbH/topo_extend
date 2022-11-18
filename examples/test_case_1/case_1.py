from topoExtend.Topo_extend_v2 import topology
import pathlib
import numpy as np

from sklearn import preprocessing


input_path = pathlib.Path('test_case_1.stl')
output_path = pathlib.Path('TOPOLOGY_EXTENSION.stl')

mesh_clean = topology(origin=[99.5, 101, 0],
                      resolution=0.5)

mesh_clean.extend_stl(input_path,
                      output_path,
                      extension_radius=100,
                      inclusion_radius=50)

no_triangles = mesh_clean.get_no_triangles()
'''
mesh_clean.reduce_mesh_target_size(50) #Megabytes
no_triangles = mesh_clean.get_no_triangles()
print('After mesh reduction: {}'.format(no_triangles))
'''