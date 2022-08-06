import pathlib

from topo_extend.Topo_extend_v2 import topology

origin = [99.5, 101, 0]
extension = 500
resolution = 1

case = 1

input_path = pathlib.Path('test_case_{}.stl'.format(case))
output_path = pathlib.Path('cleaned_v2.stl')

mesh_clean = topology(origin=origin,
                      resolution=resolution)

mesh_clean.extend_stl(input_path,
                      output_path,
                      origin=origin,
                      extension_radius=extension,
                      inclusion_radius=300)