from topo_extend.Topo_extend_v2 import topology
import pathlib

origin = [99.5, 101, 0]
extension = 500
resolution = 1

input_path = pathlib.Path('test_case_1.stl')
output_path = pathlib.Path('TOPOLOGY_EXTENSION.stl')

mesh_clean = topology(origin=origin,
                      resolution=resolution)

mesh_clean.extend_stl(input_path,
                      output_path,
                      extension_radius=extension,
                      inclusion_radius=50)