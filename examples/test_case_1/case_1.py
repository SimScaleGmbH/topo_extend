from topoExtend.Topo_extend_v2 import topology
import pathlib


input_path = pathlib.Path('test_case_1.stl')
output_path = pathlib.Path('TOPOLOGY_EXTENSION.stl')

mesh_clean = topology(origin=[99.5, 101, 0],
                      resolution=0.5)

mesh_clean.extend_stl(input_path,
                      extension_radius=500,
                      inclusion_radius=50,
                      region_of_interest_radius=25)

mesh_clean.export_mesh(output_path)