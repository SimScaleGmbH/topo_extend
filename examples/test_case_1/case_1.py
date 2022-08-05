import pathlib

from topo_extend.Topo_extend_v2 import topology

origin = [99.5, 101, 0]
extension = 500
resolution = 1

case = 1

input_path = pathlib.Path('test_case_{}.stl'.format(case))
output_path = pathlib.Path('cleaned_v2.stl')

mesh_clean = topology(origin=origin,
                      extension=extension,
                      resolution=resolution)

mesh_clean.import_mesh(input_path)

mesh_clean.create_matrix()
mesh_clean._remove_outside_roi(inclusion_radius=50)
mesh_clean.plot_topology()

#internal_functions
mesh_clean._create_polar_matrix(extension_radius=2*extension, angular_resolution=1)
mesh_clean._interpolate_to_polar()
mesh_clean._fill_radius()
mesh_clean._interpolate_missing_from_polar()
mesh_clean._create_smoothed_matrix()
mesh_clean._blend_matricies()

mesh_clean.plot_topology()
mesh_clean.plot_polar_topology()

mesh_clean.export_mesh(output_path)