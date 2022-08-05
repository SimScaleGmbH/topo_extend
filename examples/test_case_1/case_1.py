import pathlib

from topo_extend import topology

#origin = [99.5, 101, 0]
origin = [-5790, -35000, 0]
extension = 2000
resolution = 1

case = 6

input_path = pathlib.Path('test_case_{}/test_case_{}.stl'.format(case, case))
output_path = pathlib.Path('test_case_{}/cleaned_v2.stl'.format(case))

mesh_clean = topology(origin=origin,
                      extension=extension,
                      resolution=resolution)

mesh_clean.import_mesh(input_path)

mesh_clean.create_matrix()
mesh_clean._remove_outside_roi(inclusion_radius=1000)
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