import open3d as o3d
import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.ndimage import gaussian_filter

from stl import mesh
import stl
import pyvista as pv

import blend_function as bf

class topology():
    
    def __init__(self, origin=np.array([0,0,0]), resolution=1):
        
        self.input_path = None
        self.output_path = None
        
        self.mesh_min = None
        self.mesh_max = None
        
        self.origin = origin
        self.resolution = resolution 
        
        self.roi_radius = None
        self.inclusion_radius = None
        self.extension_radius = None
        self.angle_resolution = None
        
        self.mesh = None
        self.output_mesh = None
        
        self.matrix = None
        self.inclusion_matrix = None
        
        self.polar_matrix = None
        self.grid = None
        
    def import_mesh(self, input_path):
        '''
        Take path, import mesh to open3d

        Parameters
        ----------
        input_path : pathlib.Path
            A path to the .stl file to import.

        Returns
        -------
        None.

        '''
        self.input_path = input_path
        
        input_mesh = o3d.io.read_triangle_mesh(
            input_path.as_posix())
        
        input_mesh.translate(-np.array(self.origin))
        
        self.mesh = input_mesh
        
    def create_matrix(self):
        '''
        Take the mesh, create a matrix
        
        The matrix is of nx11:
            x, y, z, nx, ny, nz, theta, radius, distance, blured distance, combined distance
            
        The last column is unused but reserved for gradient.
        
        Returns
        -------
        None.

        '''
        self.mesh_min = self.mesh.get_min_bound().astype(int)
        self.mesh_max = self.mesh.get_max_bound().astype(int)

        #Geometry point matrix will have dimensions:
        #   x, y, z, nx, ny, nz, theta, radius, distance, blured distance, 
        #   combined distance, gradient, 2nd Gradient, probability, point inclusion
        #
        #Therefore, matrix will be n x 9
        
        x = np.arange(-self.extension_radius, self.extension_radius, self.resolution)
        y = np.arange(-self.extension_radius, self.extension_radius, self.resolution)
        
        z = self.mesh_max[2] + 1 #We add 1 to ensure the initial 
                                 #point is higher than the max

        xx, yy, zz = np.meshgrid(x, y, z)

        x_col = xx.reshape(-1)
        y_col = yy.reshape(-1)
        z_col = zz.reshape(-1)

        matrix = np.zeros([x_col.shape[0], 14])

        matrix[:, 0] = x_col
        matrix[:, 1] = y_col
        matrix[:, 2] = z_col
        matrix[:, 5] = -np.ones(x_col.shape)
        matrix[:, 6] = np.degrees(np.arctan2(y_col, x_col))
        matrix[:, 7] = (x_col**2 + y_col**2)**0.5


        tensor = matrix[:, 0:6]
        rays = o3d.core.Tensor(tensor,
                               dtype=o3d.core.Dtype.Float32)

        #triangles = input_mesh.triangles
        input_mesh = o3d.t.geometry.TriangleMesh.from_legacy(self.mesh)

        scene = o3d.t.geometry.RaycastingScene()
        scene.add_triangles(input_mesh)


        ans = scene.cast_rays(rays)
        t_hit = ans['t_hit'].numpy()
        
        #We subtract 1 to counter the 1 we added earlier
        height = self.mesh_max[2] - t_hit.reshape(xx.shape[0:2]) + 1

        self.grid = height.reshape(xx.shape)

        matrix[:, 8] = height.reshape(-1)
        
        self.matrix = matrix
        
    def _remove_outside_inclusion(self, inclusion_radius):
        '''
        Takes an inclusion radius and keeps anything within that radius.

        Parameters
        ----------
        inclusion_radius : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        index = np.where(self.matrix[:, 7] >= inclusion_radius)[0]
        
        matrix = self.matrix[:, 8]
        
        matrix[index] = -np.inf
        
        self.inclusion_radius = inclusion_radius
        self.matrix[:, 8] = matrix
        
    def _create_polar_matrix(self, angular_resolution=1):
        '''
        Create a matrix in polar coordinate space, theta for angle, r for radius

        Parameters
        ----------
        angular_resolution : int or float, optional
            This is the resolution in the angular space, typically 1 is used,
            where one would expect 360 as the matrix dimension. 
            
            The default is 1.

        Returns
        -------
        None.

        '''
        
        self.angle_resolution = angular_resolution
        
        a = np.arange(-180, 180 + angular_resolution, angular_resolution)
        r = np.arange(0, self.extension_radius*2, self.resolution)
        
        aa, rr = np.meshgrid(a, r)
        
        a_shape = aa.shape
        aa = aa.reshape(a_shape[0], a_shape[1], 1)
        
        r_shape = aa.shape
        rr = rr.reshape(r_shape[0], r_shape[1], 1)
        
        zz = np.zeros(rr.shape)
        dz = np.zeros(rr.shape)
        
        matrix = np.concatenate([aa, rr, zz, dz], axis=2)
        self.polar_matrix = matrix
        
    def _interpolate_to_polar(self):
        '''
        Interpoalate from cartesian to polar space

        Returns
        -------
        None.

        '''
        points = self.matrix[:, 6:8]
        #aa = self.matrix[:, 6].reshape(self.grid.shape[0:2])
        #rr = self.matrix[:, 7].reshape(self.grid.shape[0:2])
        zz = self.matrix[:, 8]#.reshape(self.grid.shape[0:2])
        
        new_a = self.polar_matrix[:,:,0]
        new_r = self.polar_matrix[:,:,1]
        
        zz1 = interpolate.griddata(points, zz, (new_a, new_r), 
                                   method='nearest', 
                                   fill_value=-np.inf)
        
        self.polar_matrix[:, :, 2] = zz1
        
    def _fill_radius(self):
        '''
        Fill the missing data with the parimeter height

        Returns
        -------
        None.

        '''
        matrix = self.polar_matrix
        
        index = np.where(np.isfinite(matrix[:,:, 2])==True)
        matrix[index[0],index[1], 3] = index[0]
        
        edge = np.argmax(matrix[:,:,3], axis=0)
        
        missing_data_index = np.where(np.isfinite(matrix[:,:, 2])==False)
        edge_index = edge[missing_data_index[1]]
        
        matrix[missing_data_index[0],missing_data_index[1], 2] \
            = matrix[edge_index[0],missing_data_index[1], 2]
        
        self.polar_matrix = matrix
        
    def _interpolate_missing_from_polar(self):
        '''
        interpolate anything missing in the cartesian space from the polar space

        Returns
        -------
        None.

        '''
        points = np.vstack(self.polar_matrix[:, :, 0:2])
        zz = self.polar_matrix[:, :, 2].reshape(-1)
        
        new_a = self.matrix[:, 6].reshape(self.grid.shape[0:2])
        new_r = self.matrix[:, 7].reshape(self.grid.shape[0:2])
        
        zz1 = interpolate.griddata(points, zz, (new_a, new_r), method='linear')
        zz1 = zz1.reshape(-1)
        index = np.where(np.isfinite(self.matrix[:, 8]) == False)[0]
        
        matrix = self.matrix[:, 8]
        
        matrix[index] = zz1[index]
        
    def _gradient(self):
        '''
        Create matrix of gradient magnitude, and second gradient magnitude
        
        Second gradients are used to identify where harsh changes in topology
        are present, and therefore, will need more resolution

        Returns
        -------
        None.

        '''
        zz = self.matrix[:, 9].reshape(self.grid.shape[0], self.grid.shape[1])
        gradient_1 = np.gradient(zz)
        gradient_1_mag = (gradient_1[0]**2 + gradient_1[1]**2)**0.5
        self.matrix[:, 10] = gradient_1_mag.reshape(-1)
        
        zz = self.matrix[:, 10].reshape(self.grid.shape[0], self.grid.shape[1])
        gradient_2 = np.gradient(zz)
        gradient_2_mag = (gradient_2[0]**2 + gradient_2[1]**2)**0.5
        self.matrix[:, 11] = gradient_2_mag.reshape(-1)
        
    def _inclusion(self):
        '''
        Determin is a pixel is included in the final mesh
        
        Probabilistic way of determining if a pixel should be included,
        if a high second gradient, a higher probability is given, and more 
        resolution is therefore present in the end mesh.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        
        def randProb(prob):
            rand = np.random.rand(*prob.shape)
                
            return np.where(rand < prob, 1, 0)
        
        self.matrix[:, 12] = self._createProabilityMatrix()
        
        self.matrix[:, 13] = randProb(self.matrix[:, 12])
        
        self._cut_circle()
        
        boolean_matrix = np.where(self.matrix[:, 13] == 1, True, False)
        
        self.inclusion_matrix = self.matrix[boolean_matrix, :]
    
    def _createProabilityMatrix(self):
        '''
        Create a probability matrix from second gradient
        
        As the previous method states, a higher probability is given to high
        second gradient values.

        Returns
        -------
        normalised : np.array
            An array of 1 or 0.

        '''
        absolute_matrix = np.abs(self.matrix[:, 11])
        
        data = absolute_matrix
        normalised_outer = bf.get_probability_from_graient2(data, 0.9, 0.1)
        inner_lower_bound = 0.25
        
        normalised_inner = bf.get_probability_from_graient2(
            data, 0.9, inner_lower_bound)
        
        normalised = np.where(
            self.matrix[:, 7] < self.inclusion_radius, 
            normalised_inner, normalised_outer)
        
        return normalised
        
        
    def _cut_circle(self):
        '''
        Cuts the grid, into a circle of same diameter

        Returns
        -------
        None.

        '''
        self.matrix[:, 13] = np.where(self.matrix[:, 7] > self.extension_radius,
                                      0, self.matrix[:, 13])
        
    def _create_smoothed_matrix(self):
        '''
        create a smoothed cartesian space using gausian blur

        Returns
        -------
        None.

        '''
        zz = self.matrix[:, 8].reshape(self.grid.shape[0:2])
        
        zz_smoothed = gaussian_filter(zz, sigma=5)
        
        self.matrix[:, 9] = zz_smoothed.reshape(-1)
        
    def _blend_matricies(self):
        '''
        Blend the original and smoothed cartesian space
        
        This blend is done at the inclusion radius, using a sigmoid function

        Returns
        -------
        None.

        '''
        zz_original = self.matrix[:, 8]
        zz_smoothed = self.matrix[:, 9]
        
        radius = self.matrix[:, 7]
        
        reference_radius = radius - self.inclusion_radius + 10
        
        #clip to 1
        blend_function = 1/(1 + np.exp(-reference_radius))
        
        blended_matrix = (zz_original*(1-blend_function)) + \
            (zz_smoothed*(blend_function))
            
        self.matrix[:, 8] = blended_matrix
        
    def plot_topology(self):
        '''
        Simply plot the height map in cartesian space

        Returns
        -------
        None.

        '''
        zz = self.matrix[:, 8].reshape(self.grid.shape)
        plt.imshow(zz)
        plt.show()
    
    def plot_topology_gradient1(self):
        '''
        Simply plot the height map in cartesian space

        Returns
        -------
        None.

        '''
        
        
        zz = self.matrix[:, 10].reshape(self.grid.shape)
        abs_zz = np.abs(zz)
        nintyith_percentile = np.percentile(abs_zz, 99)
        
        
        plt.imshow(abs_zz, vmin=0, vmax=nintyith_percentile)
        plt.show()    
    
    def plot_topology_gradient2(self):
        '''
        Simply plot the height map in cartesian space

        Returns
        -------
        None.

        '''
        
        
        zz = self.matrix[:, 11].reshape(self.grid.shape)
        abs_zz = np.abs(zz)
        nintyith_percentile = np.percentile(abs_zz, 99)
        
        
        plt.imshow(abs_zz, vmin=0, vmax=nintyith_percentile)
        plt.show()
        
    def plot_topology_probability(self):
        '''
        Simply plot the height map in cartesian space

        Returns
        -------
        None.

        '''
        zz = self.matrix[:, 12].reshape(self.grid.shape)
        plt.imshow(zz)
        plt.show()
        
    def plot_topology_points(self):
        '''
        Simply plot the height map in cartesian space

        Returns
        -------
        None.

        '''
        zz = self.matrix[:, 13].reshape(self.grid.shape)
        plt.imshow(zz, cmap='binary')
        plt.show()
        
    def plot_polar_topology(self):
        '''
        simply plot the height map in polar space

        Returns
        -------
        None.

        '''
        zz = self.polar_matrix[:, :, 2]
        plt.imshow(zz)
        plt.show()
        
    def export_mesh(self, output_path):
        
        self.output_path = output_path
        
        xyz = self.inclusion_matrix[:, [0, 1, 8]]
        points = pv.PolyData(xyz)
        
        remesh = points.delaunay_2d()
        
        cell_centres = remesh.cell_centers().points
        
        cell_radius = (cell_centres[:, 0]**2 + cell_centres[:, 1]**2)**0.5
        
        farfield_cell_idx = np.transpose(np.nonzero(cell_radius>=self.inclusion_radius))
        farfield_cells = remesh.extract_cells(farfield_cell_idx).extract_surface()
        
        if self.region_of_interest_radius == None:
            nearfield_cell_idx = np.transpose(np.nonzero(
                (cell_radius<=self.inclusion_radius)))
        elif self.region_of_interest_radius >= self.inclusion_radius:
            raise Exception("Cannot have a region of interest radius larger \
                            than inclusion radius")
        else:
            nearfield_cell_idx = np.transpose(np.nonzero(
                (cell_radius<=self.inclusion_radius) & \
                (cell_radius>self.region_of_interest_radius)))
                
            roi_cell_idx = np.transpose(np.nonzero(
                cell_radius<=self.region_of_interest_radius))
            roi_cells = remesh.extract_cells(roi_cell_idx).extract_surface()
            recentered_roi = roi_cells.translate(self.origin, inplace=True)
            
            roi_path = output_path
            roi_path = roi_path.with_stem('TOPOLOGY_ROI')
            
            recentered_roi.save(roi_path,
                                binary=False,
                                texture=None)
            
        nearfield_cells = remesh.extract_cells(nearfield_cell_idx).extract_surface()
        
        recentered_farfield = farfield_cells.translate(self.origin, inplace=True)
        recentered_nearfield = nearfield_cells.translate(self.origin, inplace=True)
        
        farfield_path = output_path
        farfield_path = farfield_path.with_stem('TOPOLOGY_EXTENSION')
        
        recentered_farfield.save(farfield_path,
                                 binary=False,
                                 texture=None)
        
        nearfield_path = output_path
        nearfield_path = nearfield_path.with_stem('TOPOLOGY_Inclusion')
        
        recentered_nearfield.save(nearfield_path, 
                                 binary=False,
                                 texture=None)
        
    def extend_stl(self, 
                   input_path,
                   extension_radius=2000,
                   inclusion_radius=400,
                   region_of_interest_radius=300,
                   debug=False
                   ):
        '''
        A function that wraps a typical process into a workflow

        Parameters
        ----------
        input_path : pathlib.Path
            The .stl file to import.
        output_path : pathlib.Path
            the .stl file to output.
        extension_radius : int or float, optional
            The distance in meters to extend the topology by. 
            
            The default is 1000.
        inclusion_radius : int or float, optional
            A distance in metres from the origin in which we want to keep. 
            
            The default is 300, the same as a default ROI in SimScale.

        Returns
        -------
        None.

        '''
        self.region_of_interest_radius= region_of_interest_radius
        self.extension_radius = extension_radius
        self.import_mesh(input_path)
        
        self.create_matrix()
        self._remove_outside_inclusion(inclusion_radius=inclusion_radius)

        self._create_polar_matrix(angular_resolution=1)
        
        self._interpolate_to_polar()
        self._fill_radius()
        self._interpolate_missing_from_polar()
        self._create_smoothed_matrix()
        self._blend_matricies()
        
        print("Number of points in original mesh: {}".format(len(self.matrix[:, 10])))
        
        self._gradient()
        self._inclusion()
        
        print("Number of points in final mesh: {}".format(np.sum(self.matrix[:, 13])))
        
        if debug:
            self.plot_topology()
            self.plot_topology_gradient1()
            self.plot_topology_gradient2()
            
            self.plot_topology_probability()
            self.plot_topology_points()
            
if __name__ == "__main__":
    import pathlib
    
    input_path = pathlib.Path.cwd().parents[0] / (
        'examples/test_case_1/test_case_1.stl')
    output_path = pathlib.Path.cwd().parents[0] / (
        'examples/test_case_1/TOPOLOGY_EXTENSION.stl')

    mesh_clean = topology(origin=[99.5, 101, 0],
                          resolution=0.5)

    mesh_clean.extend_stl(input_path,
                          extension_radius=500,
                          inclusion_radius=50, 
                          region_of_interest_radius=25)

    mesh_clean.export_mesh(output_path)