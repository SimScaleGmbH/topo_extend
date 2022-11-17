import open3d as o3d
import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.ndimage import gaussian_filter

from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler

from stl import mesh
import stl
import pyvista



class topology():
    
    def __init__(self, origin=np.array([0,0,0]), resolution=1):
        
        self.input_path = None
        self.output_path = None
        
        self.mesh_min = None
        self.mesh_max = None
        
        self.origin = origin
        self.resolution = resolution 
        
        self.disc_radius = None
        self.extension_radius = None
        self.angle_resolution = None
        
        self.mesh = None
        self.output_mesh = None
        
        self.matrix = None
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
        Take the mesh, create aa matrix
        
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
        #   combined distance, gradient, point inclusion
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

        matrix = np.zeros([x_col.shape[0], 12])

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
        
    def _remove_outside_roi(self, inclusion_radius):
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
        
        self.disc_radius = inclusion_radius
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
        self.matrix[:, 10] = np.gradient(self.matrix[:, 9])
        
    def _inclusion(self):
        
        def randProb(prob):
            rand = np.random.rand(*prob.shape)
                
            return np.where(rand < prob, 1, 0)
        
        def createProabilityMatrix(matrix):
            print("Number of points in original mesh: {}".format(len(matrix)))
            absolute_matrix = np.abs(matrix.reshape(-1, 1))
            
            nintyith_percentile = np.percentile(absolute_matrix, 90)
            absolute_matrix = np.where(absolute_matrix > nintyith_percentile, 
                                       nintyith_percentile, absolute_matrix)
            
            scaler = MinMaxScaler()
            data = np.log(1*absolute_matrix)
            scaler.fit(data)
            
            normalised = scaler.transform(data)[:,0]
            points = randProb(1.1*normalised)
            print("Number of points in final mesh: {}".format(np.sum(points)))
            return points
        
        self.matrix[:, 11] = createProabilityMatrix(self.matrix[:, 10])
        
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
        
        reference_radius = radius - self.disc_radius + 10
        
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
        
    def plot_topology_gradient(self):
        '''
        Simply plot the height map in cartesian space

        Returns
        -------
        None.

        '''
        zz = self.matrix[:, 10].reshape(self.grid.shape)
        plt.imshow(np.abs(zz), vmin=0, vmax=10)
        plt.show()
        
    def plot_topology_points(self):
        '''
        Simply plot the height map in cartesian space

        Returns
        -------
        None.

        '''
        zz = self.matrix[:, 11].reshape(self.grid.shape)
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
        '''
        Take a path for an STL, export the height map into an stl mesh

        Parameters
        ----------
        output_path : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        
        def return_vector_List(xyz_array):
            #Legacy function
            x, y, z = xyz_array[1, :, :], xyz_array[0, :, :], xyz_array[2, :, :]
            
            x, y, z = x.reshape(-1), y.reshape(-1), z.reshape(-1)
            
            vector_list = np.array([x, y, z]).T
            return vector_list

        def return_face_vector_list(list_coordinate_points):
            #Legacy function
            list_coordinate_list = []
            for _list in list_coordinate_points:
                list_coordinate_list.append(return_vector_List(_list))
            
            return np.stack(list_coordinate_list, axis=1)

        def base_array_to_face_vector(base_array):
            #Legacy function
            coordinate_p1 = base_array[:, :-1, :-1]
            coordinate_p2 = base_array[:, 1:, :-1]
            coordinate_p3 = base_array[:, 1:, 1:]
            coordinate_p4 = base_array[:, :-1, 1:]
            
            #We need two faces per pixel
            vector_list1 = return_face_vector_list([coordinate_p1, 
                                                          coordinate_p2, 
                                                          coordinate_p3])
            
            vector_list2 = return_face_vector_list([coordinate_p1, 
                                                          coordinate_p4, 
                                                          coordinate_p3])

            vector_list = np.concatenate([vector_list1, vector_list2], axis=0)
            
            return vector_list
        
        self.output_path = output_path
        
        matrix = self.matrix
        mask = np.isinf(matrix[:, 8])
        matrix[:, 8][mask] = (np.ones(matrix[:, 8].shape)*self.mesh_min[2])[mask]
        
        xx = matrix[:, 0].reshape(self.grid.shape[0:2])
        yy = matrix[:, 1].reshape(self.grid.shape[0:2])
        zz = matrix[:, 8].reshape(self.grid.shape[0:2])
        
        base_array = np.array([yy, xx, zz])
        vector_list = base_array_to_face_vector(base_array)
        
        no_faces = (xx.shape[0]-1) * (yy.shape[1]-1) * 2

        final_mesh = mesh.Mesh(np.zeros(no_faces, dtype=mesh.Mesh.dtype))
        final_mesh.vectors = vector_list
        
        final_mesh.translate(self.origin)
        
        self.output_mesh=final_mesh
        
        final_mesh.save(output_path, mode=stl.stl.Mode.ASCII)
        
    def get_no_triangles(self):
        eval_mesh = o3d.io.read_triangle_mesh(
            self.output_path.as_posix())

        triangles = len(eval_mesh.triangles)

        return triangles
    
    def reduce_mesh_percentage(self, percentage):
        mesh = pyvista.read(self.output_path.as_posix())
        decimated = mesh.decimate(percentage)
        decimated.save(self.output_path.as_posix(), 
                  binary=False,
                  texture=None)
        
    def reduce_mesh_target_number(self, number_of_triangles=200000.0):
        percentage = 1 - (number_of_triangles / self.get_no_triangles())
        
        self.reduce_mesh_percentage(percentage)
        
    def reduce_mesh_target_size(self, no_mega_bytes=50):
        x1 = 200000
        y1 = 50.9
        
        x2 = 1000000
        y2 = 243.4
        
        m = (y2-y1)/(x2-x1)
        
        x_t = no_mega_bytes/m
        
        self.reduce_mesh_target_number(number_of_triangles=x_t)
        
    def extend_stl(self, 
                   input_path,
                   output_path,
                   extension_radius=1000,
                   inclusion_radius=300,
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
        self.extension_radius = extension_radius
        self.import_mesh(input_path)
        
        self.create_matrix()
        self._remove_outside_roi(inclusion_radius=inclusion_radius)

        self._create_polar_matrix(angular_resolution=1)
        
        self._interpolate_to_polar()
        self._fill_radius()
        self._interpolate_missing_from_polar()
        self._create_smoothed_matrix()
        self._blend_matricies()
        self._gradient()
        self._inclusion()
        
        self.plot_topology()
        self.plot_topology_gradient()
        self.plot_topology_points()

        self.export_mesh(output_path)