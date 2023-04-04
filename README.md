# Short Summary
A module that takes an .stl file of a topology that is not prepared for a Pedestrian Wind Comfort Analysis on the SimScale platform, and return an stl file that is.

### Workflow
1. Export your layers into .stl file format, a list of typical layers are as follows:
    * BUILDINGS_OF_INTEREST.stl
    * CONTEXT.stl
    * TOPOLOGY.stl
2. Use the example script pointing the input files at the TOPOLOGY.stl 
3. Zip the TOPOLOGY_EXTENSION.stl, BUILDINGS_OF_INTEREST.stl and CONTEXT.stl into a .zip and upload to SimScale

### Installation
```bash
pip install git+https://github.com/DHLynch/topo_extend.git
```

### Quick Start
Import this package and pathlib
```bash
from topo_extend.Topo_extend_v2 import topology
import pathlib
```

Set the stl to import using a path, and define the output stl using a path
```bash
input_path = pathlib.Path('test_case_1.stl')
output_path = pathlib.Path('TOPOLOGY_EXTENSION.stl')
```

Initialise the topology object
```bash
mesh_clean = topology(origin=[0,0,0],
                      resolution=1)
```

Perform the extension
```bash
mesh_clean.extend_stl(input_path,
                      extension_radius=500,
                      inclusion_radius=50)
```

Export the meshes as STL's'
```bash
mesh_clean.export_mesh(output_path)
```

Output
Original Mesh             |  Mesh after extension
:-------------------------:|:-------------------------:
![](graphics/original_mesh.png)  |  ![](graphics/extended_mesh.png)

# Introduction
SimScale's LBM/PWC solution is incredibly robust to CAD/Geometry/3D model quality issues, however some preperation is still required to ensure results are reliable.

One of these steps involves ensuring that a topology or terrain (if pressent) surface is hole free (within reason) and extends past the boundary of a wind tunnel. More information as to why this is required can be found directly on the SimScale documentation pages.

Typically the user would need to open their model in a 3D moddeling tool such as Rhino to do these modifications, however, usually due to time constraints or knowledge gaps this process produces its own problems in terms of artifacts or still has issues.

# A perfect solution
There is no perfect solution, but lets talk briefly about the properties of the topology that would make it perfect.

1. The topology is a single surface that is big enough to intersect all edge faces of the wind tunnel.
2. The topology has enough resolution and no artifacts that would artificially influence the flow.
3. The topology would be completelly flat at the inlet so as not to influence the sape of the ABL.

Now, from the above we can see that to meet those points we would need to satisfy contradictary requirments and possible pay alot to obtain the size of topology required.

# This solution
To be able to obtain a reasonable topology surface from a typical users topology we need to do the following:

1. Take a complex collection of meshes (as an STL) and return the top surface mesh.
2. Isolate a circle that is smaller than some inclusion radius from a centre (origin) location.
3. Extend the parimeter out from a centre point.
4. Smooth the extension so that it does not have any artificial influence on results and that the ABL is not over fitting the topology or influenced dramatically at the boundaries.
5. Export a new STL file.

An experienced user will know that more original topology is better, but will be able to weigh up how much is needed.

# Topology holes and other imperfections
A hole in the topology is defined by an area withing the main disk, which is surrounded by topology, i.e. its missing topology not at the edge of the circle. 

Holes are not desired in a topology since flow might drastically effect CFD results. To combat this, we advise users to at least roughly patch the hole prior to running this script. However, to further safegard, if a hole is present this script automatically fills the hole with a height of the lowest point in the topology. This is the most simple form of handling and prevents the worst flow interactions with holes.

# Future development
1. Fill gaps in the mesh that lie in the inclusion circle
2. Create some reporting where the height map, and explaination of the topology elements are given. For example, what holes were filled automatically or where was the region of inclusion and region of interest defined.
3. A CLI function to do some default processing
4. integrate the CLI into a grasshopper component.

# Change log
17.08.22: Added functionality to reduce a mesh by %, also added a fix for some topology extensions were infinity was reported for the highest points.

19.08.22: Added functionality to reduce a mesh to a target number of triangles or to a target file size in megabytes.

04.04.23: Topology extension now exports two STL files, one for extension, one for the inner disk. New method implemented for mesh size reduction, what use to take hours now takes under 5 mins.
