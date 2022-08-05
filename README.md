# Short Summary
A module that takes an .stl file of a topology that is not prepared for a Pedestrian Wind Comfort Analysis on the SimScale platform, and return an stl file that is.

### Installation
```bash
pip install git+https://github.com/DHLynch/topo_extend.git
```

# Introduction
SimScale's LBM/PWC solution is incredibly robust to CAD/Geometry/3D model quality issues, however some preperation is still required to ensure results are reliable.

One of these steps involves ensuring that a topology or terrain (if pressent) surface is hole free (within reason) and extends past the boundary of a wind tunnel. More information as to why this is required can be found directly on the SimScale documentation pages.

Typically the user would need to open their model in a 3D moddeling tool such as Rhino to do these modifications, however, usually due to time constraints or knoledge gaps this process produces its own problems in terms of artifacts or still has issues.

