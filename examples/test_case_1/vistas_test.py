#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 18:24:46 2022

@author: darrenlynch
"""

import pyvista as pv
points = pv.Polygon(n_sides=30).points
circle = pv.PolyData(points)
filled_circle = circle.delaunay_2d()
filled_circle.plot(show_edges=True, line_width=5)