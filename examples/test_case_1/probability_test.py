#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:32:18 2022

@author: darrenlynch
"""

import numpy as np
from sklearn import preprocessing

#Generate random sample grid
grid = np.random.rand(10, 10)*10

def createProabilityMatrix(matrix, _max=1):
    absolute_matrix = np.abs(matrix)
    
    return preprocessing.normalize(np.log(absolute_matrix))

def randProb(prob):
    rand = np.random.rand(*prob.shape)
        
    return np.where(rand < prob, 1, 0)


prob_array = createProabilityMatrix(grid, 1)
print(randProb(prob_array))