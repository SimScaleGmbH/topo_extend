#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 09:28:31 2022

@author: darrenlynch
"""

import numpy as np
import scipy.interpolate as si


def bspline(cv, n=100, degree=3, periodic=False):
    """ 
        Example from:
        https://stackoverflow.com/questions/34803197/fast-b-spline-algorithm-with-numpy-scipy
        
        Calculate n samples on a bspline

        cv :      Array ov control vertices
        n  :      Number of samples to return
        degree:   Curve degree
        periodic: True - Curve is closed
                  False - Curve is open
    """

    # If periodic, extend the point array by count+degree+1
    cv = np.asarray(cv)
    count = len(cv)

    if periodic:
        factor, fraction = divmod(count+degree+1, count)
        cv = np.concatenate((cv,) * factor + (cv[:fraction],))
        count = len(cv)
        degree = np.clip(degree,1,degree)

    # If opened, prevent degree from exceeding count-1
    else:
        degree = np.clip(degree,1,count-1)


    # Calculate knot vector
    kv = None
    if periodic:
        kv = np.arange(0-degree,count+degree+degree-1)
    else:
        kv = np.clip(np.arange(count+degree+1)-degree,0,count-degree)

    # Calculate query range
    u = np.linspace(0,(count-degree),n)


    # Calculate result
    return np.array(si.splev(u, (kv,cv.T,degree))).T

def p_sigmoid_from_data(data, ratio, min_prob=0.1):
    maximum = np.max(data)*1.01
    minimum = np.min(data)*0.99
    
    
    mid_x = (maximum - minimum)*(1-ratio)
    
    maximum_probability = 1
    minimum_probability = min_prob
    
    cv = np.array([[ maximum,  maximum_probability],
                   [ mid_x,  maximum_probability],
                   [mid_x, minimum_probability],
                   [ minimum,  minimum_probability]])
    
    p = bspline(cv,n=100,degree=3,periodic=False)
    
    return p, cv

def interp_from_p(p, data):
    #output = np.interp(data, p[:, 0], p[:, 1])
    f = si.interp1d(p[:, 0], p[:, 1], kind='nearest')
    
    output = f(data)
    
    return output

def get_probability_from_graient2(gradient2, ratio, min_prob=0.1):
    print("Calculating the sigmoid...")
    p, cv = p_sigmoid_from_data(gradient2, 0.5, min_prob)
    print("Interpolating data to the sigmoid...")
    interp = interp_from_p(p, gradient2)
    print("...Finished interpolating")
    
    return interp

'''
import matplotlib.pyplot as plt

data = np.random.rand(100000)*10
data = np.sort(data)

p, cv = p_sigmoid_from_data(data, 0.5, 0.25)

interp = interp_from_p(p, data)

plt.plot(cv[:,0],cv[:,1], 'o-', label='Control Points')

x,y = p.T
plt.plot(x,y,'r-',label='Degree %s')
plt.plot(data,interp,'k-',label='Degree %s')

plt.minorticks_on()
plt.legend()
plt.xlabel('Second Gradient')
plt.ylabel('Probability')
plt.xlim(0, 10)
plt.show()
'''