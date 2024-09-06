# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 10:24:48 2024

@author: josue
"""

import numpy as np
from Functions.Auxiliary_functions import ary
## INPUT: #Nx: Numb_spacial_points/ c: courant number/ a: wave velocity

def space_time (c, a, xmin, xmax, Nx, tmin, tmax): #Nx: Numb_spacial_points/ c: courant number
    # Discretize
    x = np.linspace(xmin, xmax, Nx + 1)  # space discretization(attention with Nx+1)
    dx = float((xmax - xmin) / Nx)  # spatial step size
    dt = np.round( (c / a) * dx, 15) # stable time step calculated from the stability requirement
    Nt = int((tmax - tmin) / dt)  # number of time steps
    time = ary(tmin, tmax, dt) # time discretization
    return x, dx, dt, Nt, time


def space_time_v2 (c, a, xmin, xmax, Nt, tmin, tmax): #Nx: Numb_spacial_points/ c: courant number
    # Discretize
    t = np.linspace(tmin, tmax, Nt + 1)  # space discretization(attention with Nx+1)
    dt = float((tmax - tmin) / Nt)  # spatial step size
    dx = round( a / c * dt, 10) # stable time step calculated from the stability requirement
    Nx = int((xmax - xmin) / dx)  # number of time steps
    x = ary(xmin, xmax, dx) # time discretization
    return t, dx, dt, Nx, x