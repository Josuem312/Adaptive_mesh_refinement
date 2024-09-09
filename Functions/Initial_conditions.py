# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 13:07:13 2024

@author: josue
"""

import numpy as np


def init(type_initial_condition, x):
    y = np.zeros_like(x)
    if type_initial_condition==1:
        y=  np.exp(-200*(x+0.5)**2) + 0.1*np.sin(2*np.pi*(x+0.5))
    elif type_initial_condition==2:
        #print(x)
        y[np.where(x < 1.4)] = 1.0
        #print(y)
    elif type_initial_condition== 3:
        y = 0.2*np.exp(-(x-4)**2)
    elif type_initial_condition==4:
        y= np.maximum(0.5 - np.abs(x - 2) / 0.5, 0)
    elif type_initial_condition==5:    
        #tv= (x < 0.4) & (x > 0.2)
        y = 50 * (0.01 - (x - 0.3)**2) * np.exp(0.01 / ((x - 0.3)**2 - 0.01))
        y[np.where((x > 1.4) | (x < 1.2))] = 0.0
    return y
