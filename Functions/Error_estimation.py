# 
"""
Created on Thu Jul 18 12:16:56 2024

@author: josue
"""
from Functions.Flagging import flagging
import numpy as np

############################################################################
# Here we simply we calculate the formula of error estimator
############################################################################
def estimator(U_two_times, U_two_steps, order):
    error_estim= np.abs(U_two_steps-U_two_times)/(2**(order+1)-2)
    return error_estim


def error_estimation_flagging (Z2, Z3, order, threshold):
    error= estimator(Z2, Z3, order)
    flag, count= flagging(error, threshold)
    return error, flag, count