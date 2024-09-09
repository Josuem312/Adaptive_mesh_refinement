# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 15:34:15 2024
@author: josue
"""
import numpy as np
from Functions.scheme import scheme_using_node
from Functions.scheme import scheme
from Functions.Error_estimation import error_estimation_flagging
from Functions.Initial_conditions import init
from Functions.Subgridding import subgridding_process
from Functions.Interpolation import u_finer_init_new_node_version_2_2
from Functions.Updating import updating_new_node
from Functions.Updating import updating_new_node_version2
import matplotlib.pyplot as plt


########## DEFINITION OF ADAPTIVE MESH REFINEMENT ################
def AMR_cas_1(  ratio_ref, c, 
               node, type_of_scheme, threshold, type_initial_condition):
    ## Storage of initial condition ##
    if node.level==0:
        power=0
    else:
        power=1
   
    for i in np.arange(ratio_ref**power):
### First, we calcule the needed vector for calculating the error estimator    
        Z1, Z2, Z3, order, stencil = scheme(type_of_scheme, node, c, 
                                        ratio_ref, type_initial_condition)

### Second Regridding process ? ###
        # Here we ask if we need a regriding process 
        if node.level< node.max_level :
            error, flag, count = error_estimation_flagging(Z2, Z3, order, threshold)
       
            if count > 0:          
### Third SUBGRIDDING  ###

                l = subgridding_process(node, flag, stencil,
                        ratio_ref, type_initial_condition) 
        
### Fourth STEPPING step
        node.time= np.round(node.time + node.dt, 10) 
        node.value_q= node.value
        node.value= Z1
        node.value_q1= Z1
        

### Fifth, we ask if there exists a new level? ###
        if node.level< node.max_level:
            if count>=1:
### Sixth, we call again our recursive function ###    
                node.children[0].value_q1= u_finer_init_new_node_version_2_2(node, node.children[0].x_value, flag, ratio_ref, stencil)                              
                AMR_cas_1( ratio_ref, c, 
                          node.children[0], type_of_scheme, threshold, 
                          type_initial_condition)           
                        

### Seventh, Updating part ###               
                node.value= updating_new_node_version2(node, l)

    return node.value

