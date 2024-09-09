# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:21:48 2024

@author: josue
"""
import numpy as np
from Functions.Auxiliary_functions import compare_vectors_with_indices

"""
# This function creates updates the coarser grid solution:
# 1. Create a vector of the same dimensio of the vector to update
# 2. In our list l, we are going to consider its indexes and each component
# associated to his respectevely index.
# 3. We ask if the size of the object is equal to 1 because in that case 
# we are treating with a coarser point.
# 4. We print the object and his index
# 5. We print the past value of the vector to update and the new value taken from
# the fine grid solution.
# 6. Finally we replace.
"""

def updating_new_node(node):
    new_vector= np.copy(node.value)
    matches= compare_vectors_with_indices(node.x_value, node.children[0].x_value)
    for obj in matches:
        new_vector[obj[1]]=node.children[0].value[obj[2]]
    return new_vector    

def updating_new_node_version2(node, l):
    inicial_flag= node.children[0].flag[0]
    for i, obj in enumerate(l):
        node.value[inicial_flag+i]= node.children[0].value[obj]
    return node.value
