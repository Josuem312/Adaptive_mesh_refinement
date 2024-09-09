# Test.py

import numpy as np


"""
This function ary creates a new array from given minimum, maximum values 
 and the separation space between two grid points.
 See functioning of the no.lindspace 
"""
def ary(xmin, xmax, step):
    num_elements = int((xmax - xmin) / step) + 1
    v = np.linspace(xmin, xmax, num_elements)
    return v 


"""
This functions merge the values of two vectors into a single vector
and then the resulting vector sorts all its values.
"""
def merge_and_sort_vectors(space_vector1, vector1 ,  space_vector2, vector2):
    
    # Merge the two vectors
    merged_domain = np.concatenate((space_vector1 , space_vector2))
    merged_solution= np.concatenate((vector1 , vector2))
     
    sorted_indices = np.argsort(merged_domain)
    sorted_domain = merged_domain[sorted_indices]
    sorted_images = merged_solution[sorted_indices]
    
    return sorted_domain, sorted_images


"""
This function compare if there are repated values between two vectors and then
it save in a list the repeated value, and the index of the two vectors where
this repeated value is found.
"""
def compare_vectors_with_indices(vector1, vector2):
    vector1= np.round(vector1, 15)   
    vector2= np.round(vector2, 15)   
    # Create a dictionary to store the indices of the components of vector1
    indices_vector1 = {value: idx for idx, value in enumerate(vector1)}
    
    # Create a list to store the matches along with their indices
    matches = []
    
    # Iterate over vector2 and look for matches in the dictionary of vector1
    for idx2, value in enumerate(vector2):
        if value in indices_vector1:
            idx1 = indices_vector1[value]
            matches.append((value, idx1, idx2))
    
    return matches
