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

# # Example usage
# vector1 = np.array([1, 2, 3, 4, 5])
# vector2 = np.array([1, 3, 3, 0, 5, 6, 7])

# matches = compare_vectors_with_indices(vector1, vector2)
# if matches:
#     print("Matching components in the vectors with their indices:")
#     for val, idx1, idx2 in matches:
#         print(f"Value {val}: Index in vector1 = {idx1}, Index in vector2 = {idx2}")
# else:
#     print("There are no matching components in the vectors.")




# # Definir los dominios y las imágenes
# domain1 = np.array([1, 3, 5, 7])
# images1 = np.array([10, 30, 50, 70])

# domain2 = np.array([2, 4, 6, 8])
# images2 = np.array([20, 40, 60, 80])

# # # Concatenar los dominios y las imágenes
# # combined_domain = np.concatenate((domain1, domain2))
# # combined_images = np.concatenate((images1, images2))

# # # Ordenar los dominios y las imágenes de acuerdo al dominio ordenado
# # sorted_indices = np.argsort(combined_domain)
# # sorted_domain = combined_domain[sorted_indices]
# # sorted_images = combined_images[sorted_indices]

# # print(sorted_indices)
# # print("Combined Domain:", sorted_domain)
# # print("Combined Images:", sorted_images)

# print(merge_and_sort_vectors(domain1, images1, domain2, images2))