# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:29:04 2024

@author: josue
"""

import numpy as np

# Simpson's rule implementation for discrete data
def simpson_integral_discrete(y_values, x_values):
    n = len(x_values) - 1  # Number of intervals; must be even for Simpson's rule
    if n % 2 == 1:
        raise ValueError("The number of intervals must be even for Simpson's rule.")
    
    h = (x_values[-1] - x_values[0]) / n  # Step size
    integral = y_values[0] + y_values[-1]  # First and last terms
    
    # Simpson's rule summation
    for i in range(1, n, 2):  # Odd indices
        integral += 4 * y_values[i]
    for i in range(2, n-1, 2):  # Even indices
        integral += 2 * y_values[i]
    
    integral *= h / 3  # Apply Simpson's rule formula
    return integral

# Calculate L2 norm of the error between two vectors using Simpson's rule
def l2_norm_error_simpson(numerical_solution, exact_solution, x_values):
    # Compute the error vector (element-wise difference squared)
    error_squared = (numerical_solution - exact_solution) ** 2
    
    # Integrate the squared error using Simpson's rule
    integral_error = simpson_integral_discrete(error_squared, x_values)
    
    # L2 norm is the square root of the integral
    return np.sqrt(integral_error)

# # Example vectors (you can replace these with your actual data)
# x_values = np.linspace(0, np.pi, 101)  # 101 points from 0 to pi
# exact_solution = np.sin(x_values)  # Exact solution (e.g., sin(x))
# numerical_solution = np.sin(x_values) + np.random.normal(0, 0.1, len(x_values))  # Numerical solution with some noise

# # Ensure number of intervals is even (101 points -> 100 intervals)
# if (len(x_values) - 1) % 2 != 0:
#     raise ValueError("The number of intervals must be even for Simpson's rule.")

# # Calculate the L2 norm of the error
# l2_error_norm = l2_norm_error_simpson(numerical_solution, exact_solution, x_values)
# print(f"L2 norm of the error using Simpson's rule: {l2_error_norm}")
