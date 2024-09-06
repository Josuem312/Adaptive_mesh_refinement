# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:21:12 2024

@author: josue
"""

import numpy as np
import matplotlib.pyplot as plt
from Functions.scheme import scheme


def plotting_solution_flag(node, t, plot_flag,f, level):
    plt.clf()  # Limpiar la figura
    plt.plot(node.x_value, node.value)#u_solution_delta_t)
    plt.plot(node.x_value, f(node.x_value-t) )
    Z1= scheme(1, f(node.x_value), 0.4)[0] 
    plt.plot(node.x_value, Z1 , linestyle='--')
    plt.plot(node.x_value, np.zeros_like(node.x_value), linestyle='None', marker='o')

    if plot_flag==1:    
        #if level==1:
        #     for i, obj in enumerate(node.children[0].flag):
        #         if i==0 or i==(np.size(node.children[0].flag)-1):
        #             plt.vlines(x=node.children[0].x_value[obj], ymin=0, ymax=1.25, colors=['black'], linestyles=['--'], linewidth=0.5)
        #         else:            
        #             plt.vlines(x=node.children[0].x_value[obj], ymin=0, ymax=1.25, colors=['yellow'], linestyles=['--'], linewidth=0.5)        
        #elif level==2:
            for i, obj in enumerate(node.children[0].children[0].flag):
                if i==0 or i==(np.size(node.children[0].children[0].flag)-1):
                    plt.vlines(x=node.children[0].children[0].x_value[obj], ymin=0, ymax=1.25, colors=['black'], linestyles=['--'], linewidth=0.5)
                else:            
                    plt.vlines(x=node.children[0].children[0].x_value[obj], ymin=0, ymax=1.25, colors=['r'], linestyles=['--'], linewidth=0.5)
    
    plt.title(f"Solution at {t}")
    plt.xlabel("x")
    plt.ylabel("Value u")
    
    plt.draw()  
    plt.pause(1)  

    plt.ioff()  
    plt.show()  
    return












# u_final=u_solution_delta_t


# # Graficar las funciones
# plt.figure(figsize=(10, 6))

# plt.plot(x, u_final, label='Values of u coarser', color='blue')
# plt.plot(x, f(x- 1), label= 'Real solution at t', color='red')


# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('Graphic')
# plt.legend()
# plt.grid(True)

# plt.show()



################### PLOTTING 2024-08-28 ###################"

#u_final=u_solution_delta_t


# if node1.hasChildren():  
#     if node1.children[0].hasChildren():
#         sol_dom, sol_image= merge_and_sort_vectors(node1.children[0].children[0].x_value,
#                               node1.children[0].children[0].value,
#                                 node1.children[0].x_value,  node1.children[0].value) 
#         # Final_dom, Final_sol= merge_and_sort_vectors(sol_dom, sol_image, 
#         #                         node1.children[0].children[0].x_value,
#         #                         node1.children[0].children[0].value)
#         Final_dom, Final_sol= merge_and_sort_vectors(node1.x_value, node1.value, 
#                                 sol_dom,
#                                 sol_image)
#     else: 
#         Final_dom, Final_sol= merge_and_sort_vectors(node1.x_value,
#                               node1.value,
#                                 node1.children[0].x_value,  node1.children[0].value)  
# else:   
#     Final_dom = node1.x_value
#     Final_sol = node1.value
# Graficar las funciones
# plt.figure(figsize=(10, 6))

# #plt.plot(node1.x_value, u_final, label='Values of u coarser', color='blue')
# plt.plot(node1.x_value, init(type_initial_condition, node1.x_value - a*tmax), label= 'Real solution at t', color='red')
# #plt.plot(node1.x_value, f(node1.x_value), label= 'Initial sol', color='green')

# #plt.plot(node1.x_value, f(node1.x_value - 1), label= 'Real solution at t', color='red')
# #plt.plot(node1.x_value, u_final, label= 'FINAL AMR SOL', color='black')
# plt.plot(node1.x_value, 
#           node1.value
#           , label= 'FINAL AMR SOL', color='black')
# # plt.plot(node1.children[0].x_value, 
# #           node1.children[0].value
# #           , label= 'FINAL AMR SOL', color='black')


# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('Graphic')
# plt.legend()
# plt.grid(True)

# plt.show()
