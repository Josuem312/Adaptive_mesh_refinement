
from timeit import default_timer as timer

import numpy as np
import matplotlib.pyplot as plt
from Functions.Tree import Node
from Functions.Space_time_discretization import space_time
from Functions.AMR import AMR_cas_1
from Functions.scheme import scheme
from Functions.Initial_conditions import init
from Functions.Plotting import plotting_solution_flag
from Functions.Auxiliary_functions import merge_and_sort_vectors

###############################################
#     DISCRETIZATION OF TIME AND SPACE
################################################
a = 1  # wave velocity
tmin, tmax = 0.0, 1.5 # initial and final time
xmin, xmax = 0.0, 10.0  # 1D domain
Nx = 80  # "number of spacial points
c = 0.1 # Courant number, we need c<=1 for stability

x, dx, dt, Nt, time = space_time(
    c=c, a=a, xmin=xmin, xmax=xmax, Nx=Nx, tmin=tmin, tmax=tmax)
#################################################
#           OTHER PARAMETERS
##################################################
type_initial_condition = 3
U0 = init(type_initial_condition, x)
ratio_ref = 32
max_level = 4
threshold = 0.0001
type_of_scheme = 2
# We creat this variable in order to initialize solution vector in the subgrids
############################# RUNNUNG AMR CODE ##########################
#################################################################################
#               PLOTTING
#################################################################################
node1 = Node(center=0.5*(x[0] + x[-1]), dx=dx, dt=dt,
             time=0,
             value=U0, x_value=x, max_level=max_level,
             level=0)
##################################### FINAL SOLUTION ##################################



# errors=[]
# Courant= [0.1, 0.2, 0.4, 0.8]
# #inicio_parte1 = timer()  # Usar timeit.default_timer
# for co in Courant:
#     t = 0
#     x, dx, dt, Nt, time = space_time(
#         c=co, a=a, xmin=xmin, xmax=xmax, Nx=Nx, tmin=tmin, tmax=tmax)
#     node1 = Node(center=0.5*(x[0] + x[-1]), dx=dx, dt=dt,
#                  time=0,
#                  value=U0, x_value=x, max_level=max_level,
#                  level=0)
#     for _ in range(Nt):
#         # for i, t in enumerate(time[1:]):
#         # while  t<= tmax:
#         # for i in time[1:]:
#         # print (f"iter: {t//dt+1}")
#         u_solution_delta_t = AMR_cas_1(ratio_ref, co,
#                                        node1, type_of_scheme, threshold, type_initial_condition)
#         t = np.round(t + node1.dt, 10)
#         # plotting_solution_flag(node1, t, 1, f, 1)
#         # Medir tiempo de fin de la secció
#     rest = u_solution_delta_t- init(type_initial_condition, node1.x_value - a*node1.dt*Nt)
#     rest= rest**2
#     error= np.round(np.sqrt(node1.dx)*np.sqrt(sum(rest)), 15)
#     errors.append(error)
    
#     plt.figure(figsize=(10, 6))

#     plt.plot(node1.x_value, init(type_initial_condition,
#                                   node1.x_value - a*node1.dt*Nt), marker='o', label=f'Real solution at t={node1.dt*Nt}', color='red')
                                 
#                                   #a*node1.dt*Nt), marker='o', label='Real solution at t', color='red')

#     plt.plot(node1.x_value,
#               node1.value, marker='o', label='FINAL AMR sol', color='black')
#     # plt.axvline(x= node1.x_value[node1.children[0].flag[0]], color="red", linestyle='--', linewidth=1, label=f'x ={node1.x_value[node1.children[0].flag[0]]} ')  # Línea vertical en x = -5
#     # plt.axvline(x= node1.x_value[node1.children[0].flag[-1]], color="red", linestyle='--', linewidth=1, label=f'x ={node1.x_value[node1.children[0].flag[-1]]} ')  # Línea vertical en x = -5


#     plt.xlabel('x')
#     plt.ylabel('u')
#     plt.title(f'Solution at Courant number {c}')
#     plt.legend()
#     plt.grid(True)

#     plt.show()
    
#     print("ya paso")    
#     node1= None
        
t=0
for _ in range(Nt):
    # for i, t in enumerate(time[1:]):
    # while  t<= tmax:
    # for i in time[1:]:
    # print (f"iter: {t//dt+1}")
    u_solution_delta_t = AMR_cas_1(ratio_ref, c,
                                    node1, type_of_scheme, threshold, type_initial_condition)
    t = np.round(t + node1.dt, 10)
    # plotting_solution_flag(node1, t, 1, f, 1)
    # Medir tiempo de fin de la secció
                             
                              #a*node1.dt*Nt), marker='o', label='Real solution at t', color='red')



plt.figure(figsize=(10, 6))

plt.plot(node1.x_value, init(type_initial_condition,
                              node1.x_value - a*node1.dt*Nt), marker='o', label=f'Real solution at t={node1.dt*Nt}', color='red')
                             
                              #a*node1.dt*Nt), marker='o', label='Real solution at t', color='red')
plt.plot(node1.x_value,
          node1.value, marker='o', label='FINAL AMR sol', color='black')
# plt.axvline(x= node1.x_value[node1.children[0].flag[0]], color="red", linestyle='--', linewidth=1, label=f'x ={node1.x_value[node1.children[0].flag[0]]} ')  # Línea vertical en x = -5
# plt.axvline(x= node1.x_value[node1.children[0].flag[-1]], color="red", linestyle='--', linewidth=1, label=f'x ={node1.x_value[node1.children[0].flag[-1]]} ')  # Línea vertical en x = -5


plt.xlabel('x')
plt.ylabel('u')
plt.title(f'Solution at Courant number {c}')
plt.legend()
plt.grid(True)

plt.show()



