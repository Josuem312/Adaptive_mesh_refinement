
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
import matplotlib.animation as animation


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
type_initial_condition = 2
U0 = init(type_initial_condition, x)
ratio_ref = 4
max_level = 2
threshold = 0.0001
type_of_scheme = 1
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
  
t=0
for _ in range(Nt):
    u_solution_delta_t = AMR_cas_1(ratio_ref, c,
                                    node1, type_of_scheme, threshold, type_initial_condition)
    t = np.round(t + node1.dt, 10)


plt.figure(figsize=(10, 6))

plt.plot(node1.x_value, init(type_initial_condition,
                              node1.x_value - a*node1.dt*Nt), marker='o', label=f'Real solution at t={node1.dt*Nt}', color='red')
                             
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



