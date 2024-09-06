import numpy as np
import matplotlib.pyplot as plt
from Functions.Space_time_discretization import space_time
import math
from Simpson_rule import l2_norm_error_simpson


##############################################
    #INITIAL CONDITION
##############################################
def f(X):
    """Asigna un valor de 1.0 para valores menores que 0.1"""
    f = np.zeros_like(X)
    #○f[np.where(X < 1.4)] = 1.0
    # if X <= 1.4 :
    #     f=1
    # else:
    #     f=0
    f = np.exp(-2*(X-4)**2)
    #f = np.exp(-200*(X+0.5)**2) + 0.1*np.sin(2*np.pi*(X+0.5))
    #print(xmax)
    return f

###############################################

# Functions for finite difference schemes for the transport equation
# dt u + a*dx u = 0

def euler_backward(u, c):
    Z1= np.copy(u)
    Z1[1:-1]=np.round( (1 - c) * u[1:-1] + c * u[:-2],15)

    # for i in np.arange(1, np.size(u)-2):
    #     Z1[i]= (1 - c) * u[i] + c * u[i-1]
    #     print(i) 
    ### NEULANN BOUNDARY CONDITIONS ###
    Z1[0]=Z1[1]
    Z1[-1]=Z1[-2]
    return  Z1

def Lax_wendroff(u, c):
    Z1= np.copy(u)
    #Z1[1:-1]= u[1:-1] - 0.5 * c * (u[2:] - u[:-2]) + 0.5 * c**2 * (u[2:] - 2*u[1:-1] + u[:-2])
    #Z1[1:-1]= u[1:-1] - 0.5 * c * (u[2:] - u[:-2]) + 0.5 * (c) ** 2 * (u[2:] - 2 * u[1:-1] + u[:-2])
    #Z1[1:-1]=  (1 - c**2)*u[1:-1] + 0.5*c*((c-1)*u[2:] + (c+1)*u[:-2])
    Z1[1:-1]= np.round(0.5*c*(1+c)*u[:-2] + (1-c**2)*u[1:-1] - 0.5*c*(1-c)*u[2:], 15)
    ### NEULANN BOUNDARY CONDITIONS ###
    Z1[0]=Z1[1]
    Z1[-1]=Z1[-2]
    return  Z1

def Lax_friedrich(u, c):
    Z1= np.copy(u)
    Z1[1:-1] =  np.round( 0.5*(u[:-2] + u[2:]) - 0.5*c * (u[2:] - u[:-2]), 15)
    ### NEULANN BOUNDARY CONDITIONS ###
    Z1[0]=Z1[1]
    Z1[-1]=Z1[-2]
    return  Z1



# Example usage
    
a = 1  # wave velocity
tmin, tmax = 0.0, 1.   #3.125e-05  # initial and final time
xmin, xmax = 0.0, 100.0  # 1D domain
c = .3 #0.2# Courant number, we need c<=1 for stability


Nx_values = np.array([2**k for k in range(5, 15)])   # Different spatial resolutions
#Nx_values = np.array([20])
errors_euler_backward = []
errors_lax_wendroff =[]
errors_lax_friedrich =[]
euler_solutions= []
lw_solutions=[]
lf_solutions=[]
dx_values =[]
dt_values =[]
x_vectors=[]



MSE_error=[]
for Nx in Nx_values:
    # Discretisaiton of time and space
    #x, dx, dt, Nt, time = space_time(c=c, a=a, xmin=xmin, xmax=xmax, Nx=Nx, tmin=tmin, tmax=tmax)
    x = np.linspace(xmin, xmax, Nx + 1) # discretization of space
    dx = float((xmax-xmin)/Nx) # spatial step size
    dt = c/a*dx # stable time step calculated from stability requirement
    time = np.arange(tmin, tmax + dt, dt) # discretization of time
    Nt = int((tmax - tmin) / dt)
    dx_values.append(dx)
    dt_values.append(dt)
    x_vectors.append(x)
    # Initialization
    #u_initial =  np.exp(-2*(x-4)**2)
    u_initial = f(x)
    # Solution with different schemes
    u_euler_backward = u_initial.copy()
    u_lw= u_initial.copy()
    u_lf= u_initial.copy()
    
    error=0
    samplePoints=0

    #for i, t in enumerate(time[1:]): # np.arange(Nt+1): 
    for _ in range(Nt):
        u_euler_backward = euler_backward(u_euler_backward, c)
        u_lw = Lax_wendroff(u_lw, c)
        u_lf = Lax_friedrich(u_lf, c)
        #error += np.sum((u_lw - f(x-a*t))**2) # add error from this timestep
        #samplePoints += len(u_lw)

        """
        Boundary conditions
        """
        #u_lw[0]= np.exp(-200*(x[0]+0.5)**2) + 0.1*np.sin(2*np.pi*(x[0]+0.5))
        #u_lw[-1]= np.exp(-200*(x[-1]+0.5)**2) + 0.1*np.sin(2*np.pi*(x[-1]+0.5))
        


    #○error = np.sqrt(error/samplePoints) # calculate rms-error
    #MSE_error.append(error)   
        
    
    euler_solutions.append(u_euler_backward)
    lw_solutions.append(u_lw)
    lf_solutions.append(u_lf)
        
    # Calculate de real solution
    #real_sol= np.sin(2 * np.pi * (x - a*tmax) )
    #real_sol= np.exp(-2*(x-4-a*tmax)**2)
    real_sol= f(x-a*tmax)
    
    # Calculate the error with respect to the initial condition
    # Methode 2
    rest_euler = u_euler_backward- real_sol
    rest_euler = rest_euler**2
    error_euler= np.round(np.sqrt(dx)*np.sqrt(sum(rest_euler)), 15)
    # Simpson Rule
    #error_euler = l2_norm_error_simpson(u_euler_backward, real_sol, x)

    
    

    rest_lw = u_lw - real_sol
    rest_lw = rest_lw**2
    error_lw= np.round(np.sqrt(dx)*np.sqrt(sum(rest_lw)), 15)
    # Simpson rule
    # error_lw = l2_norm_error_simpson(u_lw, real_sol, x)


    rest_lf = u_lf- real_sol
    rest_lf = rest_lf**2
    error_lf= np.round(np.sqrt(dx)*np.sqrt(sum(rest_lf)), 15)
    # Simpson rule 
    #error_lf = l2_norm_error_simpson(u_lf, real_sol, x)

    
    errors_euler_backward.append(error_euler)
    errors_lax_wendroff.append(error_lw)
    errors_lax_friedrich.append(error_lf)

"""
Graphing the slope 
"""

slope_ftbs=[]
slope_lw=[]
slope_lf=[]
counting=[]


for i in np.arange(np.size(Nx_values)-1):
    slope=  math.log(errors_euler_backward[i]/errors_euler_backward[i+1] , 2)
    slope_ftbs.append(slope)
    
    slope =  math.log(errors_lax_wendroff[i]/errors_lax_wendroff[i+1] , 2)
    slope_lw.append(slope)
    
    slope =  math.log(errors_lax_friedrich[i]/errors_lax_friedrich[i+1] , 2)
    slope_lf.append(slope)    
    
    counting.append(i+1)
    

# Plotting the order of convergence
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(counting, slope_ftbs, marker='o', label='FTBS')
plt.plot(counting, slope_lw, marker='o', label='LW')
#plt.plot(counting, slope_lf, marker='o', label='LF')


plt.xlabel('index')
plt.title('Order of convergence for different advection schemes')
plt.legend()
plt.grid(True)
plt.show()




# # Plotting the MSE ERROR
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# plt.plot(dx_values[:-1], ERRORS[:-1], marker='o', label='lax WENDROFF')
# plt.yscale('log')
# plt.xlabel('dx')
# plt.ylabel('Error')
# plt.title('MSE ERROR')
# plt.legend()
# plt.grid(True)
# plt.show()


def slope(x):
    return  2*x

# Plotting the order of convergence
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(dx_values[:-1], errors_euler_backward[:-1], marker='o', label='FTBS')
plt.plot(dx_values[:-1], errors_lax_wendroff[:-1], marker='o', label='Lax Wendroff')
plt.plot(dx_values[:-1], dx_values[:-1], color= "black",  label='Straight line slope=1')
plt.plot(dx_values[:-1], [ x**2 for x in dx_values[:-1]], color="red", label='Straight line slope=2')
#plt.plot(dx_values[:-1], errors_lax_friedrich[:-1], marker='o', label='Lax Friedrich')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('dx')
plt.ylabel('Error')
plt.title(f'Error vs space resolution for Courant number {c}')
plt.legend()
plt.grid(True)
plt.show()

slope_euler = (np.log(errors_euler_backward[-1])-np.log(errors_euler_backward[-2]))/(np.log(dx_values[-1])-np.log(dx_values[-2]))
print("slope FTBS in space", slope_euler)
slope_lw = (np.log(errors_lax_wendroff[-2]/errors_lax_wendroff[-3]))/(np.log(dx_values[-2]/dx_values[-3]))
print("slope L-W in space", slope_lw)
#slope_lf = (np.log(errors_lax_friedrich[-1])-np.log(errors_lax_friedrich[-2]))/(np.log(dx_values[-1])-np.log(dx_values[-2]))
#print("slope L-F in space", slope_lf)



# # Plotting the order of convergence
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# plt.plot(dt_values[:-1], errors_euler_backward[:-1], marker='o', label='FTBS')
# plt.plot(dt_values[:-1], errors_lax_wendroff[:-1], marker='o', label='Lax Wendroff')
# #plt.plot(dt_values[:-1], errors_lax_friedrich[:-1], marker='o', label='Lax Friedrich')


# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('dt')
# plt.ylabel('Error')
# plt.title('Order of Convergence for Different Schemes for TIME')
# plt.legend()
# plt.grid(True)
# plt.show()

slope_euler = (np.log(errors_euler_backward[-1])-np.log(errors_euler_backward[-2]))/(np.log(dt_values[-1])-np.log(dt_values[-2]))
print("slope FTBS in time", slope_euler)
#slope_lw = (np.log(errors_lax_wendroff[-1])-np.log(errors_lax_wendroff[-2]))/(np.log(dt_values[-1])-np.log(dt_values[-2]))
slope_lw = (np.log(errors_lax_wendroff[-2]/ errors_lax_wendroff[-3]))/(np.log(dt_values[-2]/dt_values[-3]))
print("slope L- Win time", slope_lw)
slope_lf = (np.log(errors_lax_friedrich[-1])-np.log(errors_lax_friedrich[-2]))/(np.log(dt_values[-1])-np.log(dt_values[-2]))
print("slope L-F in time", slope_lf)


plt.figure(figsize=(10, 6))
plt.plot(x_vectors[0],  f(x_vectors[0]-a*tmax) , marker='o', label='FTBS exact')
plt.plot(x_vectors[0],  euler_solutions[0], marker='o', label='FTBS scheme')
#plt.plot(x_vectors[0],  lw_solutions[0], marker='o', label='Lax-Wendroff scheme')
#plt.plot(x_vectors[0],  lf_solutions[0], marker='o', label='Lax-Friedrich scheme')

# plt.plot(dx_values[:-1], errors_lax_wendroff[:-1], marker='o', label='Lax Wendroff')
# plt.plot(dx_values[:-1], errors_lax_friedrich[:-1], marker='o', label='Lax Friedrich')

plt.xlabel('x')
plt.ylabel('Error')
plt.title('Order of Convergence for Different Schemes for SPACE')
plt.legend()
plt.grid(True)
plt.show()


# plt.figure(figsize=(10, 6))
# #plt.plot(x_vectors[2],  np.exp(-2*(x_vectors[2]-4-a*tmax)**2), marker='o', label='FTBS exact')
# plt.plot(x_vectors[2], f(x_vectors[2]-a*tmax)  , marker='o', label='FTBS exact')
# plt.plot(x_vectors[2],  euler_solutions[2], marker='o', label='FTBS scheme')
# plt.plot(x_vectors[2],  lw_solutions[2], marker='o', label='Lax-Wendroff scheme')


# # plt.plot(dx_values[:-1], errors_lax_wendroff[:-1], marker='o', label='Lax Wendroff')
# # plt.plot(dx_values[:-1], errors_lax_friedrich[:-1], marker='o', label='Lax Friedrich')

# plt.xlabel('dx')
# plt.ylabel('Error')
# plt.title('Order of Convergence for Different Schemes for SPACE')
# plt.legend()
# plt.grid(True)
# plt.show()

