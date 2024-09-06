
import numpy as np
import matplotlib.pyplot as plt
from Functions.Space_time_discretization import space_time
from Functions.Space_time_discretization import space_time_v2
from Simpson_rule import l2_norm_error_simpson



from decimal import Decimal, getcontext

# Configurar la precisión deseada

##############################################
    #INITIAL CONDITION
##############################################
def f(X):
    """Asigna un valor de 1.0 para valores menores que 0.1"""
    f = np.zeros_like(X)
    f[np.where(X <= 1.4)] = 1.0
    # if X <= 1.4 :
    #     f=1
    # else:
    #     f=0
    #f = np.exp(-2*(X-4)**2)
    # print(xmax)
    return f
###############################################
# Functions for finite difference schemes for the transport equation
###############################################
def euler_backward(u, c):
    Z1= np.copy(u)
    Z1[1:-1]= (1 - c) * u[1:-1] + c * u[:-2]
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
    Z1[1:-1]= c/2.0*(1+c)*u[:-2] + (1-c**2)*u[1:-1] - c/2.0*(1-c)*u[2:]
    ### NEULANN BOUNDARY CONDITIONS ###
    Z1[0]=Z1[1]
    Z1[-1]=Z1[-2]
    return  Z1

def Lax_friedrich(u, c):
    Z1= np.copy(u)
    Z1[1:-1] = 0.5*(u[:-2] + u[2:])- c *0.5* (u[2:] - u[:-2]) 
    ### NEULANN BOUNDARY CONDITIONS ###
    Z1[0]=Z1[1]
    Z1[-1]=Z1[-2]
    return  Z1



# Example usage
    
a = 1  # wave velocity
tmin, tmax = 0.0, 1.0   #3.125e-05  # initial and final time
xmin, xmax = 0.0, 10.0  # 1D domain
c =np.array([0.01, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])# before was 1   #np.array([0.2, 0.4, 0.6, 0.8, 1])  # Courant number, we need c<=1 for stability


#Nx_values = np.array([2**k for k in range(8, 16)]) # Different spatial resolutions
Nx_values = 4000

errors_euler_backward_fxs = []
errors_lax_wendroff_fxs =[]
errors_lax_friedrich_fxs =[]
euler_solutions_fxs= []
lw_solutions_fxs=[]
lf_solutions_fxs=[]
dx_values_fxs =[]
dt_values_fxs =[]
x_vectors_fxs =[]



MSE_error_lw=[]
for c in c:
    # Discretisaiton of time and space
    x_fxs, dx_fxs, dt_fxs, Nt_fxs, time_fxs = space_time(c = c, a=a, xmin=xmin, xmax=xmax, Nx=Nx_values, tmin=tmin, tmax=tmax)
    dx_values_fxs.append(dx_fxs)
    dt_values_fxs.append(dt_fxs)
    x_vectors_fxs.append(x_fxs)
    # Initialization
    u_initial_fxs = f(x_fxs)
    # Solution with different schemes
    u_euler_backward_fxs = u_initial_fxs.copy()
    u_lw_fxs= u_initial_fxs.copy()
    u_lf_fxs= u_initial_fxs.copy()
    
    
    
    error_lw=0
    samplePoints_lw=0
    for i, t in enumerate(time_fxs[1:]):     #(Nt_fxs+1):
        u_euler_backward_fxs = euler_backward(u_euler_backward_fxs, c)
        u_lw_fxs = Lax_wendroff(u_lw_fxs, c)
        u_lf_fxs = Lax_friedrich(u_lf_fxs, c)
        error_lw += np.sum((u_lw_fxs - f(x_fxs-a*t))**2) # add error from this timestep
        samplePoints_lw += len(u_lw_fxs)
    
    
    euler_solutions_fxs.append(u_euler_backward_fxs)
    lw_solutions_fxs.append(u_lw_fxs)
    lf_solutions_fxs.append(u_lf_fxs)
        
    # Calculate de real solution
    real_sol_fxs= f(x_fxs-a*tmax)
    # for i, obj in enumerate(real_sol_fxs):
    #     if abs(obj)<=1e-15:
    #         real_sol_fxs[i]=0.
    
    # Calculate the error with respect to the initial condition
    """
    ERRORS Fixing the spacial step
    """
    #### L2 norm error ###
    # rest_euler_fxs = u_euler_backward_fxs- real_sol_fxs
    # rest_euler_fxs = rest_euler_fxs**2
    # error_euler_fxs= np.round(np.sqrt(dx_fxs)*np.sqrt(sum(rest_euler_fxs)), 15)
    ### SIMPSON ###
    error_euler_fxs = l2_norm_error_simpson(u_euler_backward_fxs, real_sol_fxs, x_fxs)
    ### MSE ###
    

    # rest_lw_fxs = u_lw_fxs - real_sol_fxs
    # rest_lw_fxs = rest_lw_fxs**2
    # error_lw_fxs= np.round(np.sqrt(dx_fxs)*np.sqrt(sum(rest_lw_fxs)), 15)
    ### SIMPSON ###
    error_lw_fxs = l2_norm_error_simpson(u_lw_fxs, real_sol_fxs, x_fxs)
    ### MSE ###
    error_lw = np.sqrt(error_lw/samplePoints_lw) # calculate rms-error
    MSE_error_lw.append(error_lw)  



    # rest_lf_fxs = u_lf_fxs- real_sol_fxs
    # rest_lf_fxs = rest_lf_fxs**2
    # error_lf_fxs= np.round(np.sqrt(dx_fxs)*np.sqrt(sum(rest_lf_fxs)), 15)  
    ### SIMPSON ###
    error_lf_fxs = l2_norm_error_simpson(u_lf_fxs, real_sol_fxs, x_fxs)
    ### MSE ###
 
    
    
    errors_euler_backward_fxs.append(error_euler_fxs)
    errors_lax_wendroff_fxs.append(error_lw_fxs)
    errors_lax_friedrich_fxs.append(error_lf_fxs)

    
    
    
    

"""
Calculating the SLOPE
"""
slopes_fixing_space_euler= []
slopes_fixing_space_lw= []
slopes_fixing_space_lf= []
counting=[]

for i in np.arange(np.size(dt_values_fxs)-1):
    slope_fxs_euler = (np.log(errors_euler_backward_fxs[i+1])-np.log(errors_euler_backward_fxs[i]))/(np.log(dt_values_fxs[i+1])-np.log(dt_values_fxs[i]))
    slope_fxs_lw = (np.log(errors_lax_wendroff_fxs[i+1])-np.log(errors_lax_wendroff_fxs[i]))/(np.log(dt_values_fxs[i+1])-np.log(dt_values_fxs[i]))
    slope_fxs_lf = (np.log(errors_lax_friedrich_fxs[i+1])-np.log(errors_lax_friedrich_fxs[i]))/(np.log(dt_values_fxs[i+1])-np.log(dt_values_fxs[i]))

    slopes_fixing_space_euler.append(slope_fxs_euler)
    slopes_fixing_space_lw.append(slope_fxs_lw)
    slopes_fixing_space_lf.append(slope_fxs_lf)
    counting.append(i+1)
    
SLOPE_fixing_space_euler= np.mean(slopes_fixing_space_euler)
SLOPE_fixing_space_lw= np.mean(slopes_fixing_space_lw)
SLOPE_fixing_space_lf= np.mean(slopes_fixing_space_lf)


"""
GRAPH 
"""

# Plotting the order of convergence
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(dt_values_fxs, errors_euler_backward_fxs, marker='o', label='FTBS')
plt.plot(dt_values_fxs, errors_lax_wendroff_fxs, marker='o', label='Lax Wendroff')
plt.plot(dt_values_fxs, errors_lax_friedrich_fxs, marker='o', label='Lax Friedrich')


plt.xscale('log')
plt.yscale('log')
plt.xlabel('dt')
plt.ylabel('Error')
plt.title('Order of Convergence for Different Schemes for SPACE')
plt.legend()
plt.grid(True)
plt.show()

print("slope FTBS in time", SLOPE_fixing_space_euler)
print("slope L-W in time", SLOPE_fixing_space_lw)
print("slope L-F in time", SLOPE_fixing_space_lf)





plt.figure(figsize=(10, 6))
#plt.plot(x_vectors_fxs[0],  f(x_vectors_fxs[0]) , marker='o', label='FTBS exact')
#plt.plot(x_vectors_fxs[0],  f(x_vectors_fxs[0]-a*tmax) , marker='o', label='FTBS exact')
plt.plot(x_vectors_fxs[0], np.round( abs( euler_solutions_fxs[0]-f(x_vectors_fxs[0]-a*tmax)), 15), marker='o', label='FTBS scheme')
plt.plot(x_vectors_fxs[0], np.round( abs(lw_solutions_fxs[0]-f(x_vectors_fxs[0]-a*tmax))), marker='o', label='Lax-Wendroff scheme')
#plt.plot(x_vectors_fxs[0],  lf_solutions_fxs[0]-f(x_vectors_fxs[0]-a*tmax), marker='o', label='Lax-Friedrich scheme')

# plt.xscale('log')
plt.yscale('log')
plt.xlabel('x')
plt.ylabel('Sol')
plt.title('Sol for SPACE')
plt.legend()
plt.grid(True)
plt.show()

REAL_sol= f(x_vectors_fxs[0]-a*tmax)
# for i, obj in enumerate(REAL_sol):
#     if abs(obj)<=1e-15:
#         REAL_sol[i]=0.

print("Error Euler:", abs(euler_solutions_fxs[0]- REAL_sol))
Error_euler=  np.round(abs(euler_solutions_fxs[0]- REAL_sol), 15)
#print("Error LW:", abs(lw_solutions_fxs[0]-f(x_vectors_fxs[0]-a*tmax)))



# Plotting the MSE ERROR
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(dt_values_fxs, MSE_error_lw, marker='o', label='lax WENDROFF')
plt.yscale('log')
plt.xlabel('dt')
plt.ylabel('Error')
plt.title('MSE ERROR')
plt.legend()
plt.grid(True)
plt.show()


"""
Graphing the slope 
"""
    

# Plotting the order of convergence
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(counting, slopes_fixing_space_euler, marker='o', label='FTBS')
plt.plot(counting, slopes_fixing_space_lw, marker='o', label='LW')
plt.plot(counting, slopes_fixing_space_lf, marker='o', label='LF')

plt.xlabel('index')
plt.ylabel('log_2 of ration')
plt.title('Order of convergence for different advection schemes')
plt.legend()
plt.grid(True)
plt.show()
