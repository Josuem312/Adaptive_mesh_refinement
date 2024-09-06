
import numpy as np
import matplotlib.pyplot as plt
from Functions.Space_time_discretization import space_time
from Functions.Space_time_discretization import space_time_v2

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

Nt_values = 81



errors_euler_backward_fxt = []
errors_lax_wendroff_fxt =[]
errors_lax_friedrich_fxt =[]
euler_solutions_fxt= []
lw_solutions_fxt=[]
lf_solutions_fxt=[]
dx_values_fxt =[]
dt_values_fxt =[] ## ONLY ONE 
t_vectors_fxt =[] 
x_vectors_fxt =[]


for c in c:
    
    """
    TIME
    """
    # Discretisaiton of time and space
    t_fxt, dx_fxt, dt_fxt, Nx_fxt, x_fxt  = space_time_v2(c = c, a=a, xmin=xmin, xmax=xmax, Nt=Nt_values, tmin=tmin, tmax=tmax)
    dx_values_fxt.append(dx_fxt)
    dt_values_fxt.append(dt_fxt)
    x_vectors_fxt.append(x_fxt)
    t_vectors_fxt.append(t_fxt)
    # Initialization
    u_initial_fxt = f(x_fxt)
    # Solution with different schemes
    u_euler_backward_fxt = u_initial_fxt.copy()
    u_lw_fxt= u_initial_fxt.copy()
    u_lf_fxt= u_initial_fxt.copy()
    #print (Nt_fxs,Nt_values)

    for _ in range(Nt_values+1):
        """
        TIME
        """
        u_euler_backward_fxt = euler_backward(u_euler_backward_fxt, c)
        u_lw_fxt = Lax_wendroff(u_lw_fxt, c)
        u_lf_fxt = Lax_friedrich(u_lf_fxt, c)

    """
    TIME
    """
    euler_solutions_fxt.append(u_euler_backward_fxt)
    lw_solutions_fxt.append(u_lw_fxt)
    lf_solutions_fxt.append(u_lf_fxt)
        
    # Calculate de real solution
    real_sol_fxt= f(x_fxt-a*tmax)   
    
    # Calculate the error with respect to the initial condition
    """
    ERRORS Fixing the time step
    """
    rest_euler_fxt = u_euler_backward_fxt- real_sol_fxt
    rest_euler_fxt = rest_euler_fxt**2
    error_euler_fxt= np.round(np.sqrt(dx_fxt)*np.sqrt(sum(rest_euler_fxt)), 15)

    rest_lw_fxt = u_lw_fxt - real_sol_fxt
    rest_lw_fxt = rest_lw_fxt**2
    error_lw_fxt= np.round(np.sqrt(dx_fxt)*np.sqrt(sum(rest_lw_fxt)), 15)

    rest_lf_fxt = u_lf_fxt- real_sol_fxt
    rest_lf_fxt = rest_lf_fxt**2
    error_lf_fxt= np.round(np.sqrt(dx_fxt)*np.sqrt(sum(rest_lf_fxt)), 15)  
    
    errors_euler_backward_fxt.append(error_euler_fxt)
    errors_lax_wendroff_fxt.append(error_lw_fxt)
    errors_lax_friedrich_fxt.append(error_lf_fxt)

"""
Calculating the SLOPE TIME
"""
slopes_fixing_time_euler= []
slopes_fixing_time_lw= []
slopes_fixing_time_lf= []

for i in np.arange(np.size(dt_values_fxt)-1):
    slope_fxt_euler = (np.log(errors_euler_backward_fxt[i])-np.log(errors_euler_backward_fxt[i+1]))/(np.log(dx_values_fxt[i])-np.log(dx_values_fxt[i+1]))
    slope_fxt_lw = (np.log(errors_lax_wendroff_fxt[i])-np.log(errors_lax_wendroff_fxt[i+1]))/(np.log(dx_values_fxt[i])-np.log(dx_values_fxt[i+1]))
    slope_fxt_lf = (np.log(errors_lax_friedrich_fxt[i])-np.log(errors_lax_friedrich_fxt[i+1]))/(np.log(dx_values_fxt[i])-np.log(dx_values_fxt[i+1]))

    slopes_fixing_time_euler.append(slope_fxt_euler)
    slopes_fixing_time_lw.append(slope_fxt_lw)
    slopes_fixing_time_lf.append(slope_fxt_lf)
    
SLOPE__fixing_time_euler= np.mean(slopes_fixing_time_euler)
SLOPE__fixing_time_lw= np.mean(slopes_fixing_time_lw)
SLOPE__fixing_time_lf= np.mean(slopes_fixing_time_lf)



"""
GRAPH TIME 
"""
# Plotting the order of convergence
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(dx_values_fxt[:-1], errors_euler_backward_fxt[:-1], marker='o', label='FTBS')
plt.plot(dx_values_fxt[:-1], errors_lax_wendroff_fxt[:-1], marker='o', label='Lax Wendroff')
plt.plot(dx_values_fxt[:-1], errors_lax_friedrich_fxt[:-1], marker='o', label='Lax Friedrich')


plt.xscale('log')
plt.yscale('log')
plt.xlabel('dx')
plt.ylabel('Error')
plt.title('Order of Convergence for Different Schemes for TIME')
plt.legend()
plt.grid(True)
plt.show()

slope_euler = (np.log(errors_euler_backward_fxt[-1])-np.log(errors_euler_backward_fxt[-2]))/(np.log(dx_values_fxt[-1])-np.log(dx_values_fxt[-2]))
print("slope FTBS in space", SLOPE__fixing_time_euler)
#slope_lw = (np.log(errors_lax_wendroff[-1])-np.log(errors_lax_wendroff[-2]))/(np.log(dx_values[-1])-np.log(dx_values[-2]))
slope_lw = (np.log(errors_lax_wendroff_fxt[-1]/ errors_lax_wendroff_fxt[-2]))/(np.log(dx_values_fxt[-1]/dx_values_fxt[-2]))
print("slope L-W in space", SLOPE__fixing_time_lw)
slope_lf = (np.log(errors_lax_friedrich_fxt[-1])-np.log(errors_lax_friedrich_fxt[-2]))/(np.log(dx_values_fxt[-1])-np.log(dx_values_fxt[-2]))
print("slope L-F in space", SLOPE__fixing_time_lf)




# plt.figure(figsize=(10, 6))
# plt.plot(x_vectors_fxt[0],  f(x_vectors_fxt[0]-a*tmax) , marker='o', label='FTBS exact')
# plt.plot(x_vectors_fxt[0],  euler_solutions_fxt[0], marker='o', label='FTBS scheme')
# plt.plot(x_vectors_fxt[0],  lw_solutions_fxt[0], marker='o', label='Lax-Wendroff scheme')
# plt.plot(x_vectors_fxt[0],  lf_solutions_fxt[0], marker='o', label='Lax-Friedrich scheme')

# plt.xlabel('x')
# plt.ylabel('Sol')
# plt.title('Sol for TIME')
# plt.legend()
# plt.grid(True)
# plt.show()
