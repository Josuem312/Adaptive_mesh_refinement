# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 13:36:31 2024

@author: josue
"""

import numpy as np
import matplotlib.pyplot as plt
from Functions.Space_time_discretization import space_time

# Parámetros
Nx = 160  # número de puntos espaciales
c = 1.0   # velocidad de advección
a = 1  # wave velocity
tmin, tmax = 0.0, 1.5  #3.125e-05  # initial and final time
xmin, xmax = 0.0, 10.0  # 1D domain
c = 0.4  # Courant number, we need c<=1 for stability

x, dx, dt, Nt, time = space_time(c=c, a=a, xmin=xmin, xmax=xmax, Nx=Nx, tmin=tmin, tmax=tmax)
time=time[1:]



### INITIAL CONDITION #####
###########################
def f(x):
    """Asigna un valor de 1.0 para valores menores que 0.1"""
    f=  np.exp(-200*(x + 0.5)**2) + 0.1*np.sin(2*np.pi*(x+0.5))
    return f
u = f(x)


for i, obj in enumerate(time):
    u[1:-1]  = (u[:-2] + u[2:]) / 2.0 - c * (u[2:] - u[:-2]) / 2.0 
    u[0]=f( xmin - a*obj) ## We suppose that a=1 and xmax=4
    u[-1]= f(xmax - a*obj)   

# Array para almacenar la solución en cada paso temporal


# Visualización del resultado
plt.figure(figsize=(10, 6))

plt.plot(x, u, label=f'Aprixumated solution at Time = {tmax}')
plt.plot(x, f(x-a*tmax), label=f'real solution at Time = {tmax}')
#plt.plot(x, f(x-2), label='Proof')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Graphic')
plt.legend()
plt.grid(True)

plt.show()
