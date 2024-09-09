#scheme.py
"""
Created on Thu Jul 18 12:00:01 2024

@author: josue
"""
import numpy as np
"""
In this function we apply the differen finite difference schemes###
The parameters are the type of scheme, the vector of previous solution u and
 the courant number c
""" 
def f(type_initial_condition, x):
    if type_initial_condition ==1:
        y= np.exp(-200*(x+0.5)**2) + 0.1*np.sin(2*np.pi*(x+0.5))
    elif type_initial_condition ==2:
        if x < 1.4:
            y=1
        else:
            y=0
    elif type_initial_condition==3:
        y = np.exp(-(x-4)**2)
    elif type_initial_condition==4:
        y= np.maximum(0.5 - np.abs(x - 2) / 0.5, 0)
    elif type_initial_condition==5:
        if x < 0.4 and x > 0.2:
            y= 50 * (0.01 - (x - 0.3)**2) * np.exp(0.01 / ((x - 0.3)**2 - 0.01))
        else:
            y=0
    else:
        print("No available value")
    return y

def scheme(type_of_scheme, node, c, ratio_ref,
           type_initial_condition): #type_of_scheme, node, c
    u = node.value

    u=node.value
    Z1= np.zeros_like(u)
    Z2= np.zeros_like(u) # Forward time backward space (two steps)
    Z3= np.zeros_like(u)
    #Z4= np.zeros_like(u)
    if type_of_scheme==1:
        # FORWARD TIME BACKWARD SPACE
            #Formula of Q
            Z1[1:-1]= np.round((1 - c) * u[1:-1] + c * u[:-2],15)
            
            """
            Boundary conditions
            """
            if node.level==0:
                if type_initial_condition==1:
                    Z1[0]=f(type_initial_condition, node.x_value[0] - node.time) ## We suppose that a=1 and xmax=4
                    Z1[-1]=f(type_initial_condition, node.x_value[-1] - node.time)   
                else:
                    Z1[0]=Z1[1]
                    Z1[-1]=Z1[-2]
            else: 
                #### Time INTERPOLATION ###
                Z1[0]= np.round((1 - node.counter_time/ratio_ref)*node.father.value_q[0] + node.counter_time/ratio_ref*node.value_q1[0], 15)  
                Z1[-1]=np.round((1 - node.counter_time/ratio_ref)*node.father.value_q[-1] + node.counter_time/ratio_ref*node.value_q1[-1],15) 
                node.counter_time=node.counter_time + 1
                if node.counter_time == ratio_ref:
                    node.counter_time = 0
            # Test G0 or not
            # If G0: use the usual boundary condition
            # If NOT G0: 
            
            #Another way to calculate Q^2 
            Z2[1:-1]= (1 - c) * Z1[1:-1] + c * Z1[:-2]
            Z2[1]=0
            Z2[-2]=0
            
            # Formula of Q_2h
            Z3[2:-2]= (1 - c) * u[2:-2] + c * u[:-4]
            order=1
            stencil=1


    elif type_of_scheme==2:
    #Lax Wendrof scheme
            #Formula of Q
            Z1[1:-1]= np.round(u[1:-1] - 0.5 * c * (u[2:] - u[:-2]) + 0.5 * c**2 * (u[2:] - 2*u[1:-1] + u[:-2]),15)
            
            """
            Boundary conditions
            """
            if node.level==0:
                if type_initial_condition==1:
                    Z1[0]=f(type_initial_condition, node.x_value[0] - node.time) ## We suppose that a=1 and xmax=4
                    Z1[-1]=f(type_initial_condition, node.x_value[-1] - node.time)   
                else:
                    Z1[0]=Z1[1]
                    Z1[-1]=Z1[-2]          
            else: 
                #### Time INTERPOLATION ###
                Z1[0]= np.round((1 - node.counter_time/ratio_ref)*node.father.value_q[0] + node.counter_time/ratio_ref*node.value_q1[0], 10)  
                Z1[-1]= np.round((1 - node.counter_time/ratio_ref)*node.father.value_q[-1] + node.counter_time/ratio_ref*node.value_q1[-1],10) 
                node.counter_time=node.counter_time + 1
                if node.counter_time == ratio_ref:
                    node.counter_time = 0
            # Test G0 or not
            # If G0: use the usual boundary condition
            # If NOT G0: 
            
           
            #Another way to calculate Q^2 
            Z2[1:-1]= Z1[1:-1] - 0.5 * c * (Z1[2:] - Z1[:-2]) + 0.5 * c**2 * (  Z1[2:] - 2*Z1[1:-1] + Z1[:-2])
            Z2[1]=0
            Z2[-2]=0
            #A better way to calculate Q^2
            #Z2[2:-2]=(1 - c) * Z1[2:-2] + c * Z1[1:-3]
            
            # Formula of Q_2h
            Z3[2:-2]= u[2:-2] - 0.5 * c * (u[4:] - u[:-4]) + 0.5 * c**2 * (u[4:] - 2*u[2:-2] + u[:-4])
            order=2
            stencil=1
       
            
    elif type_of_scheme==3:
        # Lax-Friedrich Advection
            #Formula of Q
            Z1[1:-1] = (u[:-2] + u[2:]) / 2.0 - c * (u[2:] - u[:-2]) / 2.0
            #Neumann boundary conditions
            Z1[0]=Z1[1]
            Z1[-1]=Z1[-2]
            #Formula of Q^2
            #Z2[1:-1] = (Z1[:-2] + Z1[2:]) / 2.0 - c * (Z1[2:] - Z1[:-2]) / 2.0
            
            #Formula of Q_2h
            Z3[2:-2] = (u[:-4] + u[4:]) / 2.0 - c * (u[4:] - u[:-4]) / 2.0

    ##  HERE we need to add the other finite difference schemes 
    else:
            print("invalid number")
            return 0
    return Z1, Z2, Z3, order, stencil
        
 # Here the only difference is that we use the node as an INPUT       
def scheme_using_node(node, c, type_of_scheme):
    Q1, Q_power_2, Q_2h, order, stencil = scheme(type_of_scheme, node.value, c)
    return Q1, Q_power_2, Q_2h, order, stencil  

    
