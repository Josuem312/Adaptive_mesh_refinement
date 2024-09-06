# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 16:29:55 2024

@author: josue
"""

def scheme(type_of_scheme, node, c, ratio_ref): #type_of_scheme, node, c
    u = node.value
    def f(x):
        """Asigna un valor de 1.0 para valores menores que 0.1"""
        if x < 1.5:
            f=1
        else:
            f=0
        return f
    #u=node.value
    Z1= np.zeros_like(u)
    Z2= np.zeros_like(u) # Forward time backward space (two steps)
    Z3= np.zeros_like(u)
    #Z4= np.zeros_like(u)
    if type_of_scheme==1:
        # FORWARD TIME BACKWARD SPACE
            #Formula of Q
            Z1[1:-1]= (1 - c) * u[1:-1] + c * u[:-2]
            #BOUNDARY CONDITIONS
            # if node.level>0:        
            # # else:
            # Z1[0]= f(node.x_value[0]- node.time)
            # Z1[-1]= f(node.x_value[-1]- node.time)
            # Z1[0]=u[1]
            # Z1[-1]=u[-2]
            """
            Boundary conditions
            """
            if node.level==0:
                Z1[0]=u[1]
                Z1[-1]=u[-2]
            else: 
                #### Time INTERPOLATION ###
                Z1[0]= np.round((1 - node.counter_time/ratio_ref)*node.father.value_q[0] + node.counter_time/ratio_ref*node.value_q1[0], 10)  
                Z1[-1]=np.round((1 - node.counter_time/ratio_ref)*node.father.value_q[-1] + node.counter_time/ratio_ref*node.value_q1[-1],10) 
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
            #A better way to calculate Q^2
            #Z2[2:-2]=(1 - c) * Z1[2:-2] + c * Z1[1:-3]
            
            # Formula of Q_2h
            Z3[2:-2]= (1 - c) * u[2:-2] + c * u[:-4]
            return Z1, Z2, Z3, 1, 1
    elif type_of_scheme==2:
        # Lax-Friedrich Advection
            #Formula of Q
            Z1[1:-1] = (u[:-2] + u[2:]) / 2.0 - c * (u[2:] - u[:-2]) / 2.0
            #Neumann boundary conditions
            Z1[0]=Z1[1]
            Z1[-1]=Z1[-2]
            #Formula of Q^2
            #Z2[1:-1] = (Z1[:-2] + Z1[2:]) / 2.0 - c * (Z1[2:] - Z1[:-2]) / 2.0
            Z2[2:-2] = (Z1[1:-3] + Z1[3:1]) / 2.0 - c * (Z1[3:1] - Z1[1:-3]) / 2.0
            
            #Formula of Q_2h
            Z3[2:-2] = (u[:-4] + u[4:]) / 2.0 - c * (u[4:] - u[:-4]) / 2.0
            return Z1, Z2, Z3, 1, 1
    ##  HERE we need to add the other finite difference schemes 
    else:
            print("invalid number")
            return 0
        