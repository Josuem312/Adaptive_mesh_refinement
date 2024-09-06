"""
Created on Fri Jul 19 10:23:14 2024

@author: josue
"""
import numpy as np
from Functions.Auxiliary_functions import compare_vectors_with_indices

""""
In this function we initialize into the new subgrid: ###
 1. Creates a vector of same dim of subgrid and a list which stores the indexes
 which are taken into account for the new vector initialization
 2. Creates a number equal to the initial index of flag vector
 3. By using a for sentence and considering the number of stencil we perform
# an interpolation.
"""

def u_finer_init(u, subgrid, flag):#    Note that here we consider stencil=1
    u0_fine=np.zeros_like(subgrid)
    counter=flag[0]
    l=[]
    for i in np.arange(np.size(u0_fine)):
        if i % 2==0:
            u0_fine[i]=0.5*(u[counter-1] + u[counter])
            #print(u[counter-1], u[counter], 'flags are', counter-1, counter)   
            l.append([counter-1, counter])
        else:
            u0_fine[i]=u[counter]
            #print(i, "flag is" , counter)
            l.append(counter)
            counter=counter+1
    return u0_fine, l



""" THIS IS THE PRINCIPAL WE ARE USING
## In this function we initialize into the new subgrid: ###
1. Creates a vector of same dim of subgrid and a list which stores the indexes
which are taken into account for the new vector initialization
2. Creates a number equal to the initial index of flag vector
3. Creates a vector of coefficients for interpolation
4. Depending on the strencil we will start by filling the boundary values
separately
5. By using a count variable we will fill out the rest N-2 components:
    a) if the count=0 means we are in the coarse grid point so we increase
    the variable counter which takes into account flags (indexes) values
    b) else we are in the interpolation part so we use the vector
    coefficients by using the coordinates count 
    c) Finally we restart count when we get to de ratio_ref
"""


#def u_finer_init_new(u, subgrid, flag, ratio_ref, stencil):#    Note that here we consider stencil=1
def u_finer_init_new(u, subgrid, flag, ratio_ref, stencil):
    u0_fine=np.zeros_like(subgrid)
    counter=flag[0]
    l=[]
    # Create the coefficients for the interpolation
    coeficients= np.linspace(0,1, ratio_ref+1)[1:-1]
    if stencil==1:
        ## Interpolation of boundary ###
        u0_fine[0]= (1-coeficients[-1])*u[counter-1]+coeficients[-1]*u[counter]
        l.append([counter-1, counter])
        count=0
        #for i, coef in enumerate(coeficients):
        for i in np.arange(np.size(u0_fine[1:-1])):  
            if count==0:
                u0_fine[i+1]=u[counter]
                l.append(counter)
                counter=counter+1
            else:
                u0_fine[i+1]=(1-coeficients[count-1])*(u[counter-1]) + coeficients[count-1]*u[counter]
                #print(u[counter-1], u[counter], 'flags are', counter-1, counter)   
                l.append([counter-1, counter])
            count=count+1
            if count==ratio_ref:
                count=0
        u0_fine[-1]= (1-coeficients[0])*u[flag[-1]]+coeficients[0]*u[flag[-1]+1] 
        l.append([flag[-1], flag[-1]+1])
                
    else:
        print("ERROR: Invalid stencil number")
    return u0_fine, l



def u_finer_init_new_node(node, subgrid, flag, ratio_ref, stencil):
    u_finer_initial=np.zeros_like(subgrid)
    counter=flag[0]
    l=[]
    # Create the coefficients for the interpolation
    coeficients= np.linspace(0,1, ratio_ref+1)[1:-1]
    if stencil==1:
        ## Interpolation of boundary ###
        u_finer_initial[0]= (1-coeficients[-1])*node.value[counter-1]+coeficients[-1]*node.value[counter]
        l.append([counter-1, counter])
        count=0
        #for i, coef in enumerate(coeficients):
        for i in np.arange(np.size(u_finer_initial[1:-1])):  
            if count==0:
                u_finer_initial[i+1]=node.value[counter]
                l.append(counter)
                counter=counter+1
            else:
                u_finer_initial[i+1]=(1-coeficients[count-1])*(node.value[counter-1]) + coeficients[count-1]*node.value[counter]
                #print(u[counter-1], u[counter], 'flags are', counter-1, counter)   
                l.append([counter-1, counter])
            count=count+1
            if count==ratio_ref:
                count=0
        u_finer_initial[-1]= (1-coeficients[0])*node.value[flag[-1]]+coeficients[0]*node.value[flag[-1]+1] 
        l.append([flag[-1], flag[-1]+1])
                
    else:
        print("ERROR: Invalid stencil number")
    
    '''
    Here we plan to construct the regridding process
    1 First we add an algorithm to keep information from old nodes
    '''
    if node.time!=0 and node.hasChildren():
        x_value_children=node.children[0].x_value
        if np.array_equal(subgrid , x_value_children):
            u_finer_initial= node.children[0].value
        else:
            matches= compare_vectors_with_indices(subgrid, x_value_children)
            if len(matches)!=0:
                for obj in matches:
                    u_finer_initial[obj[1]]= node.children[0].value[obj[2]]
     
    return u_finer_initial, l


def u_finer_init_new_node_version_2(node, subgrid, flag, ratio_ref, stencil):
    u_finer_initial=np.zeros_like(subgrid)
    counter=flag[0]
    l=[]
    # Create the coefficients for the interpolation
    coeficients1=np.linspace(0,1, ratio_ref+1)
    coeficients= coeficients1[1:-1]
    if stencil==1:
        ## Interpolation of boundary ###
        # u0_fine[0]= (1-coeficients[-1])*node.value[counter-1]+coeficients[-1]*node.value[counter]
        # l.append([counter-1, counter])
        u_finer_initial[0]= (1-coeficients1[-3])*node.value[counter-1]+coeficients1[-3]*node.value[counter]
       # l.append([counter-1, counter])
        u_finer_initial[1]= (1-coeficients1[-2])*node.value[counter-1]+coeficients1[-2]*node.value[counter]
       # l.append([counter-1, counter])
        count=0
        #for i, coef in enumerate(coeficients):
        for i in np.arange(np.size(u_finer_initial[2:-2])):  
            if count==0:
                u_finer_initial[i+2]=node.value[counter]
                l.append(i+2)
                counter=counter+1
            else:
                u_finer_initial[i+2]=(1-coeficients[count-1])*(node.value[counter-1]) + coeficients[count-1]*node.value[counter]
                #print(u[counter-1], u[counter], 'flags are', counter-1, counter)   
                #l.append([counter-1, counter])
            count=count+1
            if count==ratio_ref:
                count=0
        # u0_fine[-1]= (1-coeficients[0])*node.value[flag[-1]]+coeficients[0]*node.value[flag[-1]+1] 
        # l.append([flag[-1], flag[-1]+1])
        u_finer_initial[-2]= (1-coeficients1[1])*node.value[flag[-1]]+coeficients1[1]*node.value[flag[-1]+1]
       # l.append([flag[-1], flag[-1]+1])    
        u_finer_initial[-1]= (1-coeficients1[2])*node.value[flag[-1]]+coeficients1[2]*node.value[flag[-1]+1]
       # l.append([flag[-1], flag[-1]+1])    
    else:
        print("ERROR: Invalid stencil number")
        '''
        Here we plan to construct the regridding process
        1 First we add an algorithm to keep information from old nodes
        '''
                    
    if node.time!=0 and node.hasChildren():
        x_value_children=node.children[0].x_value
        if np.array_equal(subgrid , x_value_children):
            u_finer_initial= node.children[0].value
        else:
            matches= compare_vectors_with_indices(subgrid, x_value_children)
            if len(matches)!=0:
                for obj in matches:
                    u_finer_initial[obj[1]]= node.children[0].value[obj[2]]
                    
     
    return u_finer_initial, l


def u_finer_init_new_node_version_2_2(node, subgrid, flag, ratio_ref, stencil):
    u_finer_initial=np.zeros_like(subgrid)
    counter=flag[0]
    #l=[]
    # Create the coefficients for the interpolation
    coeficients1=np.linspace(0,1, ratio_ref+1)
    coeficients= coeficients1[1:-1]
    if stencil==1:
        ## Interpolation of boundary ###
        # u0_fine[0]= (1-coeficients[-1])*node.value[counter-1]+coeficients[-1]*node.value[counter]
        # l.append([counter-1, counter])
        u_finer_initial[0]= (1-coeficients1[-3])*node.value_q1[counter-1]+coeficients1[-3]*node.value_q1[counter]
        #l.append([counter-1, counter])
        u_finer_initial[1]= (1-coeficients1[-2])*node.value_q1[counter-1]+coeficients1[-2]*node.value_q1[counter]
        #l.append([counter-1, counter])
        count=0
        #for i, coef in enumerate(coeficients):
        for i in np.arange(np.size(u_finer_initial[2:-2])):  
            if count==0:
                u_finer_initial[i+2]=node.value_q1[counter]
                #l.append(i+2)
                counter=counter+1
            else:
                u_finer_initial[i+2]=(1-coeficients[count-1])*(node.value_q1[counter-1]) + coeficients[count-1]*node.value_q1[counter]
                #print(u[counter-1], u[counter], 'flags are', counter-1, counter)   
                #l.append([counter-1, counter])
            count=count+1
            if count==ratio_ref:
                count=0
        # u0_fine[-1]= (1-coeficients[0])*node.value[flag[-1]]+coeficients[0]*node.value[flag[-1]+1] 
        # l.append([flag[-1], flag[-1]+1])
        u_finer_initial[-2]= (1-coeficients1[1])*node.value_q1[flag[-1]]+coeficients1[1]*node.value_q1[flag[-1]+1]
        #l.append([flag[-1], flag[-1]+1])    
        u_finer_initial[-1]= (1-coeficients1[2])*node.value_q1[flag[-1]]+coeficients1[2]*node.value_q1[flag[-1]+1]
        #l.append([flag[-1], flag[-1]+1])    
    else:
        print("ERROR: Invalid stencil number")
        '''
        Here we plan to construct the regridding process
        1 First we add an algorithm to keep information from old nodes
        '''
                    
    if node.time!=0 and node.hasChildren():
        x_value_children=node.children[0].x_value
        if np.array_equal(subgrid , x_value_children):
            u_finer_initial= node.children[0].value
        else:
            matches= compare_vectors_with_indices(subgrid, x_value_children)
            if len(matches)!=0:
                for obj in matches:
                    u_finer_initial[obj[1]]= node.children[0].value[obj[2]]
                    
     
    return u_finer_initial
"""
 In this function we simply solve over the subgrid but we do not use this function
"""
def solving_subgrid(scheme, u_finer_initial,ratio_ref, subgrid, c):
    u_finer=np.zeros((ratio_ref, np.size(subgrid)))  
    if scheme==1:
        Z= np.copy(u_finer_initial)
        for i in np.arange(ratio_ref):
            # FORWARD TIME BACKWARD SPACE
                Z[1:-1]= (1 - c) * Z[1:-1] + c * Z[:-2]
                #BOUNDARY CONDITIONS (Neumann)
                Z[0]=Z[1]
                Z[-1]=Z[-2]
                u_finer[i, :]=Z
        return u_finer
    else:
        print("incorrect value")
        return 0

'''
Here are going to create a funcion which allow us to initialize solution values
to new subgrids
'''    

