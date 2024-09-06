from Functions.Auxiliary_functions import ary
"""
Created on Thu Jul 18 13:44:44 2024

@author: josue
"""
import numpy as np
from Functions.Flagging import flagging_with_subgriddings
from Functions.Interpolation import u_finer_init_new
from Functions.Interpolation import u_finer_init_new_node
from Functions.Interpolation import u_finer_init_new_node_version_2
from Functions.Interpolation import u_finer_init_new_node_version_2_2
from Functions.Initial_conditions import init
from Functions.Auxiliary_functions import compare_vectors_with_indices
from Functions.Tree import Node

"""
This functions creates a new finer subgrid, i.e. creates a vector
whose components are the spacial discretisation of the new subgrid.
NOTE: we need to consider the case there different subgrids are needed
Morever here we add the so called "GHOST CELLS"
"""


def new_finer_subgrid_node(node,  flag, stencil, ratio_ref):
    #print("node", node.x_value, "flag", node.x_value[flag])
    if node.father != None:
        if node.x_value[flag[0]]==node.x_value[0]:# or node.x_value[flag[0]]==node.x_value[0]:
            print("We arrive to the left limit at", node.time)
        if node.x_value[flag[-1]]==node.x_value[-1]: # or node.x_value[flag[-1]]==node.x_value[-2]:
            print("We arrive to the right limit at", node.time)
    new_flag=flag
    new_step= node.dx/ratio_ref  
    x_min=node.x_value[new_flag[0]]- 2*stencil*new_step# before 2
    x_max=node.x_value[new_flag[-1]]+ 2*stencil*new_step# before 2
    subgrid= ary(x_min, x_max, new_step)## CHANGE NAME ARY
    
    return subgrid
"""
Creating of subgrid and initialization of solution vector U
"""
    
def subgridding_process(node, flag, stencil, ratio_ref, 
                type_initial_condition): 
    
#Creation of the subgrid
        subgrid=new_finer_subgrid_node(node, flag, stencil, ratio_ref)
        
#Initialization of the new vector solution
        # if node.time==0:
        #     u_finer_initial= initial_conditions(type_initial_condition, subgrid)
        # else:
        #     u_finer_initial, l= u_finer_init_new_node_version_2(node, subgrid, flag, ratio_ref, stencil)
        
        u_finer_initial, l= u_finer_init_new_node_version_2(node, subgrid, flag, ratio_ref, stencil)
        if node.time==0:
            u_finer_initial= init(type_initial_condition, subgrid)
            
# Updating or creating a new node data structure
        ## Careful about having more than one child
        if node.hasChildren():
            node.children[0].value= u_finer_initial
            node.children[0].center= 0.5*(subgrid[0]+ subgrid[-1])
            node.children[0].flag= flag
            node.children[0].x_value= subgrid            
        else:
            new_level_node= Node(center= 0.5*(subgrid[0]+ subgrid[-1]), dx=node.dx/ratio_ref,
                                 dt= node.dt/ratio_ref,
                                 time=0,
                                 counter_time=0,
                                 flag=flag ,value=u_finer_initial, 
                                 value_q= None,
                                 value_q1= None,
                                 x_value=subgrid, 
                                 max_level=node.max_level, level= node.level+1, father= node)
            node.children= [new_level_node]
        return l
