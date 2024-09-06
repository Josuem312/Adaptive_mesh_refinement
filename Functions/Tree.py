# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 14:48:00 2024

@author: josue
"""

import numpy as np


"""
Here we create our tree data structure
"""
class Node():
    def __init__(self, center, dx, value, x_value, max_level, dt, #delta_t,
                 time, #new_time, 
                 counter_time=None,
                 flag=None, level=0, father=None,
             child_id=0, value_q=None, value_q1=None):
      self.center = center
      self.dx = dx
      self.xf = np.array([center-dx/2, center+dx/2])
      self.children = []
      self.value = value # scalar or vector data stored in this cell
      self.x_value= x_value
      self.level = level
      self.father= father
      self.is_destroyed = False
      self.child_id = child_id
      self.child_deletion_agreement = [False, False]
      self.max_level = max_level
      self.dt = dt # Time step 
      #self.delta_t= delta_t
      self.flag= flag
      self.time= time
      self.value_q= value_q
      self.value_q1=value_q1
      self.counter_time= counter_time
      #self.new_time=new_time
#      self.value_with_ghost= value_with_ghost
#      self.x_value_with_ghost= x_value_with_ghost
      #self.min_level = min_level
      if father is None:
          assert level==0, 'only 0-th level nodes may have no father'
          
    def find_node(self, target_level, target_index):
        if self.level == target_level and self.child_id == target_index:
            return self
        for child in self.children:
            result = child.find_node(target_level, target_index)
        if result is not None:
            return result
        return None
      
    def hasChildren(self):
        return len(self.children)>0
    

  # def destroy(self):
  #   if self.level==self.min_level:
  #     return # we cannot delete level-0 nodes
  #   # if not self.is_destroyed:
  #   self.father.require_lineage_destruction(self.child_id)
  #     # self.is_destroyed = True
      
  # def __del__(self):
  #   # print('being destroyed')
  #   pass
  #   # assert self.is_destroyed
  #   # assert not self.father.hasChildren()
    
  # def reset_lineage_destruction(self):
  #   """ Reset destruction counter, i.e. if not all children have been asked to be
  #   destroyed, we keep them and reset the requirement for destruction """
  #   self.child_deletion_agreement = [False, False]
    
  
  # def recursive_destruction_reset(self):
  #   assert not self.is_destroyed
  #   self.reset_lineage_destruction()
  #   for c in self.children:
  #     c.recursive_destruction_reset()
    
  # def require_lineage_destruction(self, child_id):
  #   """ Destroy children leafs and project values conservatively """
  #   if not self.hasChildren():
  #     raise Exception('weird')
  #   self.child_deletion_agreement[child_id] = True
    
  #   if all(self.child_deletion_agreement):
  #     # both children have asked for, we can destroy them
  #     # we recover their values first
  #     self.value = 0.5*( self.children[0].value + self.children[1].value )
  #     for c in self.children:
  #       assert not c.is_destroyed
  #       c.is_destroyed = True
  #     self.children = []
    
  # def refine(self):
  #   """ Creates children """
  #   # if self.center== -9.6875:
  #     # print('here')
  #   if self.hasChildren():
  #     raise Exception('node is already refined')
      
  #   if self.level==self.max_level:
  #     # raise Exception('maximum refinment level reached')
  #     return
      
  #   # TODO: better interpolation
  #   self.children = [Node(center=self.center-self.dx/4, dx=self.dx/2, value=self.value,
  #                         level=self.level+1, father=self, child_id=0,
  #                         min_level=self.min_level, max_level=self.max_level),
  #                    Node(center=self.center+self.dx/4, dx=self.dx/2, value=self.value,
  #                         level=self.level+1, father=self, child_id=1,
  #                         min_level=self.min_level, max_level=self.max_level)]
    
  # def recursive_refine(self, target_level):
  #   if self.level>=target_level:
  #     return
  #   if not self.hasChildren():
  #     self.refine()
  #   for c in self.children:
  #     c.recursive_refine(target_level)
      
  # # def getMesh(self):
  # #   """ Recursive function to gather the outer leafs and construct the mesh """
  # #   if not self.hasChildren():
  # #     return self.xf, self.center, self.level, self.value, self
  # #   else:
  # #     left_xfs,  left_centers,  left_levels,  left_values,  left_nodes  = self.children[0].getMesh()
  # #     right_xfs, right_centers, right_levels, right_values, right_nodes = self.children[1].getMesh()
  # #     faces = np.hstack((left_xfs[:-1], right_xfs)) # remove duplicate face between children
  # #     if BDEBUG:
  # #       assert np.unique(faces).size == faces.size
  # #     return faces, \
  # #            np.hstack((left_centers, right_centers)), \
  # #            np.hstack((left_levels, right_levels)), \
  # #            np.vstack((left_values, right_values)), \
  # #            np.hstack((left_nodes, right_nodes))
             
    
  # def getLeafs(self):
  #   """ Returns the leaf nodes of its lineage """
  #   if not self.hasChildren():
  #     return [self]
  #   else:
  #     return [c.getLeafs() for c in self.children]
    
  # def getMaxSubLevel(self):
  #   """ Returns the maximum level among its lineage (useful for grading) """
  #   if not self.hasChildren():
  #     return self.level
  #   else:
  #     return max([c.getMaxSubLevel() for c in self.children])
    
  # def getMaxLeftSublevel(self):
  #   """ Returns the maximum level among its lineage
  #   at the left boundary of its domain (useful for grading) """
  #   if not self.hasChildren():
  #     return self.level
  #   else:
  #     return self.children[0].getMaxLeftSublevel()
    
  # def getMaxRightSublevel(self):
  #   """ Returns the maximum level among its lineage
  #   at the left boundary of its domain (useful for grading) """
  #   if not self.hasChildren():
  #     return self.level
  #   else:
  #     return self.children[1].getMaxRightSublevel()
      
  # def gradeTree(self):
  #   if not self.hasChildren():
  #     return
    
  #   max_left_level  =  self.children[0].getMaxRightSublevel()
  #   max_right_level =  self.children[1].getMaxLeftSublevel()
  #   if max_left_level < max_right_level-1: # refine left child
  #     self.children[0].recursiveRightRefine(max_right_level-1)
  #   elif max_right_level < max_left_level-1: # refine right child
  #     self.children[1].recursiveLeftRefine(max_left_level-1)
      
  #   for c in self.children:
  #     c.gradeTree()
      
      
  # def recursiveRightRefine(self, target_level):
  #   if self.level>=target_level:
  #     # raise Exception('why refine ?')
  #     return
  #   if not self.hasChildren():
  #     self.refine()
  #   self.children[1].recursiveRightRefine(target_level)
    
  # def recursiveLeftRefine(self, target_level):
  #   if self.level>=target_level:
  #     # raise Exception('why refine ?')
  #     return
  #   if not self.hasChildren():
  #     self.refine()
  #   self.children[0].recursiveLeftRefine( target_level)
      
  # def recursive_set_value(self,fun):
  #   self.value = fun(self.center)
  #   if self.hasChildren():
  #     for c in self.children:
  #       c.recursive_set_value(fun)