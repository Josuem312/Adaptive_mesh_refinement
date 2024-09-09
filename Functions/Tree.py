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
      self.flag= flag
      self.time= time
      self.value_q= value_q
      self.value_q1=value_q1
      self.counter_time= counter_time
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
    

