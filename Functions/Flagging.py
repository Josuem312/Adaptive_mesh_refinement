# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 13:47:17 2024

@author: josue
"""
import numpy as np

"""
This function flags all the points where the estimator errror is bigger 
 than a threshold. Here WE ADD THE BUFFER ZONE .
 Count allows us to know if there exists any point for regriding.
"""
def flagging(error, epsilon):
    flags= []
    count=0
    for i, u in enumerate(error):
        if u > epsilon:
            count=count+1
            if count==1:
                    #flags.append(i-3)
                    #flags.append(i-2)
                    flags.append(i-1)
            flags.append(i)          
    if len(flags)!=0:
        max_flag= max(flags)+1
      #  max_flag_2= max_flag +1
       # max_flag_3= max_flag_2 +1
        flags.append(max_flag)
        #flags.append(max_flag_2)
        #flags.append(max_flag_3)
    return flags, count



########################################################################
#Here we have other versions of flagging functions
########################################################################
def flaggin_without_buffering(error, epsilon):
    flags= []
    count=0
    for i, u in enumerate(error):
        if u > epsilon:
            flags.append(i)          
    return flags, count


def flagging_with_subgriddings(flag):
    subgridings=[]
    for i in np.arange(len(flag)-1):
        if flag[i+1]-flag[i]>2:
            subgridings.append([flag[i], flag[i+1]])
    return subgridings