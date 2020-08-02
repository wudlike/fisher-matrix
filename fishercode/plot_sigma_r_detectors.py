# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 20:43:56 2019

@author: wudl

plot of sigma_r versus number of detectors (here at 95 GHz).
"""
import numpy as np
import matplotlib.pyplot as plt

#data = np.loadtxt(r'I:\work1\mypaper\code\data1\large\results\two\r_0.01\95_150_list.txt')

# data = np.loadtxt(r'../data1/large/results/two/r_0.01/95_150_list.txt')  #large
data = np.loadtxt(r'../data1/small/results/two/r_0.01/95_150_list.txt')     #small
x = data[:,0]
y = data[:,2]
y_min = min(y)
x_index = list(y).index(y_min)
x_value = x[x_index]

plt.plot(x,y)
plt.scatter(x_value,y_min,color='r')
plt.xlabel('$N_{95}$',fontsize=15)
plt.ylabel('$\sigma_r$',fontsize=15)
plt.tick_params(labelsize=15)

plt.show()
