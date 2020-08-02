# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 21:31:18 2019

@author: wudl
"""

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-10,10,100)
y = np.sin(x)

plt.plot(x,y)
plt.xlabel('$\\theta$')
plt.ylabel('$\\beta$')
plt.title(r'$\beta=sin(\theta)$')

plt.show()