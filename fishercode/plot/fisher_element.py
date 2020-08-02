#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 15:45:47 2018

@author: wudl
"""

#plot fisher element
import numpy as np
import os
os.chdir('D:/wudl/fisher_matrix/myLSPE/LSPE3') 
import fisher
import matplotlib.pyplot as plt
from settings import Settings
st = Settings()

def fisher_BB_element():
     if len(st.nu) == 2:
          fisher_BB_element = fisher.Fisher.calculation()[5]
     elif len(st.nu) == 3:
          fisher_BB_element = fisher.Fisher.calculation()[6]
     elif len(st.nu) == 4:
          fisher_BB_element = fisher.Fisher.calculation()[7]
     return fisher_BB_element
l = np.linspace(st.l_min,st.l_max-1,st.l_max-st.l_min)
fisher_BB_element = fisher_BB_element()
plt.plot(l,fisher_BB_element,label='Fisher BB element')
plt.ylabel('$F_{rr}$')
plt.xlabel('$\ell$')
plt.legend()
plt.show()




