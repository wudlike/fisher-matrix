#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 16:47:21 2018

@author: wudl
"""

#plot Scott Dodelson paper:doi 10.1103/PhysRevLett.112.191301 
#red line , 使用了l之后的求和
from fisher_element import fisher_BB_element as fBe  
import numpy as np
import os
os.chdir(r'/media/sf_Share/fisher_matrix/myLSPE/LSPE4') 
import matplotlib.pyplot as plt
from settings import Settings
st = Settings()

def Scott1(obj):
    i = 0
    j = 0
    sum = 0
    cc = np.zeros(len(obj))
    for i in range(len(obj)):
        while j< len(obj):
            sum += obj[j]
            j += 1
        cc[i] = sum
        j = i+1
        sum = 0
    fsum = np.array(cc)
    delt_r = np.sqrt(1/fsum)
    r_delt_r = st.r/delt_r
    return r_delt_r

def Scott2(kk):
    i = 0
    j = 0
    sum = 0
    hh = np.zeros(len(kk))
    for i in range(len(kk)):
        while j<= i:
            sum += kk[j]
            j += 1
        hh[i] = sum
        j = 0
        sum = 0
    fsumh = np.array(hh)
    delt_rh = np.sqrt(1/fsumh)
    r_delt_rh = st.r/delt_rh
    return r_delt_rh

if __name__=='__main__':
     fisher_BB_element = fBe
     r_delt_r = Scott1(fisher_BB_element)
     r_delt_rh = Scott2(fisher_BB_element)
           
     l = np.linspace(st.l_min,st.l_max-1,st.l_max-st.l_min)
     plt.plot(l,r_delt_r*1,'r',label="$\ell'>\ell_{min}$")           #放大4.5倍
     #blue line  使用了l之前的求和
     plt.plot(l,r_delt_rh*1,'b',label="$\ell'<\ell_{max}$")        #放大4.5倍
     #plt.axis([0,700,0,15])
#     plt.legend(title='Scott sum BB')
     plt.ylabel('$r/(\Delta r)$')
     plt.xlabel("$\ell'$")
     plt.legend()
     plt.show()
    

    


