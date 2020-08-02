#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 22:35:08 2018

@author: wudl
"""
import numpy as np
from sympy import symbols,exp,Integral
from scipy import integrate
import matplotlib.pyplot as plt

#v=symbols('v')
eta = 0.5
k = 1.3806505*10**(-23)  # J/K
#nu = 28   #unit: GHz
nu = float(input('Please enter frequency: '))
#TRJ = 5
TRJ = float(input('Please enter TRJ: '))
#epsilon = 0.05
epsilon = float(input('Please enter epsilon: '))
deltanu = 0.3*nu*10**9  # 1/s
pinternal = 1*10**(-12)    #W
tau = 0.5
h = 6.6260693*10**(-34+9) # J/GHz
T = 2.73             # k
c = 299792458
band = 0.3/2             #(band range has a great impact on NET)

pload = 2*eta*k*TRJ*deltanu+pinternal   # unit:W or J/s
print('Pload =>>>',pload)
NEP = np.sqrt(2*h*nu*pload+2*pload**2/deltanu)  # unit:J*s**(-1/2)
#change NEP unit to aW/sqrt(Hz)
NEPunit = NEP*10**(18)
print('NEP =>>> %.f'%NEPunit)
def f(v):
    return h**2*v**2/(T**2*k)*np.exp(h*v/(k*T))*tau*(1-epsilon)/(np.exp(h*v/(k*T))-1)**2 

w, err = integrate.quad(f,nu-nu*band,nu+nu*band)
# times 10**6 from k to uk; divided by 10**9 from GHz to s
NET = 10**6*NEP/w/(10**9)      #uk*sqrt(s)  
print('NET =>>> %.f'%NET)

#NET_list = []
#for nu in range(30,350,5):
##    if nu <= 200:
##        epsilon = 0.05
##    elif nu <= 250:
##        epsilon = 0.1
##    elif nu <= 300:
##        epsilon = 0.2
#    pload = 2*eta*k*TRJ*deltanu+pinternal   # unit:W or J/s
#    NEP = np.sqrt(2*h*nu*pload+2*pload**2/deltanu)  # unit:J*s**(-1/2)
#    w, err = integrate.quad(f,nu-nu*band,nu+nu*band)
#    # times 10**6 from k to uk; divided by 10**9 from GHz to s
#    NET = NEP/w*10**6/(10**9)      #uk*sqrt(s)  
#    NET_list.append(NET)
#
#nu = range(30,350,5)
#plt.plot(nu,NET_list)
#plt.xlabel('$frequency$')
#plt.ylabel('NET')
#plt.show()















