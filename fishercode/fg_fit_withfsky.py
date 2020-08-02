# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 21:24:36 2019

@author: wudl
"""

# the paper from Planck 2018 results. XI. Polarized dust foregrounds

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def polynomial_fit(x,y,N):
    return np.polyfit(x,y,N)

def exponential(p):
    a,b = p
    return Ad-(a*np.exp(b*fsky))

def power_law(p):
    a,b = p
    return Ad-(a*fsky**b)

def test_cos(p):
    a,b = p
    return Ad-(a/(1-fsky+b))

fsky = np.array([0.24,0.33,0.42,0.52,0.62,0.71])
Ad = np.array([0.48*34.3,0.45*47.3,0.50*74.7,0.53*120.1,0.53*190.7,0.53*315.4])
# x range
x = np.linspace(0,1,100)

# polynomial_fit
popt_poly = polynomial_fit(fsky,Ad,3)
p1 = np.poly1d(popt_poly)
Ad_poly = p1(x)
print(p1)
#test cos
popt_cos = leastsq(test_cos,[1,0])
a_cos,b_cos = popt_cos[0]
#exponential form
popt_exp = leastsq(exponential,[1,0])
a_exp,b_exp = popt_exp[0]
#power law
popt_power_law = leastsq(power_law,[1,0])
a_pow,b_pow = popt_power_law[0]

Ad_cos = a_cos/(1-x+b_cos)
Ad_fit_exp = a_exp*np.exp(b_exp*x)
Ad_fit_pow = a_pow*x**b_pow

if __name__=='__main__':
     #plt.figure()
     plt.scatter(fsky,Ad)
     plt.loglog(x,Ad_poly,label='polynomial fit')
     plt.loglog(x,Ad_cos,label='cos fit')
     plt.loglog(x,Ad_fit_exp,label='exponential fit')
     plt.loglog(x,Ad_fit_pow,label='power_law fit')
     plt.xlabel('fsky')
     plt.ylabel('Ad')
     
     plt.legend()
     plt.show()
