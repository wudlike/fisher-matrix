#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:44:08 2018

@author: wudl
"""

import numpy as np
import matplotlib.pyplot as plt


kB = 1.3806505*10**(-23)  # J/K
h = 6.6260693*10**(-34+9) # J·s, dimension becomes J/GHz by plus 9
c = 299792458

def x(v):
     kB = 1.3806505*10**(-23)  # J/K
     h = 6.6260693*10**(-34+9) # J·s, dimension becomes J/GHz by plus 9
     Tcmb = 2.73
     x = h*v/(kB*Tcmb)
     return x

def g(v):
     g = (np.exp(x(v))-1)**2/(x(v)**2*np.exp(x(v)))
     return g

def sync1(v):    #1502.01983  %72
     betas = -2.9; As = 2.1*10**(-5); vs = 65; 
     wvs = 1/g(vs); wv = 1/g(v)
#RJ unit
     ws = (v/vs)**betas
#cmb unit
#     ws = wvs/wv*(v/vs)**betas
     cl = (ws)**2*As
     return cl
 
def sync2(v):    #1502.01983  #53
     betas = -2.9; As = 2.1*10**(-5); vs = 65;
     wvs = 1/g(vs); wv = 1/g(v)
#RJ unit
     ws = (v/vs)**betas
#cmb unit
#     ws = wvs/wv*(v/vs)**betas
     cl = (ws)**2*As
     return cl
 
def sync3(v):    #1502.01983  #24
     betas = -2.9; As = 2.1*10**(-5); vs = 65; 
     wvs = 1/g(vs); wv = 1/g(v)
#RJ unit
     ws = (v/vs)**betas
#cmb unit
#     ws = wvs/wv*(v/vs)**betas
     cl = (ws)**2*As
     return cl
 
def sync4(v):    #1502.01983 %11
     betas = -2.9; As = 4.2*10**(-6); vs = 65; 
     wvs = 1/g(vs); wv = 1/g(v)
#RJ unit
     ws = (v/vs)**betas
#cmb unit
#     ws = wvs/wv*(v/vs)**betas
     cl = (ws)**2*As
     return cl
 
def sync5(v):    #1502.01983 %1
     betas = -2.9; As = 4.2*10**(-6); vs = 65; 
     wvs = 1/g(vs); wv = 1/g(v)
#RJ unit
     ws = (v/vs)**betas
#cmb unit
#     ws = wvs/wv*(v/vs)**betas
     cl = (ws)**2*As
     return cl
    
def dust1(v):   #1502.01983  %72
     Td = 19.6; vd = 353; betad = 1.59; Ad = 0.169
     wvd = 1/g(vd); wv = 1/g(v)
#RJ unit
     Wd = (v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
#cmb unit
#     Wd = wvd/wv*(v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
     cl = Ad*(Wd)**2
     return cl

def dust2(v):   #1502.01983  %53
     Td = 19.6; vd = 353; betad = 1.59; Ad = 0.065
     wvd = 1/g(vd); wv = 1/g(v)
#RJ unit
     Wd = (v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
#cmb unit
#     Wd = wvd/wv*(v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
     cl = Ad*(Wd)**2
     return cl
 
def dust3(v):   #1502.01983  %24
     Td = 19.6; vd = 353; betad = 1.59; Ad = 0.019
     wvd = 1/g(vd); wv = 1/g(v)
#RJ unit
     Wd = (v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
#cmb unit
#     Wd = wvd/wv*(v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
     cl = Ad*(Wd)**2
     return cl
 
def dust4(v):   #1502.01983  %11
     Td = 19.6; vd = 353; betad = 1.59; Ad = 0.013
     wvd = 1/g(vd); wv = 1/g(v)
#RJ unit
     Wd = (v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
#cmb unit
#     Wd = wvd/wv*(v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
     cl = Ad*(Wd)**2
     return cl
 
def dust5(v):   #1502.01983  %1
     Td = 19.6; vd = 353; betad = 1.59; Ad = 0.006
     wvd = 1/g(vd); wv = 1/g(v)
#RJ unit
     Wd = (v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
#cmb unit
#     Wd = wvd/wv*(v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
     cl = Ad*(Wd)**2
     return cl
 
def cmb(v):
     Acmb = 1
     cl = Acmb/g(v)
     return cl


startf = 10; endf = 1000; 
v = np.linspace(startf,endf,endf-startf+1)
s1 = sync1(v)
s2 = sync2(v)
s3 = sync3(v)
s4 = sync4(v)
s5 = sync5(v)
d1 = dust1(v)
d2 = dust2(v) 
d3 = dust3(v)
d4 = dust4(v)
d5 = dust5(v)
cmb = cmb(v)
plt.loglog(v,s1,'r',label='71%')
plt.loglog(v,s2,'k',label='53%')
plt.loglog(v,s3,'b',label='24%')
plt.loglog(v,s4,'y',label='11%')
plt.loglog(v,s5,'g',label="1%")
plt.loglog(v,d1,'r',label='71%')
plt.loglog(v,d2,'k',label='53%')
plt.loglog(v,d3,'b',label='24%')
plt.loglog(v,d4,'y',label='11%')
plt.loglog(v,d5,'g',label='1%')
plt.loglog(v,cmb,'r',label='cmb')
plt.legend()
plt.show()

if __name__ == '__main__':
     startl = 10; endl = 1000; nu = 150; ls = 80
     ld = 80; alpha_s = -2.6; alpha_d = -2.42
     beta_s = -2.9; beta_d = 1.59
     l = np.linspace(startl,endl,endl-startl+1)
     cl_s = sync3(nu)*(l/ls)**alpha_s
     cl_d = dust3(nu)*(l/ld)**(alpha_d+2)
     plt.loglog(l,cl_s)
     plt.loglog(l,cl_d)
     
     plt.show()