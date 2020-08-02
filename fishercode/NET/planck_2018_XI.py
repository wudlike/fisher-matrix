#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 16:35:00 2018

@author: wudl
"""

import numpy as np
import matplotlib.pyplot as plt


kB = 1.3806505*10**(-23)  # J/K
h = 6.6260693*10**(-34+9) # J·s, dimension becomes J/GHz by plus 9
c = 299792458
def sync1(v):     #1509.05419
     cl = 4.2*10**12*(v/0.408)**(-6)
     return cl

def x(v):
     kB = 1.3806505*10**(-23)  # J/K
     h = 6.6260693*10**(-34+9) # J·s, dimension becomes J/GHz by plus 9
     Tcmb = 2.73
     x = h*v/(kB*Tcmb)
     return x

def g(v):
     g = (np.exp(x(v))-1)**2/(x(v)**2*np.exp(x(v)))
     return g
def sync2(v):             #LSPE
     v0 = 30; betas = -3.14; As = 0.0031
#RJ unit
     cl = 0.35*As*(v/v0)**(2*betas)
#cmb unit
#     cl = 0.35*As*(v/v0)**(2*betas)*(g(v)/g(v0))**2
     return cl

def sync3(v):    #1502.01983
     betas = -2.9; As = 2.1*10**(-5); vs = 65; ls = 80; alpha_s = -2.6
     wvs = 1/g(vs); wv = 1/g(v)
#RJ unit
     ws = (v/vs)**betas
#cmb unit
#     ws = wvs/wv*(v/vs)**betas
     cl = (ws)**2*As
     return cl
def sync4(v):  #planck 2018 XI dust
    vs = 30; betas = -3.15; As = 51
    wvs = 1/g(vs); wv = 1/g(v)
    ws = wvs/wv
#    Dl = As*(v/vs)**(2*betas)*(ws)**2
    Dl = As*(v/vs)**(2*betas)
    return Dl

def dust1(obj):       #1509.05419
     Ad = 247; vd = 353; betad = 1.59; td = 19.6
     dust_BB = Ad*(v/vd)**(2*betad-4)*\
     (2*h*v**3*c**(-2)/(np.exp(h*v/(kB*td))-1)/(2*h*vd**3*c**(-2)/(np.exp(h*vd/(kB*td))-1)))**2 
     return dust_BB

def dust2(v):      #LSPE
     betad = 1.59; vd = 353; Ad = 0.2765; Td = 19.6
#RJ units
     cl = 0.56*Ad*(v/vd)**(betad+1)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
     
#cmb units
#     cl = 0.56*Ad*(v/vd)**(betad+1)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))*\
#     (g(v)/g(vd))**2
     return cl

def dust3(v):
     Td = 19.6; vd = 353; betad = 1.59; Ad = 0.169
     wvd = 1/g(vd); wv = 1/g(v)
#RJ unit
     Wd = (v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
#cmb unit
#     Wd = wvd/wv*(v/vd)**(1+betad)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
     cl = Ad*(Wd)**2
     return cl

def dust4(v): #planck
    alpha = -2.16; Aee = 34.3; Abb = 0.48*Aee; Td = 19.6; betad = 1.50; vd = 353
    Dl = Abb*(v/353)**(betad+1)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
#    Dl = Abb*(v/353)**(betad+1)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))*(g(v)/g(vd))**2
    return Dl

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
d1 = dust1(v)
d2 = dust2(v) 
d3 = dust3(v)
d4 = dust4(v)
cmb = cmb(v)
#plt.loglog(v,s1,'r',label='1509,sync1')
plt.loglog(v,s2,'k',label='LSPE,sync2')
#plt.loglog(v,s3,'b',label='1502,sync3')
#plt.loglog(v,d1,'r',label='1509,dust1')
plt.loglog(v,d2,'k',label='LSPE,dust2')
#plt.loglog(v,d3,'b',label='1502,dust3')
#plt.loglog(v,cmb,'g',label='cmb')
plt.legend()
plt.show()
plt.loglog(v,s4,'k',label='planck')
plt.loglog(v,d4,'k',label='planck')
plt.xlabel('$\mu$')
plt.legend()
plt.show()

if __name__ == '__main__':
     startl = 2; endl = 600; nu = 150; ls = 80
     ld = 80; alpha_s = -2.6; alpha_d = -2.42
     beta_s = -2.9; beta_d = 1.59
     l = np.linspace(startl,endl,endl-startl+1)
     cl_s = sync3(nu)*(l/ls)**alpha_s
     cl_d = dust3(nu)*(l/ld)**(alpha_d+2)
     plt.loglog(l,cl_s)
     plt.loglog(l,cl_d)
     
     plt.show()
     Dl_planck = dust4(95)*(l/ld)**(-2.16+2)
     Sl_planck = sync4(95)*(l/ld)**(-2.32)
     plt.loglog(l,Dl_planck)
     plt.loglog(l,Sl_planck)
     plt.xlabel('$\ell$')
     plt.axis([startl,endl,10**(-5),10**2])
     plt.show()
     
     
     