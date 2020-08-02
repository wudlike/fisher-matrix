# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 12:00:50 2018

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
def sync2(v):             #1801.04945 PLanck 2018  71%
     v0 = 30; betas = -3.14; As = 0.91;
#RJ unit
     cl = As*(v/v0)**(2*betas)
#cmb unit
#     cl = As*(v/v0)**(2*betas)*(g(v)/g(v0))**2
     return cl


def dust2(v):      #1801.04945 PLanck 2018  71%
     betad = 1.53; vd = 353; Ad_ee = 315; Td = 19.6; alpha_d = -2.54
#RJ units
     cl = 0.53*Ad_ee*(v/vd)**(betad+1)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))
     
#cmb units
#     cl = 0.56*Ad_ee*(v/vd)**(betad+1)*((np.exp(h*vd/(kB*Td))-1)/(np.exp(h*v/(kB*Td))-1))*\
#     (g(v)/g(vd))**2
     return cl


def cmb(v):
     Acmb = 1
     cl = Acmb/g(v)
     return cl


startf = 10; endf = 1000; 
v = np.linspace(startf,endf,endf-startf+1)
#s1 = sync1(v)
s2 = sync2(v)
#s3 = sync3(v)
#d1 = dust1(v)
d2 = dust2(v) 
#d3 = dust3(v)
fg2 = s2+d2
cmb = cmb(v)

#plt.loglog(v,s1,'r',label='1509,sync1')
plt.loglog(v,s2,'k',label='sync_Planck,2018,71%')
#plt.loglog(v,s3,'b',label='1502,sync3')
#plt.loglog(v,d1,'r',label='1509,dust1')
plt.loglog(v,d2,'k',label='sync_Planck,2018,71%')
plt.loglog(v,fg2,'r',label='fg')
#plt.loglog(v,d3,'b',label='1502,dust3')
#plt.loglog(v,cmb,'g',label='cmb')
plt.xlabel('$frequency$')
plt.ylabel('$Brightness temperature$')
plt.legend()
plt.show()

#if __name__ == '__main__':
#     startl = 10; endl = 1000; nu = 95; ls = 80
#     ld = 80; alpha_s = -2.6; alpha_d = -2.54
#     beta_s = -2.9; beta_d = 1.59
#     l = np.linspace(startl,endl,endl-startl+1)
#     cl_s = sync2(nu)*(l/ls)**alpha_s
#     cl_d = dust2(nu)*(l/ld)**(alpha_d+2)
#     plt.xlabel('$\ell$')
#     plt.ylabel('$l(l+1)C_l/2\pi$')
#     plt.loglog(l,cl_s)
#     plt.loglog(l,cl_d)
#     plt.legend()
#     plt.show()
















