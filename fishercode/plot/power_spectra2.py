#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 20:31:16 2018

@author: wudl
"""
from __future__ import division
#import os
#os.chdir(r'/media/sf_Share/fisher_matrix/myLSPE/LSPE4')  #ATTENTION:change work dir
import numpy as np
import matplotlib.pyplot as plt
import foreground_models as fm
import settings
import CMB

st = settings.Settings()

l = np.linspace(st.l_min,st.l_max-1,st.l_max-st.l_min)
cmb_BB = CMB.cmb_r()[2][st.l_min:st.l_max]
cmb_tensor_BB = CMB.cmb_tensor()[st.l_min:st.l_max]
cmb_BB = cmb_BB - cmb_tensor_BB
#dust BB
dust_BB_list = []
sync_BB_list = []
noise_BB_list = []
noise_Dl_BB = []
for i in range(st.l_min,st.l_max):
    dust_BB = fm.Fg_matrix().Dust_BB(i)[0]   #first frequency dust BB
    sync_BB = fm.Fg_matrix().Sync_BB(i)[0]   #first frequency sync BB
    noise_bb,Dl_BB = fm.Instrumental_noise().noise(i,5000,2500,2500)[0],\
    fm.Instrumental_noise().noise(i,10000,10000,10000)[1]    #first noise BB, the latters are number of detector
    dust_BB_list.append(dust_BB[0])
    sync_BB_list.append(sync_BB[0])
    noise_BB_list.append(noise_bb[0])
    noise_Dl_BB.append(Dl_BB[0,0])
#    noise_Dl_BB.append(Dl_BB[0,0]/4)   #here over 4 because when compute noise_bb, i forgot over 4.
fg = np.array(dust_BB_list)+np.array(sync_BB_list)

sum_all_component = fg+np.array(noise_Dl_BB)+cmb_BB



#plt.loglog(l,cmb_BB,'r',label='cmb_BB')
#plt.loglog(l,cmb_tensor_BB,'g',label='cmb_tensor_BB')
#plt.loglog(l,dust_BB_list,'g',label='dust_BB')
#plt.loglog(l,sync_BB_list,'b',label='sync_BB')
#plt.loglog(l,noise_Dl_BB,'y',label='noise_Dl_BB')
plt.loglog(l,cmb_tensor_BB,'r',label='Tensor B-mode')
plt.loglog(l,dust_BB_list,'g',label='Dust')
plt.loglog(l,sync_BB_list,'b',label='Synchrotron')
plt.loglog(l,noise_Dl_BB,'y',label='Noise')
plt.loglog(l,cmb_BB,'r',label='lensing B-mode',linestyle='-.')
#plt.loglog(l,fg,'k',label='fg')
#plt.xlabel('$\ell$')
#plt.ylabel('$\ell(\ell+1)/(2\pi)C^{BB}_l$')
#plt.legend()

#plt.show()

#plt.loglog(l,cmb_BB,'r',label='cmb_lensed_BB')
#plt.loglog(l,cmb_tensor_BB,'g',label='cmb_tensor_BB')
#plt.loglog(l,noise_Dl_BB,'y',label='noise_Dl_BB')
#plt.loglog(l,sum_all_component)
#plt.loglog(l,fg)
plt.xlabel('$\ell$')
plt.ylabel('$\ell(\ell+1)C^{BB}_\ell/(2\pi)(\mu K)^2$')
plt.axis([30,1000,10**(-7),10**(0)])  #axis for 150
#plt.axis([30,1000,10**(-6),10**(1)])  #axis for 95 
#plt.axis([30,1000,10**(-5),10**(3)])   #axis for 41
plt.legend()
plt.show()

###plot power with respect frequency
#def x(v):        
#    kB = st.kB  # J/K
#    h = st.h # J·s, dimension becomes J/GHz by plus 9
#    Tcmb = 2.73
#    x = h*v/(kB*Tcmb)
#    return x        
#def g(v):
#    g = (np.exp(x(v))-1)**2/(x(v)**2*np.exp(x(v)))
#    return g
#
#def dust_BB(nu):   #1801.04945 PLANCK 5%
#    Td = 19.6; vd = 353; betad = 1.53; Ad = 20*0.53; ld = 80; alpha_d = -2.54
#    wvd = 1/g(vd); wv = 1/g(nu)
#    Wd = wvd/wv*(nu/vd)**(1+betad)*((np.exp(st.h*vd/(st.kB*Td))-1)/(np.exp(st.h*nu/(st.kB*Td))-1))
##    Wd = (nu/vd)**(2*betad-4)*(2*st.h*nu**3*st.c**(-2)/(np.exp(st.h*nu/(st.kB*Td)-1)))\
##    /(2*st.h*vd**3*st.c**(-2)/(np.exp(st.h*vd/(st.kB*Td)-1)))
##    dust_BB = Ad*(l/ld)**(alpha_d+2)*(Wd)**2
#    dust_BB = Ad*(Wd)**2
#    return dust_BB
#
#
##    def synchrotron(self,l):   #1801.04945 PLANCK    71%
##        betas = -3.14; As = 0.91; vs = 30; ls = 80; alpha_s = -2.6
##        wvs = 1/g(vs); wv = 1/g(self.nu)
##        ws = wvs/wv*(self.nu/vs)**betas
##        n_sync_bb = (ws)**2*As*(l/ls)**alpha_s
##        return n_sync_bb
#
#def synchrotron(nu):   #1801.04945 PLANCK    5%
#    betas = -3.14; As = 0.081; vs = 30; ls = 80; alpha_s = -2.6
#    wvs = 1/g(vs); wv = 1/g(nu)
#    ws = wvs/wv*(nu/vs)**betas
##    ws = (nu/vs)**betas
##    n_sync_bb = (ws)**2*As*(l/ls)**alpha_s
#    
#    n_sync_bb = (ws)**2*As
#    return n_sync_bb
#
#nu = np.linspace(10,1000,1000-10+1)
#dust = dust_BB(nu)
#sync = synchrotron(nu)
#fg = dust+sync
#plt.loglog(nu,dust,nu,sync)
#plt.loglog(nu,fg)
















