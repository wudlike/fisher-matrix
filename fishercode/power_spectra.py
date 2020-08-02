#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 20:31:16 2018

@author: wudl
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import foreground_models as fm
import settings2
import CMB

st = settings2.Settings()

l = np.linspace(st.l_min,st.l_max-1,st.l_max-st.l_min)
cmb_BB = CMB.cmb_r()[2][st.l_min:st.l_max]
cmb_tensor_BB = CMB.cmb_tensor()[st.l_min:st.l_max]
#dust BB
dust_BB_list = []
sync_BB_list = []
noise_BB_list = []
noise_Dl_BB = []
for i in range(st.l_min,st.l_max):
    dust_BB = fm.Fg_matrix().Dust_BB(i)[0]   #first frequency dust BB
    sync_BB = fm.Fg_matrix().Sync_BB(i)[0]   #first frequency sync BB
    noise_bb,Dl_BB = fm.Instrumental_noise().noise(i,5000,2500,2500)[0],\
    fm.Instrumental_noise().noise(i,5000,2500,2500)[1]    #first noise BB, the latters are number of detector
    dust_BB_list.append(dust_BB[0])
    sync_BB_list.append(sync_BB[0])
    noise_BB_list.append(noise_bb[0])
    noise_Dl_BB.append(Dl_BB[0])
#fg = np.array(dust_BB_list)+np.array(sync_BB_list)
##plt.loglog(l,cmb_BB,'r',label='cmb_BB')
##plt.loglog(l,cmb_tensor_BB,'g',label='cmb_tensor_BB')
#plt.loglog(l,dust_BB_list,'g',label='dust_BB')
#plt.loglog(l,sync_BB_list,'b',label='sync_BB')
#plt.loglog(l,noise_Dl_BB,'y',label='noise_Dl_BB')
##plt.loglog(l,fg,'k',label='fg')
#plt.xlabel('$\ell$')
#plt.ylabel('$\ell(\ell+1)/(2\pi)C^{BB}_l$')
#plt.legend()
#plt.show()

plt.loglog(l,cmb_BB,'r',label='cmb_lensed_BB')
plt.loglog(l,cmb_tensor_BB,'g',label='cmb_tensor_BB')
plt.loglog(l,noise_Dl_BB,'y',label='noise_Dl_BB')
plt.xlabel('$\ell$')
plt.ylabel('$\ell(\ell+1)/(2\pi)C^{BB}_l$')
plt.legend()
plt.show()