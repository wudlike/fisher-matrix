#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import settings as sts 
st = sts.Settings()
#from settings import Settings

tot_unlens1 = np.loadtxt(r'../data/CMB_r/r_1/test_tensCls.dat')


class Foreground():

    def __init__(self,nu):
        self.nu = nu
        
    def x(self,v):        
        kB = st.kB  # J/K
        h = st.h # J·s, dimension becomes J/GHz by plus 9
        Tcmb = st.Tcmb
        x = h*v/(kB*Tcmb)
        return x 
       
    def g(self,v):
        g = (np.exp(self.x(v))-1)**2/(self.x(v)**2*np.exp(self.x(v)))
        return g        
        
    def cmb_r(self):
        totcl = np.loadtxt(r'../data/CMB_r/r_'+str(st.r)+'/test_lensedtotCls.dat')
#        totcl = np.loadtxt(r'D:/Share/LSPE/r_0.01_7/test_lensedtotCls.dat')
        cmb_TT = totcl[:,1]
        cmb_EE = totcl[:,2]
        cmb_BB = totcl[:,3]
        cmb_TE = totcl[:,4]
        return cmb_TT,cmb_EE,cmb_BB,cmb_TE
    
    def boltzmann_func(self):
        bolt = 2*st.h*self.nu**3*st.c**(-2)/(np.exp(st.h*self.nu/(st.kB*self.td))-1)
        return bolt    
        
    def dust_BB(self,l):   #1801.04945 PLANCK        5% and 71%
        Td = st.params_dict['Td']; vd = st.params_dict['vd']; betad = st.params_dict['betad']
        Ad = st.params_dict['Ad']; ld = st.params_dict['ld']; alpha_d = st.params_dict['alpha_d']
        wvd = 1/self.g(vd); wv = 1/self.g(self.nu)
        Wd = wvd/wv*(self.nu/vd)**(1+betad)*((np.exp(st.h*vd/(st.kB*Td))-1)/(np.exp(st.h*self.nu/(st.kB*Td))-1))
        dust_BB = Ad*(l/ld)**(alpha_d+2)*(Wd)**2
        return dust_BB

    def synchrotron(self,l):   #1801.04945 PLANCK    5% and 71%
        betas = st.params_dict['betas']; As = st.params_dict['As']; vs = st.params_dict['vs']
        ls = st.params_dict['ls']; alpha_s = st.params_dict['alpha_s']        
        wvs = 1/self.g(vs); wv = 1/self.g(self.nu)
        ws = wvs/wv*(self.nu/vs)**betas
        n_sync_bb = (ws)**2*As*(l/ls)**alpha_s
        return n_sync_bb
    
class Fg_matrix():
    '''
    The matrix contains all frequency, e.g. Dust_bb[0,0] and Dust_bb[1,1] correspond to 1st frequency and 2nd frequency,respectively.
    '''
    def Dust_BB(self,l):
        Dust_bb = np.zeros([len(st.nu),len(st.nu)])
        for i in range(len(st.nu)):
            Dust_bb[i,i] = Foreground(st.nu[i]).dust_BB(l)
        return Dust_bb
    def Dust_EE(self,l):
        Dust_ee = np.zeros([len(st.nu),len(st.nu)])
        for i in range(len(st.nu)):
            Dust_ee[i,i] = Foreground(st.nu[i]).dust_EE(l)
        return Dust_ee
    def Sync_BB(self,l):
        Sync_bb = np.zeros([len(st.nu),len(st.nu)])
        for i in range(len(st.nu)):
            Sync_bb[i,i] = Foreground(st.nu[i]).synchrotron(l)
        return Sync_bb
    
class Instrumental_noise():
    
    def __init__(self):
        pass
    
    def noise(self,l,x,y,z,m):
        Yield = 1   #NET :uk.sqrt(s)
        EffectiveDetectorSeconds = st.Tobs * Yield
        Noise_EE = np.zeros([len(st.nu),len(st.nu)])
        Dl_BB = np.zeros([len(st.nu),len(st.nu)])
        for i in range(len(st.nu)):
            if i == 0:
                Ndet = x
            elif i == 1:
                Ndet = y
            elif i == 2:
                Ndet = z
            elif i == 3:
                Ndet = m
            elif i == 4:
                Ndet = st.tot_det-x-y-z-m
            FWHMs = st.freq_dict[str(st.nu[i])+' GHz']['FWHMs']
            NET = st.freq_dict[str(st.nu[i])+' GHz']['NET']   
            w_efft = (np.sqrt(Ndet*EffectiveDetectorSeconds))/(NET*np.sqrt(st.skyam))    
            sigma_b = (FWHMs*st.arcmin_to_radian)/(np.sqrt(8.0*np.log(2.0)))    
            w_effb = w_efft/np.sqrt(2); w_effe = w_effb
            n_ee_inst_e = w_effe**2*np.exp(-l*(l+1)*sigma_b**2)
            n_ee_instrument = 1.0/n_ee_inst_e*st.arcmin_to_radian**2
            Noise_EE[i,i] = n_ee_instrument
            Dl_BB[i,i] = l*(l+1)/(2*np.pi)*Noise_EE[i,i]
        Noise_BB = Noise_EE
        Dl_BB = Dl_BB
        return Noise_BB,Dl_BB














