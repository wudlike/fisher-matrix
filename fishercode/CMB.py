#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 11:27:32 2018

@author: wudl
"""

#CMB module
import numpy as np

def cmb(config):
#    st = Settings()
    totcl = np.loadtxt(r'../data/CMB_r/r_'+str(config['r'])+'/test_lensedtotCls.dat')
    cmb_TT = totcl[:,1]
    cmb_EE = totcl[:,2]
    cmb_BB = totcl[:,3]
    cmb_TE = totcl[:,4]
    return cmb_TT,cmb_EE,cmb_BB,cmb_TE

def cmb_par_r():
    tot_unlens1 = np.loadtxt(r'../data/CMB_r/r_1/test_tensCls.dat')
    par_bb = tot_unlens1[:,3];par_ee = tot_unlens1[:,2] 
    return par_bb,par_ee

#def cmb_r():
#    st = Settings()
#    totcl = np.loadtxt(r'D:/Share/LSPE/r_'+str(st.r)+'/test_lensedtotCls.dat')
##    totcl = np.loadtxt(r'/media/sf_Share/LSPE/r_0.01_7/test_lensedtotCls.dat')
#    cmb_TT = totcl[:,1]
#    cmb_EE = totcl[:,2]
#    cmb_BB = totcl[:,3]
#    cmb_TE = totcl[:,4]
#    ell = np.linspace(0,len(cmb_TT)-1,len(cmb_TT))
#    cmb_TT = cmb_TT/(ell*(ell+1)/(2*np.pi))
#    cmb_EE = cmb_EE/(ell*(ell+1)/(2*np.pi))
#    cmb_BB = cmb_BB/(ell*(ell+1)/(2*np.pi))
#    return cmb_TT,cmb_EE,cmb_BB,cmb_TE
#    
#def cmb_par_r():
#    tot_unlens1 = np.loadtxt(r'D:/Share/LSPE/r_1/test_tensCls.dat')
#    par_bb = tot_unlens1[:,3];par_ee = tot_unlens1[:,2]
#    ell = np.linspace(0,len(par_bb)-1,len(par_bb))
#    par_bb = par_bb/(ell*(ell+1)/(2*np.pi))
#    par_ee = par_ee/(ell*(ell+1)/(2*np.pi))
#    return par_bb,par_ee

def cmb_tensor(config):
     spec = np.loadtxt(r'../data/CMB_r/r_'+str(config['r'])+'/test_tensCls.dat')
     cmb_tensor_BB = spec[:,3]
     return cmb_tensor_BB
     
#if __name__=='__main__':
#    #plot cmb power spectra
#    cmb_TT,cmb_EE,cmb_BB,cmb_TE = cmb()
#    print(cmb_BB.shape)
#    par_bb,par_ee = cmb_par_r()
#    print(par_bb.shape)
#    l = np.linspace(st.l_min,st.l_max-1,st.l_max-st.l_min)
#    plt.plot(l,cmb_EE[st.l_min:st.l_max],label='CMB_EE')   
#    plt.loglog(l,cmb_BB[st.l_min:st.l_max],label='CMB_BB')  
#    plt.legend()
#    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
