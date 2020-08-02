#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Jan 26 18:48:52 2018

@author: wudl
The aim of this code is to plot the figure which can interpret the connection of 
sigma_r with respect to fsky.
the data 

"""
import numpy as np
import settings2 as sts 
import time
import matplotlib.pyplot as plt
from fg_fit_withfsky import a_exp,b_exp

st = sts.Settings()

tot_unlens1 = np.loadtxt(r'../data/CMB_r/r_1/test_tensCls.dat')

class CMB():
    
    def __init__(self,r):
        self.r =r
        
    def cmb(self):
    #    st = Settings()
        totcl = np.loadtxt(r'../data/CMB_r/r_'+str(self.r)+'/test_lensedtotCls.dat')
        cmb_TT = totcl[:,1]
        cmb_EE = totcl[:,2]
        cmb_BB = totcl[:,3]
        cmb_TE = totcl[:,4]
        return cmb_TT,cmb_EE,cmb_BB,cmb_TE
    
    def cmb_par_r(self):
        tot_unlens1 = np.loadtxt(r'../data/CMB_r/r_1/test_tensCls.dat')
        par_bb = tot_unlens1[:,3];par_ee = tot_unlens1[:,2] 
        return par_bb,par_ee

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
        
    def dust_BB(self,l,fsky):   #1801.04945 PLANCK        5% and 71%
#        Ad = st.params_dict['Ad'];
        Ad = a_exp*np.exp(b_exp*fsky)
        Td = st.params_dict['Td']; vd = st.params_dict['vd']; betad = st.params_dict['betad']
        ld = st.params_dict['ld']; alpha_d = st.params_dict['alpha_d']
        wvd = 1/self.g(vd); wv = 1/self.g(self.nu)
        Wd = wvd/wv*(self.nu/vd)**(1+betad)*((np.exp(st.h*vd/(st.kB*Td))-1)/(np.exp(st.h*self.nu/(st.kB*Td))-1))
        dust_BB = Ad*(l/ld)**(alpha_d+2)*(Wd)**2
        return dust_BB

    def synchrotron(self,l,fsky):   #1801.04945 PLANCK    5% and 71%
#        As = st.params_dict['As']
        As_78per = 0.9 # fit As from the data of fsky=78%
        As = As_78per*np.exp(b_exp*(fsky-0.78))
        betas = st.params_dict['betas']; vs = st.params_dict['vs']
        ls = st.params_dict['ls']; alpha_s = st.params_dict['alpha_s']        
        wvs = 1/self.g(vs); wv = 1/self.g(self.nu)
        ws = wvs/wv*(self.nu/vs)**betas
        n_sync_bb = (ws)**2*As*(l/ls)**alpha_s
        return n_sync_bb
    
class Fg_matrix():
    '''
    The matrix contains all frequency, e.g. Dust_bb[0,0] and Dust_bb[1,1] correspond to 1st frequency and 2nd frequency,respectively.
    '''
    def Dust_BB(self,l,fsky):
        Dust_bb = np.zeros([len(st.nu),len(st.nu)])
        for i in range(len(st.nu)):
            Dust_bb[i,i] = Foreground(st.nu[i]).dust_BB(l,fsky)
        return Dust_bb
    def Dust_EE(self,l):
        Dust_ee = np.zeros([len(st.nu),len(st.nu)])
        for i in range(len(st.nu)):
            Dust_ee[i,i] = Foreground(st.nu[i]).dust_EE(l)
        return Dust_ee
    def Sync_BB(self,l,fsky):
        Sync_bb = np.zeros([len(st.nu),len(st.nu)])
        for i in range(len(st.nu)):
            Sync_bb[i,i] = Foreground(st.nu[i]).synchrotron(l,fsky)
        return Sync_bb
    
class Instrumental_noise():
    
    def __init__(self,fsky):
        self.fsky = fsky
        self.skyam = 4.0*np.pi*self.fsky*((180.0*60.0)/np.pi)**2 
        
    def noise(self,l,x,y,z):
        Yield = 1;   #NET :uk.sqrt(s)
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
                Ndet = st.tot_det-x-y-z
            FWHMs = st.freq_dict[str(st.nu[i])+' GHz']['FWHMs']
            NET = st.freq_dict[str(st.nu[i])+' GHz']['NET']            
            w_efft = (np.sqrt(Ndet*EffectiveDetectorSeconds))/(NET*np.sqrt(self.skyam))    
            sigma_b = (FWHMs*st.arcmin_to_radian)/(np.sqrt(8.0*np.log(2.0)));    
            w_effb = w_efft/np.sqrt(2); w_effe = w_effb
            n_ee_inst_e = w_effe**2*np.exp(-l*(l+1)*sigma_b**2)
            n_ee_instrument = 1.0/n_ee_inst_e*st.arcmin_to_radian**2
            Noise_EE[i,i] = n_ee_instrument
            Dl_BB[i,i] = l*(l+1)/(2*np.pi)*Noise_EE[i,i]
        Noise_BB = Noise_EE
        Dl_BB = Dl_BB
        return Noise_BB,Dl_BB

class Fisher(object):
    def __init__(self,Add_prior=[],params=st.params,param_values=st.param_values,prior_percents=st.prior_percents):
        self.Add_prior = Add_prior
        self.params = params
        self.param_values = param_values

    def find_min_index(self,matrix_list):
        '''
        input a matrix_list like, e.g., [mat[a],mat[b],mat[c]...]
        find minimum value of parameters and its corresponding index from a matrix-list
        '''
        Frr = []
        for i in range(len(matrix_list)):
            Frr.append(matrix_list[i][0,0])         # [0,0] for r, [1,1] for As...
        Frr_min=min(Frr)
        index_min = Frr.index(Frr_min)
        Fmatrix_min=matrix_list[index_min]
        return index_min,Frr,Frr_min,Fmatrix_min
        
    def add_priors(self,param_values,param_percents):
        
        priors = np.array(param_values)*np.array(param_percents)
        prior_matrix = np.zeros((st.nparas,st.nparas))
        for i in range(st.nparas):
            if priors[i] != 0:
                prior_matrix[i,i]=1./priors[i]**2.0
        return prior_matrix
            
    
    def sqrt(self,sigmaMat):
        '''
        This function is used to sqrt the diagonal elements of sigmaMat
        '''
        for i in range(len(sigmaMat)):
            sigmaMat[i,i] = np.sqrt(sigmaMat[i,i])
        return sigmaMat
    
    def mat_derivative(self,l,x,y,z,C_l,C_l_r,C_l_As,C_l_alps,C_l_betas,C_l_Ad,C_l_alpd,C_l_betad,cmb_BB,\
                          ls,vs,ld,vd,Ad,As,par_bb,fsky):
        '''
        Input matrix Cl, output partial Cl/ partial theta, where theta refers to parameter.
        '''
        for i in range(st.n_nu):
            '''
            Cl is a symmetric matrix, so antidiagonal elements are equal, just need to compute one time.
            Thus, for j in range(i,len(nu))
            '''
            for j in range(st.n_nu):
#            for j in range(i,st.n_nu):    
                if i == j:
                    C_l[i][j] = cmb_BB[l]+Fg_matrix().Sync_BB(l,fsky)[i,i]*st.sigma_F\
                    +Fg_matrix().Dust_BB(l,fsky)[i,i]*st.sigma_F+Instrumental_noise(fsky).noise(l,x,y,z)[1][i,i] # [1] for Dl_BB,[0] for Cl_BB
                    # partial Cl/ partial r
                    C_l_r[i][j] = par_bb[l]  
                    # partial Cl/ partial As
                    C_l_As[i][j] = Fg_matrix().Sync_BB(l,fsky)[i,i]/As*st.sigma_F
                    # partial Cl/ partial alpha_s
                    C_l_alps[i][j] = Fg_matrix().Sync_BB(l,fsky)[i,i]*np.log(l/ls)*st.sigma_F
                    # partial Cl/ partial beta_s
                    C_l_betas[i][j] = Fg_matrix().Sync_BB(l,fsky)[i,i]*st.sigma_F*2*np.log(st.nu[i]/vs)
                    # for dust
                    C_l_Ad[i][j] = Fg_matrix().Dust_BB(l,fsky)[i,i]/Ad*st.sigma_F
                    # partial Cl/ partial alpha_d
                    C_l_alpd[i][j] = Fg_matrix().Dust_BB(l,fsky)[i,i]*np.log(l/ld)*st.sigma_F
                    # partial Cl/ partial beta_d
                    C_l_betad[i][j] = Fg_matrix().Dust_BB(l,fsky)[i,i]*st.sigma_F*2*np.log(st.nu[i]/vd)                                
                else:
                    C_l[i][j] = cmb_BB[l]+np.sqrt(Fg_matrix().Sync_BB(l,fsky)[i,i]*st.sigma_F*Fg_matrix().Sync_BB(l,fsky)[j,j]*st.sigma_F)\
                    +np.sqrt(Fg_matrix().Dust_BB(l,fsky)[i,i]*st.sigma_F*Fg_matrix().Dust_BB(l,fsky)[j,j]*st.sigma_F)   
#                    C_l[j][i] = C_l[i][j]
                    # partial Cl/ partial r
                    C_l_r[i][j] = par_bb[l]   
#                    C_l_r[j][i] = C_l_r[i][j]
                    # partial Cl/ partial As
                    C_l_As[i][j] = np.sqrt(Fg_matrix().Sync_BB(l,fsky)[i,i]*st.sigma_F*Fg_matrix().Sync_BB(l,fsky)[j,j]*st.sigma_F)/As                                
#                    C_l_As[j][i] = C_l_As[i][j]
                    # partial Cl/ partial alpha_s
                    C_l_alps[i][j] = np.sqrt(Fg_matrix().Sync_BB(l,fsky)[i,i]*st.sigma_F*Fg_matrix().Sync_BB(l,fsky)[j,j]*st.sigma_F)*np.log(l/ls)
#                    C_l_alps[j][i] = C_l_alps[i][j]
                    # partial Cl/ partial beta_s
                    C_l_betas[i][j] = np.sqrt(Fg_matrix().Sync_BB(l,fsky)[i,i]*st.sigma_F*Fg_matrix().Sync_BB(l,fsky)[j,j]*st.sigma_F)*np.log(st.nu[i]*st.nu[j]/(vs*vs))     
#                    C_l_betas[j][i] =C_l_betas[i][j]
                    # for dust
                    C_l_Ad[i][j] = np.sqrt(Fg_matrix().Dust_BB(l,fsky)[i,i]*st.sigma_F*Fg_matrix().Dust_BB(l,fsky)[j,j]*st.sigma_F)/Ad                                
#                    C_l_Ad[j][i] = C_l_Ad[i][j]
                    # partial Cl/ partial alpha_d
                    C_l_alpd[i][j] = np.sqrt(Fg_matrix().Dust_BB(l,fsky)[i,i]*st.sigma_F*Fg_matrix().Dust_BB(l,fsky)[j,j]*st.sigma_F)*np.log(l/ld)
#                    C_l_alpd[j][i] = C_l_alpd[i][j]
                    # partial Cl/ partial beta_d
                    C_l_betad[i][j] = np.sqrt(Fg_matrix().Dust_BB(l,fsky)[i,i]*st.sigma_F*Fg_matrix().Dust_BB(l,fsky)[j,j]*st.sigma_F)*np.log(st.nu[i]*st.nu[j]/(vd*vd))  
#                    C_l_betad[j][i] = C_l_betad[i][j]
        return C_l,C_l_r, C_l_As, C_l_alps, C_l_betas, C_l_Ad, C_l_alpd, C_l_betad
    
    def cov_matrix(self,x,y,z,ls,vs,ld,vd,Ad,As,fmatrix,cmb_BB,par_bb,sigmaMat_list,sigma_r_list,fsky):
        for l in range(st.l_min,st.l_max): 
            C_l = np.eye(len(st.nu)); C_l_r = np.eye(len(st.nu))
            C_l_As = np.eye(len(st.nu)); C_l_alps = np.eye(len(st.nu)); C_l_betas = np.eye(len(st.nu))
            C_l_Ad = np.eye(len(st.nu)); C_l_alpd = np.eye(len(st.nu)); C_l_betad = np.eye(len(st.nu))
                
            C_l,C_l_r, C_l_As, C_l_alps, C_l_betas, C_l_Ad, C_l_alpd, C_l_betad = \
            Fisher(self.Add_prior).mat_derivative(l,x,y,z,C_l,C_l_r,C_l_As,C_l_alps,C_l_betas,C_l_Ad,C_l_alpd,C_l_betad,cmb_BB,\
                  ls,vs,ld,vd,Ad,As,par_bb,fsky)
            # A list of matrix where the order of element is r As alphas betas Ad alphad betad                    
            C_l_inv = np.linalg.inv(C_l)
            par_C_l_list=[C_l_r, C_l_As, C_l_alps, C_l_betas, C_l_Ad, C_l_alpd, C_l_betad]
            for i in range(len(par_C_l_list)):
                for j in range(len(par_C_l_list)):
                    fmatrix[i][j] += (2*l+1)/2*fsky*np.trace(np.dot(np.dot(C_l_inv,par_C_l_list[i]),np.dot(C_l_inv,par_C_l_list[j])))
        '''              
        sigma_ij = np.sqrt(F^-1_ij) but here since covariance matrix can be ninus(non-diagonal elements), so 
        we set sqrt after get min_sigma, and just for diagonal elements
        Therefore, here below we give the squre of sigma_matrix
        '''
        sigma_mat_square = np.linalg.inv(fmatrix) # the truth is sigma_mat = np.sqrt(np.linalg.inv(fmatrix))
        if self.Add_prior=='True':
            prior_matrix = self.add_priors(st.param_values,st.prior_percents)
            sigma_mat_square = np.linalg.inv(fmatrix+prior_matrix)                
        sigmaMat_list.append(sigma_mat_square)  
        sigma_r_list.append(np.sqrt(sigma_mat_square[0,0]))   # [0,0] refers to sigma_rr
        return sigmaMat_list,sigma_r_list
    
    def optimize(self,r,xnumber,ynumber,fsky):
        self.r = r
        ls = st.params_dict['ls']; vs = st.params_dict['vs']
        ld = st.params_dict['ld']; vd = st.params_dict['vd']
#        Ad = st.params_dict['Ad']
        Ad = a_exp*np.exp(b_exp*fsky)
#        As = st.params_dict['As']
        As_78per = 0.9
        As = As_78per*np.exp(b_exp*(fsky-0.78))
        sigma_r_list=[]
        cmb_TT,cmb_EE,cmb_BB,cmb_TE = CMB(self.r).cmb() 
        par_bb,par_ee = CMB(self.r).cmb_par_r()
        sigmaMat_list = []
        if st.n_nu == 2:
            x = xnumber 
            y = ynumber
            z=1
            fmatrix = np.zeros((st.nparas,st.nparas))
            sigmaMat_list,sigma_r_list=self.cov_matrix(x,y,z,ls,vs,ld,vd,Ad,As,fmatrix,cmb_BB,par_bb,sigmaMat_list,sigma_r_list,fsky)
            return x,y,sigma_r_list
            
def main():
    fsky_list = []
    sigma_r = []
    fsky_range = np.linspace(0.0001,1,50)
    for fsky in fsky_range:
        x,y,sigma_r_list = Fisher(Add_prior='True').optimize(r,xnumber,ynumber,fsky) 
        fsky_list.append(fsky)
        sigma_r.append(sigma_r_list)
    return fsky_list,sigma_r

def plot_figure(x,y):
    plt.xlabel('fsky')
    plt.ylabel('$\sigma_r$')
    plt.loglog(x,y)
    plt.show()

if __name__ == '__main__':
    '''
    for Alicpt, r is set to 0.01. Besides, the number of detector at 95GHz and 150 GHz channels are
    both set to 3424
    '''
    r = 0.01
    xnumber = 3424
    ynumber = 3424
    start = time.time()
    print("running...")
    start = time.time()
    fsky_list,sigma_r = main()
    sigma_r_min = min(sigma_r)
    sigmar_min_index = np.argmin(np.array(sigma_r))
    fsky_min = fsky_list[sigmar_min_index]
    print('minimum sigma_r %.5f : '%np.array(sigma_r_min)[0])
    print('mininum sigma_r with fsky %.3f : '%fsky_min)
    plot_figure(fsky_list,sigma_r)
    end = time.time()
    running_time = (end - start)/60    #minutes
    print('running time: %.5f mins' %running_time) 
#    print('sigma_min_rr= %.9f '%sigma_min_rr)
#    print('################### answer ################')
#    print('sigma_min_rr= %.5f '%sigma_min_rr)
#    print(str(st.nu[0])+' GHz = %d '%num_x)
#    print(str(st.nu[1])+' GHz = %d '%num_y)




