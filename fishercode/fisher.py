#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Jan 26 18:48:52 2018

@author: wudl

"""
import numpy as np
import settings as sts 
import time
import CMB
import foreground_models as fm
#import test4 as fm

st = sts.Settings()
fn_mat = fm.Instrumental_noise()


class Fisher(object):
    def __init__(self,Add_prior,params=st.params,param_values=st.param_values,prior_percents=st.prior_percents):
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
                          ls,vs,ld,vd,Ad,As,par_bb):
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
                    C_l[i][j] = cmb_BB[l]+fm.Fg_matrix().Sync_BB(l)[i,i]*st.sigma_F\
                    +fm.Fg_matrix().Dust_BB(l)[i,i]*st.sigma_F+fn_mat.noise(l,x,y,z)[1][i,i] # [1] for Dl_BB,[0] for Cl_BB
                    # partial Cl/ partial r
                    C_l_r[i][j] = par_bb[l]  
                    # partial Cl/ partial As
                    C_l_As[i][j] = fm.Fg_matrix().Sync_BB(l)[i,i]/As*st.sigma_F
                    # partial Cl/ partial alpha_s
                    C_l_alps[i][j] = fm.Fg_matrix().Sync_BB(l)[i,i]*np.log(l/ls)*st.sigma_F
                    # partial Cl/ partial beta_s
                    C_l_betas[i][j] = fm.Fg_matrix().Sync_BB(l)[i,i]*st.sigma_F*2*np.log(st.nu[i]/vs)
                    # for dust
                    C_l_Ad[i][j] = fm.Fg_matrix().Dust_BB(l)[i,i]/Ad*st.sigma_F
                    # partial Cl/ partial alpha_d
                    C_l_alpd[i][j] = fm.Fg_matrix().Dust_BB(l)[i,i]*np.log(l/ld)*st.sigma_F
                    # partial Cl/ partial beta_d
                    C_l_betad[i][j] = fm.Fg_matrix().Dust_BB(l)[i,i]*st.sigma_F*2*np.log(st.nu[i]/vd)                                
                else:
                    C_l[i][j] = cmb_BB[l]+np.sqrt(fm.Fg_matrix().Sync_BB(l)[i,i]*st.sigma_F*fm.Fg_matrix().Sync_BB(l)[j,j]*st.sigma_F)\
                    +np.sqrt(fm.Fg_matrix().Dust_BB(l)[i,i]*st.sigma_F*fm.Fg_matrix().Dust_BB(l)[j,j]*st.sigma_F)   
#                    C_l[j][i] = C_l[i][j]
                    # partial Cl/ partial r
                    C_l_r[i][j] = par_bb[l]   
#                    C_l_r[j][i] = C_l_r[i][j]
                    # partial Cl/ partial As
                    C_l_As[i][j] = np.sqrt(fm.Fg_matrix().Sync_BB(l)[i,i]*st.sigma_F*fm.Fg_matrix().Sync_BB(l)[j,j]*st.sigma_F)/As                                
#                    C_l_As[j][i] = C_l_As[i][j]
                    # partial Cl/ partial alpha_s
                    C_l_alps[i][j] = np.sqrt(fm.Fg_matrix().Sync_BB(l)[i,i]*st.sigma_F*fm.Fg_matrix().Sync_BB(l)[j,j]*st.sigma_F)*np.log(l/ls)
#                    C_l_alps[j][i] = C_l_alps[i][j]
                    # partial Cl/ partial beta_s
                    C_l_betas[i][j] = np.sqrt(fm.Fg_matrix().Sync_BB(l)[i,i]*st.sigma_F*fm.Fg_matrix().Sync_BB(l)[j,j]*st.sigma_F)*np.log(st.nu[i]*st.nu[j]/(vs*vs))     
#                    C_l_betas[j][i] =C_l_betas[i][j]
                    # for dust
                    C_l_Ad[i][j] = np.sqrt(fm.Fg_matrix().Dust_BB(l)[i,i]*st.sigma_F*fm.Fg_matrix().Dust_BB(l)[j,j]*st.sigma_F)/Ad                                
#                    C_l_Ad[j][i] = C_l_Ad[i][j]
                    # partial Cl/ partial alpha_d
                    C_l_alpd[i][j] = np.sqrt(fm.Fg_matrix().Dust_BB(l)[i,i]*st.sigma_F*fm.Fg_matrix().Dust_BB(l)[j,j]*st.sigma_F)*np.log(l/ld)
#                    C_l_alpd[j][i] = C_l_alpd[i][j]
                    # partial Cl/ partial beta_d
                    C_l_betad[i][j] = np.sqrt(fm.Fg_matrix().Dust_BB(l)[i,i]*st.sigma_F*fm.Fg_matrix().Dust_BB(l)[j,j]*st.sigma_F)*np.log(st.nu[i]*st.nu[j]/(vd*vd))  
#                    C_l_betad[j][i] = C_l_betad[i][j]
        return C_l,C_l_r, C_l_As, C_l_alps, C_l_betas, C_l_Ad, C_l_alpd, C_l_betad
    
    def cov_matrix(self,x,y,z,ls,vs,ld,vd,Ad,As,fmatrix,cmb_BB,par_bb,sigmaMat_list,sigma_r_list):
        for l in range(st.l_min,st.l_max): 
            C_l = np.eye(len(st.nu)); C_l_r = np.eye(len(st.nu))
            C_l_As = np.eye(len(st.nu)); C_l_alps = np.eye(len(st.nu)); C_l_betas = np.eye(len(st.nu))
            C_l_Ad = np.eye(len(st.nu)); C_l_alpd = np.eye(len(st.nu)); C_l_betad = np.eye(len(st.nu))
                
            C_l,C_l_r, C_l_As, C_l_alps, C_l_betas, C_l_Ad, C_l_alpd, C_l_betad = \
            Fisher(self.Add_prior).mat_derivative(l,x,y,z,C_l,C_l_r,C_l_As,C_l_alps,C_l_betas,C_l_Ad,C_l_alpd,C_l_betad,cmb_BB,\
                  ls,vs,ld,vd,Ad,As,par_bb)
            # A list of matrix where the order of element is r As alphas betas Ad alphad betad                    
            C_l_inv = np.linalg.inv(C_l)
            par_C_l_list=[C_l_r, C_l_As, C_l_alps, C_l_betas, C_l_Ad, C_l_alpd, C_l_betad]
            for i in range(len(par_C_l_list)):
                for j in range(len(par_C_l_list)):
                    fmatrix[i][j] += (2*l+1)/2*st.fsky*np.trace(np.dot(np.dot(C_l_inv,par_C_l_list[i]),np.dot(C_l_inv,par_C_l_list[j])))
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
    
    def optimize(self,config):
        ls = st.params_dict['ls']; vs = st.params_dict['vs']
        ld = st.params_dict['ld']; vd = st.params_dict['vd']
        Ad = st.params_dict['Ad']
        As = st.params_dict['As']
        x_list = []
        y_list = []
        z_list = []
        m_list = []
        sigma_r_list=[]
        cmb_TT,cmb_EE,cmb_BB,cmb_TE = CMB.cmb(config)
        par_bb,par_ee = CMB.cmb_par_r()
        sigmaMat_list = []
        if st.n_nu == 2:
            x = input('Number of detectors at first frequency')             
            y = input('Number of detectors at second frequency')
            z = 1        
            fmatrix = np.zeros((st.nparas,st.nparas))
            sigmaMat_list,sigma_r_list=self.cov_matrix(x,y,z,ls,vs,ld,vd,Ad,As,fmatrix,cmb_BB,par_bb,sigmaMat_list,sigma_r_list)       
        elif st.n_nu == 3:
            for x in range(1,st.tot_det-1,st.step):
                for y in range(1,st.tot_det-x,st.step):
                    z = st.tot_det-x-y
                    x_list.append(x)
                    y_list.append(y)
                    z_list.append(z)               
                    fmatrix = np.zeros((st.nparas,st.nparas))
                    sigmaMat_list,sigma_r_list=self.cov_matrix(x,y,z,ls,vs,ld,vd,Ad,As,fmatrix,cmb_BB,par_bb,sigmaMat_list,sigma_r_list) 
        elif st.n_nu == 4:
            for x in range(1,st.tot_det-(st.n_nu-1)+1,st.step):
                for y in range(1,st.tot_det-(st.n_nu-2)-x+1,st.step):
                    for z in range(1,st.tot_det-(st.n_nu-3)-x-y+1,st.step):
                        m = st.tot_det-x-y-z
                        x_list.append(x)
                        y_list.append(y)
                        z_list.append(z)
                        m_list.append(m)                        
                        fmatrix = np.zeros((st.nparas,st.nparas)) 
                        sigmaMat_list,sigma_r_list=self.cov_matrix(x,y,z,ls,vs,ld,vd,Ad,As,fmatrix,cmb_BB,par_bb,sigmaMat_list,sigma_r_list)
        index_min,sigma_rr_list,sigma_min_rr,sigmaMat_min=self.find_min_index(sigmaMat_list)
        num_x=x_list[index_min]
        num_y=y_list[index_min]
        num_z=z_list[index_min]
        # print np.sqrt for sigma and sigmaMat
        sigmaMat_sqrt_diag = self.sqrt(sigmaMat_min)
        if st.n_nu == 4:
            num_m = st.tot_det-num_x-num_y-num_z
            return sigma_r_list,num_x,num_y,num_z,num_m,x_list,y_list,z_list,m_list,np.sqrt(sigma_rr_list), np.sqrt(sigma_min_rr), sigmaMat_min, sigmaMat_sqrt_diag
        else:
            return sigma_r_list,num_x,num_y,num_z,x_list,y_list,z_list,np.sqrt(sigma_rr_list), np.sqrt(sigma_min_rr), sigmaMat_min, sigmaMat_sqrt_diag

if __name__ == '__main__':
    st = sts.Settings()
    start = time.time()
    print("running...")
    start = time.time()
    num_x,num_y,sigma_rr_list, sigma_min_rr, sigmaMat_min,sigmaMat_sqrt_diag = Fisher().optimize()                    
    end = time.time()
    running_time = (end - start)/60    #minutes
    print('running time: %.5f mins' %running_time) 
    print('sigma_min_rr= %.9f '%sigma_min_rr)
    print('################### answer ################')
    print('sigma_min_rr= %.5f '%sigma_min_rr)
    print(str(st.nu[0])+' GHz = %d '%num_x)
    print(str(st.nu[1])+' GHz = %d '%num_y)




