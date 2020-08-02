#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 17:21:11 2018

@author: wudl

If you run this file, it will give you an ImportError--cannot import name 'Switch'.
The reason is you 'from run import Switch as sw' here while you 'import settings as sts '
in run.py. 
So you cross call the file between settings.py and run.py
Therefore, there is a conflict here.
But anyway, you still can implement run.py file, the answer are correct.
"""
import numpy as np
#from run import Switch as sw
from run import config

class Settings(object):
    
    def __init__(self): 
#        self.NET = 350    #uk.sqrt(second)
        self.l_min = 50
        self.l_max = 600
        self.tot_det = 10000  #total number of detector
        self.step = 200   #step for optimization detector number
#        self.Tobs = 2*6*30*12*60*60      #2 years
        
        self.Tobs = 1*12*30*24*60*60          
#        self.r = 0
#        self.r = 0.1
#        self.r = 0.01
        self.r = config['r']
        self.arcmin_to_radian = np.pi/(180.0*60.0)
        self.deg_to_radian = np.pi/(180)
        self.arcmin_to_deg = 1/60 
        self.kB = 1.3806505*10**(-23)  # J/K
        self.h = 6.6260693*10**(-34+9) # J·s, dimension becomes J/GHz by plus 9
        self.c = 299792458
        self.Tcmb = 2.73
        '''
        priors for the parameters r, As, alphas, betas, Ad, alhpad, betad in percentage, 
        where 0 below refers to no-prior.
        '''
        self.params = ['r','As','alphas','betas','Ad','alphad','betad']      
        self.params_dict = {'r': self.r,
                    'alpha_s': -2.6,
                    'betas': -3.14,
                    'alpha_d': -2.54,
                    'betad': 1.53,
                    'Td': 19.6,
                    'vd': 353,
                    'ld': 80,
                    'vs': 30,
                    'ls': 80}
        if config['sky_fraction']=='small':
            self.fsky = 0.05
            self.params_dict['As']=0.081
            self.params_dict['Ad']=20*0.53
        elif config['sky_fraction']=='large':
            self.fsky = 0.71
            self.params_dict['As']=0.91
            self.params_dict['Ad']=315*0.53                      
        self.param_values = [float(self.r), self.params_dict['As'],self.params_dict['alpha_s'],self.params_dict['betas'],\
                self.params_dict['Ad'],self.params_dict['alpha_d'],self.params_dict['betad']]
        self.prior_percents = [0, 0.17, 0, 0.054, 0.015, 0, 0.013] 
#        self.prior_percents = [0, 0., 0, 0., 0., 0, 0.] 
        self.nparas = len(self.param_values)  # r, As, alphas, betas, Ad, alphad, betad
        self.skyam = 4.0*np.pi*self.fsky*((180.0*60.0)/np.pi)**2  #arcmin**2
##        self.sigma_F = 0
##        self.sigma_F = 0.01
        self.sigma_F = 1
        self.freq_dict = {
                        '30 GHz':{
                                'FWHMs': 59.9,     # arcmin
                                'NET': 474,           #
                                 },                 
                        '41 GHz':{
                                'FWHMs': 43.8,     # arcmin
                                'NET': 388,           #
                                 },                
                        '85 GHz':{
                                'FWHMs': 22.2,     # arcmin
                                'NET': 308,           #279
                                 },
                        '95 GHz':{
                                'FWHMs': 19.2,
                                'NET': 279,             #315
                                  },                
                        '105 GHz':{
                                'FWHMs': 17.4,     # arcmin
                                'NET': 301,           #
                                 },                
                        '135 GHz':{
                                'FWHMs': 14.4,     # arcmin
                                'NET': 320,           #279
                                 },
                        '145 GHz':{
                                'FWHMs': 12.6,
                                'NET': 309,             #315
                                  },
                        '150 GHz':{
                                'FWHMs': 12.0,     # arcmin
                                'NET': 315,           #
                                 },                
                        '155 GHz':{
                                'FWHMs': 11.4,     # arcmin
                                'NET': 326,           #279
                                 },
                        '220 GHz':{
                                'FWHMs': 8.4,
                                'NET': 575,             #315
                                  },                     
                        '270 GHz':{
                                'FWHMs': 6.65,
                                'NET': 1318,             
                                  },                                   
                         }

        self.nu1 = config['nu1']
        self.nu2 = config['nu2']
        self.nu3 = config['nu3']
        self.nu4 = config['nu4']
        self.nu5 = config['nu5']
        if self.nu2 == 0:
            self.nu = np.array([self.nu1])
        elif self.nu3 == 0:
            self.nu = np.array([self.nu1,self.nu2])
        elif self.nu4 == 0:
            self.nu = np.array([self.nu1,self.nu2,self.nu3])
        elif self.nu5 == 0:
            self.nu = np.array([self.nu1,self.nu2,self.nu3,self.nu4])
        else:
            self.nu = np.array([self.nu1,self.nu2,self.nu3,self.nu4,self.nu5])
        self.n_nu = len(self.nu)   #number of parameters
        
def test(config):
    '''
    Testing for whether changing to small or large sky
    Testing for frequency channel
    '''
    test_fsky = Settings().fsky
    print('Testing for sky fraction and frequencies with\nfsky= %s : %.2f'%(config['sky_fraction'],test_fsky))
    print(Settings().nu)
    print('r = : %.2f'%float(config['r']))
    print('##########################################') 



        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
