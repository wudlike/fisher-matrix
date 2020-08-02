#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 16:27:59 2018

@author: wudl
"""
import numpy as np
def FWHM():
#    fre = 41   #GHz
    fre = float(input('Please enter frequency: '))
    c = 299792458
    lam = c/(fre*10**9) #unit: m
    D = 0.70  #unit: m
    fwhm = 1.22*lam/D #unit: rad 
    fwhm = fwhm*180/np.pi
    return fwhm

if __name__ =='__main__':
    fwhm = FWHM()
    print('fwhm=>>> %.2f(unit:degree)'%fwhm)
    print('fwhm=>>> %.2f(unit:arcmin)'%(fwhm*60))