#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 23:38:35 2018

@author: wudl
"""
from __future__ import division
import os
os.chdir('/media/sf_Share/fisher_matrix/myLSPE/LSPE4')  #ATTENTION:change work dir
import matplotlib.pyplot as plt
import numpy as np
import settings
from matplotlib import colors

st = settings.Settings()

data = np.loadtxt(r'/media/sf_Share/fisher_matrix/myLSPE/data/threefre.txt')

#data = np.loadtxt(r'/media/sf_Share/fisher_matrix/myLSPE/data/threefre.txt')
z = data[:,3]   #sigma_r

an = []
i = 0
for x in range(1,st.tot_det-1,st.step):     #x_min=1,x_max=161
    for y in range(1,st.tot_det-1,st.step):
        if x+y >= st.tot_det:
            an.append(0)
        else:
            an.append(z[i])
            i += 1

x = np.arange(1,st.tot_det-1,st.step)
y = np.arange(1,st.tot_det-1,st.step)
X,Y = np.meshgrid(x,y)
Z = np.mat(an)
Z.shape = (X.shape[0],X.shape[0])
Z = Z.T

colorslist = ['w','gainsboro','gray','aqua']
#将颜色条命名为mylist，一共插值颜色条50个
#cmaps = colors.LinearSegmentedColormap.from_list('mylist',colorslist,N=200)
#cset = plt.contourf(X,Y,Z,80,cmap = cmaps) 
#cset = plt.contourf(X,Y,Z,100,cmap = plt.cm.hot)
contour = plt.contour(X,Y,Z,100,colors='k')
plt.clabel(contour,fontsize=10,colors='k',fmt='%.4f')
#plt.colorbar(cset)
#plt.xlabel(str(st.nu[0])+ ' frequency')
#plt.ylabel(str(st.nu[1])+' frequency')
plt.xlabel('$N_{95}$')
plt.ylabel('$N_{150}$')
plt.show()





















