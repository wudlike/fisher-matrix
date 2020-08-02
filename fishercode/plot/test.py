from __future__ import print_function
import sys, os
sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'..')))
from getdist import plots, MCSamples
import getdist, IPython
import pylab as plt
print('GetDist Version: %s, Matplotlib version: %s'%(getdist.__version__, plt.matplotlib.__version__))
import numpy as np
#ndim = 4

def get_diagonal_element(matrix):
    # get diagonal elements of a given matrix
    mat = matrix
    for i in range(np.shape(mat)[0]):
        mat[i,i] = mat[i,i]**2
    return mat

nsamp = 10000
np.random.seed(10)
#A = np.random.rand(ndim,ndim)
#cov = np.dot(A, A.T)
#data = np.loadtxt(r'C:\Users\wudl\Desktop\fisher\data\small\results\two\r_0.01\95_150mat.txt')
data = np.loadtxt(r'..\..\data\large\results\two\r_0.01\95_150mat.txt')
#data = np.loadtxt(r'C:\Users\Lenovo\Desktop\recent work\revise paper\data\small\results\two\r_0.01\95_150cov.txt')
data = get_diagonal_element(data)
#data.reshape(7,7)
cov = data
ndim = data.shape[0]
mean_ = [0.01,0.082,-2.6,-3.14,20*0.53,-2.54,1.53]
samps = np.random.multivariate_normal(mean_, cov, size=nsamp)

#samps = np.random.multivariate_normal([0]*ndim, cov, size=nsamp)
#A = np.random.rand(ndim,ndim)
#cov = np.dot(A, A.T)
#samps2 = np.random.multivariate_normal([0]*ndim, cov, size=nsamp)
#names = ['r','As','$\alpha_s$','$\beta_s$','Ad','$\alpha_d$','$\beta_d$']
#labels = ['r','As','$\alpha_s$','$\beta_s$','Ad','$\alpha_d$','$\beta_d$']
names = ['r','As','alpha_s','beta_s','Ad','alpha_d','beta_d']
#labels = ['r','As','alpha_s','beta_s','Ad','alpha_d','beta_d']
labels = ['r','A_s','\\alpha_s','beta_s','Ad','alpha_d','beta_d']
samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'r':(0,None),
                                                                         'As':(0,None)})
#samples2 = MCSamples(samples=samps2,names = names, labels = labels, label='Second set')

g = plots.getSubplotPlotter()
samples.updateSettings({'contours':[0.68,0.95,0.99]})
g.settings.num_plot_contours = 3
g.triangle_plot(samples,filled=True)
#g.triangle_plot([samples, samples2], filled=True)
g.export('output_file.pdf')
