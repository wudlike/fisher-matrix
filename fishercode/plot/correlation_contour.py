from __future__ import print_function
import sys, os
sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'..')))
from getdist import plots, MCSamples
import getdist, IPython
import pylab as plt
print('GetDist Version: %s, Matplotlib version: %s'%(getdist.__version__, plt.matplotlib.__version__))
import numpy as np

def get_diagonal_element(matrix):
    # get diagonal elements of a given matrix
    mat = matrix
    for i in range(np.shape(mat)[0]):
        print("===r,As,alpha_s,beta_s,Ad,alpha_d,beta_d===")
        print(mat[i,i])
        mat[i,i] = mat[i,i]**2
    return mat

nsamp = 10000
np.random.seed(10)
# data_path = r'../../data/large/results/three/r_0.01/41_95_220mat.txt'
data_path = r'../../data/large/results/two/r_0.01/95_150mat.txt'
data = np.loadtxt(data_path)
data = get_diagonal_element(data)
cov = data
ndim = data.shape[0]
if 'large' in data_path:
    mean_ = [0.01,0.91,-2.6,-3.14,315*0.53,-2.54,1.53] # large sky with fsky = 0.71
    print('fsky:>>>> large')
else:
    mean_ = [0.01,0.081,-2.6,-3.14,20*0.53,-2.54,1.53] # small sky with fsky = 0.05
    print('fsky:>>> small')
samps = np.random.multivariate_normal(mean_, cov, size=nsamp)

names = ['r','As','alpha_s','beta_s','Ad','alpha_d','beta_d']
labels = ['r','A_s',r'\alpha_s','\\beta_s','A_d','\\alpha_d','\\beta_d']
samples = MCSamples(samples=samps,names = names, labels = labels,ranges={'r':(0,None),
                                                                         'As':(0,None),
                                                                         'Ad':(0,None)})
                                                                        
#g = plots.getSubplotPlotter()
g = plots.get_subplot_plotter()
# g.settings.axes_fontsize=19
g.settings.axes_fontsize = 18
g.settings.lab_fontsize=30
samples.updateSettings({'contours':[0.68,0.95,0.99]})
g.settings.num_plot_contours = 3
g.settings.title_limit_fontsize = 18
# samples.updateSettings({'contours': [0.68, 0.95]})
# g.settings.num_plot_contours = 2
g.triangle_plot(samples,filled=True,fontsize=20)
g.rotate_xticklabels(ax=g.subplots[6, 0], rotation=45)
g.rotate_xticklabels(ax=g.subplots[6, 1], rotation=45)
g.rotate_xticklabels(ax=g.subplots[6, 2], rotation=45)
g.rotate_xticklabels(ax=g.subplots[6, 3], rotation=45)
g.rotate_xticklabels(ax=g.subplots[6, 4], rotation=45)
g.rotate_xticklabels(ax=g.subplots[6, 5], rotation=45)
g.rotate_xticklabels(ax=g.subplots[6, 6], rotation=45)
# g.rotate_yticklabels(ax=x, rotation=45)
r_fid = data_path.split('/')[-2]
bands_str = data_path.split('/')[-1].split('m')[0]
fsky_scale = data_path.split('/')[3]
# g.export(bands_str+'_'+fsky_scale+r_fid+'_corre'+'.eps',adir='./fig')
g.export(bands_str+'_'+fsky_scale+'_corre1'+'.eps', adir='./fig')





