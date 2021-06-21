import numpy as np
import gvar as gv
import math
import os

from tools import Rebin_columns, Rebin_vector
from tools import Jackknife
from tools import Jack_deviation
from tools import Filter_and_interpolate
from tools import chi2
from tools import ddt_4th
from tools import PLOT

import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 
from operator import itemgetter

from analysis import TEK_flow



wflow17 = TEK_flow('HYDRA/nfh/n289b0350k5hf1775')
# wflow13 = TEK_flow('HYDRA/nfh/n169b0350k5hf1775')
time = wflow17.time


correction = -0.6111526*(8*time/289)**2 + 0.490510*np.exp(-1.10238*289**2/(8*time)**2)

plt.plot(wflow17.t2E)
plt.plot(wflow17.t2E-correction,label='corrected')

plt.legend()
plt.show()




'''
deltaphi = gv.gvar(wflow17.t2E,wflow17.jkk_sigma(binsize=10)) - gv.gvar(wflow13.t2E,wflow13.jkk_sigma(binsize=10))

N1 = 289.
N2 = 169.
def function(t,C,A,B):
    factor = 1/N1**2 -  1/N2**2

    x = 64.*t**2

    return x*C*factor - A*(np.exp(-(N2**2)*B/x) - np.exp(-(N1**2)*B/x))

def function2(t,C):
    factor = 1/N1**2 -  1/N2**2

    x = 64.*t**2

    return x*C*factor


tmax = 150
time = time[1:tmax]
Y = np.array([i.mean for i in deltaphi])[1:tmax]
errY = np.array([i.sdev for i in deltaphi])[1:tmax]


fit = curve_fit(function2,time,Y,sigma=errY)
coeff = gv.gvar(fit[0],np.diag(fit[1]))
print(coeff)


C = fit[0]
plt.fill_between(time,Y+errY,Y-errY,alpha=0.5)
plt.plot(time,function2(time,*C))
plt.show()
'''