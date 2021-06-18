import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv

from tools import Jackknife, Rebin_vector, std_dev


def chin(nsm,a,b,c):
    return a*(1.-np.exp(-b/(nsm+c)))


def chir(r,a,b):
    cf = 0.05
    return a*(1.-np.exp(-b/(6.*cf/8.*(r*r))))

ciao = np.loadtxt('DATA/PROCESSED/n361b0340k7hf1850_creutz.dat').T


rsm = ciao[1]
for R in [int(sys.argv[1])]:
    chi = ciao[2*R]
    chi_err = ciao[2*R+1]

    minfit = int(sys.argv[2])
    maxfit = int(sys.argv[3])

    x = ciao[0][minfit:maxfit]
    y = np.array(ciao[2*R][minfit:maxfit])
    erry = np.array(ciao[2*R+1][minfit:maxfit])

    yy = gv.gvar(y,erry)
    # yy = 1/gv.log(1-yy/0.1668)


    fit = curve_fit(chin,x,y,sigma=erry,p0=[1.,1.,0.])[0]
    print(fit)


    xfit = np.arange(.01,x[-1],.01)
    yfit = chin(xfit,*fit)
    plt.errorbar(x,[i.mean for i in yy],yerr=[i.sdev for i in yy],fmt='.')
    plt.plot(xfit,yfit)

    plt.show()