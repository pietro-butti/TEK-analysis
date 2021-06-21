import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv

from tools import Jackknife, Rebin_vector, std_dev

cf = 0.05
def chin(x,a,b,c):
    return a*(   1.-np.exp( -b/(x+c) )   )


ciao = np.loadtxt('DATA/PROCESSED/n361b0340k7hf1850_creutz.dat').T


rsm = np.array(ciao[1])
for R in [int(sys.argv[1])]:
    chi = ciao[2*R]
    chi_err = ciao[2*R+1]

    if len(sys.argv)==4:
        minfit = int(sys.argv[2])
        maxfit = int(sys.argv[3])

        x = ciao[0][minfit:maxfit]/R
        y = np.array(ciao[2*R][minfit:maxfit])
        erry = np.array(ciao[2*R+1][minfit:maxfit])


        fit,errfit = curve_fit(chin,x,y,sigma=erry,p0=[0.1668,40,0.1])
        print(minfit,maxfit,gv.gvar(fit,np.diag(errfit)))
        xfit = np.arange(.01,x[-1],.01)
        yfit = chin(xfit,*fit)
        plt.plot(xfit,yfit,label=str(gv.gvar(fit,np.diag(errfit))))
        plt.errorbar(ciao[0][:maxfit],ciao[2*R][:maxfit],yerr=ciao[2*R+1][:maxfit],fmt='.',label='$\chi$('+str(R)+'.5)')
        plt.axvline(minfit,color='gray',alpha=0.7)
        plt.axvline(maxfit,color='gray',alpha=0.7)

        plt.grid(zorder=0, linewidth=0.5, color='gainsboro')
        plt.legend()
        plt.savefig('PLOTS/'+'chi('+str(R)+',5)')
        plt.show()

    else:
        plt.errorbar(ciao[0],ciao[2*R],yerr=ciao[2*R+1],fmt='.',label='$\chi$('+str(R)+'.5)')
        plt.show()
