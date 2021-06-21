import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv

from tools import chi2
from tools import color
from tools import line
from tools import erratio
from tools import Jackknife
from tools import jstd_dev
from tools import color

########## GEARBOX ##########


KA = ['1775','1800','1825','1875','1910']
# KA = ['1910']
Op = [8]

from common import basement,path,pathj,kappaf,Njack,colorii
#############################

plt.axhline(y=0,linewidth=0.5,color='grey')
plt.xlabel(r'$\frac{1}{2\kappa_{fund}}$')
plt.ylabel(r'$am_{pcac}$')
plt.grid(zorder=0, linewidth=0.5, color='gainsboro')
#plt.figure(figsize=(3.5,7))



for ka in KA:
    print(30*'=',ka,30*'=')
                
    KF = [0.15625 if el=='1562' else float(el)/10000 for el in kappaf[ka]]
    KF = np.array(KF)

    for op in Op:
        mpcac = np.loadtxt(path[ka]+basement[ka]+ka+'_mpcac_'+str(op)+'op.dat')

        mpcacj, mpcacj_err = [], []
        for kf in kappaf[ka]:
            dataj = np.loadtxt(pathj[ka]+'coeffj_'+ka+'_'+kf+'_mpcac_'+str(op)+'op.dat').T
            mpcacj.append( dataj[0] )
            mpcacj_err.append( dataj[1] )
        mpcacj, mpcacj_err = np.array(mpcacj),np.array(mpcacj_err)

        if not len(kappaf[ka])==1:
            # Fit principal value
            fit = curve_fit(line,.5/KF,mpcac[:,1],sigma=mpcac[:,2])
            coeffs = gv.gvar(fit[0],np.diag(fit[1]))
            mc = - coeffs[1]/coeffs[0]
            kc = .5/mc

            # print(coeffs[0])


            # Fit every bin
            kcj = []
            for jj in range(len(mpcacj[0])):
                fitj = curve_fit(line,.5/KF,mpcacj[:,jj],sigma=mpcacj_err[:,jj])[0]
                kcj.append(-fitj[1]/fitj[0])
            kc = gv.gvar(kc.mean,jstd_dev(kcj))

            # Calculate chi2
            chi = chi2(mpcac[:,1],line(.5/KF,*fit[0]),mpcac[:,2])


            # Plot
            x = np.arange(3.15,.5/(KF[-1]-.01),.01)
            y = line(x,*fit[0])
            plt.plot(x,y,label=r'$\kappa_a=$0.'+str(ka)+r', $\kappa_c=$'+str(kc)+r', $\chi$='+str( '{:.2f}'.format(chi) ),color=color[colorii[ka]],alpha=0.7)

            plt.errorbar(mc.mean,0,xerr=mc.sdev,fmt='.',color='red')
            plt.errorbar(.5/KF,mpcac[:,1],yerr=mpcac[:,2],fmt='.',color=color[colorii[ka]])


        else:

            # Calculate fictitious kc
            x = np.arange(3.15,.5/(KF-.01),.01)            
            y = line(x,0.927,mpcac[1]-0.927*.5/KF)

            plt.plot(x,y,color=color[colorii[ka]],alpha=0.2)
            plt.errorbar(.5/KF,mpcac[1],yerr=mpcac[2],fmt='.',color=color[colorii[ka]],label=r'$\kappa_a=$0.'+str(ka))






plt.legend()
plt.show()