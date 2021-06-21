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


KA = ['1775','1800','1825','1910']
# KA = ['1800']
Op = [8]

basement = {
    '1775':'n289b0350k5hf',
    '1800':'n289b0350k5hf',
    '1825':'n289b0350k5hf',
    '1910':'n289b0340k5hf'
}

path = {
    '1775':'OUTPUT/1775_17bin/',
    '1800':'OUTPUT/1800_17bin/',
    '1825':'OUTPUT/1825_17bin/',
    '1910':'OUTPUT/1910_half/'
}
pathj = {
    '1775':'OUTPUT/JACK/1775_17bin/',
    '1800':'OUTPUT/JACK/1800_17bin/',
    '1825':'OUTPUT/JACK/1825_17bin/',
    '1910':'OUTPUT/JACK/1910_half/'
}
kappaf = {
    '1775':['1500','1525','1550','1562'],
    '1800':['1470','1500','1525','1550','1562'],
    '1825':['1470','1500','1525','1550','1558'],
    '1910':['1570']
}
Njack =  {
    '1775':17,
    '1800':17,
    '1825':17,
    '1910':20
}

colorii = {
    '1775':0,
    '1800':4,
    '1825':10,
    '1910':11
}
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
            plt.errorbar(.5/KF,mpcac[1],yerr=mpcac[2],fmt='.',color=color[colorii[ka]],label=r'b=0.34, $\kappa_a=$0.'+str(ka))






plt.legend()
plt.show()