import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv

from tools import Jackknife
from tools import jstd_dev
from tools import chi2
from tools import line
from tools import parabola
from tools import color



# ========================================================
KA = ['1775','1800','1825','1910']
# KA = [sys.argv[1]]
Op = [8]

from common import basement,path,pathj,kappaf,Njack,colorii


minpoint = 0
maxpoint = 3
# ========================================================


cc = 0
for ka in KA:
    print(30*'=',ka,30*'=')
    # basement = basements[ka]+ka

    for op in Op:    
        # Gather data ----------------------------------------------------------------------
        mpcac  = np.loadtxt(path[ka]+basement[ka]+ka+'_mpcac_'+str(op)+'op.dat')[-len(kappaf[ka]):]
        vi     = np.loadtxt(path[ka]+basement[ka]+ka+'_vi_'+str(op)+'op.dat')[-len(kappaf[ka]):]

        mpcacj = []
        vij    = []
        for kf in kappaf[ka]:
            mpcacj.append( np.loadtxt(pathj[ka]+'coeffj_'+ka+'_'+kf+'_mpcac_'+str(op)+'op.dat').T )
            vij.append(    np.loadtxt(pathj[ka]+'coeffj_'+ka+'_'+kf+'_vi_'+str(op)+'op.dat').T )
        mpcacj = np.array(mpcacj)
        vij    = np.array(vij)




        # Fit -------------- ------------------------------------------------------------------
        fit = curve_fit(line,mpcac[minpoint:maxpoint,1],vi[minpoint:maxpoint,3],sigma=vi[minpoint:maxpoint,4])[0]
        Q  = fit[1]


        Qj = [] 
        for jj in range(Njack[ka]):
            fitj = curve_fit(line,mpcacj[minpoint:maxpoint,0,jj],vij[minpoint:maxpoint,1,jj],sigma=vij[minpoint:maxpoint,2,jj])[0]
            Qj.append( fitj[1] )
        Q_sigmaj = jstd_dev(Qj)




        erry = abs(line(mpcac[:,1],*fit))*np.sqrt( (mpcac[:,2]/mpcac[:,1])**2 + (vi[:,4]/vi[:,3])**2 )
        chi = chi2(vi[:,3],line(mpcac[:,1],*fit),erry )
        print('m=',gv.gvar(Q,Q_sigmaj),'chi=',str('{:.4f}'.format(chi)))









        # Plot -------------------------------------------------------------------------------
        plt.errorbar(mpcac[:,1],vi[:,3],xerr=mpcac[:,2],yerr=vi[:,4],fmt='.',color=color[cc])#,label=r'$\kappa_a$=0.'+ka)

        x = np.arange(0,mpcac[:,1][0]+.01,0.0001)
        y = line(x,*fit)
        label = r'$\kappa_{adj}=0$.'+ka+r', $m_\rho^\chi=$'+str(gv.gvar(Q,Q_sigmaj))+r', $\chi=$'+str( '{:.3f}'.format(chi) )
        plt.plot(x,y,label=label,alpha=0.7,color=color[cc])

        cc += 2  



plt.grid(zorder=0, linewidth=0.5, color='gainsboro')
plt.axvline(x=0,linewidth=0.5,color='gray')
plt.xlabel(r'$am_{PCAC}$')
plt.ylabel(r'$am_\rho$')
plt.legend()
plt.show()