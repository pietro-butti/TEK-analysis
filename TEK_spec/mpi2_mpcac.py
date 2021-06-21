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
from tools import parabola
from tools import erratio
from tools import Jackknife
from tools import jstd_dev
from tools import color



########## GEARBOX ##########


# KA = ['1775','1800','1825']
KA = ['1910']
Op = [8]

from common import basement,path,pathj,kappaf,Njack,colorii


minpoint = 0
maxpoint = -1
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
        # Extract data
        mpcac = np.loadtxt(path[ka]+basement[ka]+ka+'_mpcac_'+str(op)+'op.dat')
        pp    = np.loadtxt(path[ka]+basement[ka]+ka+'_pp_'+str(op)+'op.dat')
        mpi2 = gv.gvar(pp[:,3],pp[:,4])**2

        mpcacj,ppj=[],[]
        for kf in kappaf[ka]:
            mpcacj.append(np.loadtxt(pathj[ka]+'coeffj_'+ka+'_'+kf+'_mpcac_'+str(op)+'op.dat').T)
            ppj.append(   np.loadtxt(pathj[ka]+'coeffj_'+ka+'_'+kf+'_pp_'+str(op)+'op.dat').T )
        mpcacj = np.array(mpcacj)
        ppj    = np.array(ppj)  




        # Fit
        if maxpoint<0: maxpoint=len(mpi2)+1


        y    = [i.mean for i in mpi2[minpoint:maxpoint]]
        erry = [i.sdev for i in mpi2[minpoint:maxpoint]]
        fit = curve_fit(line,mpcac[minpoint:maxpoint,1],y,sigma=erry)
        fit2 = curve_fit(parabola,mpcac[minpoint:maxpoint,1],y,sigma=erry)


        mj,qj = [],[]
        for jj in range(Njack[ka]):
            mpi2j = gv.gvar(
                ppj[minpoint:maxpoint,1,jj],
                ppj[minpoint:maxpoint,2,jj]
            )
            yj = [i.mean for i in mpi2j**2]
            erryj = [i.sdev for i in mpi2j**2]
            fitj = curve_fit(line,mpcacj[minpoint:maxpoint,0,jj],yj,sigma=erryj)
            
            mj.append( fitj[0][0] )
            qj.append( fitj[0][1] )        
        errm = jstd_dev(mj) 
        errq = jstd_dev(qj) 

        chi = chi2(pp[minpoint:maxpoint,3]**2,line(np.array(mpcac[minpoint:maxpoint,1]),*fit[0]),erry)
        chi_par = chi2(pp[minpoint:maxpoint,3]**2,parabola(np.array(mpcac[minpoint:maxpoint,1]),*fit2[0]),erry,Npar=3)
        print('m=',gv.gvar(fit[0][0],errm),', q=',gv.gvar(fit[0][1],errq),'chi=',str('{:.4f}'.format(chi)))




        plt.errorbar(mpcac[minpoint:maxpoint,1],y,yerr=erry,fmt='.',color=color[colorii[ka]])
        x = np.arange(0,mpcac[0,1]+0.02,0.001)
        
        plt.plot(x,line(x,*fit[0]),label=r'$\kappa_a=$0.'+str(ka)+r', $\chi$='+str( '{:.2f}'.format(chi) ),color=color[colorii[ka]+2],alpha=0.7)
        plt.plot(x,parabola(x,*fit2[0]),label=r'$\kappa_a=$0.'+str(ka)+r', $\chi$='+str( '{:.2f}'.format(chi_par) ),color=color[colorii[ka]+1],alpha=0.2)



plt.grid(zorder=0, linewidth=0.5, color='gainsboro')
plt.axvline(x=0,linewidth=0.5,color='gray')
plt.xlabel(r'$am_{PCAC}$')
plt.ylabel(r'$am_\pi^2$')
plt.legend()
plt.show()    