import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv

from tools import index
from tools import Jackknife
from tools import jstd_dev
from tools import std_dev
from tools import chi2
from tools import color

def corr(x,A,B):
    Nhalf = 17
    return 2.*A*np.exp(-Nhalf*B)*np.cosh((Nhalf-x)*B)





#######################################################################
tag = sys.argv[1]

kadj = tag
basement = 'n289b0340k5hf'+kadj
path = 'DATA/b34half'
outtag = tag+'_half'
Njack = 20
#+++++++ IMPORTANT ++++++++++++
kappa = {
    '1775':['1500','1525','1550','1562'],
    '1800':['1470','1500','1525','1550','1562'],
    '1825':['1470','1500','1525','1550','1558'],
    '1910':['1570']
}

# Nop = [4,5,6,7,8,9]
Nop = [8]
#++++++++++++++++++++++++++++++

obs = 'pp'

tmin = 4
tmax = 8
# --------------------------------------------
Njack_jump = 16

Ntimes = 34

minfit = tmin
maxfit = tmax+1
# --------------------------------------------


# FLAGS
debug_mode = False

plot     = True
export   = True
exportjj = True
append   = False

if debug_mode:
    plot = False
    export = False
    append = False
    exportjj = False
#######################################################################


c = 0
MASSES = {}
for kf in kappa[tag]:
    for op in Nop:
        key = (obs,str(op),kf)
        print(30*'=',key,30*'=')





        # GATHER DATA ------------------------------------------------------------------
        file  = path+'/Op'+str(op)+'_t22_t11/'+basement+'_'+kf+'_'+obs+'_cxdl'
        filej = path+'/Op'+str(op)+'_t22_t11/'+basement+'_'+kf+'_'+obs+'_jcxdl'  

        data  = np.loadtxt(file,comments='#').T
        ddataj = np.loadtxt(filej,comments='#').T

        # These are the average data
        time = data[2,:Ntimes+1]
        Ct  = data[3,:Ntimes+1]
        err = data[4,:Ntimes+1]

        '''
        # Calculate effective mass
        CC = gv.gvar(Ct[:15],err[:15])
        meff = []
        for ii in range(1,len(CC)):
            meff.append( CC[ii-1]/CC[ii] )
        meff = gv.log(meff)
        # plt.errorbar(range(len(meff)),[i.mean for i in meff],yerr=[i.sdev for i in meff],fmt='.')
        # plt.show()
        '''

        # Ctj[ii] is correlator in each bin
        Ctj, Ctj_err = [], []
        for jj in range(Njack):
            off = jj*Njack_jump*(Ntimes+1)
            Ctj.append( ddataj[3][off:off+Ntimes+1] )
            Ctj_err.append( ddataj[4][off:off+Ntimes+1] )
        Ctj, Ctj_err = np.array(Ctj), np.array(Ctj_err)





        # FIT -------------------------------------------------------------------------
        # Principal value
        mass = curve_fit(corr,time[minfit:maxfit],Ct[minfit:maxfit],sigma=err[minfit:maxfit],p0=[.5,.5])[0]

        # Jackknife bins
        massj, massj_err = [], []
        for jj in range(Njack):
            fitt = curve_fit(corr,time[minfit:maxfit],Ctj[jj,minfit:maxfit],sigma=Ctj_err[jj,minfit:maxfit],p0=[.5,.5])
            massj.append(fitt[0])
            massj_err.append(np.diag(fitt[1]))
        massj, massj_err = np.array(massj), np.array(massj_err)
        mass_sigmaj = [jstd_dev(massj[:,0]),jstd_dev(massj[:,1])]

        # # Calculate chi2 in each bin
        # for jj in range(Njack):
        #     print(
        #         chi2( Ctj[jj,minfit:maxfit],corr(time[minfit:maxfit],*massj[jj]),err[minfit:maxfit] )
        #     )



        # Calculate unbiased estimator
        MASSES[key] = [float(kf)/10000,mass[0],mass_sigmaj[0],mass[1],mass_sigmaj[1]]
        chi = chi2(Ct[minfit:maxfit],corr(time[minfit:maxfit],*mass),err[minfit:maxfit])
        # mass = [mass[i] - (Njack-1)*(mass[i]-massj[:,i].mean()) for i in [0,1]]

        print(kf,op,'C = '+str(gv.gvar(mass[0],mass_sigmaj[0])),', m = '+str(gv.gvar(mass[1],mass_sigmaj[1])),str(chi)[:4],'    ')
        # print('  &  ',str(gv.gvar(mass[1],mass_sigmaj[1])),'  &  ',str(chi)[:4])










        

        # EXPORT BIN SEPARATELY --------------------------------------------------------------
        if exportjj:
            if not os.path.isdir('OUTPUT/JACK/'+outtag): os.mkdir('OUTPUT/JACK/'+outtag)
            f = open('OUTPUT/JACK/'+outtag+'/coeffj_'+str(kadj)+'_'+str(kf)+'_'+obs+'_'+str(op)+'op.dat','w+')
            print("#%10s    %10s    %10s"%('C','m','sigmaj_m'),file=f)
            for jj in range(Njack):
                print(" %10.9f   %10.9f   %10.9f"%(massj[jj,0],massj[jj,1],massj_err[jj,1]),file=f)
            f.close()
        # ------------------------------------------------------------------------------------

        # PLOT ------------------------------------------------------------------------
        if plot:
            label = r'$\kappa_f=0$.'+kf + ', op = '+str(op)+', m = '+str(gv.gvar(mass[1],mass_sigmaj[1]))+r', $\chi^2$ = '+str(chi)[:5]
            # label = r'$\kappa_f=0$.'+kf + ', op = '+str(op)+', m = '+str('{:.4f}'.format(mass[1])+'('+'{:.4f}'.format(mass_sigmaj[1]))+r'), $\chi^2$ = '+str(chi)[:5]
            plt.errorbar(time[:18],Ct[:18],yerr=err[:18],fmt='.',color=color[c])

            timeplot = np.arange(0,18,0.1)
            yplot = corr(timeplot,*mass)
            plt.plot(timeplot,yplot,label=label,alpha=0.5,color=color[c])

            plt.axvline(x=minfit,linewidth=0.5,color='gray')
            plt.axvline(x=maxfit-1,linewidth=0.5,color='gray')

            c+=2
        # -----------------------------------------------------------------------------

if plot:
    plt.title(r'$\kappa_a=$0.'+str(kadj))
    plt.yscale('log')
    plt.legend()
    #plt.legend(fontsize='small')
    plt.savefig('PLOTS/comparison/vi_'+basement[1:4])
    plt.show()

if export:
    if not os.path.isdir('OUTPUT/'+outtag): os.mkdir('OUTPUT/'+outtag)
    for op in Nop:
        fileout = 'OUTPUT/'+outtag+'/'+basement+'_'+obs+'_'+str(op)+'op.dat'
        print('#',30*'-',fileout,30*'-')

        access_mode='w+'
        if append: access_mode='a'

        f = open(fileout,access_mode)
        if not append:
            print("# %5s      %20s   %20s         %20s   %20s"%("kappa","C","sigmaj","m","sigmaj"),file=f)
        
        for k in kappa[tag]:
            key = (obs,str(op),k)
            
            C = MASSES[key][1]
            C_err = MASSES[key][2]
            M = MASSES[key][3]
            M_err = MASSES[key][4]
            # Export
            print("%1.4f      %20.19f   %20.19f         %20.19f   %20.19f"%(float(k)/10000,C,C_err,M,M_err),file=f)
            # print("%1.4f      %20.19f   %20.19f         %20.19f   %20.19f"%(float(k)/10000,C,C_err,M,M_err))
            print(float(k)/10000,gv.gvar(M,M_err))
    f.close()