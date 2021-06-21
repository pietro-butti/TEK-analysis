import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv

from tools import Jackknife
from tools import jstd_dev


#######################################################################
from common import tag,kadj,basement,path,outtag,Njack,kappa,Nop
#+++++++++++++++++++++++++++++++++
pa_mod = -1

tmin = 4
tmax = 8

# --------------------------------------------

offset = 2
Ntimes = 16 - offset

minfit = tmin-offset
maxfit = tmax-offset+1
# --------------------------------------------
# FLAGS
debug_mode= False

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



MPCAC = {}
for kf in kappa[tag]:
    for op in Nop:
        key = (str(op),kf)
        print(30*'=',key,30*'=')

        # GATHER OUTPUT ----------------------------------------------------
        filempcac  = os.path.join(path,'Op'+str(op)+'_t22_t11/'+basement+'_'+kf+'_mpcacdl')
        filempcacj = os.path.join(path,'Op'+str(op)+'_t22_t11/'+basement+'_'+kf+'_jmpcacdl')

        # data   = np.loadtxt(filempcac,comments='#').T
        # ddataj = np.loadtxt(filempcacj,comments='#').T

        data   = np.load(filempcac+'.npy').T
        ddataj = np.load(filempcacj+'.npy').T


        # Change the sign to the effective mass (only true for the new version MesonFundamental-r880)
        data[1]   *= pa_mod
        ddataj[1] *= pa_mod

        # Metti in un array i diversi bin dataj[jj]
        dataj,errj = [],[]
        for jj in range(Njack):
            dataj.append(ddataj[1][jj*(Ntimes+1):jj*(Ntimes+1)+Ntimes+1])
            errj.append(ddataj[2][jj*(Ntimes+1):jj*(Ntimes+1)+Ntimes+1])
        dataj,errj = np.array(dataj), np.array(errj)



        # FIT (weigthed mean over 1/sigma^2) ------------------------------
        err2 = data[2][minfit:maxfit]**2
        mpcac = sum(data[1][minfit:maxfit]/err2)/sum(1/err2)
        mpcac_sigma = math.sqrt(1/sum(1/err2))

        mpcacj, mpcacj_err = [], []
        for jj in range(Njack):
            errj2 = errj[jj][minfit:maxfit]**2
            mpcacj.append(
                sum(dataj[jj][minfit:maxfit]/errj2)/sum(1/errj2)
            )
            mpcacj_err.append(  math.sqrt(1/sum(1./errj2))   )
        mpcac_sigmaj = jstd_dev(np.array(mpcacj))

        
        # Calculate with unbiased estimator
        # mpcac = mpcac - (Njack-1)*(mpcac - np.array(mpcacj).mean())
        MPCAC[key] = [float(kf)/10000,mpcac,mpcac_sigmaj,mpcac_sigma]




        # Export bins separately --------------------------------
        if exportjj:
            if not os.path.isdir('OUTPUT/JACK/'+outtag): os.mkdir('OUTPUT/JACK/'+outtag)

            f = open('OUTPUT/JACK/'+outtag+'/coeffj_'+str(kadj)+'_'+str(kf)+'_mpcac_'+str(op)+'op.dat','w+')
            print("#%10s    %10s"%('mpcac','err'),file=f)
            for jj in range(Njack):
                print(" %10.9f    %10.9f"%(mpcacj[jj],mpcacj_err[jj]),file=f)
            f.close()
        # -------------------------------------------------------

        # Plot --------------------------------------------------
        if plot:
            label=r'$\kappa_f=0$.'+kf+r', op='+str(op)+r'$, m_{PCAC}=$'+str(gv.gvar(MPCAC[key][1],MPCAC[key][2])) 
            # label=r'$\kappa_f=0$.'+kf+r', op='+str(op)+r'$, m_{PCAC}=$'+str('{:.4f}'.format(MPCAC[key][1]))+'('+str('{:.4f}'.format(MPCAC[key][2])+')')

            plt.errorbar(data[0][:-7],data[1][:-7],yerr=data[2][:-7],fmt='.')
            plt.axvline(x=minfit+offset,linewidth=0.5,color='gray')
            plt.axvline(x=maxfit+offset-1,linewidth=0.5,color='gray')

            x = np.arange(minfit+offset,maxfit+offset-1,.1)
            #label=r'$m_{PCAC}=$'+str(gv.gvar(MPCAC[key][1],MPCAC[key][2]))
            plt.fill_between(x,mpcac-mpcac_sigma,mpcac+mpcac_sigma,alpha=0.5,label=label)

            plt.title(r'$\kappa_a=$0.'+kadj)
        # -------------------------------------------------------

if plot:
    plt.legend()
    #plt.savefig('PLOTS/mpcac_'+kappa_adj+'_'+k)
    # plt.savefi('PLOTS/comparison/mpcac_'+basement[1:4])
    plt.show()

if export:
    for op in Nop:
        if not os.path.isdir('OUTPUT/'+outtag): os.mkdir('OUTPUT/'+outtag)
        fileout = 'OUTPUT/'+outtag+'/'+basement+'_mpcac_'+str(op)+'op.dat'
        print(30*'-',fileout,30*'-')

        access_mode='w+'
        if append: access_mode='a'

        f = open(fileout,access_mode)
        if not append:
            print("# %5s    %20s     %20s     %20s"%("kappa","mpcac","sigma_jack","sigma"),file=f)
        
        for k in kappa[tag]:
            key = (str(op),k)
   
            Mpcac = MPCAC[key][1]
            Mpcac_errj = MPCAC[key][2]
            Mpcac_err = MPCAC[key][3]
            # Export
            print("%1.4f    %20.19f    %20.19f    %20.19f"%(float(k)/10000,Mpcac,Mpcac_errj,Mpcac_err),file=f)
            print("%1.4f    %20.19f    %20.19f    %20.19f"%(float(k)/10000,Mpcac,Mpcac_errj,Mpcac_err))
    f.close()

