import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv

from tools import Jackknife, Rebin_vector, std_dev


def function(nsm,a,b,c):
    return a*(1.-np.exp(-b/(nsm+c)))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
basement = 'n361b0340k7hf1850'
cf = 0.05

binsize = 20

Nwilson = 9
RR = [0,1,2,3,4,5]
# RR = [int(sys.argv[1])]
# ---------------------------------------------------------
Wilson = np.load('DATA/PROCESSED/'+basement+'.npy')
Nwilson = Wilson.shape[2]
Nsmears = Wilson.shape[1]
Nconfs  = Wilson.shape[0]
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++



nbin = int(Nconfs/binsize)
CHI = np.zeros((Nsmears,Nwilson,2))



for R in RR:
    creutz,creutz_err = [],[]
    for sm in range(Nsmears):
        wMM = Rebin_vector(Wilson[:,sm,R,R]    , Nbin=nbin )
        wPP = Rebin_vector(Wilson[:,sm,R+1,R+1], Nbin=nbin )
        wPM = Rebin_vector(Wilson[:,sm,R+1,R  ], Nbin=nbin )

        ratio = wMM*wPP/(wPM**2)
        ratio = ratio[ratio>0]
      
        chi = -np.log(ratio)

        creutz.append( chi.mean() )
        creutz_err.append( chi.std()/math.sqrt(len(chi)) )

        CHI[sm,R,0] = chi.mean()
        CHI[sm,R,1] = chi.std()/math.sqrt(len(chi))
        
        # yy = gv.gvar(creutz,creutz_err)/creutz[0]
        # creutz_err.append( Jackknife(chi) )


    tau = np.sqrt(np.arange(Nsmears)*cf/6.)/(R+1)
    plt.errorbar(np.arange(Nsmears),creutz,yerr=creutz_err,fmt='.',label=r'$\chi$('+str(R+1)+'.5)')

f = open('DATA/PROCESSED/'+basement+'_creutz.dat','w')
print("# This is the outcome of meson_analysis: cretuz ratios with errors chi(1.5), chi(2.5), chi(3.5), chi(4.5), chi(5.5) ....",file=f)
print("# Wilson loops are binned (binsize=48), than ratios are calculated, ratios<0 are filtered (less than 10 per nsm), than creutz ratios are calculated",file=f)
print("# nsm  sqrt(8*t)",file=f)

for sm in range(Nsmears):
    print("%1.3i  %10.9f "%(10*sm,math.sqrt(sm/6.*cf)), end='          ',file=f)
    for ii in range(len(RR)):
        print("%10.9f  %10.9f"%(CHI[sm,:,0][ii],CHI[sm,:,1][ii]),end='        ',file=f)
    print('',file=f)
