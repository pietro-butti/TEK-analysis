import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv

from tools import Jackknife, Rebin_vector, std_dev


def function(nsm,a,b):
    return a*(1.-np.exp(-b/(nsm)))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
basement = 'n361b0340k7hf1850'
cf = 0.05

binsize = 48
# ---------------------------------------------------------
Wilson = np.load('DATA/PROCESSED/'+basement+'.npy')
Nwilson = Wilson.shape[2]
Nsmears = Wilson.shape[1]
Nconfs  = Wilson.shape[0]
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++



nbin = int(Nconfs/binsize)

for R in [1]:
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
        
        # yy = gv.gvar(creutz,creutz_err)/creutz[0]
        # creutz_err.append( Jackknife(chi) )

    tau = np.sqrt(8.*np.arange(Nsmears)*cf/6.)/R
    # plt.errorbar(tau,[y.mean for y in yy],yerr=[y.sdev for y in yy],fmt='.',label=r'$\chi$('+str(R+1)+'.5)')
    plt.errorbar(tau,creutz,yerr=creutz_err,fmt='.',label=r'$\chi$('+str(R+1)+'.5)')

    









plt.grid(zorder=0, linewidth=0.5, color='gainsboro')
plt.xlabel(r'$\sqrt{8\frac{f}{6}n_{sm}}$'+' with f=0.05')
plt.legend()
plt.show()


