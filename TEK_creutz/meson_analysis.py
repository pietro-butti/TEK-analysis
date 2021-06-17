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

for R in [int(sys.argv[1])]:
    creutz,creutz_err = [],[]
    for sm in range(Nsmears):
        wMM = Rebin_vector(Wilson[:,sm,R,R]    , Nbin=nbin )
        wPP = Rebin_vector(Wilson[:,sm,R+1,R+1], Nbin=nbin )
        wPM = Rebin_vector(Wilson[:,sm,R+1,R  ], Nbin=nbin )

        ratio = wMM*wPP/(wPM**2)
        ratio = ratio[ratio>0]

        chi = -np.log(ratio)

        creutz.append( chi.mean() )
        # creutz_err.append( chi.std()/math.sqrt(len(chi)) )
        creutz_err.append( Jackknife(chi) )

    t = np.arange(Nsmears)*cf/6.

    # TRY TO FIT
    maxfit = 50
    fit = curve_fit(function,range(Nsmears)[1:maxfit],creutz[1:maxfit],sigma=creutz_err[1:maxfit])
    
    y = function(range(Nsmears)[1:],fit[0][0],fit[0][1])




    plt.plot(range(Nsmears)[1:],y)
    plt.errorbar(range(Nsmears),creutz,yerr=creutz_err,fmt='.',label=r'$\chi$('+str(R)+'.5)')

plt.grid(zorder=0, linewidth=0.5, color='gainsboro')
plt.xlabel(r'$\frac{f}{6}n_{sm}$'+' with f=0.05')
plt.legend()
plt.show()


