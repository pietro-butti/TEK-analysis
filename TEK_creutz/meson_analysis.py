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

binsize = 48

RR = [0,1,2,3,4,5]
# ---------------------------------------------------------
Wilson = np.load('DATA/PROCESSED/'+basement+'.npy')
Nwilson = Wilson.shape[2]
Nsmears = Wilson.shape[1]
Nconfs  = Wilson.shape[0]
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++



nbin = int(Nconfs/binsize)

CHI = np.zeros((Nsmears,len(RR),2))

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

    tau = np.sqrt(8.*np.arange(Nsmears)*cf/6.)/(R+1)
    plt.errorbar(np.arange(Nsmears),creutz,yerr=creutz_err,fmt='.',label=r'$\chi$('+str(R+1)+'.5)')

f = open('DATA/PROCESSED/'+basement+'_creutz.dat','w')
for sm in range(Nsmears):
    print("%1.3i  %10.9f "%(sm,math.sqrt(8.*sm/6.*cf)), end='          ',file=f)
    for ii in range(len(RR)):
        print("%10.9f  %10.9f"%(CHI[sm,:,0][ii],CHI[sm,:,1][ii]),end='        ',file=f)
    print('',file=f)




#     # FIT
#     minfit = 0
#     maxfit = 20
    
#     fit = curve_fit(function,np.arange(Nsmears)[minfit:maxfit],creutz[minfit:maxfit],sigma=creutz_err[minfit:maxfit])
#     coeff = gv.gvar(fit[0],np.diag(fit[1]))

#     yy = function(np.arange(Nsmears)[minfit:maxfit],*fit[0])
#     plt.plot(np.arange(Nsmears)[minfit:maxfit],yy)






# plt.grid(zorder=0, linewidth=0.5, color='gainsboro')
# plt.xlabel(r'$\sqrt{8\frac{f}{6}n_{sm}}$'+' with f=0.05')
# plt.legend()
# plt.show()


