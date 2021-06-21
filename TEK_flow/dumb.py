import numpy as np
import gvar as gv
import math
import os
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

#from analysis import TEK_flow
from tools import PLOT
from tools import erratio
from tools import errsqrt

colors = [
    'maroon',
    'firebrick',
    'salmon',
    'sandybrown',
    'orange',
    'gold',
    'yellow',
    'greenyellow',
    'seagreen',
    'skyblue',
    'cornflowerblue',
    'slateblue',
    'darkblue',
    'blue' 
]


files = [
    'DATA/nfh/n169b0350k5hf1775_norm.dat',
    'DATA/nfh/n169b0350k5hf1800_norm.dat',
    'DATA/nfh/n169b0350k5hf1825_norm.dat',
    'DATA/nfh/n169b0350k5hf1850_norm.dat',
    'DATA/nfh/n169b0350k5hf1875_norm.dat',
    
    #'DATA/nfh/n169b0347k5hf1850_norm.dat',
    #'DATA/nfh/n169b0347k5hf1875_norm.dat',
    #'DATA/nfh/n169b0347k5hf1890_norm.dat',
    #'DATA/nfh/n169b0340k5hf1875_norm.dat',
    #'DATA/nfh/n169b0340k5hf1890_norm.dat',
    #'DATA/nfh/n169b0340k5hf1910_norm.dat'
]

phimax = 0.15

cc = 0
for file in files:
    N = int(file[10:13])
    kappa = '0.'+file[22:26]

    Tmax = N/128
    tmax = int(Tmax/0.003125)
    print(Tmax)

    flow = np.loadtxt(file).T


    # Filter T<N/128
    time = flow[0][:tmax]
    t2E = flow[1][:tmax]
    t2E_err = flow[2][:tmax]

    # Filter t2E<0.12
    t2E = t2E[t2E<phimax]
    time = np.take(time,np.nonzero(t2E<phimax)[0])
    t2E_err = np.take(t2E_err,np.nonzero(t2E<phimax)[0])


    label = r'$\kappa_{adj}=$'+kappa
    plt.fill_between(time,t2E+t2E_err,t2E-t2E_err,color=colors[cc],alpha=0.7,label=label)
    cc += 1

    plt.axhline(y=phimax)
    #plt.axvline(y=phimax)

plt.legend()
#plt.show()







'''
t0 = [
    4.29373719168805,
    5.01837379013009,
    4.16861617908717, 
    5.01837379013009, 
    5.86529106537974, 
    7.21389423768470, 
    12.1386342143448,
    5.48485411950805,
    7.56495536372762,
    10.9709857447099,
    3.22339546433290,
    3.86353581627580,
    5.54202949924693
]

ii=0
for file in files:
    # Import data
    name = 'DATA/'+file+'_norm.dat'

    ciao = np.loadtxt(name,comments='#').T
    time = ciao[0]/t0[ii]
    ii+=1

    plt.plot(time,ciao[1])

plt.axvline(x=1)
plt.axhline(y=0.12)


plt.grid(zorder=0, linewidth=0.5, color='gainsboro')
plt.ylabel(r'$\frac{<t^2E(t)>}{N}$')
plt.xlabel(r'$\frac{t}{a^2}$')
plt.legend(loc="upper left")
#plt.savefig(outname)
plt.show()
'''