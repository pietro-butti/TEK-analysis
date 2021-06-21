import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv



# Prendiamo il caso di n361b0350k5hf1850 e proviamo a vedere come scala sul raggio

ciao = np.loadtxt('DATA/PROCESSED/n361b0340k7hf1850_creutz.dat').T

R2 = gv.gvar(ciao[4],ciao[5])#/0.1665577
R3 = gv.gvar(ciao[6],ciao[7])#/0.1053759
R4 = gv.gvar(ciao[8],ciao[9])#/0.0807069
R5 = gv.gvar(ciao[10],ciao[11])#/0.072541155
R6 = gv.gvar(ciao[12],ciao[13])#/0.072541155

plt.errorbar(ciao[1]/2,[i.mean for i in R2],yerr=[i.sdev for i in R2],fmt='.')
plt.errorbar(ciao[1]/3,[i.mean for i in R3],yerr=[i.sdev for i in R3],fmt='.')
plt.errorbar(ciao[1]/4,[i.mean for i in R4],yerr=[i.sdev for i in R4],fmt='.')
plt.errorbar(ciao[1]/5,[i.mean for i in R5],yerr=[i.sdev for i in R5],fmt='.')
plt.errorbar(ciao[1]/6,[i.mean for i in R6],yerr=[i.sdev for i in R6],fmt='.')


plt.show()

