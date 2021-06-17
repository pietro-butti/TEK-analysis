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


ciao = np.loadtxt('DATA/PROCESSED/n361b0340k7hf1850_creutz.dat').T

plt.errorbar(ciao[0],ciao[4],yerr=ciao[5],fmt='.')

plt.show()