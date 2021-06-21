import numpy as np
import gvar as gv
import math
import os

from tools import Rebin_columns, Rebin_vector
from tools import Jackknife
from tools import Jack_deviation
from tools import Filter_and_interpolate
from tools import chi2
from tools import ddt_4th
from tools import PLOT

import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 
from operator import itemgetter

from analysis import TEK_flow


# Cexp = gv.gvar(['-0.6111526(50)', '0.490510(93)', '1.10238(24)'])
# Cexp = gv.gvar(['-0.6851051(71)'])

# Cexp = gv.gvar(['-0.4845019(37)', '0.509330(48)', '1.03854(11)'])
# Cexp = gv.gvar(['-0.5344843(80)'])

Cexp = gv.gvar(['-0.3884966(27)', '0.481108(35)', '1.067388(95)'])
# Cexp = gv.gvar(['-0.4430612(18)'])


wflow = TEK_flow('HYDRA/nfh/n289b0350k5hf1825')
wflow.fit_correction([c.mean for c in Cexp],289)
t0 = wflow.find_t0(0.1)
t1 = wflow.find_t0(0.05)
print(gv.sqrt(8*t0),gv.sqrt(8*t1))

wflow = TEK_flow('HYDRA/nfh/n289b0350k5hf1825')
wflow.norm_correction(289)
t0 = wflow.find_t0(0.1)
t1 = wflow.find_t0(0.05)
print(gv.sqrt(8*t0),gv.sqrt(8*t1))



