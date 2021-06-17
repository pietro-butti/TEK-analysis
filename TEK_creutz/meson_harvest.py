import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv





# +++++++++++++++++++++++++++++++++++++++++++++++++
basement = 'n361b0340k7hf1850'
folder = 'DATA/RAW/'+basement+'/'

Nwilson = 9
Nsmears = 120
Nconfs  = 480
# +++++++++++++++++++++++++++++++++++++++++++++++++


Wilson = np.zeros((Nconfs,Nsmears,Nwilson,Nwilson))


toharvest = sorted(os.listdir(folder))

conf = 0
for file in toharvest:
    print(file)
    filename = folder+file

    f = open(filename,'r')
    for line in f.readlines():
        if line.split()[0]=='nsm,isq,jsq=':
            sm = int(line.split()[1])
            R  = int(line.split()[2])
            L  = int(line.split()[3])
            Wilson[conf,sm-1,R-1,L-1] = float(line.split()[5])
            
    conf += 1


np.save('DATA/PROCESSED/'+basement,Wilson)

