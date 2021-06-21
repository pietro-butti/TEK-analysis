import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv





# +++++++++++++++++++++++++++++++++++++++++++++++++
basement = 'n289b0350k5hf1775'
# basement = 'n361b0340k7hf1850'
folder = 'DATA/RAW/'+basement+'/'

Nwilson = 9
Nsmears = 120
Nconfs  = 900
# +++++++++++++++++++++++++++++++++++++++++++++++++


Wilson = np.zeros((Nconfs,Nsmears,Nwilson,Nwilson))
for ii in range(Nconfs):
    for jj in range(Nsmears):
        for kk in range(Nwilson):
            for ll in range(Nwilson):
                Wilson[ii,jj,kk,ll] = None



toharvest = sorted(os.listdir(folder))

conf = 0
for file in toharvest:
    filename = folder+file

    if os.path.isdir(filename):
        print(file,'is a directory...Skipping')    
    else:
        print(file)
        f = open(filename,'r')
        for line in f.readlines():
            if line.split()[0]=='nsm,isq,jsq=':
                sm = int(int(line.split()[1]))
                R  = int(line.split()[2])
                L  = int(line.split()[3])
                # print(float(line.split()[5]))
                Wilson[conf,sm-1,R-1,L-1] = float(line.split()[5])
        conf += 1


print("Datas are going to be saved in "+'DATA/PROCESSED/'+basement+'.npy')
np.save('DATA/PROCESSED/'+basement,Wilson)

