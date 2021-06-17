import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv






def Rebin_vector(vector,binsize=None,Nbin=None):
    Nobs = len(vector)

    # Rebin at fixed binsize
    if not binsize==None:
        Nbb = int(Nobs/binsize)
        rest = Nobs%binsize
        if rest==0: Nbin=Nbb
        else: Nbin=Nbb+1

        new = []
        for bb in range(Nbin):
            bin = vector[bb*binsize:bb*binsize+binsize]#-1]
            new.append(  np.array(bin).mean()  )
        if not rest==0:
            bin = vector[bb*binsize:]
            new.append( np.array(bin).mean() )

    # Rebin at fixed Nbin
    elif not Nbin==None:
        R = Nobs/Nbin
        binsize = math.floor(R)
        rest = Nobs%Nbin
        
        if rest==0: Nbb=Nbin
        else: Nbb=Nbin-1

        if not binsize==1:
            new = []
            for bb in range(Nbb):
                bin = vector[bb*binsize:bb*binsize+binsize-1]
                new.append( np.array(bin).mean() )
        else:
            new = vector[:Nbb].tolist()

        if not rest==0:
            bin = vector[-rest:]
            new.append( np.array(bin).mean() )

    else:
        print("LOCAL error in Rebin_vector(): You have to assign Nbin or binsize")
        exit()



    return np.array(new)







# This take an array, and calculate jackknife error
def Jackknife(data,binsize=None):

    if binsize==None: vector=data
    else: vector = Rebin_vector(data,binsize=binsize)


    # 'None' filtering
    vector = [i for i in vector if i!=None]
    vector = np.array(vector)

    N = len(vector)
    if N<2: 
        print('ERROR: Few values. Returning None value...')
        return None

    mean = vector.mean()
    jk_means = []
    for ii in range(N):
        jk_means.append( np.delete(vector,ii).mean() )

    jkk = sum((jk_means-mean)**2)*float(N-1)/float(N)
    return math.sqrt(jkk)






def std_dev(vector):
    # 'None' filtering
    vector = [i for i in vector if i!=None]
    vector = np.array(vector)

    N = len(vector)
    mean = vector.mean()

    return math.sqrt( sum((vector-mean)**2)/float(N*(N-1)) )




