import numpy as np
import math
import sys
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import gvar as gv


color=[
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


N_op = 4
def index(a,b): 
    return a*N_op+b

def corr(x,A,B):
    Nhalf = 17
    return 2.*A*np.exp(-Nhalf*B)*np.cosh((Nhalf-x)*B)

def line(x,m,q):
    return m*x+q

def parabola(x,a,b,c):
    return a*x*x + b*x + c

# This take an array, and calculate jackknife error
def Jackknife(vector):
    # 'None' filtering
    vector = [i for i in vector if i!=None]
    vector = np.array(vector)

    N = len(vector)
    if N<2: 
        print('ERROR: len(vector)<2. Returning None value...')
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
    if N<2: 
        print('ERROR: len(vector)<2. Returning None value...')
        return None

    mean = vector.mean()
    return math.sqrt(sum( (vector-mean)**2 )/float(N-1))
    
def jstd_dev(vector):
    # 'None' filtering
    vector = [i for i in vector if i!=None]
    vector = np.array(vector)

    N = len(vector)
    if N<2: 
        print('ERROR: len(vector)<2. Returning None value...')
        return None

    mean = vector.mean()
    return math.sqrt(sum( (vector-mean)**2 )/float(N-1)*float(N))





def chi2(exp,th,sigma,reduce=True,Npar=2):
    exp = np.array(exp)
    th  = np.array(th)
    sigma = np.array(sigma)

    chi2 = (exp-th)/sigma
    chi2 *= chi2
    chi2 = sum(chi2)

    if reduce:
      chi2 /= (exp.size - Npar)

    return math.sqrt(chi2)

def erratio(A,B,sigmaA,sigmaB):
    sigmaAB = 0 
    if isinstance(A,list) and isinstance(B,list):
        A = np.array(A)
        B = np.array(B)

        sigmaAB = sum((A-A.mean())*(B-B.mean()))/float(len(A)-1)

    return np.sqrt( (sigmaA/A)**2 + (sigmaB/B)**2 - 2*sigmaAB/(A*B) )*abs(A/B)


def errsqrt(A,sigmaA):
    if isinstance(A, list):
        return [abs(sigmaA[ii]/(2*math.sqrt(A[ii]))) for ii in range(len(A))]
    else:
        return abs(sigmaA/(2.*np.sqrt(A)))

