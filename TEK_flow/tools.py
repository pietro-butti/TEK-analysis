import numpy as np
import gvar as gv
import math
import os
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 

# This takes an array and rebin the elements
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
            bin = vector[bb*binsize:bb*binsize+binsize-1]
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



    return new


# This take a numpy matrix (vector of columns), cycles over T and rebin elements
def Rebin_columns(data,binsize=None,Nbin=None):
    if binsize==None and Nbin==None:
        print("ERROR: Select binsize or Nbin")
        return
    elif binsize!=None and Nbin!=None:
        print("ERROR: Select just one option")
        return

    copy = []
    for col in data.T:
        if not binsize==None: copy.append(  Rebin_vector(col,binsize=binsize)  )
        if not    Nbin==None: copy.append(  Rebin_vector(col,Nbin=Nbin)  )
    copy = np.array(copy).T

    return copy

# This take an array, and calculate jackknife error
def Jackknife(vector):
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

def Jack_deviation(vector):
    # 'None' filtering
    vector = [i for i in vector if i!=None]
    vector = np.array(vector)

    N = len(vector)
    if N<2: 
        print('ERROR: Few values. Returning None value...')
        return None

    mean = vector.mean()
    jkk = float(N-1)/float(N)*sum((vector-mean)**2)
    return math.sqrt(jkk)



def line(x,m,q):
    return m*x + q

def interpolate(x,y,ref):
    coeff = curve_fit(line,x,y)

    t0 = (ref-coeff[0][1])/coeff[0][0]
    return t0

def Filter_and_interpolate(xdata,ydata,refval,delta):
    ydata = np.array(ydata)
    xdata = np.array(xdata)

    # Filter: extract data in a proper range
    ymin = refval-delta
    ymax = refval+delta

    y_filtered = ydata[ydata>ymin]
    x_filtered = np.take(xdata,np.nonzero(ydata>ymin)[0])

    y_filtered = y_filtered[y_filtered<ymax]
    x_filtered = np.take(x_filtered,np.nonzero(y_filtered<ymax)[0])

    # Check if data can be interpolated
    flag = False
    if y_filtered.size>5: 
        prima = y_filtered[0]
        for el in y_filtered[1:]:
            dt = el - prima
            if dt<0:
                flag = True
                break
            prima = el
    else: flag = True

    
    if flag==True:
        print('*',end='')
    else:
        t0 = interpolate(x_filtered,y_filtered,refval)
        return t0
    print('')
    

def ddt_4th(array,deltat):
    # Calculate derivative of t^2E(t) at 4th order following
    # f'[i] = (f[i-2] - 8*f[i-1] + 8*f[i+1] - f[i+2])/(12*h)

    factor = 1./(12.*deltat)
    ddt = []
    for t in range(2,len(array)-2):
        ddt.append(
            factor*( array[t-2]-8.*array[t-1]+8.*array[t+1]-array[t+2] )
        )
    ddt.insert(0,ddt[0])
    ddt.insert(0,ddt[0])
    ddt.append(ddt[-1])
    ddt.append(ddt[-1])

    return ddt



def erratio(A,B,sigmaA,sigmaB):
    R = abs(A/B)
    err = (sigmaA/A)**2 + (sigmaB/B)**2
    err = [math.sqrt(el) for el in err]
    err *= R

    return err


def errsqrt(A,sigmaA):
    if isinstance(A, list):
        return [abs(sigmaA[ii]/(2*math.sqrt(A[ii]))) for ii in range(len(A))]
    else:
        return abs(sigmaA/(2.*np.sqrt(A)))

def chi2(exp,th,sigma,Npar=3):
    exp = np.array(exp)
    th  = np.array(th)
    sigma = np.array(sigma)

    chi2 = (exp-th)/sigma
    chi2 *= chi2
    chi2 = sum(chi2)
    chi2 /= (exp.size - Npar)

    return chi2


def PLOT(filename,color,label=None,t0=None):
    time = np.loadtxt(filename,comments="#",usecols=0)
    if not t0==None: time=time/t0
    data = np.loadtxt(filename,comments="#",usecols=1)
    jerr = np.loadtxt(filename,comments="#",usecols=2)
    plt.fill_between(time,data+jerr,data-jerr,color=color,label=label,alpha=0.7)








'''
class TEK_flow_Ninf:
    # RETURNED
    # e.g. ---> self.flowN['289'] : is a TEK_flow class, trimmed to have same Ntimes and Nconfs
    # self.basement, self.N, self.strN
    def __init__(self,N,basement):
        self.N = N
        self.strN = [str(n) for n in self.N]
        self.basement = {}
        self.shapes = {}
        
        self.flowN = {}
        self.flowN_inf = {}
        self.flowN_inf_data = {}
        self.TEK_flow = {}
        self.deltaPhi = {}
        self.deltaPhi_sigmaj = {}
        self.deltaPhi_data = {}

        self.C = {}
        self.C_sigmaj = {}

        # Create TEK_flow classes for every N
        for ii in range(len(N)):
            self.basement[self.strN[ii]] = basement[ii]
            self.flowN[self.strN[ii]] = TEK_flow(basement=basement[ii])
            self.shapes[self.strN[ii]] = self.flowN[self.strN[ii]].data.shape

        # Reshape to have same time extent
        times = []
        Nconf = []
        for nn in self.strN:
            times.append( self.shapes[nn][1] )
            Nconf.append( self.shapes[nn][0] )
        timemax  = min(times)
        Nconfmax = min(Nconf)
        self.finalshape = (Nconfmax,timemax)
        self.Nbin_max = int(self.finalshape[0]/2)


        for nn in self.strN:
            if not self.flowN[nn].data.shape[1]==timemax:
                self.flowN[nn].data = self.flowN[nn].data[:,:timemax]
                self.flowN[nn].t2E = self.flowN[nn].t2E[:timemax]
            else:
                self.time = self.flowN[nn].time

    # RETURNED
    # e.g. ---> self.deltaPhi['289-169']
    # e.g. ---> self.deltaPhi_data['289-169']: self.deltaPhi but rebinned with jackknife for different Nbins (averaged over every element but one)
    # e.g. ---> self.deltaPhi_sigmaj['289-169']
    def deltaPhi_jk(self,N1,N2):
        n1 = str(N1)
        n2 = str(N2)
        key = n1+'-'+n2

        # Calculate flow difference for averages
        self.deltaPhi[key] = self.flowN[n1].t2E - self.flowN[n2].t2E

        # Calculate flow difference for fixed Nbin
        t2E_sigmaj = []
        self.deltaPhi_data[key] = []
        for Nbin in range(10,self.Nbin_max,10):
            # Reshape
            newdataN1 = Rebin_columns(self.flowN[n1].data,Nbin=Nbin)
            newdataN2 = Rebin_columns(self.flowN[n2].data,Nbin=Nbin)

            phi = []
            for conf in range(len(newdataN1)):
                phi.append( newdataN1[conf] - newdataN2[conf] )
            phi = np.array(phi)
        
            # Compute average difference in every bin but one
            phi_jkk = []
            for skip_this in range(len(phi)):
                newphi = np.delete(phi,skip_this,0)

                phi_b = []
                for tt in newphi.T:
                    phi_b.append( np.array(tt).mean() )
                phi_jkk.append( phi_b )
            phi_jkk = np.array(phi_jkk)
            self.deltaPhi_data[key].append( phi_jkk )

            # At fixed t find jackknife error on t^2E(t)
            t2E_sigma = []
            for tt in range(self.finalshape[1]):
                t2E_sigma.append( Jack_deviation(phi_jkk[:,tt]) )
            
            t2E_sigmaj.append(t2E_sigma)
        t2E_sigmaj = np.array(t2E_sigmaj)

        err = []
        for T in t2E_sigmaj.T:
            err.append(np.array(T).mean())
        self.deltaPhi_sigmaj[key] = np.array(err)


    def pplot(self):
        key = '169-289'

        factor = 64*(1/169**2 - 1/289**2)*0.361

        y = self.deltaPhi[key] + self.time*self.time*factor + 10
        plt.plot([math.log(i) for i in y])
        plt.show()

        
    # RETURNED
    # e.g. ---> self.C['289-169']
    # e.g. ---> self.C_sigmaj['289-169']
    def fit_jkk(self,N1,N2,Tmin=2,Tmax=10,exp=False,deltat=0.03125):
        key = str(N1)+'-'+str(N2)

        # Define function
        def function(t,*C):
            if not exp:
                return 64.*C[0]*(1/N2**2-1/N1**2)*t*t
            else:
                return 64.*C[0]*(1/N2**2-1/N1**2)*t*t + C[1]*( np.exp(-C[2]*N1/(8.*t)) - np.exp(-C[2]*N2/(8.*t))   )

        # Set initial guess for fit
        if not exp: p0=.5
        else: p0=[.5,.5,.5]

        self.delta = deltat
        self.Tmin = int(Tmin/self.delta) 
        self.Tmax = int(Tmax/self.delta)
        # Fit averaged value
        time  = self.time[self.Tmin:self.Tmax]
        Phi   = self.deltaPhi[key][self.Tmin:self.Tmax]
        meanC = curve_fit(function,time,Phi,p0=p0)[0]
        self.C[key] = meanC

        # Cycle over bin
        jkk = []
        for newdata in self.deltaPhi_data[key]:
            jk = []
            for bin in newdata:
                c = curve_fit(function,time,bin[self.Tmin:self.Tmax],p0=p0)[0]
                # Filter too big values ------------------------------
                flag = False
                for el in c: 
                    if abs(el)>10.: 
                        flag=True
                if not flag: 
                    jk.append( c )
                # ----------------------------------------------------
            jkk.append( [Jack_deviation(np.array(jk)[:,ii]) for ii in range(len(meanC))] )
        jkk = np.array(jkk)

        self.C_sigmaj[key] = []
        for col in jkk.T:
            self.C_sigmaj[key].append( col.mean() )


        # Calculate square chi
        fit = function(self.time[self.Tmin:self.Tmax],*meanC)
        chi2 = sum(  ((fit-self.deltaPhi[key][self.Tmin:self.Tmax])/self.deltaPhi_sigmaj[key][self.Tmin:self.Tmax])**2 )
        chi2 /= len(self.time[self.Tmin:self.Tmax]-len(meanC))


        final = gv.gvar(self.C[key],self.C_sigmaj[key])
        print('Fit result in [',int(self.delta*self.Tmin),'-',int(self.delta*self.Tmax),'] ---> ',final,'{ chi2 = ',chi2,'}')

    def function(self,t,N1,N2,*C):
        howmany = len(C)
        if howmany==1:
            return 64.*C[0]*(1/N2**2-1/N1**2)*t*t
        elif howmany==3:
            return 64.*C[0]*(1/N2**2-1/N1**2)*t*t + C[1]*( np.exp(-C[2]*N1/(8.*t)) - np.exp(-C[2]*N2/(8.*t))   )
    def plot(self,key,color=None,alpha=0.55,label=None):
        if color==None: color='firebrick'
        if label==None: label=key

        y = self.deltaPhi[key]
        erry = self.deltaPhi_sigmaj[key]
        plt.fill_between(self.time,y-erry,y+erry,color=color,alpha=alpha)

        x = np.arange(self.Tmin*0.03125,self.Tmax*0.03125,0.03125)
        y = self.function(x,169,289,*self.C[key])
        plt.plot(x,y,color=color)

    # RETURNED
    # e.g. ---> self.flowN_inf_data[('169-289','289')] 
    # e.g. ---> self.flowN_inf[('169-289','289')]     it is a TEK_flow initialized class
    def correct_flow(self,key,N,C=None,errC=None):     
        if C==None: 
            C = self.C[key]
            errC = self.C_sigmaj[key]

        # Compute correction and error
        correction = C[0]*(self.time*8./N)**2
        if len(C)==3:
            exp = np.insert(np.exp( -C[2]/((self.time[1:]*8./N)**2) ),0,0.)
            correction += -C[1]*exp

        # Sum correction
        kkey = (key,str(N))
        self.flowN_inf_data[kkey] = []

        for flow in self.flowN[str(N)].data:
            self.flowN_inf_data[kkey].append(flow + np.array(correction))
        self.flowN_inf_data[kkey] = np.array(self.flowN_inf_data[kkey])

        self.flowN_inf[kkey] = TEK_flow(datacopy=self.flowN_inf_data[kkey])

    # First argument is kkey = ('169-289','289')
    def do_things(self,kkey,export=False,refval=None):
        # Do things with corrected flow
        if export:
            self.flowN_inf[kkey].jkk_sigma()
            filename = self.flowN[kkey[1]].basement+'_inf'+kkey[0]+'.dat'
            self.flowN_inf[kkey].export(filename)
        
        if not refval==None:
            self.flowN_inf[kkey].find_t0(refval)
            #print(kkey,refval,self.flowN_inf[kkey].t0[str(refval)],self.flowN_inf[kkey].t0_sigmaj[str(refval)])
            
    def ciao(self,kkey):
        t0 = self.flowN_inf[kkey].t0['0.1']
        t0_err = self.flowN_inf[kkey].t0_sigmaj['0.1']
        t1 = self.flowN_inf[kkey].t0['0.05']
        t1_err = self.flowN_inf[kkey].t0_sigmaj['0.05']
        print('%15.14f      %15.14f     %15.14f      %15.14f'%(t0,t0_err,t1,t1_err))
        
            
    def log(self):
        for nn in self.strN:
            print(30*'-')
            print('N of colors --->',nn)
            print('basements --->',self.basement[nn])
            print('Matrix original shape --->',self.shapes[nn])
        print(30*'-')
        print('Unique matrix final shape is --->',self.finalshape)
'''




