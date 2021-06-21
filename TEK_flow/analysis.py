import numpy as np
import gvar as gv
import math
import os

from tools import Rebin_columns
from tools import Jackknife
from tools import Jack_deviation
from tools import Filter_and_interpolate
from tools import chi2
from tools import ddt_4th
from tools import PLOT

import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 
from operator import itemgetter


class TEK_flow:
    # RETURNED:
    # self.time
    # self.t2E (self.t2E is <t^2E>)
    # self.data (matrix of <t^2E> for every configuration)
    def __init__(self,basement=None,datacopy=None,deltat=0.03125,usecols=None):
        self.deltat = deltat
        self.t0 = {}
        self.t0_sigmaj = {}

        if usecols==None: col = 3
        else: col = usecols

        if not basement==None:
            self.basement = basement
            self.filename = basement+'.dat'
            # Gather data in a matrix
            data = []
            for file in os.listdir(basement):
                file = os.path.join(basement,file)
                data.append( np.loadtxt(file,comments='#',usecols=col).tolist() )
        elif isinstance(datacopy,np.ndarray):
            data = datacopy
            self.basement = 'unknown'

        # ----------------------------------------------------------------------
        # Count how many timesteps are there and reshape in case
        Nsteps = [len(el) for el in data]
        self.time = np.arange(0,Nsteps[0]*deltat,deltat)

        if not all(i==Nsteps[0] for i in Nsteps):
            tmax = min(Nsteps)
            self.time = self.time[:tmax]
            for ii in range(len(data)):
                data[ii] = data[ii][:tmax]
        self.data = np.array(data)

        self.Tmin = 0
        self.Tmax = len(self.time)

        # # Calculate mean
        self.t2E = []
        for rows in self.data.T:
            self.t2E.append( rows.mean() )
        self.t2E = np.array(self.t2E)


        # ----------------------------------------------------------------------
        # Function attributes
        self.find_t0_has_been_called = False
        self.jkk_sigma_has_been_called = False
        self.norm_correction_has_been_called = False

    # RETURNED:
    # self.NofT                                             (= N(t))
    # self.factor                                           (= C(0) = N(0)/N(t))
    # self.ddt_factor                                       (= d/dt(C(0)) )
    # self.data and self.t2E are multiplied by self.factor  (if correct_data=True)
    def norm_correction(self,N,correct_data=True):
        Ninf = 3./(3.1415926535897932384626433832795*3.1415926535897932384626433832795*128)
        NofT = np.loadtxt('normalizations/norm_n'+str(N)+'.dat',comments="#",usecols=1)
        factor = [Ninf/x for x in NofT if x!=0]
        factor.insert(0,0.)
            
        if len(factor)!=len(self.time):
            factor=factor[:len(self.time)]
            NofT = NofT[:len(self.time)]
        self.factor = np.array(factor)
        self.NofT = np.array(NofT)

        # Calculate derivative of the correction factor
        self.ddt_factor = np.array( ddt_4th(self.factor,self.deltat) )

        # Correct data if requested
        self.data_are_corrected = correct_data
        if correct_data:
            self.t2E *= self.factor
            for data in self.data:
                data *= self.factor
        
        self.norm_correction_has_been_called = True
            

    # RETURNED:
    # self.t2E_sigmaj
    def jkk_sigma(self):
        print('Harvesting') 
        jkk = []
        for binsize in range(10,21):
            newdata = Rebin_columns(self.data,binsize=binsize)

            jk = []
            for col in newdata.T:
                jk.append( Jackknife(col) )
            jkk.append( jk )
        jkk = np.array(jkk).T
        bin = []
        for el in jkk: bin.append(el.mean())

        # Export
        self.t2E_sigmaj = bin
        print('Jkk completed')
        self.jkk_sigma_has_been_called = True
        return self.t2E_sigmaj
    
    # RETURNED:
    # self.data_w
    # self.wt2E
    # self.wt2E_sigmaj
    def w_ddt(self):
        data_w = []
        for data in self.data:
            data_w.append( self.time *ddt_4th(data,self.deltat) )
        self.data_w = np.array(data_w)

        if self.norm_correction_has_been_called:
            if not self.data_are_corrected:
                self.t2E *= self.factor
                # Compute W(t) as:
                # t d/dt(C(0) * t^2E(t)) = t*d/dt(C(t))*(t^2E(t)) + t*C(t)*d/dt(t^2E(t))
                ddt_data = []
                for confN in range(len(self.data)):
                    ddt_data.append(
                        self.time*self.ddt_factor*self.data[confN] + self.factor*self.data_w[confN]
                    )
                self.data_w = np.array(ddt_data)


        # Calculate average
        self.wt2E = []
        for T in self.data_w.T:
            self.wt2E.append(T.mean())
        self.wt2E = np.array(self.wt2E)

        # Calculate jackknife error
        jkk = []
        for binsize in range(10,21):
            newdata = Rebin_columns(self.data_w,binsize=binsize)

            jk = []
            for col in newdata.T:
                jk.append( Jackknife(col) )
            jkk.append( jk )
        jkk = np.array(jkk).T
        bin = []
        for el in jkk: bin.append(el.mean())
        
        self.wt2E_sigmaj = bin
        '''
        '''

    def export(self,filename=None,w=False):
        if filename==None:
            f = open(self.filename,"w")
        else:
            f = open(filename,"w")

        if not w:
            print("# t    %20s    %20s   "%("t2 <E>","sigma jack"),file=f)
            for flowt in range(self.Tmin,self.Tmax):
                print(" %7.8f    %20.19f    %20.19f   "%(self.time[flowt],self.t2E[flowt],self.t2E_sigmaj[flowt]),file=f)
        else:
            print("# t    %20s    %20s     %20s    %20s"%("t2 <E>","sigma jack","t d/dt(t2 <E>)","sigma jack"),file=f)
            for flowt in range(self.Tmin,self.Tmax):
                print(" %7.8f    %20.19f    %20.19f   %20.19f    %20.19f   "%(self.time[flowt],self.t2E[flowt],self.t2E_sigmaj[flowt],self.wt2E[flowt],self.wt2E_sigmaj[flowt]),file=f)
        f.close()
    
    # This takes self.data (which has (Nconf,Ntimes)) and reshape it in (newconfN,Ntimes)
    # EXPLICIT
    # return data_reshaped
    def reshape(self,newconfN):
        return Rebin_columns( self.data,Nbin=newconfN )

    ############## t0 ANALYSIS ################
    # RETURNED
    # e.g. ---> self.t0[('t','0.1')]           # fist key can also be 'w' if w argument is set to True
    # e.g. ---> self.t0_sigmaj[('t','0.1')]
    def find_t0(self,refval,delta=0.01,w=False):
        if not w:  
            self.find_t0_has_been_called=True
            t2E = self.t2E
            selfdata = self.data
            key = ('t',str(refval))
        else:      
            t2E = self.wt2E
            selfdata = self.data_w
            key = ('w',str(refval))


        # Find t0 principal value
        T0 = Filter_and_interpolate(self.time,t2E,refval,delta=delta)
        self.t0[key] = T0

        # JACKKNIFE ERROR EVALUATION
        t0sigma = []
        for binsize in range(10,21):
            # Rebin configurations
            data = Rebin_columns(selfdata,binsize=binsize)
            Nbin = data.shape[0]

            # Compute t2E averaging all bins but one
            t2E_jkk = []
            for skip_this in range(np.size(data,0)):
                newdata = np.delete(data,skip_this,0)

                t2E_b = []
                for tt in newdata.T:
                    t2E_b.append( np.array(tt).mean() )
                t2E_jkk.append( t2E_b )


            # Compute t0 in every bin
            t0 = []
            for bin in t2E_jkk:
                t0.append( Filter_and_interpolate(self.time,bin,refval,delta) )

            # None filtering
            t0 = np.array([i for i in t0 if not i==None])
            if len(t0)==0: print('no t0')

            # Calculate dispersion
            mean = np.array(t0).mean()
            jkk = sum((t0-mean)**2)*float(len(t0)-1)/float(len(t0))
            jkk = math.sqrt(jkk)
            t0sigma.append(jkk)
        
        # Export
        self.t0_sigmaj[key] = np.array(t0sigma).mean()
        print(self.basement,': t0(',refval,') ---> ',gv.gvar(self.t0[key],self.t0_sigmaj[key]))
        return gv.gvar(self.t0[key],self.t0_sigmaj[key])
        
    def print_t0(self,filename,kappa,w=False):
        obs = 't'
        if w: obs = 'w'

        g = open(filename,'a+')
        #print('%.4f     %15.14f      %15.14f'%(kappa,self.t0['0.1'],self.t0_sigmaj['0.1']),file=g)
        #print('%4f     %15.14f      %15.14f     %15.14f      %15.14f'%(kappa,self.t0['0.1'],self.t0_sigmaj['0.1'],self.t0['0.05'],self.t0_sigmaj['0.05']),file=g)
        print('%4f     %15.14f      %15.14f     %15.14f      %15.14f'%(kappa,self.t0[(obs,'0.1')],self.t0_sigmaj[(obs,'0.1')],self.t0[(obs,'0.05')],self.t0_sigmaj[(obs,'0.05')]),file=g)
        g.close()

    # RETURNED:
    # M = cuts(...) ---> M = [time,t2E,t2E_sigmaj]
    def cuts(self,N,Tmin=1.,threshold=.25,t2Emax=0.12,print_cuts=False):
        if not self.find_t0_has_been_called:
            self.find_t0(.05)
        time = self.time/self.t0[('t','0.05')]
        if not self.jkk_sigma_has_been_called:
            self.jkk_sigma()


        # APPLY CUTS AND RESTRICTIONS
        # Restrict to data with T>1
        time = time[time>Tmin]
        t2E = np.take(self.t2E,np.nonzero(time>Tmin)[0])
        t2E_sigmaj = np.take(self.t2E_sigmaj,np.nonzero(self.time>Tmin)[0])

        # Restrict to data with sqrt(8T0/N)<0.25
        Tmax = threshold*N/8.
        time = time[time<Tmax]
        t2E = np.take(t2E,np.nonzero(time<Tmax)[0])
        t2E_sigmaj = np.take(t2E_sigmaj,np.nonzero(time<Tmax)[0])

        # Restrict to data with t2E<0.12
        t2E = t2E[t2E<t2Emax]
        t2E_sigmaj = np.take(t2E_sigmaj,np.nonzero(t2E<t2Emax)[0])
        time = np.take(time,np.nonzero(t2E<t2Emax)[0])

        toreturn = []
        toreturn.append(time)
        toreturn.append(t2E)
        toreturn.append(t2E_sigmaj)
        toreturn = np.array(toreturn)

        if print_cuts:
            plt.axvline(1.)

        return toreturn
        










# ==============================================================================================
# ==============================================================================================
files = [
    #'nfh/n289b0350k5hf100',
    #'nfh/n289b0350k5hf110',
    #'nfh/n289b0350k5hf120',
    #'nfh/n289b0350k5hf130',
    #'nfh/n289b0350k5hf140',
    #'nfh/n289b0350k5hf150',
    #'nfh/n289b0350k5hf155',
    #'nfh/n289b0350k5hf160',
    #'nfh/n289b0350k5hf165',
    #'nfh/n289b0350k5hf170',
    #'nfh/n289b0350k5hf175',
    'nfh/n289b0350k5hf1775',
    'nfh/n289b0350k5hf1800',
    'nfh/n289b0350k5hf1825',
    'nfh/n289b0350k5hf1850',
    'nfh/n289b0350k5hf1875',
    #'nfh/n169b0350k5hf1775',
    #'nfh/n169b0350k5hf1800',
    #'nfh/n169b0350k5hf1825',
    #'nfh/n169b0350k5hf1850',
    #'nfh/n361b0350k7hf1825'
    #'nfh/n169b0347k5hf1850',
    #'nfh/n169b0347k5hf1875',
    #'nfh/n169b0347k5hf1890'
]

for file in files:
    print(30*'-',file,30*'-')
    kappa = float(file[17:])/10000
    N = int(file[5:8])

    flow = TEK_flow(file)
    flow.norm_correction(N,correct_data=True)
    flow.jkk_sigma()
    
    x = flow.time()
    plt.fill_between()





# ==============================================================================================
# ==============================================================================================