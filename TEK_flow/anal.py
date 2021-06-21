import numpy as np
import gvar as gv
import math
import os
import sys

import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 
from operator import itemgetter


class TEK_flow:
    def __init__(self,basement=None,deltat=.03125,usecols=None):
        # Class flags and attributes
        self.find_t0_has_been_called = False
        self.jkk_sigma_has_been_called = False
        self.norm_correction_has_been_called = False
        self.harvest_has_been_called = False

        self.deltat = deltat
        self.t0 = {}
        self.t0_sigmaj = {}

        self.data = []

        # ----------------------------------------------------------------------

        if basement==None: self.basement = 'unknown'
        else:              self.basement = basement
        
        self.deltat = deltat

        if usecols==None: self.usecols = 3
        else            : self.usecols = usecols

    # RETURNED:
    # self.time
    # self.data[conf,tt_index] (matrix of <t^2E> for every configuration)
    # self.t2E                 (self.t2E is <t^2E>)
    def harvest(self,data,outfolder=None,save_to=None):
        self.harvest_has_been_called = True

        # HARVEST DATA ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # Detect if input data is a numpy array or a path
        if isinstance(data,str):
            if os.path.isdir(data): # --------------------------------------------------------------------
                tmaxes = []
                for file in sorted(os.listdir(data)):
                    filename = os.path.join(data,file)

                    time = np.loadtxt(filename,comments='#',usecols=0).tolist()
                    t2e =  np.loadtxt(filename,comments='#',usecols=self.usecols).tolist()

                    if time[0]==self.deltat: t2e.insert(0,0.)
                    # print(time[0],t2e[0])

                    self.data.append(t2e)
                    tmaxes.append(len(t2e))
                    
                # Count how many timesteps and reshape
                self.Tmin = 0
                self.Tmax = min(tmaxes)
                for ii in range(len(self.data)):
                    self.data[ii] = self.data[ii][:self.Tmax]
                
                self.data = np.array(self.data)
                print(self.data.shape)


                print('Data have been collected in self.data[conf,t2E] from',data)
                if not save_to==None:
                    print('I am going to save them in',save_to+'.npy')
                    flag = input('IS THAT OK? (y or n)')
                    if flag.strip()=='y':
                        np.save(save_to,self.data)
                    else: exit()


            elif os.path.isfile(data) and data.endswith('.npy'):  # -----------------------------------------
                self.data = np.load(data)
                print('Data have been collected in self.data[conf,t2E] from',data)
                print('self.data.shape = ',self.data.shape)               

            else: # ------------------------------------------------------------------------------------------
                print('LOCAL error: ',data,'is neither a directory nor a .npy file. Exiting...')
                exit()

        elif isinstance(data,np.ndarray):
            self.data = data
        else:
            print('LOCAL error:',data,'has a non-recognizable format. Exiting...')
            exit()
        
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        self.t2E = []
        

   
                




flow = TEK_flow('n289b0350k5hf1775')
# flow.harvest('/home/pietro/Desktop/DATA_LOCAL/TEK_flow/nfh/n289b0350k5hf1775/',save_to='DATA/'+flow.basement)
flow.harvest('DATA/'+flow.basement+'.npy')
