import numpy as np
import math
import sys
import os



DATA_LOCAL = "/home/pietro/Desktop/DATA_LOCAL/TEK_spec/"
folder=sys.argv[1]

newfolder = os.path.join('DATA/',folder)
if not os.path.isdir(newfolder): os.mkdir(newfolder)


for subdir, dirs, files in os.walk(os.path.join(DATA_LOCAL,folder)):
    for file in files:
        newfolder_op = os.path.join(newfolder,os.path.basename(subdir))
        if not os.path.isdir(newfolder_op): os.mkdir(newfolder_op)
        filename = os.path.join(subdir, file)

        outfile = os.path.join(newfolder_op,file)
        np.save(outfile,np.loadtxt(filename))