import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import fnmatch
from matplotlib import rc, rcParams
from numpy import sin
import csv
import os

#change font style and font size
#rc('text', usetex=True)
#rc('font', family='serif', serif='Computer Modern Roman', size=8)
#rc('legend', fontsize=10)

field = np.genfromtxt('electricFieldAtTime0' + '.txt')

#plt.ion()


#generate multiple plots "in range (0,2) loops 0 and 1, i.e. 2= #Loops including 0
for i in range(0,200):
    #openfigure
    fig = plt.figure()
    field = np.genfromtxt('electricFieldAtTime' + str(i) + '.txt')
    #y = np.sin(x)
    plt.imshow(field, aspect='auto', origin='lower', vmin=0, vmax=0.1)
    #Show Colorbar
    plt.colorbar()
    #_= raw_input("Press [enter] to continue")
    #filename = 'thisIsTheFileName'
    filename = 'thisIsTheFileName' + str(i)
    #Plot Settings
    plt.xlabel("X")
    plt.ylabel("Y")
    fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight')
    #close figure
    plt.close(fig)




