import numpy as np
import glob
import matplotlib.pyplot as plt
#from matplotlib import rc
import matplotlib.patches as patches
import fnmatch
from matplotlib import rc, rcParams
from numpy import sin
import csv
import os

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

#change font style and font size
#rc('text', usetex=True)
#rc('font', family='serif', serif='Computer Modern Roman', size=8)
#rc('legend', fontsize=10)

rect = np.loadtxt('rectangleInfo.txt')
simulationInfo = np.loadtxt('simulationInfo.txt')


# simulation info from file
tN = simulationInfo[0]
tPreCalc = simulationInfo[1]
dx = simulationInfo[2]
numberOfBoxesInX = np.float(simulationInfo[3])
numberOfGridPointsPerBoxInX = np.float(simulationInfo[4])
numberOfParticles = np.int(simulationInfo[5])
dt = 0.5 * dx

# calculated values
numberOfSimulationSteps = np.int(tN / dt)
numberOfPrecalculationSteps = np.int(tPreCalc / dt)
simulationArea = np.float(numberOfBoxesInX * numberOfGridPointsPerBoxInX * dx)
sizeOfOneBox = np.float(numberOfGridPointsPerBoxInX * dx)

X = np.zeros((numberOfParticles,1))
Y = np.zeros((numberOfParticles,1))
x = []
y = []

lineWidthRectangles = 3
lineWidthParticles = 1


#generate multiple plots "in range (0,2) loops 0 and 1, i.e. 2= #Loops including 0
#here: plot electric field values on grid and particle movement.
for i in range(numberOfPrecalculationSteps, numberOfPrecalculationSteps + numberOfSimulationSteps): #Extent Loop to calculate Loop Boundary dynamically
    #openfigure
    fig = plt.figure()
    # get current axis
    ax = plt.gca()
    # draw rects, indicating near field of each area
    for p in range(numberOfParticles):
        #ax1 = fig.add_subplot(111, aspect = 'equal')
        ax.add_patch(
        patches.Rectangle((rect[i * numberOfParticles + p, 0] - sizeOfOneBox, rect[i * numberOfParticles + p, 1] - sizeOfOneBox), 3 * sizeOfOneBox, 3 * sizeOfOneBox, fill = False, edgecolor = "black", linewidth = lineWidthRectangles
                          )
                  )

    for p in range(numberOfParticles):
        data = np.genfromtxt('Particle'+ str(p) +'_timeStep' + str(i) + '.txt')
        x.append(data[0,1])
        y.append(data[0,2])

    X = np.c_[X,x]
    Y = np.c_[Y,y]

    x = []
    y = []

    for p in range(numberOfParticles):
        plt.plot(X[p,1:i+2], Y[p,1:i+2], color = 'k', linewidth = lineWidthParticles)

    plotindex = str(i - numberOfPrecalculationSteps)

    #read in data from electricField File and save it into array
    field = np.genfromtxt('electricFieldAtTime' + str(i) + '.txt')

    #plot electric field values on grid
    plt.imshow(field, aspect='auto', origin='lower', extent=(0, simulationArea, 0, simulationArea), vmin=0, vmax=0.04) #extent: rescale xy Grid: Number of columns and rows will be the initial scale, but 160 Grid Points in each direction with resolution 0.2, does mean an effective length of 0.2*160=32!
#================================================================================================================================
    #Plot Settings

    #Show Colorbar
    plt.colorbar()

    #Set x- and y lim
    plt.xlim([sizeOfOneBox, simulationArea - sizeOfOneBox])
    plt.ylim([sizeOfOneBox, simulationArea - sizeOfOneBox])

    filename = 'electricFieldAndParticleAtTime' + plotindex
    plt.xlabel(r'X $\left[\mathregular{\frac{e^2}{m_e c^2}} = 2.82\cdot 10^{-15} m\right]$')
    plt.ylabel(r'Y $\left[\mathregular{\frac{e^2}{m_e c^2}} = 2.82\cdot 10^{-15} m\right]$')
    plt.xticks(np.arange(sizeOfOneBox, simulationArea - sizeOfOneBox, sizeOfOneBox))#size of one box
    plt.yticks(np.arange(sizeOfOneBox, simulationArea - sizeOfOneBox, sizeOfOneBox))

    plt.grid(linestyle= "-",color = "r")

    for index, labels in enumerate(ax.xaxis.get_ticklabels()):
        if index % 2 != 0:
            labels.set_visible(False)

    for index, labels in enumerate(ax.yaxis.get_ticklabels()):
        if index % 2 != 0:
            labels.set_visible(False)


    fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight', dpi = 300)
    #close figure
    plt.close(fig)
