import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from numpy import sin
import csv

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

# first particle
x = []
y = []

# second particle
x2 = []
#  time
t = []

lineWidthRectangles = 3
lineWidthParticles = 1

plotEveryNthIndex = 10

for i in range(numberOfPrecalculationSteps, numberOfPrecalculationSteps + (numberOfSimulationSteps - numberOfPrecalculationSteps) / plotEveryNthIndex): #Extent Loop to calculate Loop Boundary dynamically
    # plot information
    print "plotting step", i, "of", numberOfPrecalculationSteps + (numberOfSimulationSteps - numberOfPrecalculationSteps) / plotEveryNthIndex

    if i == numberOfPrecalculationSteps:
        y.append(float(0))
    if i > numberOfPrecalculationSteps:
        newElement = float(y[-1]) + (plotEveryNthIndex * float(dt)) # only every 5-th element
        y.append(newElement)

    fig = plt.figure()
    ax = plt.gca()

    getParticleTrajFromIndex = numberOfPrecalculationSteps + (plotEveryNthIndex * (i - numberOfPrecalculationSteps))
    getParticleTrajFromIndex = str(getParticleTrajFromIndex)

    particleData = np.genfromtxt('Particle'+ str(0) +'_timeStep' + str(getParticleTrajFromIndex) + '.txt')
    particleData2 = np.genfromtxt('Particle'+ str(1) +'_timeStep' + str(getParticleTrajFromIndex) + '.txt')

    x.append(particleData[0, 1])
    x2.append(particleData2[0, 1])

    plt.plot(x, y, 'ro', markersize = 5)
    plt.plot(x2, y, 'bo', markersize = 5)

    plotindex = str(i - numberOfPrecalculationSteps)

    #Set x- and y lim
    plt.xlim([sizeOfOneBox, simulationArea - sizeOfOneBox])
    plt.ylim([0, tN - tPreCalc])

    filename = 'Minkowski-Diagram_t=' + plotindex
    plt.xlabel('x label')
    plt.ylabel('y label')

    plt.xticks(np.arange(sizeOfOneBox, simulationArea - sizeOfOneBox, sizeOfOneBox))#size of one box
    plt.yticks(np.arange(0, tN - tPreCalc, dt))
    ax.yaxis.set_ticks_position('left')

    for index, labels in enumerate(ax.xaxis.get_ticklabels()):
        if index % 2 != 0:
            labels.set_visible(False)

    for index, labels in enumerate(ax.yaxis.get_ticklabels()):
        if index % 10 != 0:
            labels.set_visible(False)


    fig.savefig("png/" + "{}.png".format(filename), bbox_inches='tight', dpi = 100)
    plt.close(fig)
