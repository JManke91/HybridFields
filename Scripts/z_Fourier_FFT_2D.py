from scipy.fftpack import fft, fftfreq, fftshift
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

# try with real data

#values from simulation file
rect = np.loadtxt('rectangleInfo.txt')
simulationInfo = np.loadtxt('simulationInfo.txt')
tN = simulationInfo[0]
tPreCalc = simulationInfo[1]
dx = simulationInfo[2]
numberOfBoxesInX = np.float(simulationInfo[3])
numberOfGridPointsPerBoxInX = np.float(simulationInfo[4])
numberOfParticles = np.int(simulationInfo[5])
numberOfGridPointsInX = np.int(simulationInfo[6])
dt = 0.5 * dx

#calculated values
numberOfSimulationSteps = np.int(tN / dt)
numberOfPrecalculationSteps = np.int(tPreCalc / dt)
simulationArea = np.float(numberOfBoxesInX * numberOfGridPointsPerBoxInX * dx)
sizeOfOneBox = np.float(numberOfGridPointsPerBoxInX * dx)

# calculate parameters for linspace

# get the last simulated file


fq = 2.0; N = 100.0 # fq = frequency of signal
#x = np.linspace(0, 1, N, endpoint = False); y = x # timestep = 8 / (N-1)
x = np.linspace(0, simulationArea, numberOfGridPointsInX, endpoint = False); y = x # timestep = 8 / (N-1)
xx, yy = np.meshgrid(x, y)
# fnc = np.sin(2 * np.pi * fq * xx) + np.sin(2 * np.pi * fq * yy)

Ex = np.genfromtxt('E_field_x100.txt') # change time to last time step
Ey = np.genfromtxt('E_field_y100.txt')
Ez = np.genfromtxt('E_field_z100.txt')

FTEx = np.fft.fft2(Ex)

#ft = np.fft.fft2(fnc)
#ft = np.fft.fftshift(ft)

dx = x[1] - x[0]
sampleFrequency = 1 / dx
nyquisitFrequency = sampleFrequency / 2.0
#print(nyquisitFrequency)

#freq_x = np.fft.fftfreq(ft.shape[0], d = dx) # first parameter: window length in wx direction, second: sample spacing --> returns: DFT sample frequencys

# signal = np.array ([...(8 Elements)])
#fourier = np.fft.fft(signal); n = signal.size; timestep = dt = 0.1
# usual input: freq = np.fft.fftfreq(n, d = timestep) --> outputs 8 elements array([0, 1.25, 2.5, 3.75, -5, -3.75, -2.5, -1.25]) (unordered)
# apply fftshift to order those frequencys

#freq_y = np.fft.fftfreq(ft.shape[1], d = dx)

#freq_x = np.fft.fftshift(freq_x)
#freq_y = np.fft.fftshift(freq_y)

half = len(FTEx) / 2  # for even number of bins

#plt.imshow(
#    2 * abs(ft[: half,: half]) / half, # fix amplitude!!!
#    aspect = 'auto',
#    extent = (freq_x.min(), freq_x.max(), freq_y.min(), freq_y.max()),
#    #extent = (0, nyquisitFrequency, 0, nyquisitFrequency),
#    origin = 'lower',
#    interpolation = 'nearest',
    #cmap = 'viridis' # another color map
#)

fig, axarr = plt.subplots(2, 3) # create subplot array of 2 x 2
fig.suptitle('Title', fontsize = 14)

#fig = plt.figure()

#plt.subplot(231)

im1 = axarr[0, 0].imshow(Ex, origin='lower', cmap = 'jet', extent=(0, simulationArea, 0, simulationArea))
#plt.colorbar(format='%.0e')
axarr[0, 0].set_xlabel('X', fontsize = 14)
axarr[0, 0].set_ylabel('Y', fontsize = 14)
axarr[0, 0].set_title('$E_x$', fontsize = 14)
fig.colorbar(im1, ax = axarr[0, 0])

im2 = axarr[1, 0].matshow(
    2 * abs(FTEx[:half, :half]) / half, # fix amplitude issue
    aspect = 'equal',
    origin = 'lower',
    interpolation = 'nearest'
)
axarr[1, 0].set_xlabel('Frequency x')
axarr[1, 0].set_ylabel('Frequency y')
axarr[1, 0].xaxis.set_ticks_position('bottom')
#axarr[1, 0].set_colorbar()
fig.colorbar(im2, ax = axarr[1, 0])


plt.grid()
#plt.xlim([0, nyquisitFrequency])
#plt.ylim([0,nyquisitFrequency])
#plt.colorbar()
plt.show()
