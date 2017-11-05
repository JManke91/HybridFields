import numpy as np
import numpy.fft as fft
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import math
import warnings
from numpy import *

warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

simulationInfo = np.loadtxt('simulationInfo.txt') #E_field_z100.txt
im = np.loadtxt('E_field_z88.txt')
dx = simulationInfo[2]
numberOfBoxesInX = np.float(simulationInfo[3])
numberOfGridPointsPerBoxInX = np.float(simulationInfo[4])
numberOfParticles = np.int(simulationInfo[5])
numberOfGridPointsInX = np.int(simulationInfo[6])

#calculated values
simulationArea = np.float(numberOfBoxesInX * numberOfGridPointsPerBoxInX * dx)
sizeOfOneBox = np.float(numberOfGridPointsPerBoxInX * dx)


def makeSpectrum(E, dx, dy, upsample = 10):

    zeropadded = np.array(E.shape) * upsample
    F = fft.fftshift(fft.fft2(E, zeropadded)) / E.size
    xf = fft.fftshift(fft.fftfreq(zeropadded[1], d = dx)) * 2 * math.pi # 2pi is needed because input signal is real space and want to show k space instead of 1/lambda
    yf = fft.fftshift(fft.fftfreq(zeropadded[0], d = dy)) * 2 * math.pi
    return (F, xf, yf)


def extents(f):
    "Convert a vector into the 2-element extents vector imshow needs"
    delta = f[1] - f[0]
    return [f[0] - delta / 2, f[-1] + delta / 2]


def plotSignalAndSpectrum(F, xf, yf, signal, plotNumber):
    "Plot a spectrum array and vectors of x and y frequency spacings"
    #fig.tight_layout()
    F = np.array(F)
    F = F + 1

    #signalMin = min(signal)
    #signalMax = max(signal)

    fig.set_tight_layout(True)
    half = len(F) / 2

    # plot signal image
    im1 = axarr[0,].imshow(signal, origin='lower', cmap = 'jet', extent = (0,simulationArea,0,simulationArea)
                #vmin = signalMin,
                #vmax = signalMax
                )
    axarr[0,].set_title('Electric field E', y = 1.05)
    cbar1 = fig.colorbar(im1, ax = axarr[0,])

    # plot frequency space image
    #logF = np.log(abs(F))
    #myArray = [0, 1, 2, 3]
    #print(myAddedArray)

    maxValue = np.max(abs(F))
    minValue = np.min(abs(F))
    im2 = axarr[1,].imshow(abs(F),
               aspect="equal",
               interpolation="none",
               origin="lower",
               extent = extents(xf) + extents(yf),
               norm = LogNorm(vmin = minValue, vmax = maxValue) # Because many values are around 0, log() produces negative values --> maybe shift those values to positive in the future, but doesn't make a difference
               )

    cbar2 = fig.colorbar(im2, ax = axarr[1,])
    cbar2.set_label('Normalized amplitude', rotation=270, labelpad = 17) # label colorbar

    axarr[1,].set_xlabel(r'$\mathregular{f_x}$', labelpad = 5) #labelpad = 25, moves xlabel further away from ticks
    axarr[1,].set_ylabel(r'$\mathregular{f_y}$', labelpad = 5)
    axarr[1,].set_title(r'$\mathregular{FFT\left\{E\right\}}$')
    # FFT for real signals is redundant, therefore only show positive frequencies
    axarr[1,].set_xlim([0, max(extents(xf))])
    axarr[1,].set_ylim([0, max(extents(yf))])

    # do some annotation
    an1 = axarr[1,].annotate("", xy=(max(extents(xf)), 0.0), xycoords="data", # position of arrow
                  va="center", ha="center")


    an2 = axarr[1,].annotate("Nyquisit x", xy = (0., 0.), xycoords = an1,
                  xytext = (0.0, -0.25), textcoords=(an1, "axes fraction"),
                  va = "bottom", ha = "center",
                  bbox=dict(boxstyle="square", fc = "w"), # also possible for boxstyle: round, round4
                  arrowprops=dict(arrowstyle="->"))


    filename = "FFT_NoPrecalc_t=88_z"
    #plt.show()
    fig.savefig("png/" + "{}.png".format(filename), bbox_inches = 'tight', dpi = 300)
    plt.close(fig)


if __name__ == '__main__':
    # In seconds
    #x = np.linspace(0, 1, 200); y = x
    x = np.linspace(0,simulationArea, numberOfGridPointsInX); y = x
    xx, yy = np.meshgrid(x,y)

    # Sinusoid frequency, in Hz
    x0 = 3.0
    y0 = 8.0

    # complex signal
    #im = np.exp(2j * np.pi * (y[:, np.newaxis] * y0 + x[np.newaxis, :] * x0))
    #imSignal = np.exp(2j * np.pi * (yy * y0 + xx * x0))
    #imSignal = np.real(imSignal)

    #real signal
    #im = 2 * np.sin(2 * np.pi * (y[:, np.newaxis] * y0)) * 3 * np.sin(2 * np.pi * x[np.newaxis, :] * x0) + 10 * np.sin(2 * np.pi * (y[:, np.newaxis] * 10)) * 20 * np.sin(2 * np.pi * x[np.newaxis, :] * 15) # This is used for FFT
    imSignal = np.sin(2 * np.pi * yy * y0) * np.sin(2 * np.pi * xx * x0) + np.sin(2 * np.pi * yy * 10) * np.sin(2 * np.pi * xx * 15) # This is just used to plot the original signal

    #maxSignal = np.max(im)
    #myTest = [3,4,0,10]
    #array = np.array(myTest)
    #myLog = np.log(array + 1)



    # Generate spectrum and plot
    spectrum, xf, yf = makeSpectrum(im, x[1] - x[0], y[1] - y[0])
    fig, axarr = plt.subplots(2, 1) # create subplot array of 2 x 2
    plotSignalAndSpectrum(spectrum, xf, yf, im, 0)
    #plotSpectrum(im, x, y)
