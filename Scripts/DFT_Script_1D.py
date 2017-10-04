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
from scipy.fftpack import fft, fftfreq, fftshift

# Calculates the DFT using FFT for a given signal and concerting frequency and amplitude spaces, such that they show correct physical

# example: if time spacing dt = 0.00629 --> sample frequency 159.00 Hz --> Frequency spectrum is from 0 to fa/2 = Niquisit Limit = 79.5 Hz --> Any frequency above that frequency will be "folded" back and is tehrefore aliased


t = np.linspace(0, 2*np.pi, 1000, endpoint = True) # parameters: start point, end point, number of samples in between
f = 10.0 #Frequency of signal in Hz
ft = 50.0
A = 100.0 #Amplitude
At = 50.0
s = A * np.sin(2*np.pi*f*t) + At * np.sin(2 * np.pi * ft * t) #signal

Y = np.fft.fft(s) # --> This only produces a strange signal, which is mirrored

#plt.plot(t, s)
#plt.plot(Y)

# define N has HALF of the length of the output signal
N = len(Y)/2+1

# x and y axis values still need to be converted into real physics values!

# first, get spacing in time domain, assuming an equal spacing (given by linspace) --> If another signal is used, make sure it is equally spaces in the time domain
dt = t[1] - t[0]
fa = 1.0 / dt # scan Frequency

# Remember: Frequencys can only be resolved up until Niquisit limit fn.
# Therefore: Create a new (frequency domain) array, from 0 to half of the FFT signal, going up until the Niquisit Frequency

# y-axis amplitude still does not fit --> Because power of signal in time and frequency domain have to be equal, but here we just used the left half of the signal (N), now we need to multiply the amplitude with a factor of 2!
# Also, in most implementations, the output Y of the FFT is normalized with the number of samples --> We have to multiply by N to get the real physical value --> Factor of 2/N

X = np.linspace(0, fa/2, N, endpoint = True)

# Note, that there is still the "leakage-effect" going on, because we don't have a signal ending at amplitude 0.
# "Windowing" has to be introduced, in order to get rid of this effect!

hann = np.hanning(len(s))
Yhann = np.fft.fft(hann * s)

#plt.plot(X, 2.0 * np.abs(Y[:N]) / N)
plt.plot(X, 2.0 * np.abs(Yhann[:N]) / N)
plt.xlabel('Frequency ($Hz$)')
plt.show()
