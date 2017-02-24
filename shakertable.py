import alsaaudio, time, audioop
import numpy as np
from math import *
import matplotlib.pyplot as plt
import scipy.signal as signal
import struct

chirpsig = lambda t, f0, f1, t1, phi_0: cos(2*pi*(f0*((f1/f0)**(1./t1))**t)*t+phi_0)

chirp_t = 5
f0 = 40.
f1 = 500.
fs = 8000
n_samples = int(round(chirp_t*fs))
times = np.linspace(0,chirp_t,n_samples)
y = [(chirpsig(t,f0,f1,chirp_t,0.)*32767,chirpsig(t,f0,f1,chirp_t,pi)*32767) for t in times]

output = alsaaudio.PCM()
output.setrate(fs)

buf = ''
for i in range(n_samples):
    output.write(struct.pack('<hh', *y[i]))
