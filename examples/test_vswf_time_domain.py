import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

# add swf path to pythons search path
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from swf import *

Nt = 20001  # number of timesteps
dt = 1e-12 # timestep in seconds

time_shift = 0.7e-9

t = np.linspace(0, (Nt-1)*dt, Nt)

# TODO: this window should be calculated by the max radius, where results should be obtained
t_shifted = t-time_shift

theta = np.array([np.pi/2])
phi = np.array([0])

c = 1    # other stuff is not implemented yet
tau = 2  # other stuff is not implemented yet
l = 3
m = 1

### SHOW TIME DOMAIN SIGNAL

r = np.array([0.1])

f_r, f_theta, f_phi = vswf_time_domain(t_shifted, c, tau, l, m, r, theta, phi)

if tau == 1:
    f_observed = f_theta
    f_observed_str = '\\theta'
else:
    f_observed = f_phi
    f_observed_str = '\\phi'


plt.figure()
plt.plot(t, np.abs(f_observed), label="r = 0.1")

r = np.array([0.2])
f_r, f_theta, f_phi = vswf_time_domain(t_shifted, c, tau, l, m, r, theta, phi)

if tau == 1:
    f_observed = f_theta
else:
    f_observed = f_phi

plt.plot(t, np.abs(f_observed), label="r = 0.2")
plt.xlabel("t / s")
plt.title('Time Domain: $F_' + f_observed_str + "(t, r, \pi/2, 0)$")
plt.legend()

### COMPARE FOURIER TRANSFORM

fft_axis = 0

# Some FFT preparation
f = np.fft.fftfreq(Nt, dt)
f_idx = f > 0    # only positive frequencies
f_pos = f[f_idx]

# FFT of time domain signals
f_r_freq = np.fft.fft(f_r, axis=fft_axis)*dt
f_theta_freq = np.fft.fft(f_theta, axis=fft_axis)*dt
f_phi_freq = np.fft.fft(f_phi, axis=fft_axis)*dt

# VSWF in frequency domain for comparison
c0 = 3e8
k = 2*np.pi*f_pos/c0
f_r_freq_ref, f_theta_freq_ref, f_phi_freq_ref = vswf(k, c, tau, l, m, r, theta, phi)

f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

# norm_it = lambda x: x/np.max(x) if max(x) > 0 else x
norm_it = lambda x: x


ax1.plot(f_pos, norm_it(np.abs(f_r_freq[f_idx])), label="FFT of TD")
ax1.plot(f_pos, norm_it(np.abs(f_r_freq_ref)), '--', label="Ref. in FD")
ax1.set_xlim([-0.5e9, 10e9])
ax1.legend()
ax1.set_ylabel('$|F_r|$')
ax1.set_title("Frequency Domain")

ax2.plot(f_pos, norm_it(np.abs(f_theta_freq[f_idx])), label="FFT of TD")
ax2.plot(f_pos, norm_it(np.abs(f_theta_freq_ref)), '--', label="Ref. in FD")
ax2.set_xlim([-0.5e9, 10e9])
ax2.legend()
ax2.set_ylabel('$|F_\\theta|$')

ax3.plot(f_pos, norm_it(np.abs(f_phi_freq[f_idx])), label="FFT of TD")
ax3.plot(f_pos, norm_it(np.abs(f_phi_freq_ref)), '--', label="Ref. in FD")
ax3.set_xlim([-0.5e9, 10e9])
ax3.legend()
ax3.set_ylabel('$|F_\\phi|$')
ax3.set_xlabel('f / Hz')

plt.show(block=True)
