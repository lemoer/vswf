#!/usr/bin/python3

# Based on http://openmx-square.org/exx/CPC2009/
# Or also: https://dlmf.nist.gov/10.59
# Or also: https://www.int.uni-rostock.de/fileadmin/user_upload/publications/spors/2019/Hahn_et_al_DAGA_Time_Domain_PW.pdf

# For the other part, maybe use FT of 1/x = sign(x) and do a convolution?
# --> see integral 310 here: https://en.wikipedia.org/wiki/Fourier_transform

import numpy as np
from scipy.special import lpmv, spherical_jn, spherical_yn, eval_legendre

l = 1

T_MAX = 3 # seems like, this must be at least a, in order to deliver correct results (Probably aliasing)
NT = 1000 # The higher this is, the better the approximations of the P at the origin become. (Probably aliasing)

a = 3

x = np.linspace(0, 100*np.pi, 1000)
t = np.linspace(-T_MAX, T_MAX, NT) # The larger T_MAX becomes, the worse the error in the bessel function for large x
                                   # becomes. I suppose, the reason is that if T_MAX is increased, dt becomes larger.

t_row = t[np.newaxis, :]
t_col = t[:, np.newaxis]
x_col = x[:, np.newaxis]
window_col = np.ones(t_col.shape)
window_col[np.abs(t/a) > 1] = 0
dt = t[1] - t[0]

J = spherical_jn(l, a*x)
J_by_x = spherical_jn(l, a*x)/(a*x)
P = 1/(2*(1j**l))/a*eval_legendre(l, t_col/a) * window_col
# This is the analytical solution of the integral of P
P_integral_ana = -1/(2*(1j**l))*1j/a*(eval_legendre(l+1, t_col/a) - eval_legendre(l-1, t_col/a))/(2*l+1)*window_col

FT = np.matrix(np.exp(1j*t_row*x_col)) * dt # fourier transform

J_integral_repr = FT*P
J_by_x_integral_repr = FT*P_integral_ana

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.plot(x, J)
ax1.plot(x, np.real(J_integral_repr), '--')
ax1.plot(x, np.imag(J_integral_repr))

ax2.plot(x, J_by_x)
ax2.plot(x, np.real(J_by_x_integral_repr), '--')
ax2.plot(x, np.imag(J_by_x_integral_repr))

plt.figure()
plt.plot(t, np.abs(P))
plt.plot(t, np.abs(P_integral_ana))
plt.show(block=True)