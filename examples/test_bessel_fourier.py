#!/usr/bin/python3

# Based on http://openmx-square.org/exx/CPC2009/
# Or also: https://dlmf.nist.gov/10.59
# Or also: https://www.int.uni-rostock.de/fileadmin/user_upload/publications/spors/2019/Hahn_et_al_DAGA_Time_Domain_PW.pdf

# For the other part, maybe use FT of 1/x = sign(x) and do a convolution?
# --> see integral 310 here: https://en.wikipedia.org/wiki/Fourier_transform

import numpy as np
from scipy.special import lpmv, spherical_jn, spherical_yn, eval_legendre

l = 1

T_MAX = 2
NT = 1000 # The higher this time is, the better the approximations of the P at the origin become. (Probably aliasing)

x = np.linspace(0, 100*np.pi, 1000)
t = np.linspace(-T_MAX, T_MAX, NT) # The larger T_MAX becomes, the worse the error in the bessel function for large x
                                   # becomes. I suppose, the reason is that if T_MAX is increased, dt becomes larger.

t_row = t[np.newaxis, :]
t_col = t[:, np.newaxis]
x_col = x[:, np.newaxis]
window_col = np.ones(t_col.shape)
window_col[np.abs(t) > 1] = 0
dt = t[1] - t[0]

J = spherical_jn(l, x)
J_by_x = spherical_jn(l, x)/x
P = eval_legendre(l, t_col) * window_col
# This is the analytical solution of the integral of P
P_integral_ana = -1j*(eval_legendre(l+1, t_col) - eval_legendre(l-1, t_col))/(2*l+1)*window_col

FT = np.matrix(np.exp(1j*t_row*x_col)) * dt # fourier transform

J_integral_repr = 1/(2*(1j**l))*FT*P
J_by_x_integral_repr = 1/(2*(1j**l))*FT*P_integral_ana

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

plt.show(block=True)

plt.figure()
plt.plot(t, P)
plt.plot(t, np.imag(P_integral_ana))