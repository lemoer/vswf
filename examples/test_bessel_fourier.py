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

a = 0.2

x = np.linspace(0, 100*np.pi, 1000)
t = np.linspace(-T_MAX, T_MAX, NT) # The larger T_MAX becomes, the worse the error in the bessel function for large x
                                   # becomes. I suppose, the reason is that if T_MAX is increased, dt becomes larger.

t_row = t[np.newaxis, :]
t_col = t[:, np.newaxis]
x_col = x[:, np.newaxis]
boxcar = lambda x: np.heaviside(x+1, 0.5)-np.heaviside(x-1, 0.5)
dt = t[1] - t[0]
dx = x[1] - x[0]

J = spherical_jn(l, a*x)
J_by_x = spherical_jn(l, a*x)/(a*x)
J_deriv = spherical_jn(l, a*x, derivative=True)
P = 1/(2*(1j**l))/a*eval_legendre(l, t_col/a) * boxcar(t_col/a)
# This is the analytical solution of the integral of P
P_integral_ana = -1/(2*(1j**l))*1j/a*(eval_legendre(l+1, t_col/a) - eval_legendre(l-1, t_col/a))/(2*l+1)*boxcar(t_col/a)
P_times_t = 1j*(t_col/a)*1/(2*(1j**l))/a*eval_legendre(l, t_col/a) * boxcar(t_col/a)

FT = np.matrix(np.exp(1j*t_row*x_col)) * dt # fourier transform
IFT = 1/np.pi*np.matrix(np.exp(1j*x_col.T*t_row.T)) * dx # fourier transform

J_integral_repr = FT*P
J_by_x_integral_repr = FT*P_integral_ana
J_deriv_integral_repr = FT*P_times_t

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

ax1.plot(x, J)
ax1.plot(x, np.real(J_integral_repr), '--')
ax1.plot(x, np.imag(J_integral_repr))

ax2.plot(x, J_by_x)
ax2.plot(x, np.real(J_by_x_integral_repr), '--')
ax2.plot(x, np.imag(J_by_x_integral_repr))

ax3.plot(x, J_deriv)
ax3.plot(x, np.real(J_deriv_integral_repr), '--')
ax3.plot(x, np.imag(J_deriv_integral_repr))


plt.figure()

if np.mod(l, 2) == 0:
    part = lambda x: np.real(x)
    part_str = "Re "
    part_other = lambda x: np.imag(x)
    part_other_str = "Im "
else:
    part = lambda x: np.imag(x)
    part_str = "Im "
    part_other = lambda x: np.real(x)
    part_other_str = "Re "

plt.plot(t, part(P), label=part_str + "P")
plt.plot(t, part(IFT*J[:,np.newaxis]), label="FT J")
plt.plot(t, part_other(P_integral_ana), label=part_other_str + "P_{integral}|")
plt.plot(t, part_other(P_times_t), label=part_other_str + "t*P")
plt.legend()
plt.show(block=True)