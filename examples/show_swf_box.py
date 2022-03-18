import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

# add swf path to pythons search path
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from swf import *

f = 1e9
c0 = 3e8
k = 2*np.pi*f/c0
c, tau, l, m = 4, 2, 1, 0

fig = plt.figure(figsize =(14, 9))
ax = plt.axes(projection ='3d')

faces = []
faces2 = []
max_E_x = None

for dn in range(3):
    for sign in [-1, 1]:
        pos, dA = face(dn, [100, 100, 100], 0.1, sign)
        x, y, z = pos[0,:,:], pos[1,:,:], pos[2,:,:]

        em = EMField(k, pos, dA)
        em.override_by_vswf(c, tau, l, m)
        faces += [em]

        em = EMField(k, pos, dA)
        em.override_by_vswf(c, tau, l, m+1)
        faces2 += [em]

        if not max_E_x or np.abs(em.E[0,:,:]).max().max() > np.abs(max_E_x):
            max_E_x = np.abs(em.E[0,:,:]).max().max()

res_ab = 0
res_aa = 0
res_bb = 0

for i in range(6):
    a = faces[i]
    b = faces2[i]

    res_ab += np.sum((a.H.conj()*b.E - b.H.conj() * a.E)*a.dA)
    res_aa += np.sum((a.H.conj()*a.E - a.H.conj() * a.E)*a.dA)
    res_bb += np.sum((a.H.conj()*a.E - a.H.conj() * a.E)*a.dA)

print('res_aa:', res_aa)
print('res_ab:', res_ab)
print('res_bb:', res_bb)

norm = matplotlib.colors.Normalize(vmin=0, vmax=max_E_x)

for i in range(6):
    x = faces[i].pos[0, :, :]
    y = faces[i].pos[1, :, :]
    z = faces[i].pos[2, :, :]
    ax.plot_surface(x, y, z, facecolors=plt.cm.jet(norm(np.abs(faces2[i].E[0,:,:]))), shade=False)

# x, y, z = coord_s2c(1, theta, phi)
#r, theta, phi = coord_c2s(x, y, z)
#x, y, z = coord_s2c(1, theta, phi)

# dA = calc_dA(theta, phi)

# test_lpmv_diff_stuff()
# check_correlation_for_my_sph_harm(theta, phi, dA)
#test_vec_sph_harm()

# m, l = 0, 4

# # Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
# # Be cationous, sph_harm takes other args! 
# #fcolors = sph_harm(m, l, phi, theta).real
# fcolors = my_sph_harm(l, m, theta, phi).real
# fmax, fmin = fcolors.max(), fcolors.min()
# fcolors = (fcolors - fmin)/(fmax - fmin)

# # Set the aspect ratio to 1 so our sphere looks spherical
# fig = plt.figure(figsize=plt.figaspect(1.))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
# # Turn off the axis planes
# ax.set_axis_off()
# plt.show()

plt.show(block=True)

print("end")
