import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import sph_harm

phi = np.linspace(0, 2*np.pi, 100)
theta = np.linspace(0, np.pi, 100)
theta, phi = np.meshgrid(theta, phi)

def coord_s2c(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z

def coord_c2s(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y, x)
    theta = np.arctan2(np.sqrt(x**2+y**2), z)

    return r, theta, phi

def vec_s2c(f_r, f_theta, f_phi, r, theta, phi):
    # from Balanis 2015, p. 955

    cp = np.cos(phi)
    sp = np.sin(phi)
    ct = np.cos(theta)
    st = np.sin(theta)

    f_x = cp * (st*f_r + ct*f_theta) - sp * f_phi
    f_y = sp * (st*f_r + ct*f_theta) + cp * f_phi
    f_z = ct*f_r - st*f_theta

    return f_x, f_y, f_z

def vec_c2s(f_x, f_y, f_z, r, theta, phi):
    # from Balanis 2015, p. 955

    cp = np.cos(phi)
    sp = np.sin(phi)
    ct = np.cos(theta)
    st = np.sin(theta)

    f_r = st * (cp*f_x + sp*f_y) + ct*f_z
    f_theta = ct * (cp*f_x + sp*f_y) - st*f_z
    f_phi = -sp*f_x + cp*f_y

    return f_r, f_theta, f_phi

## Tests for vector transformation...
# F = np.ones_like(phi)
# nul = np.zeros_like(phi)

# F_r, F_theta, F_phi = vec_c2s(F, nul, F, 1, theta, phi)
# F_x, F_y, F_z = vec_s2c(F_r, F_theta, F_phi, 1, theta, phi)

# plt.figure()
# plt.plot(F_x)
# plt.figure()
# plt.plot(F_y)
# plt.figure()
# plt.plot(F_z)


x, y, z = coord_s2c(1, theta, phi)
#r, theta, phi = coord_c2s(x, y, z)
#x, y, z = coord_s2c(1, theta, phi)


m, l = 0, 4

# Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
# Be cationous, sph_harm takes other args! 
fcolors = sph_harm(m, l, phi, theta).real
fmax, fmin = fcolors.max(), fcolors.min()
fcolors = (fcolors - fmin)/(fmax - fmin)

# Set the aspect ratio to 1 so our sphere looks spherical
fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
# Turn off the axis planes
ax.set_axis_off()
plt.show()
