import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import sph_harm, factorial, lpmv, spherical_jn, spherical_yn

phi = np.linspace(0, 2*np.pi, 100)
theta = np.linspace(0, np.pi, 50)
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

def my_sph_harm(l, m, theta, phi):
    # This is a copy of sph_harm() with different argument order based
    # on Kristenson, p. 95. It has been validated for some examples,
    # that we get the same results as sph_harm().

    # we leave out the (-1)**m here in c_lm, because it is already
    # part of the legendre definition in lpmv(). 
    c_lm = np.sqrt((2*l+1)/(4*np.pi)*factorial(l-m)/factorial(l+m))

    exp_part = np.exp(1j*m*phi)
    # NOTE: There is also lpmn(), which could be more efficient, if
    # multiple l and m should be calculated... Optimize here in
    # future, if you need better performance.
    leg_part = lpmv(m, l, np.cos(theta))

    return c_lm*exp_part*leg_part

# Check my_sph_harm() as 1d fun...
# l, m = 4, 2
# # phi = np.pi/4
# # theta = np.linspace(0, np.pi, 100)
# phi = np.linspace(0, 2*np.pi, 100)
# theta = np.pi/4

# ref = sph_harm(m, l, phi, theta).imag
# my = my_sph_harm(l, m, theta, phi).imag

# plt.figure()
# plt.plot(phi*180/np.pi, ref)
# p = plt.plot(phi*180/np.pi, my, linestyle='--')

def sph_harm_diff_phi(l, m, theta, phi):
    # differentation of my_sph_harm() w.r.t. phi (with conventions
    # that phi is the azimuth.)
    return 1j*m*my_sph_harm(l, m, theta, phi)

def lpmv_diff_times_minus_sqrt_1_minus_x_squared(m, l, x):
    # This function calculates the derivative of the legendre function lpmv()
    # multiplied with a factor of -sqrt(1-x²), as in Kristenson 2014 p. 84
    # (the last formula on the page, with an additional minus on both sides
    # of the equation). This makes sense, because we actually want to derive 
    # lpmv(m, l, np.cos(theta)) and the factor -sqrt(1-x²) exactly corresponds
    # to the inner derivative of the np.cos(theta). So this function here 
    # can be used directly to calculate the derivative. See test_lpmv_diff_stuff().
    return -(1/2)*((l-m+1)*(l+m)*lpmv(m-1, l, x)-lpmv(m+1, l, x)) # grundmann

def test_lpmv_diff_stuff():
    plt.figure()
    theta = np.linspace(0, np.pi, 50)
    dtheta = theta[1] - theta[0]

    # theta values in between the original values ("half-step")
    theta_ = np.convolve(theta, [0.5, 0.5], 'valid')
    
    for l in range(3, 4):
        m = 1
        leg_part = lpmv(m, l, np.cos(theta))
        leg_diff_analytic = lpmv_diff_times_minus_sqrt_1_minus_x_squared(m, l, np.cos(theta))
        leg_diff_numeric = np.diff(leg_part)/dtheta

        plt.plot(theta, leg_part)
        # The following two graphs should be equal
        plt.plot(theta, leg_diff_analytic)
        plt.plot(theta_, leg_diff_numeric, ':')

def sph_harm_diff_theta(l, m, theta, phi):
    # we leave out the (-1)**m here in c_lm, because it is already
    # part of the legendre definition in lpmv(). 
    c_lm = np.sqrt((2*l+1)/(4*np.pi)*factorial(l-m)/factorial(l+m))

    exp_part = np.exp(1j*m*phi)
    # Inner derivation of the lpmv() argument:
    # x = np.cos(theta)
    # dx/dtheta = -sin(theta) = -sqrt(1-cos(theta)²) = -sqrt(1-x²)
    leg_part = lpmv_diff_times_minus_sqrt_1_minus_x_squared(m, l, np.cos(theta))

    return c_lm*exp_part*leg_part

def vec_sph_harm(tau, l, m, theta, phi):
    norm_fact = np.sqrt(l*(l+1))
    if tau == 1 or tau == 2:
        sin_theta = np.sin(theta)
        idx_pole = np.where(theta == 0)[0]
        sin_theta[idx_pole] = 1 # TODO: this is a hacky way of pole-zero
                                # cancellation

    if tau == 1:
        f_theta = 1/(norm_fact*sin_theta)*sph_harm_diff_phi(l, m, theta, phi)
        f_phi = -1/norm_fact*sph_harm_diff_theta(l, m, theta, phi)
        f_r = np.array([0])
    elif tau == 2:
        f_theta = 1/norm_fact*sph_harm_diff_theta(l, m, theta, phi)
        f_phi = 1/(norm_fact*sin_theta)*sph_harm_diff_phi(l, m, theta, phi)
        f_r = np.array([0])
    elif tau == 3:
        f_theta = np.array([0])
        f_phi = np.array([0])
        f_r = my_sph_harm(l, m, theta, phi)
    else:
        raise Exception('NIY')

    return f_r, f_theta, f_phi

def fix_dim(A, axis):
    # Repeats the last column or row of a matrix once, so the dimension
    # in that axis is increased by 1.
    if axis == 0:
        return np.append(A, np.reshape(A[-1,:], (1, A.shape[1])), 0)
    elif axis == 1:
        return np.append(A, np.reshape(A[:,-1], (A.shape[0], 1)), 1)
    else:
        raise Exception('fix_dim() is not implemented for axes > 1!')

def calc_dA(theta, phi):
    assert theta[0,0] == theta[1,0], 'Grid direction is not as expected!'
    dtheta = fix_dim(np.diff(theta, axis=1), axis=1)
    dphi = fix_dim(np.diff(phi, axis=0), axis=0)
    return dtheta * dphi * np.sin(theta)

def check_correlation_for_my_sph_harm(theta, phi, dA):
    m = 0
    L = 10
    C = np.zeros((L, L), dtype=complex)

    # make arrays 1D, so we can use the dot() properly.
    theta = np.ravel(theta)
    phi = np.ravel(phi)
    dA = np.ravel(dA)

    # do the work
    for l in range(1, L+1):
        for l2 in range(1, L+1):
            f1 = my_sph_harm(l, m, theta, phi)
            f2 = my_sph_harm(l2, m, theta, phi)
            # f1 = sph_harm(m, l, phi, theta)
            # f2 = sph_harm(m, l2, phi, theta)
            
            C[l-1, l2-1] = np.dot(f2.conj(), dA*f1)

    plt.figure()
    plt.imshow(np.abs(C))
    plt.colorbar()

    return C, f1, f2

def test_correlation_for_vec_sph_harm(theta, phi, dA, tau1, tau2, m1, m2):
    L = 10
    C = np.zeros((L, L), dtype=complex)

    # make arrays 1D, so we can use the dot() properly.
    theta = np.ravel(theta)
    phi = np.ravel(phi)
    dA = np.ravel(dA)

    # do the work
    for l in range(1, L+1):
        for l2 in range(1, L+1):
            f1_r, f1_th, f1_ph = vec_sph_harm(tau1, l, m1, theta, phi)
            f2_r, f2_th, f2_ph = vec_sph_harm(tau2, l2, m2, theta, phi)
            
            C[l-1, l2-1] = np.sum(f2_r.conj()*dA*f1_r) + np.sum(f2_th.conj()*dA*f1_th) + np.sum(f2_ph.conj()*dA*f1_ph)

    plt.figure()
    plt.imshow(np.abs(C))
    plt.colorbar()
    plt.title(f"$\\tau_1 = {tau1}, \\tau_2 = {tau2}, m_1 = {m1}, m_2 = {m2}$")
    # plt.clim([0, 1])

    return C

def test_vec_sph_harm():
    test_correlation_for_vec_sph_harm(theta, phi, dA, 1, 1, 1, 1)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 1, 1, 0, 1)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 2, 2, 1, 1)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 2, 2, 0, 0)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 2, 2, 0, 1)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 1, 2, 1, 1)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 1, 2, 1, 0)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 1, 2, 0, 0)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 1, 3, 1, 1)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 1, 3, 0, 0)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 2, 3, 1, 1)
    test_correlation_for_vec_sph_harm(theta, phi, dA, 2, 3, 0, 0)

def spherical_hn1(n,z,derivative=False):
    """ Spherical Hankel Function of the First Kind """
    return spherical_jn(n,z,derivative)+1j*spherical_yn(n,z,derivative)

def spherical_hn2(n,z,derivative=False):
    """ Spherical Hankel Function of the Second Kind """
    return spherical_jn(n,z,derivative)-1j*spherical_yn(n,z,derivative)

def vswf(f, c, tau, l, m, r, theta, phi):
    # p. 4 in scuffSpherical.pdf
    # tau = 1 -> "M_lm" -> vsh = X_lm
    # tau = 2 -> "N_lm" -> vxsh = Z_lm
    # c = 1 -> regular
    # c = 3 -> incoming
    # c = 4 -> outgoing

    c0 = 3e8
    k = 2*np.pi*f/c0

    if c == 1:
        rl = spherical_jn
    elif c == 2:
        rl = spherical_yn
    elif c == 3:
        rl = spherical_hn1
    elif c == 4:
        rl = spherical_hn2

    _, vsh_lm_theta, vsh_lm_phi = vec_sph_harm(tau, l, m, theta, phi)

    if tau == 1:
        rl_res = rl(l, k*r)

        f_theta = rl_res * vsh_lm_theta
        f_phi = rl_res * vsh_lm_phi
        f_r = 0
    elif tau == 2:
        rl_intermed = rl(l, k*r)/(k*r)
        rl_res = rl_intermed + rl(l, k*r, derivative=True)

        f_theta = 1j*rl_res * vsh_lm_theta
        f_phi = 1j*rl_res * vsh_lm_phi
        f_r = -rl_intermed*np.sqrt(l*(l+1))*my_sph_harm(l, m, theta, phi)
    else:
        raise NotImplemented('NIY')

    return f_r, f_theta, f_phi

x, y, z = coord_s2c(1, theta, phi)
#r, theta, phi = coord_c2s(x, y, z)
#x, y, z = coord_s2c(1, theta, phi)

dA = calc_dA(theta, phi)

print("s=1:")
print(vswf(1e9, 1, 1, 1, 1, np.array([1]), np.array([0]), np.array([0])))
print(vswf(1e9, 4, 1, 1, 1, np.array([1]), np.array([0]), np.array([0])))
print(vswf(1e9, 3, 1, 1, 1, np.array([1]), np.array([0]), np.array([0])))
print("s=2:")
print(vswf(1e9, 1, 2, 1, 1, np.array([1]), np.array([np.pi/2]), np.array([0])))
print(vswf(1e9, 4, 2, 1, 1, np.array([1]), np.array([np.pi/2]), np.array([0])))
print(vswf(1e9, 3, 2, 1, 1, np.array([1]), np.array([np.pi/2]), np.array([0])))

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

print("end")
