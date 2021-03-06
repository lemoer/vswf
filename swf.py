import numpy as np
from scipy.special import sph_harm, factorial, lpmv, spherical_jn, spherical_yn, eval_legendre

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

def sph_harm_diff_phi(l, m, theta, phi):
    # differentation of my_sph_harm() w.r.t. phi (with conventions
    # that phi is the azimuth.)
    return 1j*m*my_sph_harm(l, m, theta, phi)

def lpmv_diff_times_minus_sqrt_1_minus_x_squared(m, l, x):
    # This function calculates the derivative of the legendre function lpmv()
    # multiplied with a factor of -sqrt(1-x??), as in Kristenson 2014 p. 84
    # (the last formula on the page, with an additional minus on both sides
    # of the equation). This makes sense, because we actually want to derive 
    # lpmv(m, l, np.cos(theta)) and the factor -sqrt(1-x??) exactly corresponds
    # to the inner derivative of the np.cos(theta). So this function here 
    # can be used directly to calculate the derivative. See test_lpmv_diff_stuff().
    return -(1/2)*((l-m+1)*(l+m)*lpmv(m-1, l, x)-lpmv(m+1, l, x)) # grundmann

def sph_harm_diff_theta(l, m, theta, phi):
    # we leave out the (-1)**m here in c_lm, because it is already
    # part of the legendre definition in lpmv(). 
    c_lm = np.sqrt((2*l+1)/(4*np.pi)*factorial(l-m)/factorial(l+m))

    exp_part = np.exp(1j*m*phi)
    # Inner derivation of the lpmv() argument:
    # x = np.cos(theta)
    # dx/dtheta = -sin(theta) = -sqrt(1-cos(theta)??) = -sqrt(1-x??)
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

def __fix_dim(A, axis):
    # Repeats the last column or row of a matrix once, so the dimension
    # in that axis is increased by 1. This is useful in combination with
    # np.diff() as np.diff() reduces the size of an array.
    if axis == 0:
        return np.append(A, np.reshape(A[-1,:], (1, A.shape[1])), 0)
    elif axis == 1:
        return np.append(A, np.reshape(A[:,-1], (A.shape[0], 1)), 1)
    else:
        raise Exception('__fix_dim() is not implemented for axes > 1!')

def calc_dA(theta, phi):
    assert theta[0,0] == theta[1,0], 'Grid direction is not as expected!'
    dtheta = __fix_dim(np.diff(theta, axis=1), axis=1)
    dphi = __fix_dim(np.diff(phi, axis=0), axis=0)
    return dtheta * dphi * np.sin(theta)

def spherical_hn1(n,z,derivative=False):
    """ Spherical Hankel Function of the First Kind """
    return spherical_jn(n,z,derivative)+1j*spherical_yn(n,z,derivative)

def spherical_hn2(n,z,derivative=False):
    """ Spherical Hankel Function of the Second Kind """
    return spherical_jn(n,z,derivative)-1j*spherical_yn(n,z,derivative)

def vswf(k, c, tau, l, m, r, theta, phi):
    # p. 4 in scuffSpherical.pdf
    # tau = 1 -> "M_lm" -> vsh = X_lm
    # tau = 2 -> "N_lm" -> vxsh = Z_lm
    # c = 1 -> regular
    # c = 3 -> incoming
    # c = 4 -> outgoing

    if c == 1:
        rl = spherical_jn
    elif c == 2:
        rl = spherical_yn
    elif c == 3:
        rl = spherical_hn1
    elif c == 4:
        rl = spherical_hn2

    _, vsh_lm_theta, vsh_lm_phi = vec_sph_harm(tau, l, m, theta, phi) # This is either called Xlm or Zlm (depending on tau)

    if tau == 1:
        rl_res = rl(l, k*r)

        f_theta = rl_res * vsh_lm_theta
        f_phi = rl_res * vsh_lm_phi
        f_r = np.zeros(f_phi.shape)
    elif tau == 2:
        rl_intermed = rl(l, k*r)/(k*r)
        rl_res = rl_intermed + rl(l, k*r, derivative=True)

        f_theta = 1*rl_res * vsh_lm_theta
        f_phi = 1*rl_res * vsh_lm_phi
        f_r = rl_intermed*np.sqrt(l*(l+1))*my_sph_harm(l, m, theta, phi)
    else:
        raise NotImplemented('NIY')

    return f_r, f_theta, f_phi

def vswf_time_domain(t, c, tau, l, m, r, theta, phi):
    # p. 4 in scuffSpherical.pdf
    # tau = 1 -> "M_lm" -> vsh = X_lm
    # tau = 2 -> "N_lm" -> vxsh = Z_lm
    # c = 1 -> regular
    # c = 3 -> incoming
    # c = 4 -> outgoing

    if c == 1:
        rl = spherical_jn
    else:
        raise Exception('NIY')

    _, vsh_lm_theta, vsh_lm_phi = vec_sph_harm(tau, l, m, theta, phi) # This is either called Xlm or Zlm (depending on tau)

    boxcar = lambda x: np.heaviside(x+1, 0.5)-np.heaviside(x-1, 0.5)

    c0 = 299792458 # speed of light
    a = r/c0 # scaling factor
    # k*r = omega*a
    if tau == 1:
        # Frequency domain: rl_res = rl(l, omega*a)
        rl_res = 1/(2*(1j**l))/a*eval_legendre(l, t/a) * boxcar(t/a)

        f_theta = rl_res * vsh_lm_theta
        f_phi = rl_res * vsh_lm_phi
        f_r = np.zeros(f_phi.shape)
    elif tau == 2:
        # Frequency domain: rl_intermed = rl(l, omega*a)/(omega*a)
        rl_intermed = -1/(2*(1j**l))*1j/a*(eval_legendre(l+1, t/a) - eval_legendre(l-1, t/a))/(2*l+1)*boxcar(t/a)
        # Frequency domain: rl_res = rl_intermed + rl(l, omega*a, derivative=True)
        rl_res = rl_intermed + 1j*(t/a)*1/(2*(1j**l))/a*eval_legendre(l, t/a) * boxcar(t/a)

        f_theta = 1*rl_res * vsh_lm_theta
        f_phi = 1*rl_res * vsh_lm_phi
        f_r = rl_intermed*np.sqrt(l*(l+1))*my_sph_harm(l, m, theta, phi)
    else:
        raise NotImplemented('NIY')

    return f_r, f_theta, f_phi
    

def face(dn, n, da, sign = 1):
    d1 = np.mod(dn+1, 3) # dimension normal to dn
    d2 = np.mod(dn+2, 3) # other dimension normal to dn

    min1, max1 = -(n[d1]-1)/2, (n[d1]-1)/2
    min2, max2 = -(n[d2]-1)/2, (n[d2]-1)/2

    x1 = np.linspace(min1, max1, n[d1])*da
    x2 = np.linspace(min2, max2, n[d2])*da

    c1, c2 = np.meshgrid(x1, x2) # coordinates in the face

    xn = sign*np.ones_like(c1)*(n[dn]-1)/2*da

    r = np.zeros((3, n[d2], n[d1]))
    r[d1, :, :] = c1
    r[d2, :, :] = c2
    r[dn, :, :] = xn

    dA = np.zeros((3, n[d2], n[d1]))
    dA[dn, :, :] = da**2

    return r, dA

class EMField:
    
    def __init__(self, k, pos, dA):
        self.k = k
        self.eta = 120*np.pi
        self.pos = pos
        self.dA = dA
        self.E = None
        self.H = None

    def override_by_vswf(self, c, tau, l, m):
        r, theta, phi = coord_c2s(self.pos[0,:,:], self.pos[1,:,:], self.pos[2,:,:])

        f_r, f_theta, f_phi = vswf(self.k, c, tau, l, m, r, theta, phi)
        f_x, f_y, f_z = vec_s2c(f_r, f_theta, f_phi, r, theta, phi)

        self.E = np.zeros(dtype='complex128', shape=(3, np.size(f_x,0), np.size(f_x,1)))
        factor = self.k*np.sqrt(self.eta)/np.sqrt(2)
        self.E[0, :, :] = factor*f_x
        self.E[1, :, :] = factor*f_y
        self.E[2, :, :] = factor*f_z
        
        tau_tilde = 2 if tau == 1 else 1
        f_r, f_theta, f_phi = vswf(self.k, c, tau_tilde, l, m, r, theta, phi)
        f_x, f_y, f_z = vec_s2c(f_r, f_theta, f_phi, r, theta, phi)
        
        self.H = np.zeros(dtype='complex128', shape=(3, np.size(f_x,0), np.size(f_x,1)))
        factor = 1j*self.k/np.sqrt(self.eta)/np.sqrt(2)
        self.H[0, :, :] = factor*f_x
        self.H[1, :, :] = factor*f_y
        self.H[2, :, :] = factor*f_z
