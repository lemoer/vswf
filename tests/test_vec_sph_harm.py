import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

# add swf path to pythons search path
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from swf import *

def __my_test_vec_sph_harm(theta, phi, dA, tau1, tau2, m1, m2):
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

    if __name__ == '__main__':
        import matplotlib.pyplot as plt
        import matplotlib
        from matplotlib import cm, colors
        from mpl_toolkits.mplot3d import Axes3D

        plt.figure()
        plt.imshow(np.abs(C))
        plt.colorbar()
        plt.title(f"$\\tau_1 = {tau1}, \\tau_2 = {tau2}, m_1 = {m1}, m_2 = {m2}$")
        plt.show(block=True)
    
    if m1 != m2 or tau1 != tau2:
        expectation = np.zeros((L, L))
    else:
        expectation = np.eye(L)

    np.testing.assert_allclose(np.abs(C), expectation, atol=0.02)

def test_it():
    phi = np.linspace(0, 2*np.pi, 100)
    theta = np.linspace(0, np.pi, 50)
    theta, phi = np.meshgrid(theta, phi)
    dA = calc_dA(theta, phi)

    __my_test_vec_sph_harm(theta, phi, dA, 1, 1, 1, 1)
    __my_test_vec_sph_harm(theta, phi, dA, 1, 1, 0, 1)
    __my_test_vec_sph_harm(theta, phi, dA, 2, 2, 1, 1)
    __my_test_vec_sph_harm(theta, phi, dA, 2, 2, 0, 0)
    __my_test_vec_sph_harm(theta, phi, dA, 2, 2, 0, 1)
    __my_test_vec_sph_harm(theta, phi, dA, 1, 2, 1, 1)
    __my_test_vec_sph_harm(theta, phi, dA, 1, 2, 1, 0)
    __my_test_vec_sph_harm(theta, phi, dA, 1, 2, 0, 0)
    __my_test_vec_sph_harm(theta, phi, dA, 1, 3, 1, 1)
    __my_test_vec_sph_harm(theta, phi, dA, 1, 3, 0, 0)
    __my_test_vec_sph_harm(theta, phi, dA, 2, 3, 1, 1)
    __my_test_vec_sph_harm(theta, phi, dA, 2, 3, 0, 0)

test_it()