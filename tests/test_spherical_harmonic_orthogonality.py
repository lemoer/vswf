#!/usr/bin/env python3

# add swf path to pythons search path
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from swf import *

def test_it():

    phi = np.linspace(0, 2*np.pi, 100)
    theta = np.linspace(0, np.pi, 50)
    theta, phi = np.meshgrid(theta, phi)
    dA = calc_dA(theta, phi)

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
            
            C[l-1, l2-1] = np.dot(f2.conj(), dA*f1)

    expectation = np.eye(L)
    np.testing.assert_allclose(np.abs(C), expectation, atol=0.01)

    if __name__ == '__main__':
        import matplotlib.pyplot as plt
        import matplotlib
        from matplotlib import cm, colors
        from mpl_toolkits.mplot3d import Axes3D

        plt.figure()
        plt.imshow(np.abs(C))
        plt.colorbar()

test_it()