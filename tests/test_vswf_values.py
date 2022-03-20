#!/usr/bin/env python3

# add swf path to pythons search path
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from swf import *

def test_it():
    f = 1e9
    c0 = 3e8
    k = 2*np.pi*f/c0

    phi = np.linspace(0, 2*np.pi, 100)
    theta = np.linspace(0, np.pi, 50)
    theta, phi = np.meshgrid(theta, phi)
    r = 0.1*np.ones(phi.shape)

    for l in range(1, 5):
        for m in range(-l, l+1):
            for tau in [1, 2]:

                vswf_c = lambda c: np.array(vswf(k, c, tau, l, m, r, theta, phi))

                np.testing.assert_allclose(vswf_c(1)*2, vswf_c(4) + vswf_c(3), atol=0.01)

test_it()