import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

# add swf path to pythons search path
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from swf import *

# In this test, we compare the analytical derivation implemented in
# lpmv_diff_times_minus_sqrt_1_minus_x_squared(..., np.cos(theta))
# with a numerical derivation of lpmv(..., np.cos(theta)).

def test_it():

    if __name__ == '__main__':
        import matplotlib.pyplot as plt
        import matplotlib
        from matplotlib import cm, colors
        from mpl_toolkits.mplot3d import Axes3D
        plt.figure()

    theta = np.linspace(0, np.pi, 50)
    dtheta = theta[1] - theta[0]

    # theta values in between the original values ("half-step")
    theta_ = np.convolve(theta, [0.5, 0.5], 'valid')
    
    for m in range(1, 5):
        for l in range(3, 4):
            leg_part = lpmv(m, l, np.cos(theta))
            leg_diff_analytic = lpmv_diff_times_minus_sqrt_1_minus_x_squared(m, l, np.cos(theta_))
            leg_diff_numeric = np.diff(leg_part)/dtheta

            if __name__ == '__main__':

                plt.plot(theta, leg_part)
                # The following two graphs should be equal
                plt.plot(theta_, leg_diff_analytic)
                plt.plot(theta_, leg_diff_numeric, ':')
                plt.show(block=True)

            np.testing.assert_allclose(leg_diff_analytic, leg_diff_numeric, atol=0.02)

test_it()