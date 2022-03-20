
```tests/test_lpmv_diff.py```:

This test compares the implementation of ```lpmv_diff_times_minus_sqrt_1_minus_x_squared()``` with a numerical derivative of ``lpmv()``.

The function ```lpmv_diff_times_minus_sqrt_1_minus_x_squared()``` is defined by:
<img src="https://render.githubusercontent.com/render/math?math=-%5Csqrt%7B1-x%5E2%7D%5C%2C%20%5Cfrac%7B%5Cmathrm%7Bd%7D%7D%7B%5Cmathrm%7Bd%7Dx%7D%20%5C%2C%20P_%7Bl%7D%5Em%28x%29#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math=\color{white}-%5Csqrt%7B1-x%5E2%7D%5C%2C%20%5Cfrac%7B%5Cmathrm%7Bd%7D%7D%7B%5Cmathrm%7Bd%7Dx%7D%20%5C%2C%20P_%7Bl%7D%5Em%28x%29#gh-light-mode-only#gh-dark-mode-only">


```tests/test_vswf_values```:

This test validates the following equality:

<img src="https://render.githubusercontent.com/render/math?math=2%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%281%29%7D%28%5Cmathbf%7Br%7D%29%20%3D%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%284%29%7D%28%5Cmathbf%7Br%7D%29%20%2B%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%283%29%7D%28%5Cmathbf%7Br%7D%29#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}2%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%281%29%7D%28%5Cmathbf%7Br%7D%29%20%3D%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%284%29%7D%28%5Cmathbf%7Br%7D%29%20%2B%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%283%29%7D%28%5Cmathbf%7Br%7D%29}#gh-dark-mode-only">