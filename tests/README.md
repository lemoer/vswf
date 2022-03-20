
```tests/test_lpmv_diff.py```:

This test compares the implementation of ```lpmv_diff_times_minus_sqrt_1_minus_x_squared()``` with a numerical derivative of ``lpmv()``.

The function ```lpmv_diff_times_minus_sqrt_1_minus_x_squared()``` is defined by:

<img src="https://render.githubusercontent.com/render/math?math=-%5Csqrt%7B1-x%5E2%7D%5C%2C%20%5Cfrac%7B%5Cmathrm%7Bd%7D%7D%7B%5Cmathrm%7Bd%7Dx%7D%20%5C%2C%20P_%7Bl%7D%5Em%28x%29#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math=\color{white}-%5Csqrt%7B1-x%5E2%7D%5C%2C%20%5Cfrac%7B%5Cmathrm%7Bd%7D%7D%7B%5Cmathrm%7Bd%7Dx%7D%20%5C%2C%20P_%7Bl%7D%5Em%28x%29#gh-light-mode-only#gh-dark-mode-only">

```tests/test_vec_sph_harm.py```:

The test validates the following orthogonalities for the vector spherical harmonics:

<img src="https://render.githubusercontent.com/render/math?math=%5Cint_0%5E%5Cpi%5Cint_0%5E%7B2%5Cpi%7D%5Cmathbf%7BX%7D_%7B%5Ctau%20lm%7D%28%5Ctheta%2C%20%5Cphi%29%5C%2C%20%5Cmathbf%7BX%7D_%7B%5Ctilde%7B%5Ctau%7D%20%5Ctilde%7Bl%7D%20%5Ctilde%7Bm%7D%7D%28%5Ctheta%2C%20%5Cphi%29%20%5C%2C%5Csin%20%5Ctheta%5C%2C%5Cmathrm%7Bd%7D%5Cphi%5C%2C%5Cmathrm%7Bd%7D%5Ctheta%20%3D%20%5Cdelta_%7B%5Ctau%5Ctilde%7B%5Ctau%7D%7D%20%5C%2C%20%5Cdelta_%7Bl%5Ctilde%7Bl%7D%7D%20%5C%2C%20%5Cdelta_%7Bm%5Ctilde%7Bm%7D%7D#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math=\color{white}%5Cint_0%5E%5Cpi%5Cint_0%5E%7B2%5Cpi%7D%5Cmathbf%7BX%7D_%7B%5Ctau%20lm%7D%28%5Ctheta%2C%20%5Cphi%29%5C%2C%20%5Cmathbf%7BX%7D_%7B%5Ctilde%7B%5Ctau%7D%20%5Ctilde%7Bl%7D%20%5Ctilde%7Bm%7D%7D%28%5Ctheta%2C%20%5Cphi%29%20%5C%2C%5Csin%20%5Ctheta%5C%2C%5Cmathrm%7Bd%7D%5Cphi%5C%2C%5Cmathrm%7Bd%7D%5Ctheta%20%3D%20%5Cdelta_%7B%5Ctau%5Ctilde%7B%5Ctau%7D%7D%20%5C%2C%20%5Cdelta_%7Bl%5Ctilde%7Bl%7D%7D%20%5C%2C%20%5Cdelta_%7Bm%5Ctilde%7Bm%7D%7D#gh-dark-mode-only">

whereby the vector spherical harmonics are defined by:

<img src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7BX%7D_%7B1lm%7D%28%5Ctheta%2C%20%5Cphi%29%3D%5Cfrac%7Bj%7D%7B%5Csqrt%7Bl%28l%2B1%29%7D%7D%20%5Cleft%5B%5Cfrac%7Bjm%7D%7B%5Csin%20%5Ctheta%7D%20Y_%7Blm%7D%28%5Ctheta%2C%5Cphi%29%20%5Cmathbf%7Be%7D_%7B%5Ctheta%7D%20-%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20%5Ctheta%7D%20Y_%7Blm%7D%28%5Ctheta%2C%5Cphi%29%20%5Cmathbf%7Be%7D_%7B%5Cphi%7D%20%5Cright%5D#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7BX%7D_%7B2lm%7D%28%5Ctheta%2C%20%5Cphi%29%3D%5Cfrac%7Bj%7D%7B%5Csqrt%7Bl%28l%2B1%29%7D%7D%20%5Cleft%5B%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20%5Ctheta%7D%20Y_%7Blm%7D%28%5Ctheta%2C%5Cphi%29%20%5Cmathbf%7Be%7D_%7B%5Ctheta%7D%20%2B%20%5Cfrac%7Bjm%7D%7B%5Csin%20%5Ctheta%7D%20Y_%7Blm%7D%28%5Ctheta%2C%5Cphi%29%20%5Cmathbf%7Be%7D_%7B%5Cphi%7D%20%5Cright%5D#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math=\color{white}%5Cmathbf%7BX%7D_%7B1lm%7D%28%5Ctheta%2C%20%5Cphi%29%3D%5Cfrac%7Bj%7D%7B%5Csqrt%7Bl%28l%2B1%29%7D%7D%20%5Cleft%5B%5Cfrac%7Bjm%7D%7B%5Csin%20%5Ctheta%7D%20Y_%7Blm%7D%28%5Ctheta%2C%5Cphi%29%20%5Cmathbf%7Be%7D_%7B%5Ctheta%7D%20-%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20%5Ctheta%7D%20Y_%7Blm%7D%28%5Ctheta%2C%5Cphi%29%20%5Cmathbf%7Be%7D_%7B%5Cphi%7D%20%5Cright%5D#gh-dark-mode-only">
<img src="https://render.githubusercontent.com/render/math?math=\color{white}%5Cmathbf%7BX%7D_%7B2lm%7D%28%5Ctheta%2C%20%5Cphi%29%3D%5Cfrac%7Bj%7D%7B%5Csqrt%7Bl%28l%2B1%29%7D%7D%20%5Cleft%5B%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20%5Ctheta%7D%20Y_%7Blm%7D%28%5Ctheta%2C%5Cphi%29%20%5Cmathbf%7Be%7D_%7B%5Ctheta%7D%20%2B%20%5Cfrac%7Bjm%7D%7B%5Csin%20%5Ctheta%7D%20Y_%7Blm%7D%28%5Ctheta%2C%5Cphi%29%20%5Cmathbf%7Be%7D_%7B%5Cphi%7D%20%5Cright%5D#gh-dark-mode-only">

Similar as in [1].

```tests/test_vswf_values```:

This test validates the following equality:

<img src="https://render.githubusercontent.com/render/math?math=2%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%281%29%7D%28%5Cmathbf%7Br%7D%29%20%3D%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%284%29%7D%28%5Cmathbf%7Br%7D%29%20%2B%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%283%29%7D%28%5Cmathbf%7Br%7D%29#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}2%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%281%29%7D%28%5Cmathbf%7Br%7D%29%20%3D%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%284%29%7D%28%5Cmathbf%7Br%7D%29%20%2B%20%5Cmathbf%7BF%7D_%7B%5Ctau%2Cl%2Cm%7D%5E%7B%283%29%7D%28%5Cmathbf%7Br%7D%29}#gh-dark-mode-only">

[1] Reid, Homer: *Electromagnetism in the Spherical-Wave Basis* [https://homerreid.github.io/scuff-em-documentation/tex/scuffSpherical.pdf](https://homerreid.github.io/scuff-em-documentation/tex/scuffSpherical.pdf) - Accessed: 2022-03-20.