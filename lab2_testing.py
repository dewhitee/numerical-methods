from scipy.interpolate import CubicSpline
import lab2_spline_interpolation.spline_interpolation as si
import lab2_least_squares.least_squares as ls
from scipy import interpolate
from scipy.interpolate import interp1d

import numpy as np
import matplotlib.pyplot as plt

vars = ['x1', 'x2', 'x3', 'x4', 'x5']
known_vectorX = [0, 3.3, 6.6, 9.9]
known_vectorY = [2.1, 5.9, 2.4, 3.4]

cubic_spline_result = si.CubicSplineInterpolator(
    known_vectorX=known_vectorX,
    known_vectorY=known_vectorY,
    vars=vars,
    known_points=None
)

interpolatedX, interpolatedY = cubic_spline_result.get_xy(
    resolution=30,
    makeplot=True
)
#print("interpolatedX:\n", interpolatedX, "\ninterpolatedY:\n", interpolatedY)

least_squares_approx_result = ls.least_squares(
    vectorX=known_vectorX,
    vectorY=known_vectorY,
    k_approx_order=1,
    makeplot=True
)

cs = CubicSpline(known_vectorX, known_vectorY, bc_type='natural')
print("S(0)", cs(0))
print("S(3.3)", cs(3.3))

#plt.plot()
