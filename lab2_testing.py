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
    makeplot=True,
    ftype="linear"
)

cs = CubicSpline(known_vectorX, known_vectorY, bc_type='natural')
print("S(0)", cs(0))
print("S(3.3)", cs(3.3))

#plt.plot()


#def fun_rosenbrock(x):
#    return np.array([10 * (x[1] - x[0]**2), (1 - x[0])])

#from scipy.optimize import least_squares
#res_1 = least_squares(fun_rosenbrock, known_vectorX)

si.CubicSplineInterpolator(
    known_vectorX=[1.2, 3.4, 6.4, 2.3, 8.5, 9.2, 3.2, 2.2, 1.2, 0.2],
    known_vectorY=[3.2, 5.2, 8.9, 16.3, 21.2, 5.3, 55.2, 1.4, 2.1, 1.8],
    vars=None,
    known_points=None
).get_xy(
    resolution=30,
    makeplot=True
)

least_squares_approx_result_2 = ls.least_squares(
    vectorX=known_vectorX,
    vectorY=known_vectorY,
    k_approx_order=1,
    makeplot=True,
    ftype="custom",
    customfunc=lambda x, y: x*y
)
