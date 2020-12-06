import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import CubicSpline, interp1d

import lab2_least_squares.least_squares as ls
import lab2_spline_interpolation.spline_interpolation as si
import matrix_helpers as mh

vars = ['x1', 'x2', 'x3', 'x4', 'x5']
known_vectorX = [0, 3.3, 6.6, 9.9]
known_vectorY = [2.1, 5.9, 2.4, 3.4]

#cubic_spline_result = si.CubicSplineInterpolator(
#    known_vectorX=known_vectorX,
#    known_vectorY=known_vectorY,
#    vars=vars,
#    known_points=None
#)

#interpolatedX, interpolatedY = cubic_spline_result.get_xy(
#    resolution=30,
#    makeplot=True
#)
#print("interpolatedX:\n", interpolatedX, "\ninterpolatedY:\n", interpolatedY)

least_squares_approx_result = ls.LeastSquaresApproximator(
    vectorX=known_vectorX,
    vectorY=known_vectorY,
    k_approx_order=1,
    makeplot=True,
    ftype="linear"
)

known_vectorX_2 = [0, 3.3, 6.6, 9.9]
known_vectorY_2 = [12.1, 15.9, 12.4, 13.4]

csresult = si.CubicSplineInterpolator(
    known_vectorX=known_vectorX_2,
    known_vectorY=known_vectorY_2,
    vars=None,
    known_points=None
)

csresult.get_xy(
    resolution=30,
    makeplot=False
)

print("Matrix cond = ", mh.get_matrix_cond(csresult.matrixA))

least_squares_approx_result_1 = ls.LeastSquaresApproximator(
    vectorX=known_vectorX_2,
    vectorY=known_vectorY_2,
    k_approx_order=2,
    makeplot=False,
    ftype="auto",
)

least_squares_approx_result_2 = ls.LeastSquaresApproximator(
    vectorX=known_vectorX_2,
    vectorY=known_vectorY_2,
    k_approx_order=3,
    makeplot=False,
    ftype="auto",
    resolution=20
)

least_squares_approx_result_3 = ls.LeastSquaresApproximator(
    vectorX=known_vectorX_2,
    vectorY=known_vectorY_2,
    k_approx_order=4,
    makeplot=False
)

least_squares_approx_result_4 = ls.LeastSquaresApproximator(
    vectorX=known_vectorX_2,
    vectorY=known_vectorY_2,
    k_approx_order=5,
    makeplot=False
)

# ------------- LAB 2 GRAPHICAL COMPARISON

plt.figure("Comparison of Cubic-spline interpolation with Least Squares")
plt.title("Comparison of Cubic-spline interpolation with Least Squares")
plt.plot(
    csresult.get_xy()[0],
    csresult.get_xy()[1],
    "g--",
    least_squares_approx_result_1.interpolated_vectorX,
    least_squares_approx_result_1.interpolated_vectorY,
    "r--",
    least_squares_approx_result_2.interpolated_vectorX,
    least_squares_approx_result_2.interpolated_vectorY,
    "y--",
    least_squares_approx_result_3.interpolated_vectorX,
    least_squares_approx_result_3.interpolated_vectorY,
    "k--",
    least_squares_approx_result_4.interpolated_vectorX,
    least_squares_approx_result_4.interpolated_vectorY,
    "m--",
    known_vectorX_2,
    known_vectorY_2,
    "bo")
plt.legend(
    [
        "Cubic-spline interpolation",
        "Least Squares with approximation order k = 2",
        "Least Squares with approximation order k = 3",
        "Least Squares with approximation order k = 4",
        "Least Squares with approximation order k = 5",
        "Known points"
    ]
)
plt.show()
