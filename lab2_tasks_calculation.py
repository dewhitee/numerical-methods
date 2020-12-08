import lab2_spline_interpolation.spline_interpolation as si
import lab2_least_squares.least_squares as ls
from matplotlib import pyplot as plt
import math
import numpy as np

# 1 task
least_squares_approx_result_0_0 = ls.LeastSquaresApproximator(
    vectorX=[-2.2 + delta for delta in np.arange(0, 4, 0.5)],
    vectorY=[0] * 8,
    ftype="custom",
    customfunc=lambda solvec, x: math.sin(5 * x) * math.exp(x),
    customstep=None,
    makeplot=True,
    k_approx_order=2,
    print_matrix=True,
    resolution=10
)

least_squares_approx_result_0_1 = ls.LeastSquaresApproximator(
    vectorX=[-2.2 + delta for delta in np.arange(0, 5, 0.5)],
    vectorY=[0] * 10,
    ftype="custom",
    customfunc=lambda solvec, x: math.sin(5 * x) * math.exp(x),
    customstep=None,
    makeplot=True,
    k_approx_order=4,
    print_matrix=True,
    resolution=20
)

#from scipy.optimize import brentq
#print("Root of Y = 1.5 is ", brentq(lambda x: math.sin(5 * x) * math.exp(x), 0, 10))

#from scipy.optimize import root
#print("Root of Y = 1.5 is ", root(lambda x: math.sin(5 * x) * math.exp(x), 0.5).x[0])

# initial_x = 1.337572
# y_from_x_0 = least_squares_approx_result_0_0.get_y_from_x(initial_x)
# initial_y = 1.7
# x_from_y_1 = least_squares_approx_result_0_0.get_x_from_y(initial_y, -0.0001, +0.0001, 10000, 0.000001)
# print("Y of X=", initial_x,":", y_from_x_0)
# print("X of Y=",initial_y,":", x_from_y_1)
#print("Y of X=", x_from_y_0, ":", y_from_x_1)

plt.figure("Comparison Least Squares")
plt.title("Comparison of Least Squares")
plt.plot(
    least_squares_approx_result_0_0.interpolated_vectorX,
    least_squares_approx_result_0_0.interpolated_vectorY,
    "go-",
    least_squares_approx_result_0_1.interpolated_vectorX,
    least_squares_approx_result_0_1.interpolated_vectorY,
    "y--",
    )
plt.legend(
    [
        "Least Squares with 5 points",
        "Least Squares with 10 points",
        "Known points"
    ]
)
plt.xlabel("X values")
plt.ylabel("Y values")
plt.show()

# 2 task
known_vectorX = [-3.2, -2.1, 0.4, 0.7, 2, 2.5, 2.777]
known_vectorY = [10, -2, 0, -7, 7, 0, 0]

csresult = si.CubicSplineInterpolator(
    known_vectorX=known_vectorX,
    known_vectorY=known_vectorY,
    known_points=None,
    vars=None
).get_xy(
    makeplot=False,
    resolution=30,
    showpoints=False
)

least_squares_approx_result_1 = ls.LeastSquaresApproximator(
    vectorX=known_vectorX,
    vectorY=known_vectorY,
    k_approx_order=2,
    makeplot=False,
    resolution=30
)

least_squares_approx_result_2 = ls.LeastSquaresApproximator(
    vectorX=known_vectorX,
    vectorY=known_vectorY,
    k_approx_order=3,
    makeplot=False,
    resolution=30
)

least_squares_approx_result_3 = ls.LeastSquaresApproximator(
    vectorX=known_vectorX,
    vectorY=known_vectorY,
    k_approx_order=5,
    makeplot=False,
    resolution=30,
    print_matrix=True
)

# Compare cubic spline interpolation with least squares graphically

plt.figure("Comparison of Cubic-spline interpolation with Least Squares")
plt.title("Comparison of Cubic-spline interpolation with Least Squares")
plt.plot(
    csresult[0],
    csresult[1],
    "g--",
    least_squares_approx_result_1.interpolated_vectorX,
    least_squares_approx_result_1.interpolated_vectorY,
    "c--",
    least_squares_approx_result_2.interpolated_vectorX,
    least_squares_approx_result_2.interpolated_vectorY,
    "r--",
    least_squares_approx_result_3.interpolated_vectorX,
    least_squares_approx_result_3.interpolated_vectorY,
    "y--",
    known_vectorX,
    known_vectorY,
    "bo")
plt.legend(
    [
        "Cubic-spline interpolation",
        "Least Squares with approximation order k = 2",
        "Least Squares with approximation order k = 3",
        "Least Squares with approximation order k = 5",
        "Known points"
    ]
)
plt.xlabel("X values")
plt.ylabel("Y values")
plt.show()
