import lab2_spline_interpolation.spline_interpolation as si
import lab2_least_squares.least_squares as ls
from matplotlib import pyplot as plt
import math
import numpy as np

# 1 task
least_squares_approx_result_0_0 = ls.LeastSquaresApproximator(
    vectorX=[-2.2 + delta for delta in np.arange(0, 2.5, 0.5)],
    vectorY=[0] * 5,
    ftype="custom",
    customfunc=lambda solvec, x: math.sin(5 * x) * math.exp(x),
    customstep=None,
    makeplot=False,
    k_approx_order=2,
    print_matrix=True,
    resolution=5
)

least_squares_approx_result_0_1 = ls.LeastSquaresApproximator(
    vectorX=[-2.2 + delta for delta in np.arange(0, 5, 0.5)],
    vectorY=[0] * 10,
    ftype="custom",
    customfunc=lambda solvec, x: math.sin(5 * x) * math.exp(x),
    customstep=None,
    makeplot=False,
    k_approx_order=4,
    print_matrix=True,
    resolution=10
)

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
    resolution=30
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