import roots
import lab3_parabolic_simpsons_rule.parabolic as p
import matplotlib.pyplot as plt
import numpy as np

roots.RootFinder(lambda x: x ** 2 - x - 1, None) \
    .bisection(1, 2)

res = p.Parabolic(left_border_a=[1, None], right_border_b=[2, None],
                  customfunc=lambda val: val ** 3 - val - 1)
#p.Parabolic(
#    left_border_a=[-4, 14],
#    right_border_b=[0, -3],
#    customfunc=lambda val: val ** 2 - 3
#)

# Plot with spline interpolation
import lab2_spline_interpolation.spline_interpolation as si
import lab2_least_squares.least_squares as ls

#xs = [-4, 0, 2]
#ys = [35, -5, 11]
#ys = [14, -3, 1]
#si.CubicSplineInterpolator(xs, ys).get_xy(makeplot=True, showpoints=True)
#ls.LeastSquaresApproximator(xs, ys, makeplot=True,
#                            customfunc=lambda solvec, newx: newx ** 2)

cubicxs, cubicys = si.CubicSplineInterpolator(res.initial_xs, res.initial_ys).get_xy()

plt.figure("Least Squares approximation with different k values")
plt.title("Least Squares approximation with different k values")
plt.plot(
    cubicxs,
    cubicys,
    "c--",
    res.initial_xs,
    res.initial_ys,
    "yo-",
    res.xs,
    res.ys,
    "ro-",
    res.initial_approx_x,
    res.customfunc(res.initial_approx_x),
    "g^",
    res.final_estimate,
    res.customfunc(res.final_estimate),
    "go",
    res.found_root_1,
    res.customfunc(res.found_root_1),
    "bo",
    res.found_root_2,
    res.customfunc(res.found_root_2),
    "b^"
)
x = np.linspace(-10, 10, 1000)
y = res.customfunc(x)
plt.plot(x, y)
plt.legend(
    [
        "Cubic spline",
        "Parabolic (initial)",
        "Parabolic (final segment)",
        "Initial approximation",
        "Final estimate (" + str(res.final_estimate) + " ; " + str(res.customfunc(res.final_estimate)) + ")",
        "Root from parabolic 1 (" + str(res.found_root_1) + " ; " + str(res.customfunc(res.found_root_1)) + ")",
        "Root from parabolic 2 (" + str(res.found_root_2) + " ; " + str(res.customfunc(res.found_root_2)) + ")",
        "Parabola with linspace"
    ]
)
plt.axvline(x=0, color="grey")
plt.axhline(y=0, color='grey')
plt.axvline(x=res.final_estimate, color="green", linestyle="--")
plt.axhline(y=res.customfunc(res.final_estimate), color='green', linestyle="--")
plt.xlabel("X values")
plt.ylabel("Y values")
plt.show()


def calc_parabola_vertex(x1, y1, x2, y2, x3, y3):
    '''
	Adapted and modifed to get the unknowns for defining a parabola:
	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	'''

    denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
    A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
    B = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;
    C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

    return A, B, C


# Define your three known points
x1, y1 = [2, 11]
x2, y2 = [-4, 35]
x3, y3 = [0, -5]

# Calculate the unknowns of the equation y=ax^2+bx+c
a, b, c = calc_parabola_vertex(x1, y1, x2, y2, x3, y3)

# Define x range for which to calc parabola
import numpy as np

x_pos = np.arange(-30, 30, 1)
y_pos = []

# Calculate y values
for x in range(len(x_pos)):
    x_val = x_pos[x]
    y = (a * (x_val ** 2)) + (b * x_val) + c
    y_pos.append(y)

# Plot the parabola (+ the known points)
import matplotlib.pyplot as plt

plt.plot(x_pos, y_pos, linestyle='-.', color='black')  # parabola line
plt.scatter(x_pos, y_pos, color='gray')  # parabola points
plt.scatter(x1, y1, color='r', marker="D", s=50)  # 1st known xy
plt.scatter(x2, y2, color='g', marker="D", s=50)  # 2nd known xy
plt.scatter(x3, y3, color='k', marker="D", s=50)  # 3rd known xy
plt.show()

from scipy.interpolate import lagrange
x = np.array(xs)
y = np.array(ys)
poly = lagrange(x, y)

from numpy.polynomial.polynomial import Polynomial
print(Polynomial(poly).coef)
