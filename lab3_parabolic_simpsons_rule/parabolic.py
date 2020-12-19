import math
import numpy as np
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial

class Parabolic:
    """ Parabolic method is based on the calculation of an intersection point with x of the parabola that was passed
    through three given points.
    """
    def __init__(self, left_border_a, right_border_b, initial_approximation_x0=None, max_iterations=20, tolerance=0.0001):
        """
        left_border_a -- left point (x, y)
        right_border_b -- right point (x, y)
        initial_approximation_x0 -- (x0, y0)
        """

        self.a = left_border_a
        self.b = right_border_b
        self.x0 = initial_approximation_x0

        current_approximation_x = (left_border_a[0] + right_border_b[0]) / 2 \
            if initial_approximation_x0 is None else initial_approximation_x0[0]
        current_a_x = left_border_a[0]
        current_b_x = right_border_b[0]

        # Getting coefficients of quadratic equation
        xs = [left_border_a[0], current_approximation_x, right_border_b[0]]
        ys = [left_border_a[1], 0 if initial_approximation_x0 is None else initial_approximation_x0[1],
              right_border_b[1]]

        iterations = 0
        while iterations < max_iterations:
            c0, c1, c2 = self.get_coefficients(xs, ys)
            print("\n--- Parabolic ---\nc0 =", c0, ", c1 =", c1, ", c2 =", c2)

            new_roots = self.get_roots(c0, c1, c2)
            print("roots = [", new_roots[0], ",", new_roots[1], "]")

            # If root is in the interval from x0 to b
            chosen_root_x = new_roots[0] if self.a[0] < new_roots[0] < self.b[0] else new_roots[1]

            if current_approximation_x < chosen_root_x < current_b_x:
                current_a_x = current_approximation_x
            else:
                current_b_x = current_approximation_x

            print("error =", "|", chosen_root_x, "-", current_approximation_x, "|",
                  abs(chosen_root_x - current_approximation_x), "<=", tolerance)
            if chosen_root_x == 0 or abs(chosen_root_x - current_approximation_x) <= tolerance:
                self.final_estimate = current_approximation_x
                print("---\niterations =", iterations, ", final estimate =", self.final_estimate)
                print("error =", "|", chosen_root_x, "-", current_approximation_x, "|",
                      abs(chosen_root_x - current_approximation_x), "<=", tolerance)
                print("\nEnd Parabolic ---\n")
                return

            current_approximation_x = chosen_root_x
            iterations += 1

            # Getting coefficients of quadratic equation
            xs = [current_a_x, current_approximation_x, current_b_x]
            ys = [self.func(c0, c1, c2, current_a_x),
                  self.func(c0, c1, c2, current_approximation_x),
                  self.func(c0, c1, c2, current_b_x)]
            print("xs =", xs)
            print("ys =", ys)

        self.final_estimate = current_approximation_x
        self.iterations = iterations
        print("---\niterations =", iterations, ", final estimate =", self.final_estimate)
        print("\nEnd Parabolic ---\n")

    def get_coefficients(self, x, y):
        denom = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
        print("denom =", denom)

        a = (x[2] * (y[1]-y[0]) + x[1] * (y[0]-y[2]) + x[0] * (y[2] - y[1])) / denom
        b = (x[2] * x[2] * (y[0]-y[1]) + x[1] * x[1] * (y[2]-y[0]) + x[0]*x[0] * (y[1]-y[2])) / denom
        c = (x[1] * x[2] * (x[1]-x[2]) * y[0]+x[2] * x[0] * (x[2]-x[0]) * y[1]+x[0] * x[1] * (x[0]-x[1]) * y[2]) / denom;
        return a, b, c
        #npx = np.array(x)
        #npy = np.array(y)
        #poly = lagrange(npx, npy)
        #coefs = Polynomial(poly).coef
        #return coefs[0], coefs[1], coefs[2]

    def get_interpolated_xy(self, resolution=10):
        out_vectorX = list()
        out_vectorY = list()

        # Interpolating a, x0, b
        step_a_x0 = (self.x0[0] - self.a[0]) / resolution
        current_x = self.a[0]
        for i in range(0, resolution):
            out_vectorY.append(self.func(self.a[1], self.x0[1], self.b[1], current_x))
            out_vectorX.append(current_x)
            current_x += step_a_x0

        step_x0_b = (self.b[0] - self.x0[0]) / resolution
        current_x = self.x0[0]
        for i in range(0, resolution):
            out_vectorY.append(self.func(self.a[1], self.x0[1], self.b[1], current_x))
            out_vectorX.append(current_x)
            current_x += step_x0_b

        return out_vectorX, out_vectorY

    def interpolate(self, a, x0, b, current_x):
        return current_x, self.func(self.a[1], self.x0[1], self.b[1], current_x)

    def func(self, c0, c1, c2, x):
        print("c0*x^2 + c1*x + c2 =", c0, "*", x ** 2, "+", c1, "*", x, "+", c2)
        return c0 + c1 * x + c2 * (x ** 2)
        #print("c2*x^2 + c1*x + c0 =", c2, "*", x**2, "+", c1, "*", x, "+", c0)
        #return c2 * (x ** 2) + c1 * x + c0

    def get_roots(self, c0, c1, c2):
        """ Returns roots of the quadratic equation.
        (-a1 +- sqrt(a1^2 - 4*a0*a2)) / 2a2
        """
        return (-c1 + math.sqrt(c1 ** 2 - 4*c0*c2)) / 2*c2, (-c1 - math.sqrt(c1 ** 2 - 4*c0*c2)) / 2*c2

