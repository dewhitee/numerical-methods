import numpy as np
import math
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial


class RootFinder:
    def __init__(self, function, solution_vectorX=None):
        self.function = function
        self.solvec = solution_vectorX
        self.final_estimate = 0
        self.iterations = 0

    def _description_line(self, method_name, a, b, tolerance):
        cond = self.function(a) * self.function(b) <= 0
        print(f'\n--- {method_name} ---\na = {a}, b = {b}, '
              f'Tolerance =', tolerance, ", Root condition (f(a) * f(b) <= 0) =",
              cond, "(",
              f'{self.function(a)} * {self.function(b)} = '
              f'{"%0.12f" % (self.function(a) * self.function(b))}',
              "<= 0)")

        if not cond:
            print("No roots found on the [", a, ";", b, "] interval.")

    def _initial_check(self, a, b, tolerance, check_initial):
        if check_initial:
            if abs(self.function(a)) <= tolerance:
                print(f'initial y value of a was ~== 0 (with precision = {tolerance}), returning a ({a})')
                return True, a
            elif abs(self.function(b)) <= tolerance:
                print(f'initial y value of b was ~== 0 (with precision = {tolerance}), returning b ({b})')
                return True, b
        return False, 0

    def bisection(self, lower, upper, max_iterations=500, tolerance=0.0001, check_initial=True):
        """ Finds the root of the function by the Bisection Method
        lower -- a
        upper -- b
        """

        # Step 1 - Print available info, checking root condition f(a) * f(b) <= 0
        self._description_line("Bisection", lower, upper, tolerance)
        print(f'{"Iteration":<10} | {"a":<16} | {"Estimate":<20} | {"b":<16} | {"Error":<16} | {"f(a) * f(x)":<16}')
        print('{:-<100}'.format(""))

        # (Modification) check for initial value of f(a) and f(b). If one of them is ~= 0 (with tolerance) then
        # return corresponding x value
        initial_check = self._initial_check(lower, upper, tolerance, check_initial)
        if initial_check[0]:
            return initial_check[1]

        iterations = 0
        while iterations < max_iterations:
            midpoint = (lower + upper) / 2.0
            midpoint_y = self.function(midpoint)
            print(f'{iterations:<10} | {"%0.8f" % lower:<16} | {"%0.8f" % midpoint:<20} | {"%0.8f" % upper:<16} | '
                  f'{"%0.8f" % abs(upper - lower):<16} | {"%0.8f" % (self.function(lower) * midpoint_y):<16}')

            # Check if f(a) and f(x) have opposite signs
            if self.function(lower) * midpoint_y <= 0:
                upper = midpoint
            else:
                lower = midpoint

            # Checking if current error is less than tolerance
            if midpoint_y == 0 or abs(upper - lower) <= tolerance:
                print("--- End of Bisection ---\n")
                return midpoint

            iterations += 1
        self.final_estimate = (lower + upper) / 2.0
        self.iterations = iterations
        print(f'{iterations:<10} | {"%0.8f" % lower:<16} | {"%0.8f" % self.final_estimate:<20} | '
              f'{"%0.8f" % upper:<16} | {"%0.8f" % abs(upper - lower):<16}')
        print("--- End of Bisection ---\n")
        return self.final_estimate

    def parabolic(self, a, b, max_iterations=50, tolerance=0.0001, check_initial=True):
        iterations = 0

        # Step 1 - Print available info, checking root condition f(a) * f(b) <= 0
        self._description_line("Parabolic (Mueller's)", a, b, tolerance)

        # Step 2 - Initializing initial x
        x0 = (a + b) / 2

        print("Initial x's =", a, ",", x0, ",", b)
        print(f'{"Iteration":<10} | {"a0":<8} | {"a1":<8} | {"a2":<8} | {"Root 1":<12} | {"Root 2":<12} | '
              f'{"a":<8} | {"b":<8} | {"Estimate (xi)":<20} | {"Error":<16}')
        print('{:-<136}'.format(""))

        # (Modification) check for initial value of f(a) and f(b). If one of them is ~= 0 (with tolerance) then
        # return corresponding x value
        initial_check = self._initial_check(a, b, tolerance, check_initial)
        if initial_check[0]:
            return initial_check[1]

        current_x = x0
        prev_x = 0.0
        while iterations < max_iterations:
            # Step 3.1 - get Y values from X ---
            xs = [a, current_x, b]
            ys = [self.function(a), self.function(x0), self.function(b)]

            # Step 3.2 - get coefficients of a quadratic function using lagrange or other method ---
            a0, a1, a2 = self._get_coefficients(xs, ys)

            # Step 4.1 - solving quadratic equation ---
            root_1, root_2 = self._solve_quadratic_equation(a0, a1, a2, current_x)

            # Step 4.2 - Update boundaries (checking for opposite signs)
            # If the f(a) * f(x) < 0    => set right border to x
            # If the f(b) * f(x) < 0    => set left border to x
            if self.function(a) * self.function(current_x) < 0:
                a = a
                b = current_x
            elif self.function(current_x) * self.function(b) < 0:
                a = current_x
                b = b

            # Step 4.3 - Chose more appropriate root (new approximation of x) that belongs to the current [a; b] interval
            if a <= root_1 <= b or a >= root_1 >= b:
                current_x = root_1
            elif a <= root_2 <= b or a >= root_2 >= b:
                current_x = root_2

            print(f'{iterations:<10} | {"%0.4f" % a0:<8} | {"%0.4f" % a1:<8} | {"%0.4f" % a2:<8} | '
                  f'{"%0.4f" % root_1:<12} | {"%0.4f" % root_2:<12} | '
                  f'{"%0.4f" % a:<8} | {"%0.4f" % b:<8} | {"%0.8f" % current_x:<20} | '
                  f'{"%0.8f" % abs(current_x - prev_x):<16}')

            # Step 5 - checking error
            if abs(current_x - prev_x) <= tolerance:
                print("--- End of Parabolic (Mueller's) ---\n")
                return current_x

            # Set current root as the previous root
            prev_x = current_x
            iterations += 1

        print("--- End of Parabolic (Mueller's) ---\n")
        return current_x

    def _get_coefficients(self, x, y):
        #print("(", x[0], "-", x[1], ") * (", x[0], "-", x[1], ") * (", x[1], "—", x[2], ") =",
        #      (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2]))
        #denom = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
        #print("denom =", denom)

        #a = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / denom
        #b = (x[2] * x[2] * (y[0] - y[1]) + x[1] * x[1] * (y[2] - y[0]) + x[0] * x[0] * (y[1] - y[2])) / denom
        #c = (x[1] * x[2] * (x[1] - x[2]) * y[0] + x[2] * x[0] * (x[2] - x[0]) * y[1] + x[0] * x[1] * (x[0] - x[1]) * y[
        #    2]) / denom
        #return a, b, c
        #npx = np.array(x)
        #npy = np.array(y)
        #poly = lagrange(npx, npy)
        #coefs = Polynomial(poly).coef
        #return coefs[0], coefs[1], coefs[2]

        # https://atozmath.com/example/CONM/Bisection.aspx?he=e&q=mu&ex=1
        # http://simenergy.ru/math-analysis/solution-methods/43-muller-method
        h1 = x[2] - x[0]
        h2 = x[1] - x[2]

        delta1 = (self.function(x[2]) - self.function(x[0])) / h1
        delta2 = (self.function(x[1]) - self.function(x[2])) / h2

        a = (delta2 - delta1) / (h2 + h1)
        b = a * h2 + delta2
        c = self.function(x[1])
        return a, b, c


    def _solve_quadratic_equation(self, a, b, c, x):
        """ Returns roots of the quadratic equation.
        (-a1 +- sqrt(a1^2 - 4*a0*a2)) / 2a2
        """
        #return (-c1 + math.sqrt(c1 ** 2 - 4 * c0 * c2)) / 2 * c2, (-c1 - math.sqrt(c1 ** 2 - 4 * c0 * c2)) / 2 * c2
        return x + (-2*c)/(b + math.sqrt(b ** 2 - 4*a*c)), x + (-2*c)/(b - math.sqrt(b ** 2 - 4*a*c))

    def secant(self, a, b, max_iterations=500, tolerance=0.0001):
        self._description_line("Secant", a, b, tolerance)
        iterations = 0
        condition = True
        #while condition:
        #    if self.function()


    #def brent(self, )

    #def custom(self, x, y)
