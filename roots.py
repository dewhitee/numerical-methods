import numpy as np
import math
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial


class RootFinder:
    def __init__(self, function, solution_vectorX):
        self.function = function
        self.solvec = solution_vectorX
        self.final_estimate = 0
        self.iterations = 0

    def _description_line(self, method_name, a, b, tolerance):
        print(f'\n--- {method_name} ---\nTolerance =', tolerance, ", Root condition (f(a) * f(b) <= 0) =",
              self.function(a) * self.function(b) <= 0, "(",
              f'{"%0.5f" % self.function(a)} * {"%0.5f" % self.function(b)} = '
              f'{"%0.5f" % (self.function(a) * self.function(b))}',
              "<= 0)")

    def bisection(self, lower, upper, max_iterations=500, tolerance=0.0001):
        """ Finds the root of the function by the Bisection Method
        lower -- a
        upper -- b
        """
        self._description_line("Bisection", lower, upper, tolerance)
        print(f'{"Iteration":<10} | {"Estimate":<20} | {"Error":<16}')
        print('{:-<80}'.format(""))

        iterations = 0
        while iterations < max_iterations:
            midpoint = (lower + upper) / 2.0

            if midpoint == 0 or abs(upper - lower) <= tolerance:
                print(f'{iterations:<10} | {"%0.8f" % midpoint:<20} | {"%0.8f" % abs(upper - lower):<16}')
                print("--- End of Bisection ---\n")
                return midpoint
            #if self.function(self.solvec, midpoint) > 0:
            if self.function(midpoint) * self.function(lower) <= 0:
                upper = midpoint
            else:
                lower = midpoint
            iterations += 1
            print(f'{iterations:<10} | {"%0.8f" % ((lower + upper) / 2.0):<20} | {"%0.8f" % abs(upper - lower):<16}')
            #print("i =", iterations, "current estimate =", (lower + upper) / 2.0, ", error =", abs(upper - lower))
        self.final_estimate = (lower + upper) / 2.0
        self.iterations = iterations
        #print("---\niterations =", iterations, ", final estimate =", self.final_estimate)
        print(f'{iterations:<10} | {"%0.8f" % self.final_estimate:<16} | {"%0.8f" % abs(upper - lower):<16}')
        print("--- End of Bisection ---\n")
        return self.final_estimate

    def parabolic(self, a, b, max_iterations=50, tolerance=0.0001):
        iterations = 0
        x0 = (a + b) / 2
        self._description_line("Parabolic (Mueller's)", a, b, tolerance)
        print("Initial x's =", a, ",", x0, ",", b)
        print(f'{"Iteration":<10} | {"a0":<8} | {"a1":<8} | {"a2":<8} | {"Root 1":<12} | {"Root 2":<12} | '
              f'{"a":<8} | {"b":<8} | {"Estimate (xi)":<20} | {"Error":<16}')
        print('{:-<136}'.format(""))

        if self.function(a) == 0:
            return a
        elif self.function(b) == 0:
            return b

        current_x = x0
        prev_x = 0.0
        while iterations < max_iterations:
            # Step 3 - get Y values from X ---
            xs = [a, current_x, b]
            ys = [self.function(a), self.function(x0), self.function(b)]

            # Step 3 - get coefficients of a quadratic function ---
            a0, a1, a2 = self.get_coefficients_lagrange(xs, ys)

            # Step 4 - solving quadratic equation ---
            root_1, root_2 = self.solve_quadratic_equation(a0, a1, a2, current_x)

            # Update boundaries
            if self.function(a) * self.function(current_x) < 0:
                a = a
                b = current_x
            elif self.function(current_x) * self.function(b) < 0:
                a = current_x
                b = b

            # Chose more appropriate root
            if a <= root_1 <= b or a >= root_1 >= b:
                current_x = root_1
            elif a <= root_2 <= b or a >= root_2 >= b:
                current_x = root_2

            print(f'{iterations:<10} | {"%0.4f" % a0:<8} | {"%0.4f" % a1:<8} | {"%0.4f" % a2:<8} | '
                  f'{"%0.4f" % root_1:<12} | {"%0.4f" % root_2:<12} | '
                  f'{"%0.4f" % a:<8} | {"%0.4f" % b:<8} | {"%0.8f" % current_x:<20} | '
                  f'{"%0.8f" % abs(current_x - prev_x):<16}')

            if abs(current_x - prev_x) <= tolerance:
                print("--- End of Parabolic (Mueller's) ---\n")
                return current_x

            prev_x = current_x
            iterations += 1
        print("--- End of Parabolic (Mueller's) ---\n")
        return current_x

    def get_coefficients_lagrange(self, x, y):
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

        #print("delta2 =", delta2, ", delta1 =", delta1, ", h2 =", h2, ", h1 =", h1)

        a = (delta2 - delta1) / (h2 + h1)
        b = a * h2 + delta2
        c = self.function(x[1])
        return a, b, c


    def solve_quadratic_equation(self, a, b, c, x):
        """ Returns roots of the quadratic equation.
        (-a1 +- sqrt(a1^2 - 4*a0*a2)) / 2a2
        """
        #return (-c1 + math.sqrt(c1 ** 2 - 4 * c0 * c2)) / 2 * c2, (-c1 - math.sqrt(c1 ** 2 - 4 * c0 * c2)) / 2 * c2
        return x + (-2*c)/(b + math.sqrt(b ** 2 - 4*a*c)), x + (-2*c)/(b - math.sqrt(b ** 2 - 4*a*c))

    def secant(self, lower_approx_x, upper_approx_x, max_iterations=500, tolerance=0.0001):
        exit
        #iterations = 1
        #while iterations < max_iterations and abs(upper_approx_x - lower_approx_x) > tolerance:
        #    out_x = upper_approx_x - (self.function())

    #def brent(self, )

    #def custom(self, x, y)
