
class RootFinder:
    def __init__(self, function, solution_vectorX):
        self.function = function
        self.solvec = solution_vectorX

    def bisection(self, lower, upper, max_iterations=500, tolerance=0.0001):
        """ Finds the root of the function by the Bisection Method
        """
        iterations = 0
        while iterations < max_iterations:
            midpoint = (lower + upper) / 2.0

            if midpoint == 0 or abs(upper - lower) < tolerance:
                return midpoint
            if self.function(self.solvec, midpoint) > 0:
                upper = midpoint
            else:
                lower = midpoint
            iterations += 1
        self.final_estimate = (lower + upper) / 2.0
        self.iterations = iterations
        return self.final_estimate

    def secant(self, lower_approx_x, upper_approx_x, max_iterations=500, tolerance=0.0001):
        exit
        #iterations = 1
        #while iterations < max_iterations and abs(upper_approx_x - lower_approx_x) > tolerance:
        #    out_x = upper_approx_x - (self.function())

    #def brent(self, )

    #def custom(self, x, y)
