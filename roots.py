
class RootFinder:
    def __init__(self, function, solution_vectorX):
        self.function = function
        self.solvec = solution_vectorX
        self.final_estimate = 0
        self.iterations = 0

    def bisection(self, lower, upper, max_iterations=500, tolerance=0.0001):
        """ Finds the root of the function by the Bisection Method
        lower -- a
        upper -- b
        """
        iterations = 0
        while iterations < max_iterations:
            midpoint = (lower + upper) / 2.0

            if midpoint == 0 or abs(upper - lower) < tolerance:
                print("---\niterations =", iterations, ", final estimate =", midpoint)
                return midpoint
            #if self.function(self.solvec, midpoint) > 0:
            if self.function(midpoint) * self.function(lower) <= 0:
                upper = midpoint
            else:
                lower = midpoint
            iterations += 1
            print("i =", iterations, "current estimate =", (lower + upper) / 2.0, ", error =", abs(upper - lower))
        self.final_estimate = (lower + upper) / 2.0
        self.iterations = iterations
        print("---\niterations =", iterations, ", final estimate =", self.final_estimate)
        return self.final_estimate

    def parabolic(self):
        exit

    def secant(self, lower_approx_x, upper_approx_x, max_iterations=500, tolerance=0.0001):
        exit
        #iterations = 1
        #while iterations < max_iterations and abs(upper_approx_x - lower_approx_x) > tolerance:
        #    out_x = upper_approx_x - (self.function())

    #def brent(self, )

    #def custom(self, x, y)
