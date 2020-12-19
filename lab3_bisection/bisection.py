
class Bisection:
    """ Method of partial division. The slowest method from all of all methods of a linear convergence.
    However, this method has absolute resistance, meaning that the convergence is guaranteed.
    """
    def __init__(self):
        exit

    def check_roots(self, a, b, func):
        return func(a) * func(b) <= 0

    def calculate_middle_of_line_segment(self, a, b):
        x = (a + b) / 2
        return x

    def set_a_or_b(self, a, b, x, func):
        if func(a) * func(x) <= 0:
            b = x
        else:
            a = x


