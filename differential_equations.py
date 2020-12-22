import numpy as np
import copy


class DifferentialEquationsSolver:
    def __init__(self, function, step_h, exact_ys=None):
        """
        :param function: Two parameter function
        """
        self.function = function
        self.step_h = step_h
        self.method_name = ""
        self.exact_ys = exact_ys

    @staticmethod
    def compare_errors(approx_ys, exact_ys):
        return [abs(y2 - y1) for y1, y2 in zip(approx_ys[1:], exact_ys[1:])]

    @staticmethod
    def get_global_errors(errors):
        sum = 0
        global_errors = []
        for e in errors:
            sum += e
            global_errors.append(sum)
        return global_errors

    def _get_next_y(self, current_x, current_y):
        return self.function(current_x, current_y)

    def _description_line(self, initial_x, initial_y, last_x):
        print(f'\n--- {self.method_name} ---')
        print("h (step) =", self.step_h, ", initial x =", initial_x, ", initial y =", initial_y, ", last x =", last_x)

    def _table_header(self):
        if self.exact_ys is None:
            print(f'{"x":<16} | {"y":<16}')
            print('{:-<50}'.format(""))
        else:
            print(f'{"x":<16} | {"y":<16} | {"Exact value":<16} | {"Error":<16}')
            print('{:-<70}'.format(""))

    def _get_exact_y_and_error(self, current_x, current_y, i):
        exact_y = self.function(current_x, current_y)
        error = abs(current_y - exact_y)
        return exact_y, error

    def _table_line(self, current_x, current_y, i):
        if self.exact_ys is None:
            print(f'{"%0.8f" % current_x:<16} | {"%0.8f" % current_y:<16}')
        else:
            error = abs(current_y - self.exact_ys[i])
            print(f'{"%0.8f" % current_x:<16} | {"%0.8f" % current_y:<16} | {"%0.8f" % self.exact_ys[i]:<16} | '
                  f'{"%0.8f" % error:<16}')

    def get_xy(self, initial_x, initial_y, last_x):
        # Initializing vector of xs with values from initial_x until last_x + step
        xs = np.arange(initial_x, last_x + self.step_h, self.step_h)
        ys = []
        current_y = initial_y
        self._description_line(initial_x, initial_y, last_x)
        self._table_header()
        for i, current_x in enumerate(xs):
            self._table_line(current_x, current_y, i)
            ys.append(current_y)
            current_y = self._get_next_y(current_x, current_y)
        print(f'--- End of {self.method_name} ---')
        return xs, ys


class Exact(DifferentialEquationsSolver):
    def __init__(self, function, step_h, exact_ys=None):
        super().__init__(function, step_h, exact_ys)
        self.method_name = "Exact"


class Euler(DifferentialEquationsSolver):
    def __init__(self, function, step_h, exact_ys=None):
        super().__init__(function, step_h, exact_ys)
        self.method_name = "Euler"

    def _get_next_y(self, current_x, current_y):
        return current_y + self.step_h * self.function(current_x, current_y)


class RungeKutta(DifferentialEquationsSolver):
    def __init__(self, function, step_h, exact_ys=None):
        super().__init__(function, step_h, exact_ys)
        self.method_name = "Runge Kutta"

    def _get_next_y(self, current_x, current_y):
        # Using Simpson's integration formula
        adjusted_step_h = self.step_h / 2
        f1_left_point = self.function(
            current_x, current_y)
        # Calculating with Euler method
        f2_central_point_euler = self.function(
            current_x + adjusted_step_h, current_y + adjusted_step_h * f1_left_point)
        # Correction
        f3_central_point_adjusted = self.function(
            current_x + adjusted_step_h, current_y + adjusted_step_h * f2_central_point_euler)
        # Prediction to the full step
        f4_right_point = self.function(
            current_x + self.step_h, current_y + self.step_h * f3_central_point_adjusted)
        return current_y + (self.step_h / 6) * (
            f1_left_point + 2 * f2_central_point_euler + 2 * f3_central_point_adjusted + f4_right_point
        )

