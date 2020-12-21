import numpy as np


class DifferentialEquationsSolver:
    def __init__(self, function, step_h, exact_ys=None):
        """
        :param function: Two parameter function
        """
        self.function = function
        self.step_h = step_h
        self.method_name = ""
        self.exact_ys = exact_ys

    def _get_next_y(self, current_x, current_y):
        return self.function(current_x, current_y)

    def _description_line(self, initial_x, initial_y, last_x):
        print(f'\n--- {self.method_name} ---')
        print("h (step) =", self.step_h, ", initial x =", initial_x, ", initial y =", initial_y, ", last x =", last_x)

    def _table_header(self):
        if self.exact_ys is None:
            print(f'{"x":<12} | {"y":<12}')
            print('{:-<30}'.format(""))
        else:
            print(f'{"x":<12} | {"y":<12} | {"Exact value":<12}')
            print('{:-<40}'.format(""))

    def _table_line(self, current_x, current_y, i):
        if self.exact_ys is None:
            print(f'{current_x:<12} | {"%0.8f" % current_y:<12}')
        else:
            print(f'{current_x:<12} | {"%0.8f" % current_y:<12} | {self.exact_ys[i]}')

    def get_xy(self, initial_x, initial_y, last_x):
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
        adjusted_step_h = self.step_h / 2
        f1_left_point = self.function(
            current_x, current_y)
        f2_central_point_euler = self.function(
            current_x + adjusted_step_h, current_y + adjusted_step_h * f1_left_point)
        f3_central_point_adjusted = self.function(
            current_x + adjusted_step_h, current_y + adjusted_step_h * f2_central_point_euler)
        f4_right_point = self.function(
            current_x + self.step_h, current_y + self.step_h * f3_central_point_adjusted)
        return current_y + (self.step_h / 6) * (
            f1_left_point + 2 * f2_central_point_euler + 2 * f3_central_point_adjusted + f4_right_point
        )

