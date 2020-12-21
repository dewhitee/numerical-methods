import numpy as np


class RungeKutta:
    def __init__(self, function, step_h):
        """
        :param function: Two parameter function
        """
        self.function = function
        self.step_h = step_h

    def _get_next_y(self, current_x, current_y):
        return current_y + self.step_h * self.function(current_x, current_y)

    def get_xy(self, initial_x, initial_y, last_x):
        xs = np.arange(initial_x, last_x + self.step_h, self.step_h)
        ys = []
        current_y = initial_y
        print("\n--- Euler ---")
        print("h (step) =", self.step_h, ", initial x =", initial_x, ", initial y =", initial_y, ", last x =", last_x)
        print(f'{"x":<12} | {"y":<12}')
        print('{:-<30}'.format(""))
        for current_x in xs:
            print(f'{current_x:<12} | {current_y:<12}')
            ys.append(current_y)
            current_y = self._get_next_y(current_x, current_y)
        print("--- End of Euler ---\n")
        return xs, ys
