import copy


# Algorithm
# 1. Start
# 2. Arrange given system of linear equations in diagonally dominant form
# 3. Read tolerable error (e)
# 4. Convert the first equation in terms of first variable, second eq in terms of second var and so on.
# 5. Set initial guesses for x0, y0, z0 and so on
# 6. Substitute value of y0, z0 ... from step 5 in first equation obtained from step 4 to calculate
# new value of x1. Use x1, z0, u0, ... in second equation obtained from step 4 to calculate new value of y1.
# Similarly, use x1, y1, u0 ... to find new z1 and so on.
# 7. If | x0 - x1 | > e and | y0 - y1 | > e and | z0 - z1 | > e and so on then goto step 9
# 8. Set x0 = x1, y0 = y1, z0 = z1 and so on and goto step 6
# 9. Print value of x1, y1, z1 and so on
# 10. Stop
class GaussSeidel:
    def __init__(self, equations: list, vars: list, e: float, print_only_results: bool = False, matrix_name: str = "",
                 matrixAB: list = None,
                 print_results=True, without_print=False, max_iterations=1000, show_errors_list=False,
                 auto_adjust_matrix=True):
        """ 
        equations -- list of lambda equations with any count of arguments.
        Example: [lambda x, y: x + y, lambda x, y: x - y]

        vars -- list of symbols to specify the names of the equations variables.
        Example: ['x', 'y', 'z']

        e -- tolerable error
        Example: 0.001

        matrix -- the whole matrix (AB)

        Note: you can also pass the matrix instead of equations

        returns list of solution variables X
        """
        if not without_print:
            print('\n-------------------------------------Gauss Seidel - ' + matrix_name)
            print('Using equations list...') if matrixAB is None else print('Using matrix...')
            print('With accuracy of ' + str(e) + '\n')
            if not print_only_results:
                print('Iter', *[f'{v:<12}' for v in vars], sep="\t")

        condition = True

        # Initializing all variables values with zero if using equations
        # or set each X as B if using matrix, as the first approximation
        vars_values = [0] * len(equations) if matrixAB is None else [row[-1] for row in matrixAB]

        if not without_print and not print_only_results:
            print(0, *[f'{"%0.4f":<12}' % elem for elem in vars_values], sep="\t")

        iteration = 1

        # Adjusting matrix (swapping rows if necessary to put all largest elements on main diagonal)
        adjusted_matrixAB = None
        if auto_adjust_matrix:
            adjusted_matrixAB = self.adjust_matrix(matrixAB)
            if not without_print:
                print("(Auto adjustment of matrix is enabled)\nAdjusted matrix:\n", adjusted_matrixAB)
        else:
            adjusted_matrixAB = matrixAB

        # Main loop
        while condition:
            e_list = []

            # Calculating all variables on the current iteration
            if adjusted_matrixAB is None:
                self.iterate_equations(equations, e_list, vars_values)
            else:
                self.iterate_matrix(adjusted_matrixAB, e_list, vars_values)

            if not without_print and not print_only_results:
                print(iteration, *[f'{"%0.4f":<12}' % elem for elem in vars_values], sep="\t")
            iteration += 1

            if show_errors_list:
                print("Errors:", *[f'{"%0.4f":<12}' % elem for elem in e_list], sep="\t")

            # Checking if all current errors are greater than required error e
            condition = self.check_error_rate(e_list, e)

            if iteration > max_iterations:
                break

        if not without_print:
            print('\nSolution:')
            for (var, val) in zip(vars, vars_values):
                print(var, '= %0.3f' % (val))
            print('-------------------------------------\n')

        self.vars_values = vars_values
        self.solution_vectorX = vars_values

    def iterate_equations(self, equations, e_list, vars_values):
        # Calculating all variables
        for i, eq in enumerate(equations):
            # Calculating the i-th lambda of equations list
            new_value = eq(*vars_values)

            # Adding i-th error to the e_list
            e_list.append(abs(vars_values[i] - new_value))

            # Set current i-th vars_values variable to it's newly calculated new_value
            vars_values[i] = new_value

    def iterate_matrix(self, matrix, e_list, vars_values):
        # Calculating all variables
        for i, row in enumerate(matrix):
            # Calculating the i-th X variable
            # Initializing current X with B coefficient value
            new_value = row[-1]
            for j in range(0, len(row) - 1):
                if i != j:
                    # Subtracting the non-diagonal values multiplied by calculated previously values
                    new_value -= (row[j] * vars_values[j])

            # Adding the value on a diagonal of the matrix
            new_value *= 1 / row[i]

            # Adding i-th error to the e_list
            e_list.append(abs(vars_values[i] - new_value))

            # Set current i-th vars_values variable to it's newly calculated new_value
            vars_values[i] = new_value

    def check_error_rate(self, e_list: list, e):
        return all([current_e > e for current_e in e_list])

    def adjust_matrix(self, matrixAB):
        """ Returns the matrixAB with all rows placed in such way,
        that main diagonal has largest values by module
        """
        if matrixAB is None:
            return None

        out_matrix = copy.deepcopy(matrixAB)
        for i, row in enumerate(out_matrix):
            current_max = 0
            current_max_index = 0

            # Iterating over the i-th column values
            for j in range(i, len(out_matrix)):
                # print(abs(current_max), "<", abs(out_matrix[j][i]))
                if abs(current_max) < abs(out_matrix[j][i]):
                    current_max = out_matrix[j][i]
                    current_max_index = j

            # print("Current max list is: ",
            #      out_matrix[current_max_index], "with index", current_max_index)

            # Swap max row with the current row
            out_matrix[i], out_matrix[current_max_index] = out_matrix[current_max_index], out_matrix[i]

        return out_matrix


# Testing

test_equations = [
    lambda x, y, z: (17 - y + 2 * z) / 20,
    lambda x, y, z: (-18 - 3 * x + z) / 20,
    lambda x, y, z: (25 - 2 * x + 3 * y) / 20
]

# gauss_seidel(test_equations, ['x','y','z'], 0.001)

test_equations_2 = [
    lambda x1, x2, x3, x4: 1 / 20.9 * (21.70 - 1.2 * x2 - 2.1 * x3 - 0.9 * x4),
    lambda x1, x2, x3, x4: 1 / 21.2 * (27.46 - 1.2 * x1 - 1.5 * x3 - 2.5 * x4),
    lambda x1, x2, x3, x4: 1 / 19.8 * (28.76 - 2.1 * x1 - 1.5 * x2 - 1.3 * x4),
    lambda x1, x2, x3, x4: 1 / 32.1 * (49.72 - 0.9 * x1 - 2.5 * x2 - 1.3 * x3)
]

# gauss_seidel(test_equations_2, ['x1', 'x2', 'x3', 'x4'], 0.001)


# 1, 3, 4 var systems testing

var_1_equations = [
    lambda x1, x2, x3, x4: 1 / (2 - x1 - 2 * x2 - x3),
    lambda x1, x2, x3, x4: 1 / 3 * (-1 * x1 - 5 * x2 - x3),
    lambda x1, x2, x3, x4: 1 / (-7) * (-2 * x1 + x3)
]

# gauss_seidel(var_1_equations, ['x1', 'x2', 'x3', 'x4'], 0.001)
