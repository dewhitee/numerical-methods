import copy

# Algorithm
# Straight step:
# Adjusting matrix by finding the pivot element for each row
# Nulling bottom-left corner
#
# Reverse step:
# Getting the solution vector X (Calculating from bottom to up)
# ...
class GaussElimination:
    def __init__(self, matrixAB: list, vars: list, print_only_results: bool = False, matrix_name: str = "", matrixA=None, vectorB=None, print_results=True, without_print=False):
        """
        returns list of solution variables X
        """

        self.matrix_len = len(matrixAB)
        matrix_copy = copy.deepcopy(matrixAB)

        def print_matrix():
            for (var, row) in zip(vars, matrix_copy):
                print(var, *["%0.4f" % elem for elem in row], sep='\t')

        if not without_print:
            print('\n-------------------------------------Gauss Elimination - ' + matrix_name)
            print('Length of the matrix (n):', self.matrix_len)
            if not print_only_results:
                print('Matrix before pivotisation (initial):\n')
                print_matrix()

        # Find the pivot element - we need to put it as the first row in the matrix
        # This will impove calculation accuracy and make division by zero impossible
        for i in range(self.matrix_len):
            # Iterating as columns starting from the i + 1 element, pushing larger rows up
            for k in range(i + 1, self.matrix_len):
                # Checking if module of main diagonal element is less than module of k-element of i-th column
                if abs(matrix_copy[i][i]) < abs(matrix_copy[k][i]):
                    # Swapping each element of current i-th row with each element of larger k-th row
                    for j in range(0, self.matrix_len + 1):
                        # Swapping elements (current row element with larger k-th row element)
                        matrix_copy[i][j], matrix_copy[k][j] = matrix_copy[k][j], matrix_copy[i][j]

        if not without_print and not print_only_results:
            print('\nMatrix after pivotisation:\n')
            print_matrix()

        # Main Gauss Elimination loop
        # Forward elimination -- Straight step (Nulling the bottom-left corner)
        # Making Upper-Triangle (TU) matrix from the initial matrix
        for i in range(self.matrix_len - 1):
            # Iterating as columns
            for k in range(i + 1, self.matrix_len):
                # Dividing k-th row i-th element by the current main diagonal element
                coefficient = matrix_copy[k][i] / matrix_copy[i][i]  # Coefficient
                print(f'Coefficient ({coefficient}) =", {matrix_copy[k][i]} / {matrix_copy[i][i]}')
        
                # Make the elements below the pivot elements equal to zero or eliminate the variables
                # Iterating all elements of rows (including the vectorB)
                for j in range(0, self.matrix_len + 1):
                    print(f'{j}-th element of {k}-th row ({matrix_copy[k][j]}) -= {coefficient} * {matrix_copy[i][j]}')
                    matrix_copy[k][j] -= coefficient * matrix_copy[i][j]
                    print(f'{j}-th element of {k}-th row = {matrix_copy[k][j]}')

        if not without_print and not print_only_results:
            print('\nMatrix after gauss elimination:\n')
            print_matrix()

        # List of variables values (x, y, z, ...)
        # Initializing all variables values with zero
        vars_values = [0] * self.matrix_len

        # Back substitution -- Reversed step
        # Finding the solution_vectorX values
        for i in range(self.matrix_len - 1, -1, -1):
            # Make the variable (Xi) to be calculated equal to the vectorB[i] from the last equation
            vars_values[i] = matrix_copy[i][self.matrix_len]

            for j in range(i + 1, self.matrix_len):
                # Subtracting all values of matrixA except the main diagonal i-th value
                # whose value is being calculated
                if j != i:
                    vars_values[i] -= matrix_copy[i][j] * vars_values[j]

            # Dividing the Xi by the main diagonal i-th value
            vars_values[i] /= matrix_copy[i][i]

        if not without_print and print_results:
            print('\nSolution:')
            for (var, val) in zip(vars, vars_values):
                print(var, '= %0.4f' % val)
            print('-------------------------------------\n')

        self.vars_values = vars_values
        self.matrix_copy = matrix_copy
        self.solution_vectorX = vars_values

test_matrix = [
    [25, 5, 1, 106.8],
    [64, 8, 1, 177.2],
    [144, 12, 1, 279.2],
]
#gauss_elimination(test_matrix, ['x1', 'x2', 'x3'])
test_matrix_2 = [
    [2, 3, -1, 7],
    [1, -1, 6, 14],
    [6, -2, 1, 11]
]
#gauss_elimination(test_matrix_2, ['x1', 'x2', 'x3'])
