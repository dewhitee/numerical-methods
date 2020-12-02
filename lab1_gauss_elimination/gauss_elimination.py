# Algorithm
#
#
#
#
# ...

def gauss_elimination(matrix: list, vars: list) -> list:
    """
    returns list of solution variables X
    """

    matrix_len = len(matrix)
    print('Length of the matrix (n):', matrix_len)

    print('Matrix before pivotisation (initial):\n')
    for (var, row) in zip(vars, matrix):
        print(var, *row, sep='\t')

    # Find the pivot element - we need to put it as the first row in the matrix
    for i in range(matrix_len):
        for k in range(i + 1, matrix_len):
            if abs(matrix[i][i]) < abs(matrix[k][i]):
                for j in range(0, matrix_len + 1):
                    # Swapping elements
                    matrix[i][j], matrix[k][j] = matrix[k][j], matrix[i][j]

    print('\nMatrix after pivotisation:\n')
    for (var, row) in zip(vars, matrix):
        print(var, *row, sep='\t')

    # Main Gauss Elimination loop
    # Forward elimination -- Straight step (Nulling the bottom-left corner)
    for i in range(matrix_len - 1):
        for k in range(i + 1, matrix_len):
            coefficient = matrix[k][i] / matrix[i][i] # Coefficient
    
            # Make the elements below the pivot elements equal to zero 
            # or eliminate the variables
            for j in range(0, matrix_len + 1):
                matrix[k][j] -= coefficient * matrix[i][j]

    print('\nMatrix after gauss elimination:\n')
    for (var, row) in zip(vars, matrix):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')

    # List of variables values (x, y, z, ...)
    # Initializing all variables values with zero
    vars_values = [0] * matrix_len

    # Back substitution -- Reversed step (Nulling upper-right corner)
    for i in range(matrix_len - 1, -1, -1):
        # Make the variable to be calculated equal to the rhs of the
        # last equation
        vars_values[i] = matrix[i][matrix_len]

        for j in range(i + 1, matrix_len):
            # Subtracting all the lhs values except the coefficient 
            # of the variable whose value is being calculated
            if j != i:
                vars_values[i] -= matrix[i][j] * vars_values[j]

        # Finally, divide the rhs by the coefficient of the variable
        # to be calculated
        vars_values[i] /= matrix[i][i]

    print('\nSolution:')
    for (var, val) in zip(vars, vars_values):
        print(var, '= %0.4f' %(val))

    return vars_values


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
