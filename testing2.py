import matrix_helpers as mh
import numpy as np

test_1_equations = [
    lambda x1, x2, x3, x4: 1/20.9*(21.70 - 1.2*x2 - 2.1*x3 - 0.9*x4),
    lambda x1, x2, x3, x4: 1/21.2*(27.46 - 1.2*x1 - 1.5*x3 - 2.5*x4),
    lambda x1, x2, x3, x4: 1/19.8*(28.76 - 2.1*x1 - 1.5*x2 - 1.3*x4),
    lambda x1, x2, x3, x4: 1/32.1*(49.72 - 0.9*x1 - 2.5*x2 - 1.3*x3)
]

test_1_matrix = [
    [20.9, 1.2, 2.1, 0.9, 21.70],
    [1.2, 21.2, 1.5, 2.5, 27.46],
    [2.1, 1.5, 19.8, 1.3, 28.76],
    [0.9, 2.5, 1.3, 32.1, 49.72]
]

test_1_vars = ['x1', 'x2', 'x3', 'x4']

# -----
mh.solve_with_gauss_elimination(test_1_matrix, test_1_vars, True, "Initial matrix")

# -----
mh.solve_with_gauss_seidel(test_1_matrix, test_1_equations, test_1_vars, 0.0001, True, "Matrix with B + deltaB")

#print(mh.get_matrix_determinant(mh.get_matrixA(test_1_matrix)))