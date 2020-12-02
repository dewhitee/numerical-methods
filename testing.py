import lab1_gauss_elimination.gauss_elimination as ge
import matrix_helpers as mh
import numpy as np
from numpy import array

#print(ge.test_matrix_2)

#test_matrix_2_without_b = [row[:-1] for row in ge.test_matrix_2]

#print(test_matrix_2_without_b)
#print("Condition number (Cond(A)) of matrix: ", cond.matrix_cond(test_matrix_2_without_b))

matrix_0 = [
    [1, 2],
    [2, 3.999]
]
#print(matrix_0)

matrix_0_inversed = mh.get_matrix_inversed(matrix_0)
#print(matrix_0_inversed)

#print("matrix_0_norm = ", norm.matrix_norm(matrix_0))
#print("matrix_0_norm_inversed = ",norm.matrix_norm(matrix_0_inversed))
#print("matrix_0_cond = ", cond.matrix_cond(matrix_0))

# -----------------------------------------------------

matrix_1_vars = ['x1', 'x2', 'x3']
matrix_1 = [
    [25, 5, 1, 106.8],
    [64, 8, 1, 177.2],
    [144, 12, 1, 279.2],
]

matrix_1_vectorX = ge.gauss_elimination(matrix_1, matrix_1_vars)

mh.A_print(matrix_1, matrix_1_vars)
mh.B_print(matrix_1, ['b1', 'b2', 'b3'])
mh.X_print(matrix_1_vectorX, matrix_1_vars)

mh.show_A_B(matrix_1)

#print("X = ", matrix_1_X)
#print("A = ", mh.get_matrixA(matrix_1))

# -- AX = B
matrix_1_AX = mh.get_matrixAX(matrix_1, matrix_1_vectorX)
matrix_1_B = mh.get_vectorB_unpacked(matrix_1)
#print("AX = B;\n", mh.sum_matrix_to_vector(matrix_1_AX), "=", matrix_1_B) # == true

mh.show_AX_B(matrix_1, matrix_1_vectorX)

# Log(Cond(A)) == number of decimals of accuracy that are lost
# -- Theoretical condition number:
matrix_1_cond = mh.get_matrix_cond(mh.get_matrixA(matrix_1))
print("\nmatrix_1_cond =", matrix_1_cond)

# -- Experimental condition number:
# delta B -> nulling all elements of vector B except the abs(max element of B)
matrix_1_deltaB = mh.get_deltaB(matrix_1)

print("\ndelta B vector:", matrix_1_deltaB)

matrix_1_modified = mh.append_vectorB(matrix_1, array(matrix_1_B) + array(matrix_1_deltaB))
mh.full_print(matrix_1_modified, matrix_1_vars, 'Modified matrix_1')

matrix_1_deltaX = ge.gauss_elimination(matrix_1_modified, matrix_1_vars)

print("delta X vector = ", matrix_1_deltaX)

mh.check_condition_number(matrix_1, matrix_1_vectorX, matrix_1_deltaX)
