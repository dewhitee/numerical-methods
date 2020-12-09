import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import least_squares as ls
import matrix_helpers as mh

#vectorX = [2, 3, 5]
#vectorY = [1, 0, 4]

#matrix, bcoef = ls.get_power_basis_matrix(2, [2, 3, 5], [1, 0, 4])

#print(matrix)
#print(bcoef)

#print("determinant of the matrix: ", mh.get_matrix_determinant(matrix))

#whole_matrix_0 = mh.get_whole_matrix(matrix, bcoef)
#vectorF, vector_deltaF = ls.least_squares(vectorX, vectorY, k_approx_order=1)

vectorX2 = [1, 2, 3, 4, 5]
vectorY2 = [5.3, 6.3, 4.8, 3.8, 3.3]
# matrix2, bcoef2 = ls.get_power_basis_matrix(1, vectorX2, vectorY2)
# print("determinant of matrix2: ", mh.get_matrix_determinant(matrix2))
# whole_matrix = mh.only_append_vectorB(matrix2.tolist(), mh.unpack_vector(bcoef2.tolist()))
# vars = ['x1', 'x2']
# mh.full_print(whole_matrix, vars=vars, title='Matrix2')
# mh.solve_with_gauss_elimination(
#     matrix=whole_matrix, 
#     vars=vars, 
#     matrix_name="Matrix2")
ls.LeastSquaresApproximator(
    vectorX=vectorX2,
    vectorY=vectorY2,
    k_approx_order=1,
    ftype="linear",
    makeplot=True)

ls.LeastSquaresApproximator(
    vectorX=vectorX2,
    vectorY=vectorY2,
    k_approx_order=1,
    ftype="auto",
    makeplot=True
)

# print("vectorF:\n", vectorF, "\nvector_deltaF:\n", vector_deltaF)


# plt.plot(vectorX2, vectorY2, 'bs', vectorX2, vectorF2, 'b--')
# plt.xlabel("X values")
# plt.ylabel("Y values")
# plt.show()
