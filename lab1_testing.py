import lab1_gauss_elimination.gauss_elimination as ge
import matrix_helpers as mh
import condition_number as cn
import summary

matrix_1_vars = ['x1', 'x2', 'x3']
matrix_1 = [
    [25, 5, 1, 106.8],
    [64, 8, 1, 177.2],
    [144, 12, 1, 279.2],
]

matrix_1_vectorX = ge.GaussElimination(
    matrix_1, matrix_1_vars).solution_vectorX

adjusted_matrix_1 = mh.get_adjusted_matrix(
    mh.get_matrixA(matrix_1), mh.get_vectorB_unpacked(matrix_1))

matrix_1_deltaX = ge.GaussElimination(
    matrixAB=adjusted_matrix_1, 
    vars=matrix_1_vars).solution_vectorX

cn.ConditionNumber(
    matrixA=mh.get_matrixA(matrix_1),
    vectorB=mh.get_vectorB_unpacked(matrix_1)
).experimental(matrix_1_vectorX, matrix_1_deltaX, printall=True)

#summary.Summary(matrix_1, adjusted_matrix_1, matrix_1_vars, matrix_1_vectorX)