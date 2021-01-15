import lab1_gauss_elimination.gauss_elimination as ge

ng = 33
ns = 14
#ns = 15  # eugene
j = 1

a11 = (ng + 4) + j*5
a12 = -3 - j*4
a13 = 4 - j*4
b1 = 3 + j*6

a21 = -3 + j*2
a22 = 8 + j*(10 - ns)
a23 = 1 + j*2
b2 = 1 - j*(ns - 20)

a31 = j*(ng + 1)
a32 = ns - 10
a33 = ns - j * ng
b3 = j*10


matrix = [
    [a11, a12, a13, b1],
    [a21, a22, a23, b2],
    [a31, a32, a33, b3]
]

ge.GaussElimination(
    matrixAB=matrix,
    vars=['x1', 'x2', 'x3'],
    print_only_results=False,
    matrix_name="First task matrix"
)

matrix_2 = [
    [42, -7, 0, 9],
    [-1, -1, 3, 2],
    [34, 9, -14, 10]
]

ge.GaussElimination(
    matrixAB=matrix_2,
    vars=['x1', 'x2', 'x3'],
    print_only_results=False,
    matrix_name="First task Dmitry matrix"
)
