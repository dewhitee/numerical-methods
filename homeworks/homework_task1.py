import lab1_gauss_elimination.gauss_elimination as ge

matrix = [
    [18, -3, 4, 3],
    [-3, 4, -1, 7],
    [34, 4, -19, 0]
]

ge.GaussElimination(
    matrixAB=matrix,
    vars=['x1', 'x2', 'x3'],
    print_only_results=False,
    matrix_name="First task matrix"
)
