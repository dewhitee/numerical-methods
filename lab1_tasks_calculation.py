import lab1_gauss_elimination.gauss_elimination as ge
import lab1_gauss_seidel.gauss_seidel as gs
import condition_number as cn
import matrix_helpers as mh
import timeit

# 1) For both elimination and seidel
matrix_1 = [
    [1, -2, 1, 2],
    [2, -5, -1, -1],
    [-7, 0, 1, -2]
]

ge.GaussElimination(
    matrixAB=matrix_1,
    vars=['x1','x2','x3'],
    print_only_results=True,
    matrix_name="First matrix (Elimination)"
)

gs.GaussSeidel(
    matrixAB=matrix_1,
    equations=None,
    vars=['x1','x2','x3'],
    e=0.001,
    print_only_results=True,
    matrix_name="First matrix (Seidel)"
)

# 2) For elimination only, test the speed with both elimination and seidel ---------------------
matrix_2 = [
    [5, 5, -3, 4, -11],
    [1, -4, 6, -4, -10],
    [-2, -5, 4, -5, -12],
    [-3, -3, 5, -5, 8]
]

matrix_2_updated = [
    [5, 5, -3, 4, -11],
    [-2, -5, 4, -5, -12],
    [1, -4, 6, -4, -10],
    [-3, -3, 5, -5, 8]
]

def timer_elimination():
    ge.GaussElimination(
        matrixAB=matrix_2,
        vars=['x1', 'x2', 'x3', 'x4'],
        print_only_results=True,
        matrix_name="Second matrix (Elimination)",
        print_results=False,
        without_print=True
    )

def timer_seidel():
    gs.GaussSeidel(
        matrixAB=matrix_2,
        equations=None,
        vars=['x1', 'x2', 'x3', 'x4'],
        e=0.00001,
        print_only_results=True,
        matrix_name="Second matrix (Seidel)",
        max_iterations=1000,
        show_errors_list=False,
        auto_adjust_matrix=True,
        without_print=True
    )

def timer_seidel2():
    gs.GaussSeidel(
        matrixAB=matrix_2_updated,
        equations=None,
        vars=['x1', 'x2', 'x3', 'x4'],
        e=0.00001,
        print_only_results=True,
        matrix_name="Second matrix (Seidel)",
        max_iterations=1000,
        show_errors_list=False,
        auto_adjust_matrix=False,
        without_print=True
    )

def timer_seidel3():
    gs.GaussSeidel(
        matrixAB=matrix_2_updated,
        equations=None,
        vars=['x1', 'x2', 'x3', 'x4'],
        e=0.0001,
        print_only_results=True,
        matrix_name="Second matrix (Seidel)",
        max_iterations=1000,
        show_errors_list=False,
        auto_adjust_matrix=False,
        without_print=True
    )

print("Gauss Elimination - time elapsed\n(number of calls = 10000):\n", timeit.timeit("timer_elimination()", 
    setup="from __main__ import timer_elimination", number=10000), 'seconds \n')
print("Gauss Seidel with auto-adjustment of matrix and e = 0.00001 - time elapsed\n (number of calls = 10000):\n", timeit.timeit('timer_seidel()',
    setup="from __main__ import timer_seidel", number=10000), 'seconds\n')
print("Gauss Seidel without auto-adjustment of matrix and e = 0.00001 - time elapsed\n (number of calls = 10000):\n", timeit.timeit('timer_seidel2()',
    setup="from __main__ import timer_seidel2", number=10000), 'seconds\n')
print("Gauss Seidel without auto-adjustment of matrix and e = 0.0001 - time elapsed\n (number of calls = 10000):\n", timeit.timeit('timer_seidel3()',
    setup="from __main__ import timer_seidel3", number=10000), 'seconds\n')


ge.GaussElimination(
    matrixAB=matrix_2,
    vars=['x1', 'x2', 'x3', 'x4'],
    print_only_results=True,
    matrix_name="Second matrix (Elimination)",
    without_print=False
)

gs.GaussSeidel(
    matrixAB=matrix_2,
    equations=None,
    vars=['x1', 'x2', 'x3', 'x4'],
    e=0.00001,
    print_only_results=True,
    matrix_name="Second matrix (Seidel)",
    max_iterations=1000,
    show_errors_list=False,
    auto_adjust_matrix=True
)

# 3) For seidel only, check the condition number, find how big are errors of calculation -----
matrix_3 = [
    [2, -1, -1, 5],
    [1, 3, -2, 7],
    [1, 2, 3, 10]
]

gs.GaussSeidel(
    matrixAB=matrix_3,
    equations=None,
    vars=['x1', 'x2', 'x3'],
    e=0.001,
    print_only_results=True,
    matrix_name="Third matrix (Seidel)"
)

cn.ConditionNumber(
    matrixAB=matrix_3,
    print_result=True
)

# 4) For seidel only ------------------------------------------------------------------------
matrix_4 = [
    [8, 5, 3, 30],
    [-2, 8, 1, 15],
    [1, 3, -10, 42]
]

ge.GaussElimination(
    matrixAB=matrix_4,
    vars=['x1', 'x2', 'x3'],
    print_only_results=True,
    matrix_name="Fourth matrix (Elimination)"
)

gs.GaussSeidel(
    matrixAB=matrix_4,
    equations=None,
    vars=['x1', 'x2', 'x3'],
    e=0.001,
    print_only_results=True,
    matrix_name="Fourth matrix (Seidel)"
)

# For gauss only, check the condition number, find how big are erros of calculation --------
matrix_5 = [
    [0.78, 0.563, 0.217],
    [0.913, 0.659, 0.254]
]

ge.GaussElimination(
    matrixAB=matrix_5,
    vars=['x1', 'x2'],
    print_only_results=True,
    matrix_name="Fifth matrix (Elimination)"
)

cn.ConditionNumber(
    matrixAB=matrix_5,
    print_result=True
)
