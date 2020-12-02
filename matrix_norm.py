import numpy as np

#
#
#
#

def matrix_norm(matrix: list) -> float:
    """
    returns || matrix ||_1 = max || matrix-col ||
    """
    result = float()
    for i in range(0, len(matrix)):
        temp = float()
        for row in matrix:
            #print(row[i])
            temp += abs(row[i])

        result = max(result, temp)
    #print("Matrix norm = ", result)
    return result


def matrix_norm_inversed(matrix: list) -> float:
    """
    returns || matrix ||_infinity = max || matrix-row ||
    """
    result = float()
    for row in matrix:
        result = max(result, sum(np.abs(row)))

    #print("Matrix norm inversed = ", result)
    return result

test_norm_matrix = [
    [-3, 5, 7],
    [2, 6, 4],
    [0, 2, 8]
]

#matrix_norm(test_matrix)
#matrix_norm_inversed(test_matrix)
