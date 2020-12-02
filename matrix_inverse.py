import numpy as np

def matrix_inverse(matrix: list) -> list:
    return np.linalg.inv(matrix)

test_inverse_matrix = [
    [1, 3, 3],
    [1, 4, 3],
    [1, 3, 4]
]

#print(matrix_inverse(test_matrix))
