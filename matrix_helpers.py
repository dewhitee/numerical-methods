import numpy as np

def get_matrixA(matrix: list) -> list:
    return [row[:-1] for row in matrix]

def get_matrixAX(matrix: list, vectorX: list):
    return np.multiply(get_matrixA(matrix), vectorX)

def show_A_B(matrix: list):
    print("A", "B", sep="\t\t\t\t\t\t")
    for (a, b) in zip(get_matrixA(matrix), get_vectorB_unpacked(matrix)):
        print(["%0.4f" % elem for elem in a], "%0.4f" % b, sep='\t\t\t')

def show_AX_B(matrix: list, vectorX: list):
    print(f'{"AX":<16} {"B":<16}')
    matrix_AX = get_matrixAX(matrix, vectorX)
    for (ax, b) in zip(sum_matrix_to_vector(matrix_AX), get_vectorB_unpacked(matrix)):
        print(f'{"%0.4f" % ax:<16} {"%0.4f" % b:<16}')

def get_vectorB(matrix: list) -> list:
    return [row[-1:] for row in matrix]

def get_vectorB_unpacked(matrix: list) -> list:
    return [row[0] for row in get_vectorB(matrix)]

def sum_matrix_to_vector(matrix: list) -> list:
    return [sum(row) for row in matrix]

def get_deltaB(matrix: list) -> list:
    vectorB = get_vectorB_unpacked(matrix)
    deltaB = list()
    largest = max(vectorB)
    print("largest value of vectorB =", largest)
    for elem in vectorB:
        deltaB.append(elem * 0.01 if elem == largest else 0)
    return deltaB

def get_deltaX(matrix: list) -> list:
    exit

def A_print(matrix: list, vars: list):
    print("\nA matrix:\n")
    matrixA = [row[:-1] for row in matrix]
    for (var, row) in zip(vars, matrixA):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')


def B_print(matrix: list, vars: list):
    print("\nB vector:\n")
    vectorB = [row[-1:] for row in matrix]
    for (var, row) in zip(vars, vectorB):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')


def X_print(vectorX: list, vars: list):
    print("\nX vector:\n")
    for (var, val) in zip(vars, vectorX):
        print(var, '= %0.4f' % (val))
