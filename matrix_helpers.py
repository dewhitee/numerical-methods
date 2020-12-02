import numpy as np
import math
import matrix_cond as cond

def get_matrixA(matrix: list) -> list:
    return [row[:-1] for row in matrix]

def get_matrixAX(matrix: list, vectorX: list):
    return np.multiply(get_matrixA(matrix), vectorX)

def show_A_B(matrix: list):
    print("\nA", '{:>50}'.format("B"))
    for (a, b) in zip(get_matrixA(matrix), get_vectorB_unpacked(matrix)):
        print(*['{:10}'.format("%0.4f" % elem) for elem in a], " | " + '{:>2}'.format("%0.4f" % b), sep="\t")

def show_AX_B(matrix: list, vectorX: list):
    print(f'\n{"AX":<16} {"B":<16}')
    matrix_AX = get_matrixAX(matrix, vectorX)
    for (ax, b) in zip(sum_matrix_to_vector(matrix_AX), get_vectorB_unpacked(matrix)):
        print(f'{"%0.4f" % ax:<16} {"%0.4f" % b:<16}')

def get_vectorB(matrix: list) -> list:
    return [row[-1:] for row in matrix]

def get_vectorB_unpacked(matrix: list) -> list:
    return [row[0] for row in get_vectorB(matrix)]

def get_vector_norm_euc(vec: list):
    result = float()
    for elem in vec:
        result += (elem ** 2)
    return math.sqrt(result)

def get_vector_norm_man(vec: list):
    result = float()
    for elem in vec:
        elem += abs(elem)
    return result

def sum_matrix_to_vector(matrix: list) -> list:
    return [sum(row) for row in matrix]

def append_vectorB(matrix: list, vectorB: list) -> list:
    matrixA = get_matrixA(matrix)
    for i, elem in enumerate(vectorB):
        matrixA[i].append(elem)
    return matrixA

def check_condition_number(matrix: list, vectorX: list, deltaX: list) -> bool:
    condA = cond.matrix_cond(get_matrixA(matrix))
    vectorB_norm = get_vector_norm_euc(get_vectorB_unpacked(matrix))
    deltaB_norm = get_vector_norm_euc(get_deltaB(matrix))
    vectorX_norm = get_vector_norm_euc(vectorX)
    deltaX_norm = get_vector_norm_euc(deltaX)
    print("Checking Cond(A) for matrix:")
    print(condA, ">= (", deltaX_norm, "/", vectorX_norm, ") * (", vectorB_norm, "/", deltaB_norm, ")")
    condition = condA >= (deltaX_norm / vectorX_norm) * (vectorB_norm / deltaB_norm)
    print("Condition returned =", condition)
    return condition

def get_deltaB(matrix: list) -> list:
    vectorB = get_vectorB_unpacked(matrix)
    deltaB = list()
    largest = max(vectorB)
    print("Largest value of vectorB =", largest)
    for elem in vectorB:
        deltaB.append(elem * 0.01 if elem == largest else 0)
    return deltaB

def get_deltaX(vectorX: list, vector_deltaX: list) -> list:
    return vectorX - vector_deltaX

def full_print(matrix: list, vars: list, title: str):
    print("\n"+title+" matrix:")
    for (var, row) in zip(vars, matrix):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')

def A_print(matrix: list, vars: list):
    print("\nA matrix:")
    matrixA = [row[:-1] for row in matrix]
    for (var, row) in zip(vars, matrixA):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')


def B_print(matrix: list, vars: list):
    print("\nB vector:")
    vectorB = [row[-1:] for row in matrix]
    for (var, row) in zip(vars, vectorB):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')


def X_print(vectorX: list, vars: list):
    print("\nX vector:")
    for (var, val) in zip(vars, vectorX):
        print(var, '= %0.4f' % (val))
