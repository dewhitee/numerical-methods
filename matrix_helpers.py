import numpy as np
import math
import lab1_gauss_elimination.gauss_elimination as ge
import lab1_gauss_seidel.gauss_seidel as gs
from numpy import array

def get_matrix_determinant(matrix: list) -> float:
    return np.linalg.det(matrix)

def get_matrix_inversed(matrix: list) -> list:
    return np.linalg.inv(matrix)

def get_matrix_cond(matrix: list) -> float:
    """ Returns the condition number of the matrix.
    Calculated as ||A|| * ||A^-1|| or ||matrix_norm|| * ||inversed_matrix_norm||
    """
    return get_matrix_norm(matrix) * get_matrix_norm(get_matrix_inversed(matrix))

def get_matrixA(matrix: list) -> list:
    """ Returns the matrixA from the whole matrix
    """
    return [row[:-1] for row in matrix]

def get_matrixAX(matrix: list, vectorX: list):
    """ Returns the matrixA multiplied by the answers vectorX
    """
    return np.multiply(get_matrixA(matrix), vectorX)

def show_A_B(matrix: list):
    print("\nA", f'{"B":>{len(matrix) * 16 + 1}}')
    for (a, b) in zip(get_matrixA(matrix), get_vectorB_unpacked(matrix)):
        print(*['{:10}'.format("%0.4f" % elem) for elem in a], "| " + '{:>2}'.format("%0.4f" % b), sep="\t")

def show_AX_B(matrix: list, vectorX: list):
    print(f'\n{"AX":<16} {"B":<16}')
    matrix_AX = get_matrixAX(matrix, vectorX)
    for (ax, b) in zip(sum_matrix_to_vector(matrix_AX), get_vectorB_unpacked(matrix)):
        print(f'{"%0.4f" % ax:<16} {"%0.4f" % b:<16}')

def get_vectorB(matrix: list) -> list:
    """ Returns the vectorB from the whole matrix as the new matrix with each element as list
    as: [[1], [5], [3.6]]
    """
    return [row[-1:] for row in matrix]

def get_vectorB_unpacked(matrix: list) -> list:
    """ Returns the vectorB from the whole matrix as: [1, 5, 3.6]
    """
    return [row[0] for row in get_vectorB(matrix)]

def get_vector_norm_euc(vec: list) -> float:
    """ Calculates the norm of the vector by the Euclidean formula.
    """
    result = float()
    for elem in vec:
        result += (elem ** 2)
    return math.sqrt(result)

def get_vector_norm_man(vec: list) -> float:
    """ Calculates the norm of the vector by the Manhattan formula.
    """
    result = float()
    for elem in vec:
        elem += abs(elem)
    return result


def get_matrix_norm(matrix: list) -> float:
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
    return result


def get_matrix_norm_inversed(matrix: list) -> float:
    """
    returns || matrix ||_infinity = max || matrix-row ||
    """
    result = float()
    for row in matrix:
        result = max(result, sum(np.abs(row)))
    return result

def sum_matrix_to_vector(matrix: list) -> list:
    """ Function sums each row of the matrix and returns the vector of them
    """
    return [sum(row) for row in matrix]

def append_vectorB(matrix: list, vectorB: list) -> list:
    """ This function gets the whole matrix, finds the matrixA and
    returns the matrixA with the specified vectorB as the NEW whole matrix
    """
    matrixA = get_matrixA(matrix)
    for i, elem in enumerate(vectorB):
        matrixA[i].append(elem)
    return matrixA

def check_condition_number(matrix: list, vectorX: list, deltaX: list) -> bool:
    """ Finds the condition number of the matrix as well as norm of the answers vector X  
    and delta vector X, norm of the coefficient vector B and delta vector B.
    Checks condition: cond(A) >= (deltaX_norm / vectorX_norm) * (vectorB_norm / deltaB_norm)
    """
    condA = get_matrix_cond(get_matrixA(matrix))
    vectorB_norm = get_vector_norm_euc(get_vectorB_unpacked(matrix))
    deltaB_norm = get_vector_norm_euc(get_deltaB(matrix))
    vectorX_norm = get_vector_norm_euc(vectorX)
    deltaX_norm = get_vector_norm_euc(deltaX)
    print("Matrix determinant =", get_matrix_determinant(get_matrixA(matrix)))
    print("Checking Cond(A) for matrix:")
    print(condA, ">= (", deltaX_norm, "/", vectorX_norm, ") * (", vectorB_norm, "/", deltaB_norm, ")")
    print(condA, ">=", (deltaX_norm / vectorX_norm) * (vectorB_norm / deltaB_norm))
    condition = condA >= (deltaX_norm / vectorX_norm) * (vectorB_norm / deltaB_norm)
    print("Condition returned", condition)
    return condition

def get_deltaB(matrix: list) -> list:
    """ DeltaB is the B vector with artificially added error.
    Is computed by finding the max element in B vector and multipling it by 0.01 (1%)
    All other elements of the B vector will be null, meaning that there will be the only
    one variable with the none-zero value in the whole out B vector.
    """
    vectorB = get_vectorB_unpacked(matrix)
    deltaB = list()
    largest = max(vectorB)
    print("Largest value of vectorB =", largest)
    for elem in vectorB:
        deltaB.append(elem * 0.01 if elem == largest else 0)
    return deltaB

def get_deltaX(vectorX: list, vector_deltaX: list) -> list:
    """ DeltaX is the difference vector between X (vector of answers) and modified X
    (vector of answers got from calculating the same system but with B + deltaB instead of B)
    """
    return vectorX - vector_deltaX

def full_print(matrix: list, vars: list, title: str):
    """ Prints the whole matrix with the custom title
    """
    print("\n"+title+" matrix:")
    for (var, row) in zip(vars, matrix):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')

def A_print(matrix: list, vars: list):
    """ Prints the matrixA of the specified whole matrix
    """
    print("\nA matrix:")
    matrixA = [row[:-1] for row in matrix]
    for (var, row) in zip(vars, matrixA):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')


def B_print(matrix: list, vars: list):
    """ Prints the coefficient vectorB of the specified whole matrix
    """
    print("\nB vector:")
    vectorB = [row[-1:] for row in matrix]
    for (var, row) in zip(vars, vectorB):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')

def X_print(vectorX: list, vars: list):
    """ Prints the answers vectorX
    """
    print("\nX vector:")
    for (var, val) in zip(vars, vectorX):
        print(var, '= %0.4f' % (val))


def solve_with_gauss_elimination(matrix: list, vars: list):
    """ Solves the specified matrix using the Gauss Elimination
    """
    vectorX = ge.gauss_elimination(matrix, vars)
    summary(matrix, vars, vectorX)

def solve_with_gauss_seidel(matrix: list, equations: list, vars: list, e: float):
    """ Solve the specified equations using the Gauss Seidel
    """
    vectorX = gs.gauss_seidel(equations, vars, e)
    summary(matrix, vars, vectorX)

def summary(matrix: list, vars: list, vectorX: list):
    """ Summarizes the data available after solving the matrix
    """
    show_A_B(matrix)
    AX = get_matrixAX(matrix, vectorX)
    B = get_vectorB_unpacked(matrix)
    show_AX_B(matrix, vectorX)
    condA = get_matrix_cond(get_matrixA(matrix))
    print("\nCondition number Cond(A) of this matrix is =", condA)
    deltaB = get_deltaB(matrix)
    print("\nDeltaB vector:", deltaB)
    modified_matrix = append_vectorB(matrix, array(B) + array(deltaB))
    full_print(modified_matrix, vars, 'Modified matrix')
    deltaX = ge.gauss_elimination(modified_matrix, vars)
    print("DeltaX vector = ", deltaX)
    check_condition_number(matrix, vectorX, deltaX)
