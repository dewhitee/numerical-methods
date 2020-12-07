import numpy as np
import math
import lab1_gauss_elimination.gauss_elimination as ge
import lab1_gauss_seidel.gauss_seidel as gs
from numpy import array
import summary

def get_matrix_determinant(matrixA: list) -> float:
    return np.linalg.det(matrixA)

def get_matrix_inversed(matrixA: list) -> list:
    return np.linalg.inv(matrixA)

def get_matrix_cond(matrixA: list) -> float:
    """ Returns the condition number of the matrix.
    Calculated as ||A|| * ||A^-1|| or ||matrix_norm|| * ||inversed_matrix_norm||
    """
    return get_matrix_norm(matrixA) * get_matrix_norm(get_matrix_inversed(matrixA))

def get_matrixA(matrixAB: list) -> list:
    """ Returns the matrixA from the whole matrix
    """
    return [row[:-1] for row in matrixAB]

def get_matrixAX(matrixA: list, vectorX: list):
    """ Returns the matrixA multiplied by the answers vectorX
    """
    return np.multiply(get_matrixA(matrixA), vectorX)

def show_A_B(matrixAB: list):
    print("\nA", f'{"B":>{len(matrixAB) * 16 + 1}}')
    for (a, b) in zip(get_matrixA(matrixAB), get_vectorB_unpacked(matrixAB)):
        print(*['{:10}'.format("%0.4f" % elem) for elem in a], "| " + '{:>2}'.format("%0.4f" % b), sep="\t")

def show_AX_B(matrixAB: list, vectorX: list):
    print(f'\n{"AX":<16} {"B":<16}')
    matrix_AX = get_matrixAX(matrixAB, vectorX)
    for (ax, b) in zip(sum_matrix_to_vector(matrix_AX), get_vectorB_unpacked(matrixAB)):
        print(f'{"%0.4f" % ax:<16} {"%0.4f" % b:<16}')

def get_vectorB(matrixAB: list) -> list:
    """ Returns the vectorB from the whole matrix as the new matrix with each element as list
    as: [[1], [5], [3.6]]
    """
    return [row[-1:] for row in matrixAB]

def get_vectorB_unpacked(matrixAB: list) -> list:
    """ Returns the vectorB from the whole matrix as: [1, 5, 3.6]
    """
    return [row[0] for row in get_vectorB(matrixAB)]

def unpack_vector(v: list) -> list:
    """
    """
    return [row[0] for row in v]

def get_whole_matrix(matrixA: np.ndarray, vectorB: np.ndarray):
    return append_vectorB_to_matrixA(matrixA.tolist(), unpack_vector(vectorB.tolist()))

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

def get_matrix_norm(matrixA: list) -> float:
    """
    returns || matrix ||_1 = max || matrix-col ||
    """
    result = float()
    colsum = float()
    for i in range(0, len(matrixA[0])):
        temp = float()
        #for row in matrixA:
        #    #print(row[i])
        #    temp += abs(row[i])
        #result = max(result, temp)
        #result = max(result, get_vector_norm_euc(matrixA[i]))
        for row in matrixA:
            print(row[i])
            temp += abs(row[i])
            colsum += (row[i] ** 2)
        result = max(result, temp)
        colsum += temp
    #print("Max value of matrixA columns was:", result, "\nSquare root of a ** 2 was:", math.sqrt(colsum))
    return result


def get_matrix_norm_inversed(matrixA: list) -> float:
    """
    returns || matrix ||_infinity = max || matrix-row ||
    """
    result = float()
    for row in matrixA:
        result = max(result, sum(np.abs(row)))
    return result

def sum_matrix_to_vector(any_matrix: list) -> list:
    """ Function sums each row of the matrix and returns the vector of them
    """
    return [sum(row) for row in any_matrix]

def append_vectorB_to_matrixA(matrixA: list, vectorB: list) -> list:
    """ This function gets the matrixA, and returns the matrixA 
    with appended specified vectorB
    """
    for i, elem in enumerate(vectorB):
        matrixA[i].append(elem)
    return matrixA

def append_vectorB_to_whole_matrix(matrixAB: list, vectorB: list) -> list:
    """ This function gets the whole matrix, finds the matrixA and
    returns the matrixA with the specified vectorB as the NEW whole matrix
    """
    matrixA = get_matrixA(matrixAB)
    for i, elem in enumerate(vectorB):
        matrixA[i].append(elem)
    return matrixA

def get_deltaB_from_vectorB(vectorB: list):
    """ DeltaB is the B vector with artificially added error.
    Is computed by finding the max element in B vector and multipling it by 0.01 (1%)
    All other elements of the B vector will be null, meaning that there will be the only
    one variable with the none-zero value in the whole out B vector.
    """
    deltaB = list()
    largest = max(vectorB)
    #print("Largest value of vectorB =", largest)
    for elem in vectorB:
        deltaB.append(elem * 0.01 if elem == largest else 0)
    return deltaB

def get_deltaB_from_whole_matrix(matrixAB: list) -> list:
    """ DeltaB is the B vector with artificially added error.
    Is computed by finding the max element in B vector and multipling it by 0.01 (1%)
    All other elements of the B vector will be null, meaning that there will be the only
    one variable with the none-zero value in the whole out B vector.
    """
    vectorB = get_vectorB_unpacked(matrixAB)
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

def full_print(matrixAB: list, vars: list, title: str):
    """ Prints the whole matrix with the custom title
    """
    print("\n"+title+" matrix:")
    if vars is not None:
        for (var, row) in zip(vars, matrixAB):
            print(var, *["%0.4f" % elem for elem in row], sep='\t  ')
    else:
        for row in matrixAB:
            print(*["%0.4f" % elem for elem in row], sep='\t  ')

def A_print(matrix: list, vars: list):
    """ Prints the matrixA of the specified whole matrix
    """
    print("\nA matrix:")
    matrixA = [row[:-1] for row in matrix]
    for (var, row) in zip(vars, matrixA):
        print(var, *["%0.4f" % elem for elem in row], sep='\t  ')


def B_print(matrixAB: list, vars: list):
    """ Prints the coefficient vectorB of the specified whole matrix
    """
    print("\nB vector:")
    vectorB = [row[-1:] for row in matrixAB]
    if vars is not None:
        for (var, row) in zip(vars, vectorB):
            print(var, *["%0.4f" % elem for elem in row], sep='\t  ')
    else:
        for row in vectorB:
            print(*["%0.4f" % elem for elem in row], sep='\t  ')

def X_print(vectorX: list, vars: list):
    """ Prints the answers vectorX
    """
    print("\nX vector:")
    for (var, val) in zip(vars, vectorX):
        print(var, '= %0.4f' % (val))


def solve_with_gauss_elimination(matrixAB: list, vars: list, print_only_results: bool = False, matrix_name: str = ""):
    """ Solves the specified matrix using the Gauss Elimination
    """
    matrix_vectorX = ge.GaussElimination(
        matrixAB=matrixAB, 
        vars=vars,
        print_only_results=print_only_results,
        matrix_name=matrix_name).solution_vectorX

    adjusted_matrix = get_adjusted_matrix(get_matrixA(matrixAB), get_vectorB_unpacked(matrixAB))
    summary.Summary(
        matrixAB=matrixAB, 
        adjusted_matrix=adjusted_matrix, 
        vars=vars, 
        X=matrix_vectorX,
        print_only_results=print_only_results,
        matrix_name=matrix_name).compare_with_adjusted_gauss_elimination(
            printall=True
        )

def solve_with_gauss_seidel(matrixAB: list, equations: list, vars: list, e: float, print_only_results: bool = False, 
matrix_name: str = "", use_matrix: bool = False):
    """ Solve the specified equations using the Gauss Seidel
    """
    matrix_vectorX = gs.GaussSeidel(
        matrixAB=matrixAB,
        equations=equations, 
        vars=vars, 
        e=e, 
        print_only_results=print_only_results, 
        matrix_name=matrix_name).solution_vectorX

    adjusted_matrix = get_adjusted_matrix(get_matrixA(matrixAB), get_vectorB_unpacked(matrixAB))
    summary.Summary(
        matrixAB=matrixAB,
        adjusted_matrix=adjusted_matrix,
        vars=vars,
        X=matrix_vectorX,
        print_only_results=print_only_results,
        matrix_name=matrix_name).compare_with_adjusted_gauss_seidel(
            equations=equations,
            e=e,
            use_matrix=True,
            printall=True
        )

def get_adjusted_matrix(matrixA, vectorB):
    return append_vectorB_to_matrixA(
        matrixA, array(vectorB) + array(get_deltaB_from_vectorB(vectorB)))
