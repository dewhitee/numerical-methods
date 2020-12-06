import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matrix_helpers as mh
import lab1_gauss_elimination.gauss_elimination as ge

def least_squares(vectorX: list, vectorY: list, k_approx_order: int = 2, ftype: str = "linear", makeplot=False, customfunc=None, resolution=10):
    """ Calculates least squares for the choosen ftype

    !Returns:
        vectorF,
        vector_deltaF
    """
    
    print("\nLeast Squares -------------\n")

    n = len(vectorX)

    # Initializing list of a arguments and the list of approximation functions with degree of k
    A = list()
    F = list()
    k = 3

    # List of x arguments
    X = list()
    Y = list()

    # Getting matrix and b coefficients vector from the power basis
    matrixA, vectorB = get_power_basis_matrix(k_approx_order, vectorX, vectorY)

    # Making the one whole matrix from the matrixA and vectorB to pass into the Gauss Elimination solving function
    whole_matrix = mh.only_append_vectorB(
        matrixA.tolist(), mh.unpack_vector(vectorB.tolist()))

    # Solving matrix with Gauss Elimination, getting the X vector of solutions for power basis matrix
    solution_vectorX = ge.gauss_elimination(
        matrix=whole_matrix, 
        vars=['x'+str(i) for i in range(0, n)],
        matrix_name="")

    print("solution_vectorX = ", solution_vectorX)
    print("vectorX[0]=",vectorX[0])

    # Initializing vectorF and vector_deltaF as empty lists
    vectorF = []
    vector_deltaF = []

    # Then we need to calculate the approximation vector from the solution vector
    # using the choosen type of approximation function (linear by default)
    if ftype == "linear":
        print("Using linear approximation basis function")
        def linfunc(a, x): 
            print("a * x + sum(solutions) = ", a, "*", x, "+",
                  sum(solution_vectorX[:-1]), "=", a * x + sum(solution_vectorX[:-1]))
            return a * x + sum(solution_vectorX[:-1])
        vectorF = [linfunc(solution_vectorX[-1], vectorX[i]) for i in range(0, n)]
        vector_deltaF = [(vectorY[i] - vectorF[i]) ** 2 for i in range(0, n)]

    elif ftype == "exponential":
        exit

    elif ftype == "custom":
        print("Using custom approximation basis function")
        # Uses the customfunc lambda to approximate given data
        vectorF = [customfunc(solution_vectorX, vectorX[i]) for i in range(0, n)]
        vector_deltaF = [(vectorY[i] - vectorF[i]) ** 2 for i in range(0, n)]

    elif ftype == "auto":
        print("Using auto approximation basis function")
        vectorF = [autofunc(solution_vectorX, vectorX[i]) for i in range(0, n)]
        vector_deltaF = [(vectorY[i] - vectorF[i]) ** 2 for i in range(0, n)]

    print("Sum of vector_deltaF =", sum(vector_deltaF))
    print("vectorF:\n", vectorF, "\nvector_deltaF:\n", vector_deltaF)
    print("\n--------------------------- End Least Squares")

    interpolated_vectorX, interpolated_vectorY = get_interpolated_xy_vectors(vectorX, vectorF, solution_vectorX, resolution)

    make_plot(vectorX, vectorY, vectorF, interpolated_vectorX, interpolated_vectorY, makeplot)

    return vectorF, vector_deltaF

    # Linear superposition as approximation function
    #def approximation_function(x: float, k: int):
    #    result = float()
    #    # f[m] - system of the basis functions
    #    # k - is the approximation order
    #    for m in range(0, k):
    #        result += A[m] * F[m](x)
    #    return result

    #result = float()
    #for i in range(0, n):
    #    result += (approximation_function(X[i], k) - Y[i]) ** 2

def autofunc(solvec, x):
    return sum([solvec[i] * (x ** i) for i in range(0, len(solvec))])

def get_interpolated_xy_vectors(vectorX, vectorF, solution_vectorX, resolution=10):
    """ Creates the interpolated lists of X and Y with points count specified by resolution parameter.
    Set makeplot to true to construct the plt.plot for the splines
    """
    out_vectorX = list()
    out_vectorY = list()

    for i in range(0, len(vectorX) - 1):
        current_step = (vectorX[i + 1] - vectorX[i]) / resolution
        current_x = vectorX[i]
        for j in range(0, resolution):
            out_vectorY.append(autofunc(solution_vectorX, current_x))
            out_vectorX.append(current_x)
            current_x += current_step

    # Adding the last point X and Y coords to the out_vectors
    out_vectorX.append(vectorX[-1])
    out_vectorY.append(vectorF[-1])

    return out_vectorX, out_vectorY

def make_plot(vectorX, vectorY, vectorF, interpolated_vectorX, interpolated_vectorY, makeplot):
    if makeplot:
        plt.plot(vectorX, vectorY, 'bs', vectorX, vectorF, 'g--', vectorX, vectorF, 'g^')
        plt.plot(interpolated_vectorX, interpolated_vectorY, 'y-')
        plt.xlabel("X values")
        plt.ylabel("Y values")
        for i in range(0, len(vectorX) - 1):
            plt.plot([vectorX[i], vectorX[i]], [vectorY[i], vectorF[i]], 'r--')

        plt.show()

def get_power_basis_matrix(current_k: int, vectorX: list, vectorY: list):
    """ Function must return the matrix of equations.

    If current_k is >= 10, then Condition Number of the output matrix will be infinity.

    algorithm reference:
    https://e.tsi.lv/pluginfile.php/130692/mod_resource/content/2/%D0%A7%D0%B8%D1%81%D0%BB%D0%B5%D0%BD%D0%BD%D1%8B%D0%B5%20%D0%BC%D0%B5%D1%82%D0%BE%D0%B4%D1%8B_LECTURES-2017.pdf
    
    matrix construction formula:
    https://studopedia.ru/9_211036_stepennoy-bazis.html
    """
    variables_count = len(vectorX)
    out_matrix = np.ndarray(shape=(current_k + 1, current_k + 1))
    out_b_coefficients = np.ndarray(shape=(current_k + 1, 1))
    
    # Filling the out_matrix
    for i in range(0, current_k + 1):
        # Iterating over the elements of a row
        for j in range(0, current_k + 1):
            # Calculating each element of the matrix A
            out_matrix[i, j] = sum([x ** (j + i) for x in vectorX])

            # Calculating each B coefficient
            out_b_coefficients[i, 0] = sum([(vectorX[index] ** i) * vectorY[index] for index in range(0, len(vectorX))])
        
    return out_matrix, out_b_coefficients

    

