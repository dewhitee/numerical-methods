import numpy as np

def least_squares(matrix: list):
    
    
    n = len(matrix)

    # Initializing list of a arguments and the list of approximation functions with degree of k
    A = list()
    F = list()
    k = 3

    # List of x arguments
    X = list()
    Y = list()

    # Linear superposition as approximation function
    def approximation_function(x: float, k: int):
        result = float()

        # f[m] - system of the basis functions
        # k - is the approximation order
        for m in range(0, k):
            result += A[m] * F[m](x)
        return result

    result = float()
    for i in range(0, n):
        result += (approximation_function(X[i], k) - Y[i]) ** 2

def get_power_basis_matrix(current_k: int, vectorX: list, vectorY: list):
    """ Function must return the matrix of equations
    algorithm reference:
    https://e.tsi.lv/pluginfile.php/130692/mod_resource/content/2/%D0%A7%D0%B8%D1%81%D0%BB%D0%B5%D0%BD%D0%BD%D1%8B%D0%B5%20%D0%BC%D0%B5%D1%82%D0%BE%D0%B4%D1%8B_LECTURES-2017.pdf
    
    matrix construction formula:
    https://studopedia.ru/9_211036_stepennoy-bazis.html
    """
    variables_count = len(vectorX)
    out_matrix = np.ndarray(shape=(current_k + 1, current_k + 1))
    out_b_coefficients = np.ndarray(shape=(current_k + 1, 1))
    
    # Initializing first element (in the upper-left corner)
    #out_matrix[0, 0] = variables_count
    # Actually, this element can also be found in a loop if we sum the x^0 or y^0 elements
    # Iterating starting from the second element, as the first is set to variables_count already
    for i in range(0, current_k + 1):
        # Iterating over the elements of a row
        for j in range(0, current_k + 1):
            # Calculating each element of the matrix A
            out_matrix[i, j] = sum([x ** (j + i) for x in vectorX])

            # Calculating each B coefficient
            out_b_coefficients[i, 0] = sum([(vectorX[index] ** i) * vectorY[index] for index in range(0, len(vectorX))])
        
    return out_matrix, out_b_coefficients

    

