import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matrix_helpers as mh
import lab1_gauss_elimination.gauss_elimination as ge
import roots
import scipy

class LeastSquaresApproximator:
    def __init__(self, vectorX: list, vectorY: list, k_approx_order: int = 2, ftype: str = "auto", makeplot=False, customfunc=None, resolution=10,
    print_matrix=False, customstep=None):
        """ Calculates least squares for the choosen ftype. Set makeplot to true to construct the plt.plot for the splines

        customfunc: Can be specified to use instead of default autofunc

        customstep: Can be specified to be used instead of the current_step ( (next X - current X) / resolution )

        !Returns:
            vectorF,
            vector_deltaF
        """

        print("\nLeast Squares -------------\n")
        # Initialize class fields
        self.k_approx_order = k_approx_order
        self.resolution = resolution
        self.vectorX = vectorX
        self.vectorY = vectorY
        self.customfunc = customfunc
        self.customstep = customstep

        # Set n to the size of known X vector
        self.n = len(vectorX)

        # Getting matrix and b coefficients vector from the power basis
        self.matrixA, self.vectorB = self.get_power_basis_matrix()

        # Making the one whole matrix from the matrixA and vectorB to pass into the Gauss Elimination solving function
        self.matrixAB = mh.append_vectorB_to_matrixA(
            self.matrixA.tolist(), mh.unpack_vector(self.vectorB.tolist()))

        # Print the newly created power basis matrix if necessary
        if print_matrix:
            mh.full_print(
                matrixAB=self.matrixAB,
                vars=['x'+str(i) for i in range(0, self.n)],
                title="Power basis matrix with k = " + str(k_approx_order))

        # Solving matrix with Gauss Elimination, getting the X vector of solutions for power basis matrix
        self.solution_vectorX = ge.GaussElimination(
            matrixAB=self.matrixAB,
            vars=['x'+str(i) for i in range(0, self.n)],
            matrix_name="",
            print_only_results=True).solution_vectorX

        print("solution_vectorX (a coefficients) =", self.solution_vectorX)
        print("vectorX[0] =",self.vectorX[0])

        # Initializing vectorF and vector_deltaF as empty lists
        self.vectorF = []
        self.vector_deltaF = []

        # Then we need to calculate the approximation vector from the solution vector
        # using the choosen type of approximation function (auto by default)
        if ftype == "linear_deprecated":
            print("Using linear approximation basis function")
            def linfunc(a, x): 
                print("a * x + sum(solutions) = ", a, "*", x, "+",
                sum(self.solution_vectorX[:-1]), "=", a * x + sum(self.solution_vectorX[:-1]))
                return a * x + sum(self.solution_vectorX[:-1])
            self.vectorF = [linfunc(self.solution_vectorX[-1], self.vectorX[i]) for i in range(0, self.n)]
            self.vector_deltaF = [(self.vectorY[i] - self.vectorF[i]) ** 2 for i in range(0, self.n)]

        elif ftype == "exponential":
            exit

        elif ftype == "custom":
            print("Using custom approximation basis function")
            self.vectorF = [customfunc(self.solution_vectorX, self.vectorX[i]) for i in range(0, self.n)]
            self.vector_deltaF = [(self.vectorY[i] - self.vectorF[i]) ** 2 for i in range(0, self.n)]

        elif ftype in ("auto", "", None):
            print("Using auto (default) approximation basis function")
            self.vectorF = [self.autofunc(self.solution_vectorX, self.vectorX[i]) for i in range(0, self.n)]
            self.vector_deltaF = [(self.vectorY[i] - self.vectorF[i]) ** 2 for i in range(0, self.n)]

        else:
            raise ValueError("You need to choose ftype!")

        print("Sum of vector_deltaF =", sum(self.vector_deltaF))
        print("vectorF:\n", self.vectorF, "\nvector_deltaF:\n", self.vector_deltaF)
        print("\n--------------------------- End Least Squares")

        # Getting interpolated vectors of X and Y to build smooth curve (count of points depends on the resolution parameter)
        self.interpolated_vectorX, self.interpolated_vectorY = self.get_interpolated_xy_vectors()

        # Constructing the plot
        if makeplot:
            self.make_plot()

    def autofunc(self, solvec, x):
        return sum([solvec[i] * (x ** i) for i in range(0, len(solvec))])

    def get_x_from_y_estimated(self, y, lower_border = -0.5, upper_border = 0.5, max_iterations=500, tolerance=0.0001):
        """ Finds the estimated x by the Bisection Method
        """
        return roots.RootFinder(
            self.autofunc if self.customfunc is None else self.customfunc,
            solution_vectorX=self.solution_vectorX
        ).bisection(
            lower=y + lower_border,
            upper=y + upper_border,
            max_iterations=max_iterations,
            tolerance=tolerance
        )

    def get_x_from_y_closest(self, y, vectorX=None, vectorY=None):
        """ Returns the closest_x value with the index
        """
        # Get the closest known x from the vectorX and the corresponding y from the vectorY
        closest_list = []
        for v in self.interpolated_vectorY if vectorY is None else vectorY:
            closest_list.append(abs(y - v))
        closest_y_index = closest_list.index(min(closest_list))
        closest_x = self.interpolated_vectorX[closest_y_index] if vectorX is None else vectorX[closest_y_index]
        return closest_x, closest_y_index

    def get_closest_boundary_points(self, y):
        closest_x = self.get_x_from_y_closest(y)
        index = 0
        for i in range(0, len(self.vectorX)):
            if self.vectorX[i] > closest_x[0]:
                index = i
                break
        return self.vectorX[index - 1], self.vectorF[index - 1], self.vectorX[index], self.vectorF[index]

    def get_x_from_y_interpolated(self, y, resolution = 100):
        # Get the closest known x from the reinterpolated vector X and the corresponding y from the reinterpolated vector Y
        xvec, yvec = self.get_interpolated_xy_vectors(resolution)
        new_closest = self.get_x_from_y_closest(y, xvec, yvec)
        return xvec[new_closest[1]]
    
    def interpolate(self, currentX, deltaX, resolution):
        interpolated_vectorX = list()
        interpolated_vectorY = list()

        step = deltaX / resolution
        for i in range(0, resolution):
            new_y = self.autofunc(self.solution_vectorX, currentX) if self.customfunc is None else self.customfunc(self.solution_vectorX, currentX)
            interpolated_vectorY.append(new_y)
            interpolated_vectorX.append(currentX)
            currentX += step

        return interpolated_vectorX, interpolated_vectorY

    def get_y_from_x(self, x):
        if self.customfunc is None:
            return self.autofunc(self.solution_vectorX, x)
        else:
            return self.customfunc(self.solution_vectorX, x)

    def get_interpolated_xy_vectors(self, custom_resolution=None):
        """ Creates the interpolated lists of X and Y with points count specified by resolution parameter.
        """
        # Initializing empty out vectors
        out_vectorX = list()
        out_vectorY = list()

        resolution = self.resolution if custom_resolution is None else custom_resolution

        # Getting values interpolated between x[i] and x[i+1] points
        for i in range(0, len(self.vectorX) - 1):
            current_step = (self.vectorX[i + 1] - self.vectorX[i]) / resolution if self.customstep is None else self.customstep
            current_x = self.vectorX[i]
            # Adding each point to the out vector
            for j in range(0, resolution):
                # Calculating interpolated Y point and adding it to the out Y vector
                out_vectorY.append(self.autofunc(self.solution_vectorX, current_x) 
                    if self.customfunc is None else self.customfunc(self.solution_vectorX, current_x))
                # Adding the current x (which is previous X + current step) to the out X vector
                out_vectorX.append(current_x)
                # Adding the current step value to the current_x to get the next x value
                current_x += current_step

        # Adding the last point X and Y coords to the out_vectors
        out_vectorX.append(self.vectorX[-1])
        out_vectorY.append(self.vectorF[-1])

        return out_vectorX, out_vectorY

    def make_plot(self):
        plt.figure("Least Squares by dewhitee")
        plt.title("Least Squares approximation with k = " + str(self.k_approx_order))
        plt.plot(self.vectorX, self.vectorY, 'bs', self.vectorX, self.vectorF, 'g--', self.vectorX, self.vectorF, 'g^')
        plt.plot(self.interpolated_vectorX, self.interpolated_vectorY, 'yo-')
        plt.xlabel("X values")
        plt.ylabel("Y values")
        for i in range(0, len(self.vectorX) - 1):
            plt.plot([self.vectorX[i], self.vectorX[i]], [self.vectorY[i], self.vectorF[i]], 'r--')
        plt.plot([self.vectorX[-1], self.vectorX[-1]], [self.vectorY[-1], self.vectorF[-1]], 'r--')
        plt.show()

    def get_power_basis_matrix(self):
        """ Function must return the matrix of equations.

        If current_k is >= 10, then Condition Number of the output matrix will be infinity.

        algorithm reference:
        https://e.tsi.lv/pluginfile.php/130692/mod_resource/content/2/%D0%A7%D0%B8%D1%81%D0%BB%D0%B5%D0%BD%D0%BD%D1%8B%D0%B5%20%D0%BC%D0%B5%D1%82%D0%BE%D0%B4%D1%8B_LECTURES-2017.pdf

        matrix construction formula:
        https://studopedia.ru/9_211036_stepennoy-bazis.html

        !Returns:
            (Wandermonde like) matrixA, vectorB 
        """
        current_k = self.k_approx_order

        variables_count = len(self.vectorX)
        out_matrix = np.ndarray(shape=(current_k + 1, current_k + 1))
        out_b_coefficients = np.ndarray(shape=(current_k + 1, 1))

        # Filling the out_matrix
        for i in range(0, current_k + 1):
            # Iterating over the elements of a row
            for j in range(0, current_k + 1):
                # Calculating each element of the matrix A
                out_matrix[i, j] = sum([x ** (j + i) for x in self.vectorX])

                # Calculating each B coefficient
                out_b_coefficients[i, 0] = sum([(self.vectorX[index] ** i) * self.vectorY[index] for index in range(0, len(self.vectorX))])

        return out_matrix, out_b_coefficients

    

