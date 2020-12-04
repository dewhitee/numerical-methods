import numpy as np
import copy
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lab1_gauss_elimination.gauss_elimination as ge
import lab1_gauss_seidel.gauss_seidel as gs
import matrix_helpers as mh

class CubicSplineInterpolator:

    def __init__(self, known_vectorX: list, known_vectorY: list, known_points: list, vars = []):
        """ 
        known_points -- is the list of pairs of the points (nodes) that are known at the start

        Each polynomial segment of a cubic spline has degree 3
        Each cubic polynomial has 4 unknown coefficients (A, B, C, D)

            a) S(x) is a cubic polynomial, denoted Sj(x), on the subinterval [xj, xj+1] for each j = 0, 1, ..., n - 1
            b) Sj(xj) = f(xj) and Sj(x(j+1)) for each j = 0, 1, ..., n - 1
            c) Sj+1(x(j+1)) = Sj(x(j+1)) for each j = 0, 1, ..., n - 2
            d) S'j+1(x(j+1)) = S'j(x(j+1)) for each j = 0, 1, ..., n - 2
            e) S''j+1(x(j+1)) = S''j(x(j+1)) for each j = 0, 1, ..., n - 2
            f) One of the following sets of boundary conditions is satisfied:
                (I)     S''(x0) = S''(xn) = 0                   ( natural or free boundary )
                (II)    S'(x0) = f'(x0) and S'(xn) = f'(xn)     ( clamped boundary )

        Every Si(X) equals to = ai + bi(x - xi) + ci(x - xi)^2 + di(x - xi)^3

            S0(x) = a0 + b0(x - x0) + ci(x - x0)^2 + d0(x - x0)^3
            S1(x) = a1 + b1(x - x1) + ci(x - x1)^2 + d1(x - x1)^3
            S2(x) = a2 + b2(x - x2) + ci(x - x2)^2 + d1(x - x2)^3
            and so on...

        Matching condition:

            S0(x1) = S1(x1)         --- joint points of the splines must be equal.
            S0'(x1) = S1'(x1)       --- first derivative of the joint points also must be equal.
            S0''(x1) = S1''(x1)     --- as well as second derivative.
            ------------------------------------------------------------------------------------------------
            
            S0(x0) = y0
            S0(x1) = y1

        Boundary condition:

            For example, if we have four points: x0, x1, x2, x3 and, therefore, 3 splines: S0, S1, S2
            then the bounding condition will be as such:

            For nature splines:
                S0''(x0) = 0        -- as the starting point
                S2''(x3) = 0        -- as the end point

            For clamped:
                S0'(x0) = f'(x0)
                S2'(x3) = f'(x3)

        Example:
        
        if we have list of points like so: [[1, 2], [2, 3], [3, 5]] then we can get the first spline S0(X) as:

            First spline:
            1) [1, 2]: S0(x) = a0 + b0(x - 1) + c0(x - 1)^2 + d0(x - 1)^3           - where 1 - is the X of the first point in the list
            
            Second spline:
            2) [2, 3]: S1(x) = a1 + b1(x - 2) + c0(x - 2)^2 + d0(x - 2)^3           - where 2 - is the X of the second point in the list

            --- Matching conditions ---
            S0 of 1 will be 2, according to the first point:            S0(1) = 2; 
                => a0 = 2
            S0 of 2 will be 3, as the y of the second point is 3:       S0(2) = 3; a0 + b0 + c0 + d0:
                => 2 + b0(2 - 1) + c0(2 - 1) + d0(2 - 1) = 1
                => a0 = 3
                => ___first equation___

            S1 of 2 will be 3:      a1 = a1 + b1(2 - 2) + c1(2 - 2) + d1(2 - 2) = 3
            S1 of 3 = 5:            3 + b1 + c1 + d1 = 5
                => b1 + c1 + d1 = 2
                => ___second equation___

            S0'(x) = b0 + 2c0(x - 1) + 3d0(x - 1)^2
                => S0''(x) = 2c0 + 6d0(x-1)
            S1'(x) = b1 + 2c1(x - 2) + 3d1(x - 2)^2
                => S1''(x) = 2c1 + 6d1(x-2)

            S0'(2) = S1'(2) => b0 + 2c0 + 3d0 = b1
                => ___third equation___

            S0''(2) = S1''(2) => 2c0 + 6d0 = 2c1
                => ___fourth equation___
            
            --- Boundary conditions ---
            S0''(1) = 0         -- starting point (first X)
            S1''(3) = 0         -- end point      (last X)
                => c0 = 0
                => 2c1 + 6d1 = 0
                => ___5 equation___

            --- Solving augmented matrix using the Gauss Elimination ---
            ...

            --- S(x) = { 2 + (3/4)*(x-1) + (1/4)*(x-1)^3                        for x in [1, 2] (first range)
                    { 3 + (3/2)*(x-2) + (3/4)*(x-2)^2 - (1/4)*(x-2)^3        for x in [2, 3] (second range)
        """
        self.known_points = known_points
        self.points_count = len(known_points if known_points is not None else known_vectorX)
        self.spline_count = self.points_count - 1

        # Get the X and Y vectors from the known_points, or use known_vectorX and known_vectorY if provided
        self.vectorX = known_vectorX if known_vectorX is not None else [elem[0] for elem in known_points]
        self.vectorY = known_vectorY if known_vectorY is not None else [elem[1] for elem in known_points]

        # Initialize vectorH - vector of vectorX[i] - vectorX[i-1] or vectorX[i + 1] - vectorX[i]
        self.vectorH = [self.vectorX[i] - self.vectorX[i - 1] for i in range(1, self.points_count)]

        # Initializing vector of undefined variables A
        self.coefficientsA = copy.deepcopy(self.vectorY)

        # Get the dirrerences of Y, to easily iterate
        self.deltaY = np.diff(self.vectorY)

        # Construct the matrix of coefficients
        # Initializing matrix A and getting the calculated C coefficients from solving the system
        self.matrixA, self.vectorC = self.construct_tridiagonal_matrix()

        # Solving tridiagonal matrix and getting calculated unknown C coefficients vector
        self.coefficientsC = self.solve_tridiagonal_matrix(self.matrixA, self.vectorC, vars)
        self.coefficientsC = np.array(self.coefficientsC)

        # Initialize empty 1D vectors of B and D unknown coefficients of size n - 1
        self.coefficientsB = np.zeros(shape=(self.points_count - 1, 1))
        self.coefficientsD = np.zeros(shape=(self.points_count - 1, 1))

        # Iterating to calculate unknown B and D coefficients from the C coefficients vector
        for i in range(0, len(self.coefficientsD)):
            self.coefficientsD[i] = (self.coefficientsC[i + 1] - self.coefficientsC[i]) / (3 * self.vectorH[i])
            self.coefficientsB[i] = (self.deltaY[i]/self.vectorH[i]) - (self.vectorH[i]/3) * (2*self.coefficientsC[i] + self.coefficientsC[i+1])

        # Printing final table of coefficients
        self.print_results_table()

        # def get_tridiagonal_matrix(a: list, b: list, c: list, k1=-1, k2=0, k3=1):
        #     return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

    def get_xy(self, step: float = 0.1, resolution=10):
        # May be used lated to create the vectors of points for each spline
        self.step = step
        self.resolution = resolution

        for i in len(self.spline_count):
            exit

    def print_sx(self, x: float, spline_index: int):
        print(f'S{spline_index}({x}) = {self.get_sx(x, spline_index)}')

    def get_sx(self, x: float, spline_index: int):
        i = spline_index
        previous_x = self.vectorX[spline_index]
        h = x - previous_x
        print("h = ", h)
        return self.coefficientsA[i] + self.coefficientsB[i][0] * h + self.coefficientsC[i] * (h ** 2) + self.coefficientsD[i][0] * (h ** 3)

    def solve_tridiagonal_matrix(self, matrixA, vectorB, vars) -> list:
        # Solve the matrix using the Gaussian Elimination method (or tridiagonal matrix solve algorithm)
        # Calculating C unknown coefficients vector
        # Thomas algorithm
        n = len(matrixA)

        # see https://e.tsi.lv/pluginfile.php/130692/mod_resource/content/2/%D0%A7%D0%B8%D1%81%D0%BB%D0%B5%D0%BD%D0%BD%D1%8B%D0%B5%20%D0%BC%D0%B5%D1%82%D0%BE%D0%B4%D1%8B_LECTURES-2017.pdf
        # Initialize two 1D vectors to hold Ai and Bi values (or ak and bk, as in the pdf)
        alphas = [0] * n
        betas = [0] * n

        # Initializing a1 and b1
        alphas[0] = matrixA[0, 1] / matrixA[0, 0]
        betas[0] = vectorB[0][0] / matrixA[0, 0]

        # Calculate values for each alpha and beta
        for i in range(1, n - 1):
            alphas[i] = matrixA[i, i + 1] / (matrixA[i, i] - matrixA[i, i-1] * alphas[i-1])
            betas[i] = (vectorB[i][0] - matrixA[i, i-1] * betas[i-1]) / (matrixA[i, i] - matrixA[i, i-1] * alphas[i - 1])
        
        # Initialize vector of out X results, where x_n = beta_n
        X = [0] * n
        X[n - 1] = betas[n - 1]

        print("alphas=",alphas)
        print("betas=",betas)

        # Backward substitution
        for i in range(n - 2, -1, -1):
            X[i] = betas[i] - alphas[i] * X[i + 1]

        return X

        #import testing2
        #return testing2.jacobi(matrixA, vectorB, np.zeros(len(self.matrixA)), tol=1e-100, n_iterations=1000)

    def print_results_table(self):
        print("Cubic Spline Interpolation results table ---------------------------------------------")
        A = list(self.coefficientsA)
        B = list(self.coefficientsB.tolist())
        C = list(self.coefficientsC.tolist())
        D = list(self.coefficientsD.tolist())
        print(f'{"Spline number":<16} | {"ai":<16} | {"bi":<16} | {"ci":<16} | {"di":<16}')
        print('{:-<80}'.format(""))
        for i in range(0, self.spline_count):
            print(
                f'{i:<16} | {A[i]}{"":<{16 - len(str(A[i]))}} | {B[i][0]}{"":<{16 - len(str(B[i][0]))}} | {C[i]}{"":<{16 - len(str(C[i]))}} | {D[i][0]}{"":<{16 - len(str(D[i][0]))}}')
        print("--------------------------------------------------------------------------------------")

    def construct_tridiagonal_matrix(self):
        # --- Constructing the tridiagonal matrix A as hi*ci + 2*(hi-1 + hi)*ci + hi*ci+1

        # Initializing the squared A matrix with zeros
        A = np.zeros(shape=(self.points_count, self.points_count))

        # Initializing the 1D vector B with zeros
        b = np.zeros(shape=(self.points_count, 1))

        # Set the first element (the most upper-left element in the matrix) to 1
        A[0, 0] = 1

        # Set the last element (the most lower-right element in the matrix) to 1
        A[-1, -1] = 1

        # See https://e.tsi.lv/pluginfile.php/130692/mod_resource/content/2/%D0%A7%D0%B8%D1%81%D0%BB%D0%B5%D0%BD%D0%BD%D1%8B%D0%B5%20%D0%BC%D0%B5%D1%82%D0%BE%D0%B4%D1%8B_LECTURES-2017.pdf
        for i in range(1, self.points_count - 1):
            A[i, i-1] = self.vectorH[i-1]                   # Set lower-diagonal element
            A[i, i+1] = self.vectorH[i]                     # Set upper-diagonal element
            A[i, i] = 2*(self.vectorH[i-1]+self.vectorH[i]) # Set main-diagonal element

            # Get vector B (C coefficients)
            b[i, 0] = 3*(self.deltaY[i]/self.vectorH[i] - self.deltaY[i-1]/self.vectorH[i-1])

        # Returning matrixA and vectorB. 
        # Converting them to list, because ndarray from numpy don't have an append method, which is required in the call of gauss_elimination
        return A, b
