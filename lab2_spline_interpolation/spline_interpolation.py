import numpy as np
import copy
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lab1_gauss_elimination.gauss_elimination as ge
import matrix_helpers as mh

class CubicSplineInterpolator:

    def __init__(self, known_vectorX: list, known_vectorY: list, known_points: list, step: float = 0.1, resolution = 10, vars = []):
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
        self.step = step
        self.resolution = resolution

        self.points_count = len(known_points if known_points is not None else known_vectorX)
        self.spline_count = self.points_count - 1

        # Get the X and Y vectors from the known_points, or use known_vectorX and known_vectorY if provided
        self.vectorX = known_vectorX if known_vectorX is not None else [elem[0] for elem in known_points]
        self.vectorY = known_vectorY if known_vectorY is not None else [elem[1] for elem in known_points]

        # Initialize vectorH - vector of vectorX[i] - vectorX[i-1] or vectorX[i + 1] - vectorX[i]
        self.vectorH = [self.vectorX[i] - self.vectorX[i - 1] for i in range(1, self.points_count)]

        print('vectorX =', self.vectorX, '\nvectorY =', self.vectorY, '\nvectorH =', self.vectorH)

        # Initializing vectors of undefined variables A, B, C, D
        self.vectorA = copy.deepcopy(self.vectorY)
        self.vectorB = list()
        self.vectorC = list()
        self.vectorD = list()
        
        # Initializing vectors of interpolated X and Y
        self.interpolated_vectorX = list()
        self.interpolated_vectorY = list()

        # Initializing list of boundary conditions for natural cubic spline interpolation
        #self.boundary_conditions = list().append((2, np.zeros(self.vectorY.shape[1:])))

        # Set the initial value of the current_value as the x of the first known_point
        #current_value = self.vectorX[0]

        # Populate interpolated_vectorX with the specified step
        #for i in np.arange(0, self.resolution / self.step, self.step):
        #    self.interpolated_vectorX.append(current_value)
        #    current_value += self.step

        #print("Interpolated vectorX =", self.interpolated_vectorX)

        # Initializing the variable that holds the length of the newly populated interpolated_vectorX
        #self.interpolated_vectorX_size = len(self.interpolated_vectorX)

        # Get the dirrerences of X (vectorH, in our case) and Y, to easily iterate
        #delta_x = np.diff(self.vectorX)
        self.delta_y = np.diff(self.vectorY)

        # Construct the matrix of coefficients
        # Initializing matrix
        self.matrix, self.vectorC = self.construct_tridiagonal_matrix()

        #mh.full_print(matrix=self.matrix, vars=vars, title="Constructed matrix")
        #self.matrix = mh.append_vectorB(self.matrix, mh.unpack_vector(self.vectorC))

        print('matrix with appended vector:', self.matrix)
        mh.show_A_B(matrix=self.matrix)

        self.matrix = np.array(self.matrix)
        self.vectorC = np.array(self.vectorC)

        # Solving tridiagonal matrix and getting calculated unknown C coefficients vector
        self.coefficientsC = self.solve_tridiagonal_matrix(self.matrix, self.vectorC, vars)

        # Initialize empty 1D vectors of B and D unknown coefficients of size n - 1
        self.coefficientsB = np.zeros(shape=(self.points_count - 1, 1))
        self.coefficientsD = np.zeros(shape=(self.points_count - 1, 1))

        # Iterating to calculate unknown B and D coefficients from the C coefficients vector
        for i in range(0, len(self.coefficientsD)):
            self.coefficientsD[i] = (self.coefficientsC[i + 1] - self.coefficientsC[i]) / (3 * self.vectorH[i])
            self.coefficientsB[i] = (self.delta_y[i]/self.vectorH[i]) - (self.vectorH[i]/3) * (2*self.coefficientsC[i] + self.coefficientsC[i+1])

        # Printing final table of coefficients
        self.print_results_table()

        # END   --------------------

        #---------------------------
        # --- Helpers ---

        #def get_x(known_point: list) -> float:
        #    return known_point[0]
        
        #def get_y(known_point: list) -> float:
        #    return known_point[1]
        
        #def get_h(index) -> float:
        #    """ Returns the difference between current x (at specified index) and previous x
        #    """
        #    return get_x(known_points[index - 1]) - get_x(known_points[index])
        
        #def get_a(index) -> float:
        #    return get_y(index)
        
        #def get_b(index, m_current, m_next):
        #    return (get_y(index + 1) - get_y(index)) / get_h(index) - (get_h(index) / 2) * m_current - (get_h(index) / 6)*(m_next - m_current)

        #def get_c(index, m_current):
        #    return m_current/2

        #def get_d(index, m_current, m_next):
        #    return (m_next - m_current) / (6 * get_h(index))

        #def get_s_0(index, a, b, c, d):
        #    return a + b*get_h(index) + c*get_h(index) ** 2 + d*get_h(index) ** 3

        #def get_s_1(index, b, c, d):
        #    return b + 2*c*get_h(index) + 3*d*get_h(index) ** 2

        #def get_s_2(index, c, d):
        #    return 2*c + 6*d*get_x(index)

        #def solve_tridiag_matrix(sub: list, diag: list, sup: list, b: list, n: int):
        #    for i in range(2, n + 1):
        #        sub[i] /= diag[i - 1]
        #        diag[i] -= sub[i] * sup[i - 1]
        #        b[i] -= sub[i] * b[i - 1]
#
        #    b[n] /= diag[n]
#
        #    for i in range(n - 1, 0, -1):
        #        b[i] = (b[i] - sup[i] * b[i + 1]) / diag[i]

        def get_tridiagonal_matrix(a: list, b: list, c: list, k1=-1, k2=0, k3=1):
            return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

    def solve_tridiagonal_matrix(self, matrixA, vectorB, vars) -> list:

        # Solve the matrix using the Gaussian Elimination method (or tridiagonal matrix solve algorithm)
        # Calculating C unknown coefficients vector
        #self.coefficientsC = ge.gauss_elimination(
        #    matrix=self.matrix[:-1],
        #    vars=vars,
        #    print_only_results=False,
        #    matrix_name="Solved augmented matrix")

        import testing2
        return testing2.jacobi(matrixA, vectorB, np.zeros(len(self.matrix)), tol=1e-100, n_iterations=1000)

    def print_results_table(self):
        print("Cubic Spline Interpolation results table ---------------------------------------------")
        A = list(self.vectorA)
        B = list(self.coefficientsB.tolist())
        C = list(self.coefficientsC.tolist())
        D = list(self.coefficientsD.tolist())
        print(f'{"Spline number":<16} | {"ai":<16} | {"bi":<16} | {"ci":<16} | {"di":<16}')
        print('{:-<80}'.format(""))
        for i in range(0, self.spline_count):
            print(
                f'{i:<16} | {A[i]}{"":<{16 - len(str(A[i]))}} | {B[i]}{"":<{16 - len(str(B[i]))}} | {C[i]}{"":<{16 - len(str(C[i]))}} | {D[i]}{"":<{16 - len(str(D[i]))}}')
        print("--------------------------------------------------------------------------------------")

    # def interpolate(self):
    #     interp_resolution = self.interpolated_vectorX_size / self.points_count
    #     for i in range(0, len(vectorH)):
    #         for k in range(0, self.resolution):
    #             deltaX = k / self.resolution * vectorH[i]
    #             termA = vectorA[i]
    #             termB = vectorB[i]
    #             termC = vectorC[i]
    #             termD = vectorD[i]
    #             interpolated_index = i * self.resolution + k
    #             interpolated_vectorX[interpolated_index] = deltaX + self.vectorX[i]
    #             interpolated_vectorY[interpolated_index] = termA + termB + termC + termD
        # Remove uninitialized variables
        #points_to_keep = resolution * (points_count - 1)
        #interpolated_vectorX_copy = copy.deepcopy()
        # ...
        # https://github.com/swharden/ScottPlot/blob/404105e5d7ae8399b2e40e9bd64b246d3b3b80dd/src/ScottPlot/Statistics/Interpolation/SplineInterpolator.cs

    # def integrate(self):
    #     integral = 0
    #     for i in range(0, len(self.vectorH)):
    #         termA = self.vectorA[i] * self.vectorH[i]
    #         termB = self.vectorB[i] * (self.vectorH[i] ** 2) / 2.0
    #         termC = self.vectorC[i] * (self.vectorH[i] ** 3) / 3.0
    #         termD = self.vectorD[i] * (self.vectorH[i] ** 4) / 4.0
    #     return integral

    def construct_tridiagonal_matrix(self):
        # Constructing right coefficients vector C (or b, for the convenience)
        print("self.spline_count=", self.spline_count)
        print("len(self.vectorY)=", len(self.vectorY))
        print("len(self.vectorH)=", len(self.vectorH))
        vectorC = [0] * 2
        for i in range(2, self.spline_count):
            print("iteration:", i)
            vectorC.append(
                3 * (self.vectorY[i] - self.vectorY[i - 1]) / self.vectorH[i] - (self.vectorY[i - 1] - self.vectorY[i - 2]) / self.vectorH[i - 1]
            )

        #matrix = list()
        matrix = np.zeros((self.spline_count, self.spline_count))
        print('zeroed_matrix=',matrix)

        # Constructing the tridiagonal matrix as hi*ci + 2*(hi-1 + hi)*ci + hi*ci+1
        print("self.spline_count=", self.spline_count)
        print("len(self.vectorC)=", len(vectorC))
        print("vectorC=",vectorC)
        #for i in range(2, self.spline_count - 1):
        #    print("iteration:",i)
        #    matrix.append(
        #        self.vectorH[i - 1]*vectorC[i - 1] + 2*(self.vectorH[i - 1] + self.vectorH[i]) * vectorC[i] + self.vectorH[i - 1]*vectorC[i + 1]
        #    )

        #for i in range(2, self.spline_count - 1):
        #    print("iteration:",i)
        #    matrix[i][i] = self.vectorH[i - 1]*vectorC[i - 1] + 2*(self.vectorH[i - 1] + self.vectorH[i]) * vectorC[i] + self.vectorH[i - 1]*vectorC[i + 1]
        
        # upper_diagonal = np.zeros((self.spline_count, self.spline_count))
        # print('init upper_diagonal=', upper_diagonal)
        # for i in range(2, self.spline_count - 1):
        #     print("iteration:", i)
        #     print("upper_diagonal[i][i]=", upper_diagonal[i][i])
        #     print("self.vectorH[i-1]=",self.vectorH[i - 1])
        #     print("self.vectorC[i-1]=", self.vectorC[i - 1])
        #     upper_diagonal[i][i] = self.vectorH[i - 1] * self.vectorC[i - 1]
        # print('upper_diagonal=', upper_diagonal)

        # main_diagonal = np.zeros((self.spline_count, self.spline_count))
        # print('init main_diagonal=',main_diagonal)
        # for i in range(1, self.spline_count):
        #     print("iteration:", i)
        #     main_diagonal[i][i] = 2 * (self.vectorH[i - 1] * self.vectorC[i])
        # print('main_diagonal=',main_diagonal)

        # lower_diagonal = np.zeros((self.spline_count, self.spline_count))
        # print('init lower_diagonal=', lower_diagonal)
        # for i in range(1, self.spline_count - 1):
        #     print("iteration:", i)
        #     lower_diagonal[i][i] = self.vectorH[i] * self.vectorC[i + 1]
        # print('lower_diagonal=', lower_diagonal)

        # Constructing the matrix --------------------------------
        ### Get matrix A
        A = np.zeros(shape=(self.points_count, self.points_count))
        b = np.zeros(shape=(self.points_count, 1))
        A[0, 0] = 1
        A[-1, -1] = 1

        # See https://e.tsi.lv/pluginfile.php/130692/mod_resource/content/2/%D0%A7%D0%B8%D1%81%D0%BB%D0%B5%D0%BD%D0%BD%D1%8B%D0%B5%20%D0%BC%D0%B5%D1%82%D0%BE%D0%B4%D1%8B_LECTURES-2017.pdf
        for i in range(1, self.points_count - 1):
            # Set lower-diagonal element
            A[i, i-1] = self.vectorH[i-1]
            # Set upper-diagonal element
            A[i, i+1] = self.vectorH[i]
            # Set main-diagonal element
            A[i, i] = 2*(self.vectorH[i-1]+self.vectorH[i])
            ### Get matrix b (C coefficients)
            b[i, 0] = 3*(self.delta_y[i]/self.vectorH[i] - self.delta_y[i-1]/self.vectorH[i-1])

        #return matrix, vectorC
        # Returning matrixA and vectorB. 
        # Converting them to list, because ndarray from numpy don't have an append method, which is required in the call of gauss_elimination
        return A.tolist(), b.tolist()

    #def interpolate_cubic(p0, p1, p2, p3, p4, x):
    #    return (-(1/2)*p0 + (3/2)*p1 - (3/2)*p2 + (1/2)*p3)*(x**3) 
    #    + (p0 - (5/2)*p1 + 2*p2 - (1/2)*p3)*(x**2) + (-(1/2)*p0 + (1/2)*p2)*x + p1
