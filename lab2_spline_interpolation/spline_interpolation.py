import numpy as np
import copy
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lab1_gauss_elimination.gauss_elimination as ge

class CubicSplineData:
    def __init__(self, known_vectorX, known_vectorY, known_points, step, resolution, 
    vectorA, vectorB, vectorC, vectorD, vectorH, interpolated_vectorX_size, points_count):
        self.known_vectorX = known_vectorX
        self.known_vectorY = known_vectorY
        self.known_points = known_points
        self.step = step
        self.resolution = resolution
        self.vectorA = vectorA
        self.vectorB = vectorB
        self.vectorC = vectorC
        self.vectorD = vectorD
        self.vectorH = vectorH

    def interpolate(self):
        interp_resolution = self.interpolated_vectorX_size / self.points_count
        for i in range(0, len(vectorH)):
            for k in range(0, resolution):
                deltaX = k / resolution * vectorH[i]
                termA = vectorA[i]
                termB = vectorB[i]
                termC = vectorC[i]
                termD = vectorD[i]
                interpolated_index = i * resolution + k
                interpolated_vectorX[interpolated_index] = deltaX + self.vectorX[i]
                interpolated_vectorY[interpolated_index] = termA + termB + termC + termD
        # Remove uninitialized variables
        #points_to_keep = resolution * (points_count - 1)
        #interpolated_vectorX_copy = copy.deepcopy()
        # ...
        # https://github.com/swharden/ScottPlot/blob/404105e5d7ae8399b2e40e9bd64b246d3b3b80dd/src/ScottPlot/Statistics/Interpolation/SplineInterpolator.cs

    def integrate(self):
        integral = 0
        for i in range(0, len(self.vectorH)):
            termA = self.vectorA[i] * self.vectorH[i]
            termB = self.vectorB[i] * (self.vectorH[i] ** 2) / 2.0
            termC = self.vectorC[i] * (self.vectorH[i] ** 3) / 3.0
            termD = self.vectorD[i] * (self.vectorH[i] ** 4) / 4.0
        return integral

def cubic_spline_interpolation(known_vectorX: list, known_vectorY: list, known_points: list, step: float = 0.1, resolution = 10):
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
    points_count = len(known_points if known_points is not None else known_vectorX)

    # Get the X and Y vectors from the known_points, or use known_vectorX and known_vectorY if provided
    vectorX = known_vectorX if known_vectorX is not None else [elem[0] for elem in known_points]
    vectorY = known_vectorY if known_vectorY is not None else[elem[1] for elem in known_points]

    # Initialize vectorH - vector of vectorX[i] - vectorX[i-1] or vectorX[i + 1] - vectorX[i]
    vectorH = [vectorX[i + 1] - vectorX[i] for i in range(0, points_count - 1)]

    print('vectorX =', vectorX, '\nvectorY =', vectorY, '\nvectorH =', vectorH)

    # Initializing vectors of undefined variables A, B, C, D
    vectorA = copy.deepcopy(vectorY)
    vectorB = list()
    vectorC = list()
    vectorD = list()
    
    # Initializing vectors of interpolated X and Y
    interpolated_vectorX = list()
    interpolated_vectorY = list()

    # Set the initial value of the current_value as the x of the first known_point
    current_value = vectorX[0]

    # Populate interpolated_vectorX with the specified step
    for i in range(0, resolution, step):
        interpolated_vectorX.append(current_value)
        current_value += step

    print("Interpolated vectorX =", interpolated_vectorX)

    # Initializing the variable that holds the length of the newly populated interpolated_vectorX
    interpolated_vectorX_size = len(interpolated_vectorX)

    # Construct the matrix of coefficients
    matrix = exit

    # Solve the matrix using the Gaussian Elimination method (or tridiagonal matrix solve algorithm)
    #ge.gauss_elimination(matrix, ['b0', 'd0', 'b1', 'c1', 'd1'], matrix_name="Solved augmented matrix")

    #---------------------------

    # Initializing the DX and DY vectors
    #vectorDX = list()
    #vectorDY = list()

    # Populating the DX and DY vectors
    #if known_vectorX is not None:
    #    for i, (x, y) in enumerate(zip(known_vectorX, known_vectorY)):
    #        vectorDX.append(vectorX[i + 1] - x)
    #        vectorDY.append(vectorY[i + 1] - y)
    #else:
    #    for i in range(0, len(known_points) - 1):
    #        vectorDX.append(vectorX[i + 1] - vectorX[i])
    #        vectorDY.append(vectorY[i + 1] - vectorY[i])

    

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

    def solve_tridiag_matrix(sub: list, diag: list, sup: list, b: list, n: int):
        for i in range(2, n + 1):
            sub[i] /= diag[i - 1]
            diag[i] -= sub[i] * sup[i - 1]
            b[i] -= sub[i] * b[i - 1]

        b[n] /= diag[n]

        for i in range(n - 1, 0, -1):
            b[i] = (b[i] - sup[i] * b[i + 1]) / diag[i]

    def tridiag(a: list, b: list, c: list, k1=-1, k2=0, k3=1):
        return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def interpolate_cubic(p0, p1, p2, p3, p4, x):
    return (-(1/2)*p0 + (3/2)*p1 - (3/2)*p2 + (1/2)*p3)*(x**3) 
    + (p0 - (5/2)*p1 + 2*p2 - (1/2)*p3)*(x**2) + (-(1/2)*p0 + (1/2)*p2)*x + p1
