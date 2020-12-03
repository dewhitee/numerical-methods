import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
import lab2_spline_interpolation.spline_interpolation as si
import lab1_gauss_elimination.gauss_elimination as ge

x_points = [0, 1, 2, 3, 4, 5]
y_points = [12, 14, 22, 39, 58, 77]

# def f(x):
#     tck = interpolate.splrep(x_points, y_points)
#     return interpolate.splev(x, tck)


# print(f(1.25))

# x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
# y = np.sin(x)
# tck = interpolate.splrep(x, y, s=0)
# xnew = np.arange(0, 2*np.pi, np.pi/50)
# ynew = interpolate.splev(xnew, tck, der=0)

# plt.figure()
# plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
# plt.legend(['Linear', 'Cubic Spline', 'True'])
# plt.axis([-0.05, 6.33, -1.05, 1.05])
# plt.title('Cubic-spline interpolation')
# plt.show()


from scipy.interpolate import CubicSpline

cs = CubicSpline(x_points, y_points, bc_type='natural')

print("S(1.25) =", cs(1.25))

print(cs.c)

# Polynomial coefficients for 0 <= x <= 1
a0 = cs.c.item(3, 0)
b0 = cs.c.item(2, 0)
c0 = cs.c.item(1, 0)
d0 = cs.c.item(0, 0)

# Polynomial coefficients for 1 < x <= 2
a1 = cs.c.item(3, 1)
b1 = cs.c.item(2, 1)
c1 = cs.c.item(1, 1)
d1 = cs.c.item(0, 1)

# ...

# Polynomial coefficients for 4 < x <= 5
a4 = cs.c.item(3, 4)
b4 = cs.c.item(2, 4)
c4 = cs.c.item(1, 4)
d4 = cs.c.item(0, 4)

# Print polynomial equations for different x regions
print('S0(0<=x<=1) = ', a0, ' + ', b0, '(x-0) + ', c0, '(x-0)^2  + ', d0, '(x-0)^3')
print('S1(1< x<=2) = ', a1, ' + ', b1, '(x-1) + ', c1, '(x-1)^2  + ', d1, '(x-1)^3')
print('...')
print('S5(4< x<=5) = ', a4, ' + ', b4, '(x-4) + ', c4, '(x-4)^2  + ', d4, '(x-4)^3')

# So we can calculate S(1.25) by using equation S1(1< x<=2)
print('S(1.25) = ', a1 + b1*0.25 + c1*(0.25**2) + d1*(0.25**3))

print('S(3.6')

# Cubic spline interpolation calculus example
#  https://www.youtube.com/watch?v=gT7F3TWihvk

#ge.gauss_elimination([[-1, 0.5], [0, 0], [3, 3]], ['k0', 'k1', 'k2'])

si.CubicSplineInterpolator(
    known_vectorX=x_points, 
    known_vectorY=y_points,
    known_points=None, 
    vars=['x1','x2','x3','x4','x5','x6'])
