ng = 33
a = -4
b = -3.2

# lambda x: (x / (x + 3) ** 3) - 33)
# lambda x: (x ** 2) * math.sin(x ** 2) - 33

import roots
import math

roots.RootFinder(lambda x: (x / (x + 3) ** 3) - 33).bisection(a, b, tolerance=0.00001)
#roots.RootFinder(lambda x: (x ** 2) * math.sin(x ** 2) - 33).bisection(-7, 8, tolerance=0.0001)


def f(x):
    return (x / (x + 3) ** 3) - 33


def xi(a, b):
    return b - ((b - a) / (f(b) - f(a))) * f(b)


def df(x):
    return (x ** 2) * math.sin(x ** 2) - 33


def dxi(a, b):
    return b - ((b - a) / (df(b) - df(a))) * df(b)


def newton(x):
    return (3-2*x)*(x+3) ** (-4)


def newtonxi(x):
    return x - (f(x) / newton(x))

it1x = xi(-4, -3.2)
it1gx = f(it1x)
it1gagx = it1gx * f(-4)

