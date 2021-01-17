ng = 33
a = -4
b = -3.2

# lambda x: (x / (x + 3) ** 3) - 33)
# lambda x: (x ** 2) * math.sin(x ** 2) - 33

import roots
import math

roots.RootFinder(lambda x: (x / (x + 3) ** 3) - 33).bisection(a, b, tolerance=0.00001)
#roots.RootFinder(lambda x: (x ** 2) * math.sin(x ** 2) - 33).bisection(-7, 8, tolerance=0.0001)
roots.RootFinder(lambda x: math.sqrt(4 - x ** 2)/x - 33).bisection(0.05, 0.07)

def initf(x):
    return (x / (x + 3) ** 3)


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

def dsimpx(x):
    return 1 - 0.002 * x * (math.sin(x**2) + (x ** 2) * math.cos(x ** 2))


def dimf(x):
    return -0.001*((x**2)*math.sin(x**2)-33)


#roots.RootFinder(lambda x: -0.001*((x ** 2) * math.sin(x ** 2) - 33)+x).fixed_point_iteration(
#    x0=7.88,
#    e=0.001,
#    N=100,
#    gfunc=lambda x: 1 - 0.002 * x * (math.sin(x**2) + (x ** 2) * math.cos(x ** 2))
#)

def adjf(x):
    #print("Adjf:", -0.001*((x/(x+3) ** 3) - 33))
    return -0.001*((x/(x+3) ** 3) - 33)


def simpx(x):
    #print("Simpx:", 1+(0.002*x-0.003)*(x+3)**(-4))
    return 1+(0.002*x-0.003)*(x+3)**(-4)

def dimas(it, x):
    print(f"iter={it}, x={x}, fi`(x)={'%0.6f' % dsimpx(x)}, delta={dimf(x)}, new_x={x+dimf(x)}")

def myf(x):
    return -0.001*((x / (x + 3) ** 3) - 33)

def my(it, x):
    print(f"iter={it}, x={x}, fi`(x)={'%0.6f' % simpx(x)}, delta={adjf(x)}, new_x={x+adjf(x)}")

dimas(0, 7.88)
dimas(1, 7.95475257883887)
dimas(2, 7.960443480519134)
dimas(3, 7.961037064476303)
dimas(4, 7.96111254832348)
dimas(5, 7.961122343397132)
dimas(6, 7.961123617708043)
dimas(7, 7.961123783547402)
my(0, -3.5)
my(1, -3.495)
my(2, -3.4908158598535093)
my(3, -3.487339586605744)
my(4, -3.4844695891278286)
my(5, -3.482113000805781)
my(6, -3.4801869414452717)
my(7, -3.47861890679347)


print((2/2.25)*f(1))
print((5/9)*f(1.7745966692414834))
print((5/9)*f(0.2254033307585166))


def ef(x):
    return math.sqrt(4-x ** 2)/x


print("\n", (2/2.25)*ef(0.6))
print((5/9)*ef(0.9098386676965934))
print((5/9)*ef(0.29016133230340657))

