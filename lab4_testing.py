import lab4_euler.euler as eu
import differential_equations as de
import matplotlib.pyplot as plt
import numpy as np
import math

#xy1_eu = eu.Euler(lambda x, y: -3*y + 2*(x ** 2), 0.25).get_xy(
#    initial_x=0, initial_y=13, last_x=1
#)
#plt.plot(xy1_eu[0], xy1_eu[1])
#plt.show()

#xy1_rk = rk.RungeKutta()

funcs = [
    lambda x, y: -3*y + 2*(x ** 2),
    lambda x, y: 0.9*y - 0.2*(y ** 2),
    lambda x, y: y,
    lambda x, y: y**2,
    lambda x, y: math.sin(x) ** 2 * y
]
y2_exact = [13, 6.15, 2.96, 1.54, 1.01]

#ex2 = de.Exact(funcs[0], 0.25)
#xy2_ex = ex2.get_xy(initial_x=0, initial_y=13, last_x=1)

eu2 = de.Euler(funcs[0], 0.25, y2_exact)
xy2_eu = eu2.get_xy(initial_x=0, initial_y=13, last_x=1)

rk2 = de.RungeKutta(funcs[0], 0.25, y2_exact)
xy2_rk = rk2.get_xy(initial_x=0, initial_y=13, last_x=1)

plt.plot(
    xy2_eu[0], y2_exact, "ro--",
    xy2_eu[0], xy2_eu[1], "go--",
    xy2_rk[0], xy2_rk[1], "bo--"
)
plt.legend(
    [
        "Exact",
        "Euler",
        "Runge Kutta",
    ]
)
plt.show()

diff_euler = abs(np.array(xy2_eu[1]) - np.array(y2_exact))
diff_rk = abs(np.array(xy2_rk[1]) - np.array(y2_exact))

plt.plot(
    xy2_eu[0], diff_euler, "ro--",
    xy2_rk[0], diff_rk, "bo--"
)
plt.show()


eu3 = de.Euler(funcs[1], 0.05)
xy3_eu = eu3.get_xy(initial_x=0, initial_y=1, last_x=10)
eu4 = de.Euler(funcs[1], 0.01)
xy4_eu = eu4.get_xy(initial_x=0, initial_y=2, last_x=10)
eu5 = de.Euler(funcs[1], 0.5)
xy5_eu = eu5.get_xy(initial_x=0, initial_y=3, last_x=10)

rk3 = de.RungeKutta(funcs[1], 0.05)
xy3_rk = rk3.get_xy(initial_x=0, initial_y=4, last_x=10)
rk4 = de.RungeKutta(funcs[1], 0.01)
xy4_rk = rk4.get_xy(initial_x=0, initial_y=5, last_x=10)
rk5 = de.RungeKutta(funcs[1], 0.5)
xy5_rk = rk5.get_xy(initial_x=0, initial_y=6, last_x=10)

plt.plot(
    xy3_eu[0], xy3_eu[1], "g--",
    xy3_rk[0], xy3_rk[1], "b--",
    xy4_eu[0], xy4_eu[1], "y--",
    xy4_rk[0], xy4_rk[1], "r--",
    xy5_eu[0], xy5_eu[1], "m--",
    xy5_rk[0], xy5_rk[1], "--"
)
plt.legend(
    [
        "Euler (0.05)",
        "Runge Kutta (0.05)",
        "Euler (0.01)",
        "Runge Kutta (0.01)",
        "Euler (0.5)",
        "Runge Kutta (0.5)"
    ]
)
plt.show()

x6_exact = np.linspace(0, 2, 21)
y6_exact = np.exp(x6_exact)
xy6_eu = de.Euler(funcs[2], 0.5).get_xy(0, 1, 2)
xy6_rk = de.RungeKutta(funcs[2], 0.5).get_xy(0, 1, 2)
plt.plot(
    xy6_eu[0], xy6_eu[1], "g--",
    xy6_rk[0], xy6_rk[1], "b--",
    x6_exact, y6_exact
)
plt.legend(
    [
        "Euler (0.1)",
        "Runge Kutta (0.1)",
        "Exact"
    ]
)
plt.show()

def solved(x, y):
    a = 0.9
    b = 0.2
    #return (a*math.exp(a*x)) / (a - b + b*math.exp(a*x))
    #return a/b
    #return ((a**2) * math.exp(a*b)) / ((a**2) + a*b*math.exp(a*b) - math.exp(a*b) + 1)
    return (4.5 * math.exp(0.9 * x)) / (3.5 + math.exp(0.9 * x))

#x7_exact = np.linspace(0, 2, 21)
#y7_exact = [solved(x, y) for ]


xy7_ex = de.Exact(solved, 0.01).get_xy(0, 1, 10)
xy7_eu = de.Euler(funcs[1], 0.01).get_xy(0, 1, 10)
xy7_rk = de.RungeKutta(funcs[1], 0.01).get_xy(0, 1, 10)
plt.plot(
    xy7_eu[0], xy7_eu[1], "g--",
    xy7_rk[0], xy7_rk[1], "b--",
    xy7_ex[0], xy7_ex[1]
)
plt.legend(
    [
        "Euler (0.01)",
        "Runge Kutta (0.01)",
        "Exact"
    ]
)
plt.show()

xy8_eu = de.Euler(funcs[4], 0.25).get_xy(0, 1, 5)
xy8_rk = de.RungeKutta(funcs[4], 0.25).get_xy(0, 1, 5)
plt.plot(
    xy8_eu[0], xy8_eu[1], "g--",
    xy8_rk[0], xy8_rk[1], "b--",
)
plt.legend(
    [
        "Euler (0.25)",
        "Runge Kutta (0.25)",
        "Exact"
    ]
)
plt.show()

