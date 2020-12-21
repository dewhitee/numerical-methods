import lab4_euler.euler as eu
import differential_equations as de
import matplotlib.pyplot as plt

#xy1_eu = eu.Euler(lambda x, y: -3*y + 2*(x ** 2), 0.25).get_xy(
#    initial_x=0, initial_y=13, last_x=1
#)
#plt.plot(xy1_eu[0], xy1_eu[1])
#plt.show()

#xy1_rk = rk.RungeKutta()

funcs = [
    lambda x, y: -3*y + 2*(x ** 2)
]
y2_exact = [13, 6.15, 2.96, 1.54, 1.01]
xy2_eu = de.Euler(funcs[0], 0.25, y2_exact).get_xy(
    initial_x=0, initial_y=13, last_x=1
)

xy2_rk = de.RungeKutta(funcs[0], 0.25, y2_exact).get_xy(
    initial_x=0, initial_y=13, last_x=1
)

plt.plot(
    [0, 0.25, 0.5, 0.75, 1], y2_exact, "ro--",
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
