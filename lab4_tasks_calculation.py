import differential_equations as de
import matplotlib.pyplot as plt

funcs = [
    lambda x, y: 0.9*y - 0.2*(y**2)
]
xy1_euler = de.Euler(funcs[0], 0.05).get_xy(0, 1, 10)
xy1_rk = de.RungeKutta(funcs[0], 0.05).get_xy(0, 1, 10)

xy2_euler = de.Euler(funcs[0], 0.01).get_xy(0, 1, 10)
xy2_rk = de.RungeKutta(funcs[0], 0.01).get_xy(0, 1, 10)

plt.plot(
    xy1_euler[0], xy1_euler[1], "r",
    xy2_euler[0], xy2_euler[1], "b",
    xy1_rk[0], xy1_rk[1], "g",
    xy2_rk[0], xy2_rk[1], "y"
)
plt.legend(
    [
        "Euler (0.05)",
        "Euler (0.01)",
        "Runge Kutta (0.05)",
        "Runge Kutta (0.01)"
    ]
)
plt.show()

plt.plot(

)
