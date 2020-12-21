import differential_equations as de
import matplotlib.pyplot as plt
import math

funcs = [
    lambda x, y: 0.9*y - 0.2*(y ** 2),
    lambda x, y: (4.5 * math.exp(0.9 * x)) / (3.5 + math.exp(0.9 * x))
]

# 1 ----
ex1 = de.Exact(funcs[1], 0.05)
xy1_ex = ex1.get_xy(0, 1, 10)

euler1 = de.Euler(funcs[0], 0.05)
xy1_euler = euler1.get_xy(0, 1, 10)

rk1 = de.RungeKutta(funcs[0], 0.05)
xy1_rk = rk1.get_xy(0, 1, 10)

euler1_errors = de.DifferentialEquationsSolver.compare_errors(xy1_ex[1], xy1_euler[1])
rk1_errors = de.DifferentialEquationsSolver.compare_errors(xy1_ex[1], xy1_rk[1])

# 2 ----
ex2 = de.Exact(funcs[1], 0.01)
xy2_ex = ex2.get_xy(0, 1, 10)

euler2 = de.Euler(funcs[0], 0.01)
xy2_euler = euler2.get_xy(0, 1, 10)

rk2 = de.RungeKutta(funcs[0], 0.01)
xy2_rk = rk2.get_xy(0, 1, 10)

euler2_errors = de.DifferentialEquationsSolver.compare_errors(xy2_ex[1], xy2_euler[1])
rk2_errors = de.DifferentialEquationsSolver.compare_errors(xy2_ex[1], xy2_rk[1])

plt.figure("Methods comparison")
plt.title("Methods comparison")
plt.plot(
    xy1_euler[0], xy1_euler[1], "r--",
    xy2_euler[0], xy2_euler[1], "r",
    xy1_rk[0], xy1_rk[1], "g--",
    xy2_rk[0], xy2_rk[1], "g",
    xy1_ex[0], xy1_ex[1], "y--",
    xy2_ex[0], xy2_ex[1], "y",
)
plt.legend(
    [
        "Euler (0.05)",
        "Euler (0.01)",
        "Runge Kutta (0.05)",
        "Runge Kutta (0.01)",
        "Exact (0.05)",
        "Exact (0.01)"
    ]
)
plt.show()

plt.figure("Errors comparison")
plt.title("Errors comparison")
plt.plot(
    xy1_euler[0][1:], euler1_errors, "r--",
    xy1_rk[0][1:], rk1_errors, "y--",
    xy2_euler[0][1:], euler2_errors, "r",
    xy2_rk[0][1:], rk2_errors, "y",
)
plt.legend(
    [
        "Euler (0.05)",
        "Runge Kutta (0.05)",
        "Euler (0.01)",
        "Runge Kutta (0.01)",
    ]
)
plt.show()
