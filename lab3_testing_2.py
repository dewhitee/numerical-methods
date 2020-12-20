import roots
import lab3_parabolic_simpsons_rule.parabolic as p
import matplotlib.pyplot as plt
import numpy as np
import math


def f1(x):
    2 * (x ** 3) - 2*x - 5


def f2(x):
    x ** 3 - x - 1

x1, x2 = 1.5, 2

rf = roots.RootFinder(lambda x: math.sqrt(4 - x ** 2) / x, None)
rf_bis = rf.bisection(x1, x2, tolerance=0.000001)
rf_parab = rf.parabolic(x1, x2, max_iterations=200)
print("Bisection =", '%0.5f' % rf_bis)
print("Parabolic =", '%0.5f' % rf_parab)

plt.figure("Parabolic v Bisection")
plt.title("Parabolic v Bisection")
plt.plot(
    [x1, x2],
    [rf.function(x1), rf.function(x2)],
    "bo",
    rf_parab,
    rf.function(rf_parab),
    "go",
    rf_bis,
    rf.function(rf_bis),
    "ro"
)
x = np.linspace(x1, x2, 1000)
y = rf.function(x)
plt.plot(x, y)
plt.legend(
    [
        "Known points",
        "Found root (Parabolic) (" + str('%0.4f' % rf_parab) + " ; " + str('%0.4f' % rf.function(rf_parab)) + ")",
        "Found root (Bisection) (" + str('%0.4f' % rf_bis) + " ; " + str('%0.4f' % rf.function(rf_bis)) + ")",
        "Parabola with linspace"
    ]
)
plt.axvline(x=0, color="grey")
plt.axhline(y=0, color='grey')
plt.axvline(x=rf_parab, color="green", linestyle="--")
plt.axhline(y=rf.function(rf_parab), color='green', linestyle="--")
plt.axvline(x=rf_bis, color="red", linestyle="--")
plt.axhline(y=rf.function(rf_bis), color='red', linestyle="--")
plt.xlabel("X values")
plt.ylabel("Y values")
plt.show()

