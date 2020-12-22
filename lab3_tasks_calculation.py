import roots
import lab3_parabolic_simpsons_rule.parabolic as p
import matplotlib.pyplot as plt
import numpy as np
import math

funcs = [
    lambda x: x / (x + 3) ** 3,
    lambda x: 2 ** (3 * x),
    lambda x: math.sqrt(4 - x ** 2) / x
]


def testing(function, x1, x2):
    rf = roots.RootFinder(function, None)
    rf_bis = rf.bisection(x1, x2, max_iterations=200)
    rf_parab = rf.parabolic(x1, x2, max_iterations=200)
    # print("Bisection =", '%0.5f' % rf_bis)
    # print("Parabolic =", '%0.5f' % rf_parab)

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
    xs = np.linspace(x1, x2, 1000)
    #print("len(xs) =", len(xs), ", xs =", xs)
    ys = np.vectorize(rf.function)(xs)
    plt.plot(xs, ys)
    plt.legend(
        [
            "Known points",
            "Found root (Parabolic) (" + str('%0.4f' % rf_parab) + " ; " + str(
                '%0.4f' % rf.function(rf_parab)) + ")",
            "Found root (Bisection) (" + str('%0.4f' % rf_bis) + " ; " + str('%0.4f' % rf.function(rf_bis)) + ")",
            "Function graph by linspace"
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


testing(funcs[0], -1.5, 3)
testing(funcs[1], -24, 2)
testing(funcs[2], 0.05, 2)

f1 = lambda x: x**2 - x - 1
print("Bisection =", roots.RootFinder(f1, None).bisection(lower=1, upper=2))
