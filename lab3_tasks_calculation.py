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


def testing(function, x1, x2, formula, tolerance, initcheck_bis=False, initcheck_par=False):
    rf = roots.RootFinder(function, None)
    rf_bis = rf.bisection(x1, x2, max_iterations=500, tolerance=tolerance, check_initial=initcheck_bis)
    rf_parab = rf.parabolic(x1, x2, max_iterations=500, tolerance=tolerance, check_initial=initcheck_par)
    # print("Bisection =", '%0.5f' % rf_bis)
    # print("Parabolic =", '%0.5f' % rf_parab)

    plt.figure("Parabolic v Bisection " + formula)
    plt.title(formula + str(f' with tolerance = {tolerance}'))
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
    plt.grid()
    plt.axvline(x=0, color="grey")
    plt.axhline(y=0, color='grey')
    plt.axvline(x=rf_parab, color="green", linestyle="--")
    plt.axhline(y=rf.function(rf_parab), color='green', linestyle="--")
    plt.axvline(x=rf_bis, color="red", linestyle="--")
    plt.axhline(y=rf.function(rf_bis), color='red', linestyle="--")
    plt.xlabel("X values")
    plt.ylabel("Y values")
    plt.show()


testing(funcs[0], -1.5, 3, "x / (x + 3) ** 3", 0.01)
testing(funcs[1], -9, 2, "2 ** (3 * x)", 0.01, initcheck_bis=False)
testing(funcs[2], 0.05, 1.99, "math.sqrt(4 - x ** 2) / x", 0.01, initcheck_par=True)

testing(funcs[0], -1.5, 3, "x / (x + 3) ** 3", 0.0001)
testing(funcs[1], -9, 2, "2 ** (3 * x)", 0.0001, initcheck_bis=False)
testing(funcs[2], 0.05, 1.99, "math.sqrt(4 - x ** 2) / x", 0.0001, initcheck_par=True)

f1 = lambda x: x**2 - x - 1
print("Bisection =", roots.RootFinder(f1, None).bisection(lower=1, upper=2))

print("My bisection =", roots.RootFinder(funcs[2], None).bisection(
    lower=0.05, upper=2, tolerance=0.001))
print("My bisection 2 =", roots.RootFinder(lambda x: x ** 3 - x - 2, None).bisection(1, 2, tolerance=0.0001))

print("My bisection 3 =", roots.RootFinder(funcs[1], None).bisection(
    lower=-13, upper=2, tolerance=0.001))

print("My bisection 3 =", roots.RootFinder(funcs[1], None).bisection(
    lower=-13, upper=2, tolerance=0.001))

f2 = lambda x: (x ** 3) + 2 * (x ** 2) + 10*x - 20
print("My parabolic 4 =", roots.RootFinder(f2, None).parabolic(
    a=0, b=2, tolerance=0.0001))
print("My bisection 4 =", roots.RootFinder(f2, None).bisection(
    lower=0, upper=2, tolerance=0.0001))
