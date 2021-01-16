import roots
import math

funcs = [
    lambda x: (x / (x + 3) ** 3),
    lambda x: (x / (x + 3) ** 3) - 33,
    lambda x: (x ** 2) / (math.sqrt(2*x + 1)),
    lambda x: (x ** 2) * math.sin(x ** 2)
]

def cals(fid, a, b, h):
    def right_rectangle(f, x):
        return f(x) * h

    def trapezoidal(f, x_prev, x_curr):
        return ((f(x_prev) + f(x_curr)) / 2) * h

    def simpsons(f, x_prev, x_curr, x_next):
        return (f(x_prev) + 4 * f(x_curr) + f(x_next)) * (h / 6)

    def deriv(f, x_prev, x_next):
        return (f(x_next) - f(x_prev)) / (2 * h)

    def largef(f, x_prev, si):
        return f(x_prev) + si

    x_1 = a - h
    x0 = a
    x1 = a + h
    x2 = a + h + h
    x3 = a + h + h + h
    x4 = a + h + h + h + h
    x5 = a + h + h + h + h + h
    x6 = a + h + h + h + h + h + h

    fxs = [
        funcs[fid](x_1),
        funcs[fid](x0),
        funcs[fid](x1),
        funcs[fid](x2),
        funcs[fid](x3),
        funcs[fid](x4),
        funcs[fid](x5),
        funcs[fid](x6)
    ]
    print("fxs:", fxs)

    rrect_fxs = [
        right_rectangle(funcs[fid], x0),
        right_rectangle(funcs[fid], x1),
        right_rectangle(funcs[fid], x2),
        right_rectangle(funcs[fid], x3),
        right_rectangle(funcs[fid], x4),
        #right_rectangle(funcs[fid], x5)
    ]
    print("Right rect fxs:", rrect_fxs, "| rect sum =", sum(rrect_fxs))

    trap_fxs = [
        trapezoidal(funcs[fid], x0, x1),
        trapezoidal(funcs[fid], x1, x2),
        trapezoidal(funcs[fid], x2, x3),
        trapezoidal(funcs[fid], x3, x4),
        trapezoidal(funcs[fid], x4, x5)
    ]
    print("trap fxs:", trap_fxs, "| trap sum =", sum(trap_fxs))

    #simp_fxs = [
    #    simpsons(funcs[fid], x0, x1, x2),
    #    simpsons(funcs[fid], x1, x2, x3),
    #    simpsons(funcs[fid], x2, x3, x4),
    #    simpsons(funcs[fid], x3, x4, x5),
    #    simpsons(funcs[fid], x4, x5, x5 + h),
    #]
    #print("simp fxs:", simp_fxs, "| simp sum =", sum(simp_fxs))

    f_derivxs = [
        deriv(funcs[fid], x_1, x1),
        #0,
        deriv(funcs[fid], x0, x2),
        deriv(funcs[fid], x1, x3),
        deriv(funcs[fid], x2, x4),
        deriv(funcs[fid], x3, x5),
        #0
        deriv(funcs[fid], x4, x6),
    ]
    print("deriv fxs:", f_derivxs)

    f_largexs = [
        #0,
        trap_fxs[0],
        trap_fxs[0] + trap_fxs[1],
        trap_fxs[0] + trap_fxs[1] + trap_fxs[2],
        trap_fxs[0] + trap_fxs[1] + trap_fxs[2] + trap_fxs[3],
        trap_fxs[0] + trap_fxs[1] + trap_fxs[2] + trap_fxs[3] + trap_fxs[4],
    ]
    print("large f fxs:", f_largexs)

#cals(fid=0, a=0, b=2)
cals(fid=0, a=0, b=2, h=0.4)
print("\n")
cals(fid=3, a=0, b=1, h=0.2)

#rrect_fxs = [
#    right_rectangle(funcs[0], x0),
#    right_rectangle(funcs[0], x1),
#    right_rectangle(funcs[0], x2),
#    right_rectangle(funcs[0], x3),
#    right_rectangle(funcs[0], x4),
#    right_rectangle(funcs[0], x5)
#]
#print("Right rect fxs:", rrect_fxs, "| rect sum =", sum(rrect_fxs))
#
#trap_fxs = [
#    trapezoidal(funcs[0], x0, x1),
#    trapezoidal(funcs[0], x1, x2),
#    trapezoidal(funcs[0], x2, x3),
#    trapezoidal(funcs[0], x3, x4),
#    trapezoidal(funcs[0], x4, x5)
#]
#print("trap fxs:", trap_fxs, "| trap sum =", sum(trap_fxs))
#
#simp_fxs = [
#    simpsons(funcs[0], x0, x1, x2),
#    simpsons(funcs[0], x1, x2, x3),
#    simpsons(funcs[0], x2, x3, x4),
#    simpsons(funcs[0], x3, x4, x5),
#    simpsons(funcs[0], x4, x5, x5+h),
#]
#print("simp fxs:", simp_fxs, "| simp sum =", sum(simp_fxs))

roots.RootFinder(
    function=funcs[1],
).bisection(0, 2)

roots.RootFinder(
    function=lambda x: (x ** 2) * math.sin(x ** 2) - 33,
).parabolic(-7, 8, 100)

#roots.RootFinder(
#    function=lambda x: (x ** 2) * math.sin(x ** 2) - 33,
#).bisection(0, 1)
