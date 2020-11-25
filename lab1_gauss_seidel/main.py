# Algorithm
# 1. Start
# 2. Arrange given system of linear equations in diagonally dominant form
# 3. Read tolerable error (e)
# 4. Convert the first equation in terms of first variable, second eq in terms of second var and so on.
# 5. Set initial guesses for x0, y0, z0 and so on
# 6. Substitute value of y0, z0 ... from step 5 in first equation obtained from step 4 to calculate
# new value of x1. Use x1, z0, u0, ... in second equation obtained from step 4 to calculate new value of y1.
# Similarly, use x1, y1, u0 ... to find new z1 and so on.
# 7. If | x0 - x1 | > e and | y0 - y1 | > e and | z0 - z1 | > e and so on then goto step 9
# 8. Set x0 = x1, y0 = y1, z0 = z1 and so on and goto step 6
# 9. Print value of x1, y1, z1 and so on
# 10. Stop

test_equations = [
    lambda x, y, z: (17 - y + 2*z) / 20,
    lambda x, y, z: (-18 - 3*x + z) / 20,
    lambda x, y, z: (25 - 2*x + 3*y) / 20
]

def gauss_seidel(equations: list, vars: list):
    """ 
    equations -- list of lambda equations with any count of arguments.
    Example: [lambda x, y: x + y, lambda x, y: x - y]
    vars -- list of symbols to specify the names of the equations variables.
    Example: ['x', 'y', 'z']
    """
    count = 1

    # Reading tolerable error
    e = float(input('Enter tolerable error: '))

    # Implementation of Gauss Seidel Iteration
    print('\nCount', *vars, sep = "\t")

    condition = True

    # Initializing all variables values with zero
    vars_values = [0] * len(equations)

    while condition:
        e_list = []
        for i, eq in enumerate(equations):
            v1 = eq(*vars_values)
            e_list.append(abs(vars_values[i] - v1))
            vars_values[i] = v1

        print(count, *["%0.4f" % elem for elem in vars_values], sep="\t")
        count += 1
        condition = check_error_rate(e_list, e)

    print('\nSolution')
    for (var, val) in zip(vars, vars_values):
        print(var,'= %0.3f' %(val))

def testing_gauss_seidel():
    f1 = lambda x, y, z: (17 - y + 2 * z) / 20
    f2 = lambda x, y, z: (-18 - 3*x + z) / 20
    f3 = lambda x, y, z: (25 - 2*x + 3*y) / 20

    # initial setup
    x0 = 0
    y0 = 0
    z0 = 0

    count = 1

    # Reading tolerable error
    e = float(input('Enter tolerable error: '))


    # Implementation of Gauss Seidel Iteration
    print('\nCount\tx\ty\tz\n')

    condition = True

    while condition:
        x1 = f1(x0, y0, z0)
        y1 = f2(x1, y0, z0)
        z1 = f3(x1, y1, z0)
        print('%d\t%0.4f\t%0.4f\t%0.4f\n' % (count, x1, y1, z1))
        e1 = abs(x0-x1)
        e2 = abs(y0-y1)
        e3 = abs(z0-z1)

        count += 1
        x0 = x1
        y0 = y1
        z0 = z1

        condition = e1 > e and e2 > e and e3 > e

    print('\nSolution: x=%0.3f, y=%0.3f and z = %0.3f\n' % (x1, y1, z1))



def check_error_rate(e_list: list, e):
    return all([current_e > e for current_e in e_list])

gauss_seidel(test_equations, ['x','y','z'])
testing_gauss_seidel()
