
from sympy import symbols, powsimp, simplify
from sympy.solvers.polysys import solve_poly_system

def main():

    x, y, z = symbols('x y z')
    x1, y1, z1 = symbols('x1 y1 z1')
    x2, y2, z2 = symbols('x2 y2 z2')
    ct = symbols('ct')

    F = [x1*x + y1*y + z1*z - ct,
         x2*x + y2*y + z2*z - ct,
         x**2 + y**2 + z**2 - 1]

    s = solve_poly_system(F, x, y, z)
    print(simplify(s[0][0]))
    print('\n')
    print(simplify(s[0][1]))
    print('\n')
    print(simplify(s[0][2]))
    print('\n')
    print(simplify(s[1][0]))
    print('\n')
    print(simplify(s[1][1]))
    print('\n')
    print(simplify(s[1][2]))

    print(s[0][0] == s[1][0])
    print(s[0][1] == s[1][1])
    print(s[0][2] == s[1][2])

if __name__ == '__main__':
    main()

