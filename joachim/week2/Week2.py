from casadi import *
import numpy as np

# Symbols/expressions
x = MX.sym('x1')
y = MX.sym('x2')
z = MX.sym('x3')
f = 3 * x ** 2 + 2 * x * y + x * z + 2.5 * y ** 2 + 2 * y * z + 2 * z ** 2 - 8 * x - 3 * y - 3 * z
g1 = x + z - 3
g2 = y + z

nlp = {}  # NLP declaration
nlp['x'] = vertcat(x, y, z)  # decision vars
nlp['f'] = f  # objective
nlp['g'] = vertcat(g1, g2)

# Create solver instance
F = nlpsol('F', 'ipopt', nlp);

# Solve the problem using a guess
r = F(x0=[2.5, 3.0, 0.75], ubg=0, lbg=0)
print(r['x'])
