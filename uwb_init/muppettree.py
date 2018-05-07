from sympy import *
import numpy as np

p = Matrix(symbols("p1:4", real=True))
v = Matrix(symbols("v1:4", real=True))

r_i = Symbol("r_i", real=True)
t_i = Symbol("t_i", real=True)
p_i = Matrix(symbols("p_i(1:4)", real=True))


pdotp = Symbol("pdotp", real=True)
pdotv = Symbol("pdotv", real=True)
vdotv = Symbol("vdotv", real=True)

x = Matrix([p, v, Matrix([pdotp, pdotv, vdotv])])

eqns = [Eq(r_i**2, (pdotp + t_i**2*vdotv + 2.*t_i*pdotv - 2.*p_i.dot(p) - 2.*(t_i*p_i).dot(v) + p_i.dot(p_i)))]

_A, _b = linear_eq_to_matrix(eqns, list(x))

_calc_A = lambdify([t_i, r_i, p_i], _A)
_calc_b = lambdify([t_i, r_i, p_i], _b)

def muppettree_fit(observations):
    A = np.vstack([_calc_A(*obs) for obs in observations])
    b = np.vstack([_calc_b(*obs) for obs in observations])
    return np.linalg.lstsq(A,b,rcond=None)
