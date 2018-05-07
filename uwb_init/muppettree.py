from sympy import *
import numpy as np

p = Matrix(symbols("p1:4"))
v = Matrix(symbols("v1:4"))

r_i = Symbol("r_i")
t_i = Symbol("t_i")
p_i = Matrix(symbols("p_i(1:4)"))


pdotp = Symbol("pdotp")
pdotv = Symbol("pdotv")
vdotv = Symbol("vdotv")

x = Matrix([p, v, Matrix([pdotp, pdotv, vdotv])])

eqns = [Eq(r_i**2, (pdotp + t_i**2*vdotv + 2.*t_i*pdotv - 2.*p_i.dot(p) - 2.*(t_i*p_i).dot(v) + p_i.dot(p_i)))]

_A, _b = linear_eq_to_matrix(eqns, list(x))

_calc_A = lambdify([t_i, r_i, p_i], _A)
_calc_b = lambdify([t_i, r_i, p_i], _b)

def muppettree_fit(observations):
    A = np.vstack([_calc_A(*observations[i]) for i in range(len(observations))])
    b = np.vstack([_calc_b(*observations[i]) for i in range(len(observations))])
    return np.linalg.lstsq(A,b,rcond=None)
