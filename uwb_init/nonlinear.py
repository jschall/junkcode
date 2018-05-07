from scipy.optimize import least_squares
import numpy as np
from sympy import *

p = Matrix(symbols("p1:4", real=True))
v = Matrix(symbols("v1:4", real=True))

r_i = Symbol("r_i", real=True)
t_i = Symbol("t_i", real=True)
p_i = Matrix(symbols("p_i(1:4)", real=True))

x = Matrix([p, v])

_resid = r_i**2 - ((p+v*t_i) - p_i).norm()**2
_jacob = Matrix([_resid]).jacobian(x)

_calc_resid = lambdify([x, t_i, r_i, p_i], _resid)
_calc_jacob = lambdify([x, t_i, r_i, p_i], _jacob)

def nonlinear_fit(observations, x0):
    resid = lambda x: np.hstack([_calc_resid(x, *obs) for obs in observations])
    jacob = lambda x: np.vstack([_calc_jacob(x, *obs) for obs in observations])

    return least_squares(resid, x0, jacob, method='trf')
