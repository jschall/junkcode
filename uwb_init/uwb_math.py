import numpy as np
import math
from scipy.optimize import leastsq
from sympy import *
from sympy.utilities import lambdify
import sys
from helpers import CCodePrinter_float

def compressedSymmetricMatrix(prefix, N):
    ret = zeros(N,N)

    r = lambda k: int(floor((2*N+1-sqrt((2*N+1)*(2*N+1)-8*k))/2))
    c = lambda k: int(k - N*r(k) + r(k)*(r(k)-1)/2 + r(k))

    for k in range((N**2-N)/2+N):
        ret[r(k),c(k)] = ret[c(k),r(k)] = Symbol('%s[%u]' % (prefix,k))
    return ret

def vec_norm(v):
    return sqrt(sum([x**2 for x in v]))

def rk(a,b,c,x_0,x_dot,dt):
    N = a.rows
    assert a.cols == N and len(b) == N and len(c) == N

    k = []

    for i in range(N):
        x_n = x_0
        for j in range(1,i):
            x_n += dt*a[i,j]*k[j]
        k.append(x_dot.xreplace(dict(zip(x_0, x_n))))

    x_n = x_0
    for i in range(N):
        x_n += dt*b[i]*k[i]

    return x_n

def rk3(x_0, x_dot, dt):
    a = Matrix([[0, 0, 0],
                [Rational(1,2),0,0],
                [-1,2,0]])
    b = Matrix([Rational(1,6), Rational(2,3), Rational(1,6)])
    c = Matrix([0, Rational(1,2),1])

    return rk(a,b,c,x_0,x_dot,dt)

def copy_upper_to_lower_offdiagonals(M):
    assert isinstance(M,MatrixBase) and M.rows == M.cols

    ret = M[:,:]

    for r in range(ret.rows):
        for c in range(ret.cols):
            if r > c:
                ret[r,c] = ret[c,r]
    return ret

def quickinv_sym(M):
    assert isinstance(M,MatrixBase) and M.rows == M.cols
    n = M.rows
    A = Matrix(n,n,symbols('_X[0:%u][0:%u]' % (n,n)))
    A = copy_upper_to_lower_offdiagonals(A)
    B = Matrix(simplify(A.inv()))
    return B.xreplace(dict(zip(A,M)))

def run_derivations():
    global get_range_residual, get_range_residual_jacobian
    global ekf_predict_x, ekf_predict_P, ekf_update_y, ekf_update_NIS, ekf_update_x, ekf_update_P

    r_obs = Symbol('obs->range')
    obs_t = Symbol('obs->t')
    a_pos = Matrix(symbols('positions[obs.anchor_idx][0:3]'))
    t_pos = Matrix(symbols('positions[obs.tag_idx][0:3]'))
    v_pos = Matrix(symbols('state((0:3)\\,0)'))
    v_vel = Matrix(symbols('state((3:6)\\,0)'))
    dt = Symbol('dt')
    ant_del = Matrix(symbols('state((6:16)\\,0)'))

    x = Matrix([v_pos, v_vel, ant_del])
    nStates = len(x)

    # range residuals and jacobian for initialization fit
    h = vec_norm(a_pos-(v_pos+t_pos+v_vel*obs_t))-ant_del[0]-ant_del[1]
    resid = Matrix([r_obs - h])
    jacob = resid.jacobian(x)

    h = ant_del[0]
    ant_del_resid = Matrix([0 - h])
    ant_del_jacob = ant_del_resid.jacobian(x)

    print(CCodePrinter_float().doprint(resid[0]))
    print("")
    for i in range(16):
        print("%s = %s;\n" % ("ret(0,%u)" % (i,), CCodePrinter_float().doprint(jacob[i])))

    print("")
    print(CCodePrinter_float().doprint(ant_del_resid[0]))
    for i in range(16):
        print("%s = %s;\n" % ("ret(0,%u)" % (i,), CCodePrinter_float().doprint(ant_del_jacob[i])))


def fit_pos(observations, guess=(0.,0.,-1.,0.,0.,0.)):
    def range_residuals(x):
        return [get_range_residual(t, anchor_pos, tag_pos, anchor_range, x) for t, anchor_pos, tag_pos, anchor_range in observations]

    def range_residual_jacob(x):
        return np.concatenate([get_range_residual_jacobian(t, anchor_pos, tag_pos, anchor_range, x) for t, anchor_pos, tag_pos, anchor_range in observations])

    popt, pcov, infodict, errmsg, ier = leastsq(range_residuals, guess, Dfun=range_residual_jacob, full_output=True)

    print errmsg

    return popt, pcov, infodict['fvec'], range_residual_jacob(popt)

def ekf_predict(*args):
    x_n = ekf_predict_x(*args)
    P_n = ekf_predict_P(*args)
    return x_n, P_n

def ekf_update(*args):
    x_n = ekf_update_x(*args)
    P_n = ekf_update_P(*args)
    return x_n, P_n

run_derivations()
