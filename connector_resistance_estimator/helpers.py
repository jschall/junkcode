from math import sqrt, floor, fmod
from sympy import *
from sympy.printing.ccode import *
import math

class CCodePrinter_float(C99CodePrinter):
    def __init__(self,settings={}):
        C99CodePrinter.__init__(self, settings)
        self.known_functions = {
            "Abs": [(lambda x: not x.is_integer, "fabsf")],
            "gamma": "tgammaf",
            "sin": "sinf_fast",
            "cos": "cosf_fast",
            "tan": "tanf",
            "asin": "asinf",
            "acos": "acosf",
            "atan": "atanf",
            "atan2": "atan2f",
            "exp": "expf",
            "log": "logf",
            "erf": "erff",
            "sinh": "sinhf",
            "cosh": "coshf",
            "tanh": "tanhf",
            "asinh": "asinhf",
            "acosh": "acoshf",
            "atanh": "atanhf",
            "floor": "floorf",
            "ceiling": "ceilf",
        }

    def _print_Pow(self, expr):
        if "Pow" in self.known_functions:
            return self._print_Function(expr)
        PREC = precedence(expr)

        invert = expr.exp < 0

        if expr.exp == -1:
            ret = self._print(expr.base)
        else:
            if invert:
                expr = 1/(expr)

            if expr.exp == 0.5:
                ret = 'sqrtf(%s)' % self._print(expr.base)
            elif expr.exp.is_integer and expr.exp <= 4:
                ret = "(%s)" % ('*'.join(["(%s)" % (self._print(expr.base)) for _ in range(expr.exp)]),)
            else:
                ret = 'powf(%s, %s)' % (self._print(expr.base), self._print(expr.exp))

        if invert:
            return '1/(%s)' % (ret,)
        else:
            return ret

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '%d.0f/%d.0f' % (p, q)

def skew(_v):
    v = toVec(_v)
    assert v.rows == 3

    return Matrix([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def vec_norm(v):
    return sqrt(sum([x**2 for x in v]))

def toVec(*args):
    ret = Matrix(map(lambda x: Matrix([x]), args)).vec()
    return ret

def Rz(theta):
    return Matrix([[cos(theta), -sin(theta), 0],[sin(theta),cos(theta), 0], [0,0,1]])

def rot_vec_to_quat_approx(_v):
    v = toVec(_v)
    assert v.rows == 3

    return toVec(1,v*0.5)

def quat_rotate(_q, _v):
    return quat_multiply(_q,rot_vec_to_quat(_v))

def quat_rotate_approx(_q, _v):
    return quat_multiply(_q,rot_vec_to_quat_approx(_v))

def rot_vec_to_quat(_v):
    v = toVec(_v)
    assert v.rows == 3

    theta = sqrt(v[0]**2+v[1]**2+v[2]**2)
    axis = v/theta
    return toVec(cos(theta/2.), sin(theta/2.) * axis[0], sin(theta/2.) * axis[1], sin(theta/2.) * axis[2])

def quat_to_rot_vec(_q):
    q = toVec(_q)
    assert q.rows == 4

    theta = 2.*acos(q[0])
    axis = toVec(q[1],q[2],q[3])/sqrt(q[1]**2+q[2]**2+q[3]**2)

    return theta*axis

def quat_inverse(_q):
    q = toVec(_q)
    assert q.rows == 4

    q[1] = -q[1]
    q[2] = -q[2]
    q[3] = -q[3]
    return q

def quat_normalize(_q):
    q = toVec(_q)
    assert q.rows == 4

    return q/sqrt(q[0]**2+q[1]**2+q[2]**2+q[3]**2)

def quat_multiply(_q1, _q2):
    q1 = toVec(_q1)
    q2 = toVec(_q2)
    assert q1.rows == 4 and q2.rows == 4

    return toVec(q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3],
                 q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2],
                 q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1],
                 q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0])

def quat_to_matrix(_q):
    q = toVec(_q)
    assert q.rows == 4

    return (q[0]**2-(q[1:,0].T*q[1:,0])[0])*eye(3) + 2.*(q[1:,0]*q[1:,0].T) + 2.*q[0]*skew(q[1:,0])

def upperTriangularToVec(M):
    assert M.rows == M.cols

    N = M.rows
    r = lambda k: int(floor((2*N+1-sqrt((2*N+1)*(2*N+1)-8*k))/2))
    c = lambda k: int(k - N*r(k) + r(k)*(r(k)-1)/2 + r(k))
    return Matrix([M[r(k),c(k)] for k in range((N**2-N)/2+N)])

def uncompressSymMatrix(M):
    x = len(M)
    N = int(floor(sqrt(8*x + 1)/2 - 1/2))
    ret = zeros(N)
    r = lambda k: int(floor((2*N+1-sqrt((2*N+1)*(2*N+1)-8*k))/2))
    c = lambda k: int(k - N*r(k) + r(k)*(r(k)-1)/2 + r(k))
    for k in range(x):
        ret[r(k),c(k)] = ret[c(k),r(k)] = M[k]
    return ret


def copy_upper_to_lower_offdiagonals(M):
    assert isinstance(M,MatrixBase) and M.rows == M.cols

    ret = M[:,:]

    for r in range(ret.rows):
        for c in range(ret.cols):
            if r > c:
                ret[r,c] = ret[c,r]
    return ret

def compressedSymmetricMatrix(prefix, N):
    ret = zeros(N,N)

    r = lambda k: int(floor((2*N+1-sqrt((2*N+1)*(2*N+1)-8*k))/2))
    c = lambda k: int(k - N*r(k) + r(k)*(r(k)-1)/2 + r(k))

    for k in range((N**2-N)/2+N):
        ret[r(k),c(k)] = ret[c(k),r(k)] = Symbol('%s[%u]' % (prefix,k))
    return ret

def count_subexpression(subexpr, expr):
    if hasattr(expr, "__getitem__"):
        return sum(map(lambda x: count_subexpression(subexpr, x), expr))
    else:
        return expr.count(subexpr)

def extractSubexpressions(inexprs, prefix='X', threshold=0, prev_subx=[]):
    subexprs, outexprs = cse(inexprs, symbols=numbered_symbols('__TMP__'), order='none')

    subexprs = prev_subx+subexprs

    for i in reversed(range(len(subexprs))):
        from sympy.logic.boolalg import Boolean
        ops_saved = (count_subexpression(subexprs[i][0], [[x[1] for x in subexprs], outexprs])-1)*subexprs[i][1].count_ops()
        if ops_saved < threshold or isinstance(subexprs[i][1], Boolean):
            sub = dict([subexprs.pop(i)])
            subexprs = map(lambda x: (x[0],x[1].xreplace(sub)), subexprs)
            outexprs = map(lambda x: x.xreplace(sub), outexprs)

    for i in range(len(subexprs)):
        newSym = Symbol('%s%u' % (prefix,i+len(prev_subx)))
        sub = {subexprs[i][0]:newSym}
        subexprs[i] = (newSym,subexprs[i][1])
        subexprs = map(lambda x: (x[0],x[1].xreplace(sub)), subexprs)
        outexprs = map(lambda x: x.xreplace(sub), outexprs)

    outexprs = map(lambda x: Matrix(x) if type(x) is ImmutableDenseMatrix else x, outexprs)

    return tuple(outexprs+[subexprs])

def quickinv_sym(M):
    assert isinstance(M,MatrixBase) and M.rows == M.cols
    n = M.rows
    A = Matrix(n,n,symbols('_X[0:%u][0:%u]' % (n,n)))
    A = copy_upper_to_lower_offdiagonals(A)
    B = Matrix(simplify(A.inv()))
    return B.xreplace(dict(zip(A,M)))

def wrap_pi(x):
    while x >= math.pi:
        x -= 2*math.pi

    while x < -math.pi:
        x += 2*math.pi

    return x

def UDUdecomposition(M):
    assert M.rows == M.cols
    assert M.is_symmetric()

    P = M[:,:]

    n = P.rows

    U = zeros(*P.shape)
    D = zeros(*P.shape)

    for j in range(n-1, 0, -1):
        D[j,j] = P[j,j]
        alpha = 1/D[j,j]
        for k in range(j):
            beta = P[k,j]
            U[k,j] = alpha*beta
            for i in range(k+1):
                P[i,k] = P[i,k]-beta*U[i,j]
    D[0,0] = P[0,0]
    for i in range(n):
        U[i,i] = 1

    return U,D
