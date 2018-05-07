from sympy import *
import random

truth_pos = (0.,0.,-10.)

gen_obs = lambda t,pos: (t, ((pos[0]-truth_pos[0])**2.+(pos[1]-truth_pos[1])**2.+(pos[2]-truth_pos[2])**2.)**0.5, pos)

observations = [
    #t   r      pos
    gen_obs(-0.00, (-1., -1., 0.)),
    gen_obs(-0.01, (-1.,  1., 0.)),
    gen_obs(-0.02, ( 1., -1., 0.)),
    gen_obs(-0.03, ( 0.,  0., 1.)),
    
    gen_obs(-1.00, (-1., -1., 0.)),
    gen_obs(-1.01, (-1.,  1., 0.)),
    gen_obs(-1.02, ( 1., -1., 0.)),
    gen_obs(-1.03, ( 0.,  0., 1.)),
    
    gen_obs(-2.00, (-1., -1., 0.)),
    gen_obs(-2.01, (-1.,  1., 0.)),
    gen_obs(-2.02, ( 1., -1., 0.)),
    gen_obs(-2.03, ( 0.,  0., 1.)),
    ]

for obs in observations:
    print(obs)

N = len(observations)

p = Matrix(symbols("p1:4"))
v = Matrix(symbols("v1:4"))

ri = [Symbol("r%u" % i) for i in range(N)]
ti = [Symbol("t%u" % i) for i in range(N)]
pi = [Matrix(symbols("p%u(1:4)" % i)) for i in range(N)]


pdotp = Symbol("pdotp")
pdotv = Symbol("pdotv")
vdotv = Symbol("vdotv")

x = Matrix([p, v, Matrix([pdotp, pdotv, vdotv])])

eqns = [Eq(ri[i]**2, (pdotp + ti[i]**2*vdotv + 2.*ti[i]*pdotv - 2.*pi[i].dot(p) - 2.*(ti[i]*pi[i]).dot(v) + pi[i].dot(pi[i]))) for i in range(N)]

expr = pdotp + ti[i]**2*vdotv + 2.*ti[i]*pdotv - 2.*pi[i].dot(p) - 2.*(ti[i]*pi[i]).dot(v) + pi[i].dot(pi[i])

pprint(simplify(expr.subs({pdotp: p.dot(p), pdotv: p.dot(v), vdotv: v.dot(v)})))

A, b = linear_eq_to_matrix(eqns, list(x))

print(latex(Eq(Symbol("A"), A, evaluate=False), mode='inline'))

print(latex(Eq(Symbol("b"), b, evaluate=False), mode='inline'))

for i in range(N):
    A = A.xreplace({ti[i]: observations[i][0]})
    A = A.xreplace({ri[i]: observations[i][1]})
    A = A.xreplace(dict(zip(pi[i], observations[i][2])))
    
    b = b.xreplace({ti[i]: observations[i][0]})
    b = b.xreplace({ri[i]: observations[i][1]})
    b = b.xreplace(dict(zip(pi[i], observations[i][2])))

#pprint(A)
#pprint(b)

#pprint(A)

#pprint(A.T*A)

print((A.T*A).det())

pprint((A.T*A).pinv()*A.T*b)
