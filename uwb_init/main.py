import numpy as np
from muppettree import muppettree_fit
from nonlinear import nonlinear_fit
import random

np.set_printoptions(suppress=True)


truth_vel = np.array([10.,0.,-3.])
truth_pos = lambda t: np.array([1.,0.,-50.]) + t*truth_vel

gen_obs = lambda t,pos,sigma: (t, ((pos[0]-truth_pos(t)[0])**2.+(pos[1]-truth_pos(t)[1])**2.+(pos[2]-truth_pos(t)[2])**2.)**0.5 + random.gauss(0,sigma), pos)

observations = [
    #t   r      pos
    gen_obs(-0.00, (-1., -1., 0.),0.1),
    gen_obs(-0.01, (-1.,  1., 0.),0.1),
    gen_obs(-0.02, ( 1., -1., 0.),0.1),
    gen_obs(-0.03, ( 0.,  0., 1.),0.1),

    gen_obs(-1.00, (-1., -1., 0.),0.1),
    gen_obs(-1.01, (-1.,  1., 0.),0.1),
    gen_obs(-1.02, ( 1., -1., 0.),0.1),
    gen_obs(-1.03, ( 0.,  0., 1.),0.1),

    gen_obs(-2.00, (-1., -1., 0.),0.1),
    gen_obs(-2.01, (-1.,  1., 0.),0.1),
    gen_obs(-2.02, ( 1., -1., 0.),0.1),
    gen_obs(-2.03, ( 0.,  0., 1.),0.1),

    gen_obs(-0.00, (-1., -1., 0.),0.1),
    gen_obs(-0.01, (-1.,  1., 0.),0.1),
    gen_obs(-0.02, ( 1., -1., 0.),0.1),
    gen_obs(-0.03, ( 0.,  0., 1.),0.1),

    gen_obs(-1.00, (-1., -1., 0.),0.1),
    gen_obs(-1.01, (-1.,  1., 0.),0.1),
    gen_obs(-1.02, ( 1., -1., 0.),0.1),
    gen_obs(-1.03, ( 0.,  0., 1.),0.1),

    gen_obs(-2.00, (-1., -1., 0.),0.1),
    gen_obs(-2.01, (-1.,  1., 0.),0.1),
    gen_obs(-2.02, ( 1., -1., 0.),0.1),
    gen_obs(-2.03, ( 0.,  0., 1.),0.1),

    gen_obs(-0.00, (-1., -1., 0.),0.1),
    gen_obs(-0.01, (-1.,  1., 0.),0.1),
    gen_obs(-0.02, ( 1., -1., 0.),0.1),
    gen_obs(-0.03, ( 0.,  0., 1.),0.1),

    gen_obs(-1.00, (-1., -1., 0.),0.1),
    gen_obs(-1.01, (-1.,  1., 0.),0.1),
    gen_obs(-1.02, ( 1., -1., 0.),0.1),
    gen_obs(-1.03, ( 0.,  0., 1.),0.1),

    gen_obs(-2.00, (-1., -1., 0.),0.1),
    gen_obs(-2.01, (-1.,  1., 0.),0.1),
    gen_obs(-2.02, ( 1., -1., 0.),0.1),
    gen_obs(-2.03, ( 0.,  0., 1.),0.1),

    gen_obs(-0.00, (-1., -1., 0.),0.1),
    gen_obs(-0.01, (-1.,  1., 0.),0.1),
    gen_obs(-0.02, ( 1., -1., 0.),0.1),
    gen_obs(-0.03, ( 0.,  0., 1.),0.1),

    gen_obs(-1.00, (-1., -1., 0.),0.1),
    gen_obs(-1.01, (-1.,  1., 0.),0.1),
    gen_obs(-1.02, ( 1., -1., 0.),0.1),
    gen_obs(-1.03, ( 0.,  0., 1.),0.1),

    gen_obs(-2.00, (-1., -1., 0.),0.1),
    gen_obs(-2.01, (-1.,  1., 0.),0.1),
    gen_obs(-2.02, ( 1., -1., 0.),0.1),
    gen_obs(-2.03, ( 0.,  0., 1.),0.1),
    ]

#for obs in observations:
    #print(obs)

#N = len(observations)

muppettree_result = muppettree_fit(observations)[0].flatten()[0:6]
nonlinear_result = nonlinear_fit(observations, muppettree_result)
print nonlinear_result.
print muppettree_result

print(muppettree_result-nonlinear_result.x)
