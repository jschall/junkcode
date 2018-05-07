import numpy as np
from muppettree import muppettree_fit
import random

np.set_printoptions(suppress=True)

truth_pos = (1.,0.,-50.)

gen_obs = lambda t,pos,sigma: (t, ((pos[0]-truth_pos[0])**2.+(pos[1]-truth_pos[1])**2.+(pos[2]-truth_pos[2])**2.)**0.5 + random.gauss(0,sigma), pos)

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
    ]

#for obs in observations:
    #print(obs)

#N = len(observations)

print(muppettree_fit(observations))
