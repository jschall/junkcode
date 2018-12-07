#!/usr/bin/python
from sympy import *
from sympy.printing.ccode import *
from math import sqrt, floor, fmod
import math
from helpers import *
import sys

V_drop_meas = Symbol('V_drop_meas', real=True)
I_meas = Symbol('I_meas', real=True)
dt = Symbol('dt', real=True, nonnegative=True)

# Tuning parameters
V_drop_meas_noise = Symbol('V_drop_meas_noise', real=True, nonnegative=True)
R_pin_pnoise = Symbol('R_pin_pnoise', real=True, nonnegative=True)
V_drop_bias_pnoise = Symbol('V_drop_bias_pnoise', real=True, nonnegative=True)

# States
R_pin = Symbol('state[0]', real=True)
V_drop_bias = Symbol('state[1]', real=True)

x = Matrix([R_pin, V_drop_bias])
nStates = len(x)

# Covariance matrix
P = compressedSymmetricMatrix('cov', nStates)

# f: state-transtition model
f = x

assert f.shape == x.shape

# F: linearized state-transition model, AKA "A" in literature
F = f.jacobian(x)

# Q: covariance of additive noise on x
Q = diag((R_pin_pnoise*dt)**2, (V_drop_bias_pnoise*dt)**2)

# x_p: state vector at time k+1
x_p = x

# P_p: covariance matrix at time k+1
P_p = F*P*F.T + Q
assert P_p.shape == P.shape
P_p = upperTriangularToVec(P_p)

print('static void connector_ekf_predict(float* state, float* cov, const float R_pin_pnoise, const float V_drop_bias_pnoise, const float dt) {')

print('    float state_ret[2];')
print('    float cov_ret[3];')

for i in range(len(x_p)):
    print('    state_ret[%u] = %s;' % (i, CCodePrinter_float().doprint(x_p[i])))

for i in range(len(P_p)):
    print('    cov_ret[%u] = %s;' % (i, CCodePrinter_float().doprint(P_p[i])))

print('')
print('    memcpy(state, state_ret, sizeof(state_ret));')
print('    memcpy(cov, cov_ret, sizeof(cov_ret));\n')

print('}\n')


# h: predicted measurement
h = Matrix([I_meas*R_pin+V_drop_bias])

# z: observation
z = Matrix([V_drop_meas])

# R: observation covariance
R = Matrix([V_drop_meas_noise**2])

# y: innovation vector
y = z-h

# H: measurement sensitivity matrix
H = h.jacobian(x)

# S: innovation covariance
S = H*P*H.T + R

S_I = quickinv_sym(S)

# K: Kalman gain
K = P*H.T*S_I

I = eye(nStates)

NIS = y.T*S_I*y # normalized innovation squared

x_n = x + K*y
P_n = upperTriangularToVec((I-K*H)*P)

x_n,P_n,NIS,y,subx = extractSubexpressions([x_n,P_n,NIS,y],'subx',threshold=1)

print('static bool connector_ekf_update(float* state, float* cov, const float V_drop_meas, const float V_drop_meas_noise, const float I_meas, const float NIS_thresh) {')

for i in range(len(subx)):
    print('    float %s = %s;' % (subx[i][0], CCodePrinter_float().doprint(subx[i][1])))

print('')
print('    float NIS = %s;' % (CCodePrinter_float().doprint(NIS[0]),))
print('    if (NIS > NIS_thresh) return false;\n')

print('    float state_ret[2];')
print('    float cov_ret[3];\n')

for i in range(len(x_n)):
    print('    state_ret[%u] = %s;' % (i, CCodePrinter_float().doprint(x_n[i])))

for i in range(len(P_n)):
    print('    cov_ret[%u] = %s;' % (i, CCodePrinter_float().doprint(P_n[i])))

print('')
print('    memcpy(state, state_ret, sizeof(state_ret));')
print('    memcpy(cov, cov_ret, sizeof(cov_ret));\n')

print('    return true;')
print('}')

#state_names = ['QUAT_W', 'QUAT_X', 'QUAT_Y', 'QUAT_Z', 'GBIAS_X', 'GBIAS_Y', 'GBIAS_Z', 'ABIAS_X', 'ABIAS_Y', 'ABIAS_Z', 'POS_N', 'POS_E', 'POS_D', 'VEL_N', 'VEL_E', 'VEL_D', 'ADEL_0', 'ADEL_1', 'ADEL_2', 'ADEL_3', 'ADEL_4', 'ADEL_5', 'ADEL_6', 'ADEL_7', 'ADEL_8', 'ADEL_9']

#assert len(state_names) == nStates

#with open(sys.argv[1],'wb') as f:
    #f.write('#pragma once\n')
    #f.write('#include <math.h>\n\n')

    #f.write('#define EKF_NUM_STATES %u\n' % (nStates,))
    #f.write('#define EKF_NUM_COV %u\n\n' % ((nStates**2-nStates)/2+nStates,))

    #for i in range(len(state_names)):
        #f.write('#define STATE_IDX_%s %u\n' % (state_names[i],i))

    #f.write('\n')

    #f.write(generatePrediction())
    #f.write('\n')
    #f.write(generateRangeFusion())
    #f.write('\n')
    #f.write(generatePixyFusion())
    #f.write('\n')




