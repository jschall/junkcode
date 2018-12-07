#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#define SQ(x) ((x)*(x))

static void connector_ekf_predict(float* state, float* cov, const float R_pin_pnoise, const float V_drop_bias_pnoise, const float dt) {
    float state_ret[2];
    float cov_ret[3];
    state_ret[0] = state[0];
    state_ret[1] = state[1];
    cov_ret[0] = ((R_pin_pnoise)*(R_pin_pnoise))*((dt)*(dt)) + cov[0];
    cov_ret[1] = cov[1];
    cov_ret[2] = ((V_drop_bias_pnoise)*(V_drop_bias_pnoise))*((dt)*(dt)) + cov[2];

    memcpy(state, state_ret, sizeof(state_ret));
    memcpy(cov, cov_ret, sizeof(cov_ret));

}

static bool connector_ekf_update(float* state, float* cov, const float V_drop_meas, const float V_drop_meas_noise, const float I_meas, const float NIS_thresh) {
    float subx0 = -I_meas*state[0] + V_drop_meas - state[1];
    float subx1 = I_meas*cov[0] + cov[1];
    float subx2 = I_meas*subx1;
    float subx3 = I_meas*cov[1];
    float subx4 = cov[2] + subx3;
    float subx5 = 1/(((V_drop_meas_noise)*(V_drop_meas_noise)) + subx2 + subx4);
    float subx6 = subx1*subx5;
    float subx7 = subx4*subx5;
    float subx8 = -subx2*subx5 + 1;

    float NIS = ((subx0)*(subx0))*subx5;
    printf("%f\n", NIS);
    if (NIS > NIS_thresh) return false;

    float state_ret[2];
    float cov_ret[3];

    state_ret[0] = state[0] + subx0*subx6;
    state_ret[1] = state[1] + subx0*subx7;
    cov_ret[0] = cov[0]*subx8 - cov[1]*subx6;
    cov_ret[1] = cov[1]*subx8 - cov[2]*subx6;
    cov_ret[2] = cov[2]*(-subx7 + 1) - subx3*subx7;

    memcpy(state, state_ret, sizeof(state_ret));
    memcpy(cov, cov_ret, sizeof(cov_ret));

    return true;
}


int main(void) {
    // Note: re-initialize estimator when battery is initially inserted
    const float init_R_pin = 1e-3; // Initial guess of pin resistance
    const float init_R_pin_variance = SQ(20e-3); // Indicates init_R_pin has a 68% confidence interval of +/- 10e-3 Ohms.
    const float init_V_drop_bias = 0; // Initial guess of voltage drop measurement bias
    const float init_V_drop_bias_variance = SQ(1e-3); // Indicates init_V_drop_bias has a 68% confidence interval of +/- 1e-3 V.

    float x[2] = {init_R_pin, init_V_drop_bias};
    float P[3] = {init_V_drop_bias, 0, init_V_drop_bias_variance};

    for (float t=0; t<10; t += 0.001) {
        const float I = (1.0+cosf(t))*20;

        printf("\n\n%f\n", I);

        const float R_pin_pnoise = 5e-3; // Indicates that the rate of change of R_pin is 0 Ohms/sec with a 68% confidence interval of 5e-3 Ohms/sec;
        const float V_drop_bias_pnoise = 1e-4; // Indicates that the rate of change of V_drop_bias is 0 V/sec with a 68% confidence interval of 1e-4 V/sec;
        const float dt = 1e-3; // Indicates that 1e-3 seconds have passed since the last iteration of the filter

        connector_ekf_predict(x, P, R_pin_pnoise, V_drop_bias_pnoise, dt);
        printf("%f %f %f %f %f\n", x[0], x[1], P[0], P[1], P[2]);

        const float V_drop_meas = I*1e-3; // Measured voltage drop across the pin in V
        const float V_drop_meas_noise = 1e-2; // Indicates V_drop_meas has a 68% confidence interval of 1e-2 V.
        const float I_meas = I; // Measured current in A.
        const float NIS_thresh = 3; // Indicates that the consistency check threshold should be 3 sigma, or roughly 0.3% of non-outliers rejected.

        connector_ekf_update(x, P, V_drop_meas, V_drop_meas_noise, I_meas, NIS_thresh);

        printf("%f %f %f %f %f\n", x[0], x[1], sqrtf(P[0]), P[1], sqrtf(P[2]));
    }

    return 0;
}
