#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

#ifndef SQ
#define SQ(x) ((x)*(x)) // 8008135!!!
#endif

using Eigen::Matrix;
using Eigen::Vector3f;
using namespace std;

typedef Eigen::Matrix<float,16,1> Vector16f;
typedef Eigen::Matrix<float,1,16> RowVector16f;
typedef Eigen::Matrix<float,16,16> Matrix16f;
typedef Eigen::DiagonalMatrix<float,16> DiagonalMatrix16f;

typedef struct {
    float t;
    float range;
    uint8_t anchor_idx;
    uint8_t tag_idx;
} Observation;

Vector3f positions[] = {
    {0,0,0},
    {0,1,0},
    {1,0,0},
    {0,1,0},
    {1,0,0},

    {0,0,0},
    {0,1,0},
    {1,0,0},
    {0,-1,0},
    {-1,0,0},
};


static float calc_objective(const Vector16f& state, const size_t num_obs, const Observation* obs, float range_sigma, float ant_del_sigma);
static float calc_range_resid(const Vector16f& state, const Observation& obs);
static void calc_range_jacob(const Vector16f& state, const Observation& obs, RowVector16f& ret);
static void calc_JTWJ_and_b(const Vector16f& state, const size_t num_obs, const Observation* obs, float range_sigma, float ant_del_sigma, Matrix16f& JTWJ, Vector16f& b);
static bool _multilaterate(Vector16f& state, Matrix16f& covariance, const size_t num_obs, const Observation* obs, float range_sigma, float ant_del_sigma);

bool print;
Vector16f best_state;

Observation observations[] = {
    {0, 100, 0, 5},
    {0, 100, 1, 5},
    {0, 100, 2, 5},
    {0, 100.05, 3, 5},

    {0.5, 100, 0, 5},
    {0.5, 100, 1, 5},
    {0.5, 100, 2, 5},

    {0, 100, 0, 5},
    {0, 100, 1, 5},
    {0, 100, 2, 5},

    {0.5, 100, 0, 5},
    {0.5, 100, 1, 5},
    {0.5, 100, 2, 5},

    {0, 100, 0, 5},
    {0, 100, 1, 5},
    {0, 100, 2, 5},

    {0.5, 100, 0, 5},
    {0.5, 100, 1, 5},
    {0.5, 100, 2, 5},

    {0, 100, 0, 5},
    {0, 100, 1, 5},
    {0, 100, 2, 5},

    {0.5, 5, 0, 5},
    {0.5, 5, 1, 5},
    {0.5, 5, 2, 5},
};

int main(void) {

    std::default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<float> distribution(0.0,0.05);

    for (unsigned i=0; i<sizeof(observations)/sizeof(observations[0]); i++) {
        observations[i].range += distribution(generator);
    }

    Vector16f state;
    Matrix16f covariance;

    unsigned N = 1;

    using std::chrono::steady_clock;
    auto start = steady_clock::now();
    for (unsigned i=0; i<N; i++) {
//         print = false;
//         best_state << 0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0;
//         _multilaterate(best_state, covariance, sizeof(observations)/sizeof(observations[0]), observations, 0.3);
//         print = true;

        state << 0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0;
        if(_multilaterate(state, covariance, sizeof(observations)/sizeof(observations[0]), observations, 0.03, 0.1)) {

            cout << state << endl << endl;
            cout << covariance << endl << endl;
        } else {
            cout << "fail" << endl;
        }
    }
    auto end = steady_clock::now();

    double elapsedSeconds = ((end - start).count()) * steady_clock::period::num / static_cast<double>(steady_clock::period::den);
    cout << "time: " << 1e6*elapsedSeconds/(double)N << " us" << endl;

    return 0;
}

static bool _multilaterate(Vector16f& state, Matrix16f& covariance, const size_t num_obs, const Observation* obs, float range_sigma, float ant_del_sigma) {
    state << 0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0;

    const float v_factor = 10;
    const float epsilon_gradient = 0.1f;
    float lambda = 1e-4f;
    float fval = calc_objective(state, num_obs, obs, range_sigma, ant_del_sigma);
    uint8_t lambda_incr_count = 0;

    float gradient;
    Vector16f h;
    Matrix16f JTWJ;
    Vector16f b;

    for (uint32_t i=0; i<50; i++) {
        for (uint32_t j=0; j<1; j++) {
            calc_JTWJ_and_b(state, num_obs, obs, range_sigma, ant_del_sigma, JTWJ, b);

            Matrix16f A1 = JTWJ;
            A1 += (lambda * JTWJ.diagonal()).asDiagonal();
            Matrix16f A2 = JTWJ;
            A2 += ((lambda/v_factor) * JTWJ.diagonal()).asDiagonal();

            Eigen::LLT<Matrix16f> llt = A1.llt();
            Vector16f h1 = llt.solve(b);
            bool fit1_valid = (llt.info() == Eigen::Success) && h1.allFinite();

            llt = A2.llt();
            Vector16f h2 = llt.solve(b);
            bool fit2_valid = (llt.info() == Eigen::Success) && h2.allFinite();

            float fit1 = calc_objective(state-h1, num_obs, obs, range_sigma, ant_del_sigma);
            float fit2 = calc_objective(state-h2, num_obs, obs, range_sigma, ant_del_sigma);

            if (fit1_valid && fit2_valid && fit2 < fval && fit2 < fit1) {
                lambda_incr_count=0;
                lambda /= v_factor;
                fval = fit2;
                state -= h2;
            } else if (fit1_valid && fit1 < fval) {
                lambda_incr_count=0;
                fval = fit1;
                state -= h1;
            } else {
                lambda_incr_count++;
                lambda *= v_factor;
            }
        }

        if (lambda_incr_count >= 6) {
            Eigen::LDLT<Matrix16f> llt = JTWJ.ldlt();
            covariance = llt.solve(JTWJ.Identity());

            return (llt.info() == Eigen::Success) && covariance.allFinite();
        }
    }

    return false;
}

static void calc_JTWJ_and_b(const Vector16f& state, const size_t num_obs, const Observation* obs, float range_sigma, float ant_del_sigma, Matrix16f& JTWJ, Vector16f& b) {
    JTWJ = JTWJ.Zero();
    b = b.Zero();

    for (size_t k=0; k<num_obs; k++) {
        RowVector16f J;
        calc_range_jacob(state, obs[k], J);

        float resid = calc_range_resid(state, obs[k]);

        // Compute huber weight
        float k_huber = 1.345f * range_sigma;
        float w_huber = fabsf(resid) <= k_huber ? 1.0f : k_huber/fabsf(resid);

        Vector16f JTW = J.transpose() * w_huber/SQ(range_sigma);
        JTWJ.noalias() += JTW*J;
        b.noalias() += JTW * resid;
    }

    for (size_t i=6; i<16; i++) {
        RowVector16f J;
        J = J.Zero();
        J(0,i) = 1;

        float resid = state(i,0);

        Vector16f JTW = J.transpose() * 1/SQ(ant_del_sigma);
        JTWJ.noalias() += JTW*J;
        b.noalias() += JTW * resid;
    }
}

static float calc_objective(const Vector16f& state, const size_t num_obs, const Observation* obs, float range_sigma, float ant_del_sigma) {
    float ret = 0;
    for (size_t i=0; i<num_obs; i++) {
        float resid = calc_range_resid(state, obs[i]);
        float k_huber = 1.345f * range_sigma;
        float rho = (fabsf(resid) <= k_huber) ? (0.5*SQ(resid)) : (k_huber * fabsf(resid) - 0.5 * SQ(k_huber));
        ret += rho/SQ(range_sigma);
    }

    for (size_t i=6; i<16; i++) {
        float resid = state(i,0);
        float rho = 0.5*SQ(resid);
        ret += rho/SQ(ant_del_sigma);
    }
    return ret;
}


static float calc_range_resid(const Vector16f& state, const Observation& obs) {
    return obs.range + state(obs.tag_idx+6,0) + state(obs.anchor_idx+6,0) - sqrtf(((-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))*(-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))) + ((-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))*(-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))) + ((-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))*(-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))));
}

static void calc_range_jacob(const Vector16f& state, const Observation& obs, RowVector16f& ret) {
    ret(0,0) = -(obs.t*state(3,0) - positions[obs.anchor_idx][0] + positions[obs.tag_idx][0] + state(0,0))/sqrtf(((-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))*(-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))) + ((-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))*(-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))) + ((-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))*(-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))));

    ret(0,1) = -(obs.t*state(4,0) - positions[obs.anchor_idx][1] + positions[obs.tag_idx][1] + state(1,0))/sqrtf(((-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))*(-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))) + ((-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))*(-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))) + ((-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))*(-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))));

    ret(0,2) = -(obs.t*state(5,0) - positions[obs.anchor_idx][2] + positions[obs.tag_idx][2] + state(2,0))/sqrtf(((-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))*(-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))) + ((-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))*(-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))) + ((-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))*(-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))));

    ret(0,3) = obs.t*(-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))/sqrtf(((-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))*(-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))) + ((-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))*(-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))) + ((-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))*(-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))));

    ret(0,4) = obs.t*(-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))/sqrtf(((-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))*(-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))) + ((-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))*(-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))) + ((-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))*(-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))));

    ret(0,5) = obs.t*(-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))/sqrtf(((-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))*(-obs.t*state(3,0) + positions[obs.anchor_idx][0] - positions[obs.tag_idx][0] - state(0,0))) + ((-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))*(-obs.t*state(4,0) + positions[obs.anchor_idx][1] - positions[obs.tag_idx][1] - state(1,0))) + ((-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))*(-obs.t*state(5,0) + positions[obs.anchor_idx][2] - positions[obs.tag_idx][2] - state(2,0))));

    for (uint8_t i=6; i<16; i++) {
        ret(0,i) = 0;
    }

    ret(0,obs.tag_idx+6) = 1;
    ret(0,obs.anchor_idx+6) = 1;
}
