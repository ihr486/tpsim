#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct params_tag
{
    double F;   //PWM base clock frequency
    double P;   //PWM period in clock cycles
    double E;   //Power supply voltage
    double Rf, Cf;  //Filter parameters
    double R, L;    //Coil parameters
    double Rs;  //Current sensor resistance
    long N;     //Number of simulated steps
    double Kp;  //Proportional coefficient
    double Ki;  //Integral coefficient
    double Kd;  //Differential coefficient
} params_t;

typedef struct state_tag
{
    double I[2];
    double Vf[3];
    double Ir[2];
} state_t;

static const double PI = 3.1415926;

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

static void pid_control(const params_t *params, const state_t *state, int D[3])
{
    double V0 = params->E * 0.5;
    double V3 = 3.0 * V0;
    for (int i = 0; i < 2; i++)
    {
        double Is = (state->Vf[i] - state->Vf[2]) / params->Rs;
        double V = params->Kp * (state->Ir[i] - Is);
        if (V + V0 > params->E) V = params->E - V0;
        if (V + V0 < 0) V = -V0;
        D[i] = params->P * (V + V0) / params->E;
        if (D[i] > params->P) D[i] = params->P;
        if (D[i] < 0) D[i] = 0;
        V3 -= V;
    }
    D[2] = params->P * V3 / params->E;
}

static const int *sorted_index(int v[3])
{
    static const int table[8][4] = {
        {-1, -1, -1},   //0 < 1 && 1 < 2 && 2 < 0 (impossible)
        {1, 2, 0, 3},      //0 > 1 && 1 < 2 && 2 < 0
        {2, 0, 1, 3},      //0 < 1 && 1 > 2 && 2 < 0
        {2, 1, 0, 3},      //0 > 1 && 1 > 2 && 2 < 0
        {0, 1, 2, 3},      //0 < 1 && 1 < 2 && 2 > 0
        {1, 0, 2, 3},      //0 > 1 && 1 < 2 && 2 > 0
        {0, 2, 1, 3},      //0 < 1 && 1 > 2 && 2 > 0
        {-1, -1, -1}    //0 > 1 && 1 > 2 && 2 > 0 (impossible)
    };
    int mask = 0;
    if (v[0] > v[1]) mask |= 0x01;
    if (v[1] > v[2]) mask |= 0x02;
    if (v[2] > v[0]) mask |= 0x04;
    return table[mask];
}

static double simulate_current(const params_t *params, double V, double I0, double t)
{
    return V / params->R + (I0 - V / params->R) * exp(-params->R / params->L * t);
}

static double simulate_sense_filter(const params_t *params, double V, double V0, double I0, double Vf0, double t)
{
    return (V0 + params->Rs / params->R * V) * (-expm1(-t / (params->Rf * params->Cf))) + \
          (params->L * params->Rs) / (params->L - params->R * params->Rf * params->Cf) * (I0 - V / params->R) * (exp(-params->R / params->L * t) - exp(-t / (params->Rf * params->Cf))) + \
          Vf0 * exp(-t / (params->Rf * params->Cf));
}

static double simulate_centroid_filter(const params_t *params, double V0, double Vf0, double t)
{
    return V0 + (Vf0 - V0) * exp(-t / (params->Rf * params->Cf));
}

static void update_state(const params_t *params, state_t *state, double E[3], double t)
{
    double V0 = (E[0] + E[1] + E[2]) / 3.0;
    for (int i = 0; i < 2; i++)
    {
        state->Vf[i] = simulate_sense_filter(params, E[i] - V0, V0, state->I[i], state->Vf[i], t);
        state->I[i] = simulate_current(params, E[i] - V0, state->I[i], t);
    }
    state->Vf[2] = simulate_centroid_filter(params, V0, state->Vf[2], t);
}

int main(int argc, const char *argv[])
{
    //printf("Three-phase coil driver simulator\n");

    const params_t params = {
        .F = 48E+6,
        .E = 5,
        .Rf = 10E+3,
        .Cf = 1000E-12,
        .R = 1.68E-8 * 180 * PI * 15E-3 / (PI * 0.2E-3 * 0.2E-3) + 1.24,
        .L = 150E-6,
        .Rs = 0.24,
        .P = 512,
        .N = 100000,
        .Kp = 0.1,
    };

    state_t state = {
        .I = {0, 0},
        .Vf = {0, 0, 0},
        .Ir = {0, 0}
    };

    for (long i = 0; i < params.N; i++)
    {
        double t0 = 1.0 / params.F * params.P * i;

        int D[4] = {0, 0, 0, params.P};
        pid_control(&params, &state, D);
        double E[3] = {params.E, params.E, params.E};

        const int *sorted = sorted_index(D);

        double t1 = 1.0 / params.F * D[sorted[0]];
        update_state(&params, &state, E, t1);

        for (int j = 0; j < 3; j++)
        {
            E[sorted[j]] = 0;
            t1 = 1.0 / params.F * (D[sorted[j + 1]] - D[sorted[j]]);
            update_state(&params, &state, E, t1);
        }

        printf("%lf %lf %lf %lf\n", t0, state.I[0] * 1E+3, state.I[1] * 1E+3, -(state.I[0] + state.I[1]) * 1E+3);
        //printf("%lf %lf %lf %lf\n", t0, state.Vf[0], state.Vf[1], state.Vf[2]);
    }

    return 0;
}
