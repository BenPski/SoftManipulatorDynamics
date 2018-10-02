/*
Run the dynamics, mostly a test

For some reason only compiles with g++, gcc gives a linker error of not being able to find force and r functions
*/

#include <gsl/gsl_matrix.h>
#include <math.h>

#include "dynamics.h"

int main(int argc, char *argv[]) {
    
    int n = atoi(argv[1]);

    //setup for system
    //const int n = 5;
    const int N = 3;
    double L = 0.04;
    double p = 1000;
    double E = 37.8e3;
    double G = E/3;
    double D = 1e-2;
    double R = D/2;
    double A = M_PI*R*R;
    double I = M_PI*pow(D,4)/32;
    double J = 2*I;

    double dt = 0.01;
    double ds = L/(n-1);
    gsl_matrix *g0 = gsl_matrix_alloc(4,4);
    gsl_matrix_set_identity(g0);
    gsl_vector *xi0 = gsl_vector_alloc(6);
    gsl_vector_set_zero(xi0);
    gsl_vector_set(xi0,5,1);
    gsl_vector *xi_ref = gsl_vector_alloc(6);
    gsl_vector_set_zero(xi_ref);
    gsl_vector_set(xi_ref,5,1);
    gsl_vector *eta0 = gsl_vector_alloc(6);
    gsl_vector_set_zero(eta0);

    gsl_vector *q = gsl_vector_calloc(N);
    gsl_vector_set(q,0,50);
    gsl_vector_set(q,1,0);
    gsl_vector_set(q,2,0);


    gsl_matrix *extra = gsl_matrix_calloc(n,N);
    
    struct SystemParameters sys_params = { dt, ds, n, N, g0, eta0, xi0, q, false, extra};
    struct BodyParameters body_params = { E, G, p, J, I, A, R, L, xi_ref};
    struct ActuatorParameters act_params;
    act_params.rx = rx_fun;
    act_params.ry = ry_fun;
    act_params.force = cableForce;
    //act_params.force = tcaForce;

    gsl_matrix *eta = gsl_matrix_alloc(n,6);
    gsl_matrix *eta_prev = gsl_matrix_alloc(n,6);
    gsl_matrix *xi_prev = gsl_matrix_alloc(n,6);
    gsl_matrix *g = gsl_matrix_alloc(n,12);
    gsl_matrix *xi = gsl_matrix_alloc(n,6);

    gsl_matrix_set_zero(eta);
    gsl_matrix_set_zero(eta_prev);
    gsl_matrix_set_zero(xi_prev);
    gsl_matrix_set_zero(g);
    gsl_matrix_set_zero(xi);
    int j;
    for (j=0;j<n;j++) {
        gsl_matrix_set(xi_prev,j,5,1);
    }

    struct SimulationParameters params = {body_params, sys_params, act_params, g, xi, xi_prev, eta_prev, eta};

    for (j=0;j<100;j++) {
        stepDynamics(params);
        gsl_matrix_memcpy(params.xi_prev,params.xi);
        gsl_matrix_memcpy(params.eta_prev,params.eta);
    }

    gsl_matrix_free(extra);
    
    gsl_matrix_free(g0);
    gsl_vector_free(xi0);
    gsl_vector_free(xi_ref);
    gsl_vector_free(eta0);

    gsl_vector_free(q);

    gsl_matrix_free(eta);
    gsl_matrix_free(eta_prev);
    gsl_matrix_free(xi);
    gsl_matrix_free(xi_prev);
    gsl_matrix_free(g);

    return 0;
}