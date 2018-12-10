#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <gsl/gsl_matrix.h>
#include <stdbool.h>

struct BodyParameters {
    // material properties
    double E;
    double G;
    double p; //density

    //geometry (assumed constant for now)
    double J;
    double I;
    double A;
    double R;
    double L;

    //reference xi
    gsl_vector *xi_ref;
};

struct SystemParameters {
    double dt; //time step
    double ds; //spatial discretization
    int n; // the number of discretizations
    int N; //the number of actuators
    //the initial conditions
    gsl_matrix *g0;
    gsl_vector *eta0;
    gsl_vector *xi0;

    //actuation
    gsl_vector *q;
    gsl_vector *q_dot;

    //flags
    bool delta; //compute delta in the integration

    gsl_matrix *extra; //i don't know if this is the best choice, but it is what I'm using for now
};


struct ActuatorParameters {
    double (*rx)(double R, int i, int N, double s); //the function to determine the displacement of the actuator, ACT# -> s -> r
    double (*ry)(double R, int i, int N, double s);
    double (*force)(struct SystemParameters sys_params, int i, int j, double s); //the function to determine the force, input -> ACT# -> POS# -> s -> r
};

struct SimulationParameters {
    struct BodyParameters body_params;
    struct SystemParameters sys_params;
    struct ActuatorParameters act_params;

    gsl_matrix *g;
    gsl_matrix *xi;
    gsl_matrix *xi_prev;
    gsl_matrix *eta_prev;
    gsl_matrix *eta;
};


void stepDynamics(struct SimulationParameters sim_params);

double tcaForce(struct SystemParameters sys_params, int i, int j, double s);
double cableForce(struct SystemParameters sys_params, int i, int j, double s);
double rx_fun(double R, int i, int N, double s);
double ry_fun(double R, int i, int N, double s);

#endif
