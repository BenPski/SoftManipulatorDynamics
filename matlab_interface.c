#include "mex.h"
#include "dynamics.h"
#include <math.h>
#include <gsl/gsl_matrix.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//provide the interface for running the dynamics with matlab
	//this should only provide stepping the dynamics, can be considered not efficient and that is correct
	//though it seems that maintaining the state would be cumbersome without just passing the data back and forth
	//if it is just passing pointers then the only overhead is reinitializing the data, which is a minor slow down

	//this should be able to switch between tca and cable via a flag


	//the data to pass in is the current values of xi, eta, q, dt, actuator flag
	//want to return g,xi,eta since the equations are solved in c now

	//let the actuator flag be an int for simplicity
		//0 -> cable
		//1 -> tca

	//can specify more, but for now only make some defaults

    int i,j;

    //the actuator flag
    int act_flag;
    act_flag = mxGetPr(prhs[0])[0];

	double *q_in;
	const mwSize *dimq;
    int dimqx, dimqy;
    int N; /* number of actuators */
    dimq = mxGetDimensions(prhs[1]);
    dimqy = (int)dimq[0]; dimqx = (int)dimq[1];
    N = dimqy*dimqx;
    q_in = mxGetPr(prhs[1]);

    //setting the q vector
    gsl_vector *q = gsl_vector_calloc(N);
    for (i=0;i<N;i++) {
    	gsl_vector_set(q,i,q_in[i]);
    }
    double *q_dot_in = mxGetPr(prhs[2]);
    gsl_vector *q_dot = gsl_vector_calloc(N);
    for (i=0;i<N;i++) {
        gsl_vector_set(q_dot,i,q_dot_in[i]);
    }

    //set eta_prev
	int n;
	const mwSize *dimEta;
    dimEta = mxGetDimensions(prhs[3]);
    n = (int)dimEta[0];
    double *eta_prev_trans;
    eta_prev_trans = mxGetPr(prhs[3]);
	gsl_matrix *eta_prev = gsl_matrix_alloc(n,6);
    for (i=0;i<n;i++) {
        for (j=0;j<6;j++) {
        	gsl_matrix_set(eta_prev,i,j,eta_prev_trans[j*n+i]);
        }
    }

    //set xi_prev
    gsl_matrix *xi_prev = gsl_matrix_alloc(n,6);
    double *xi_prev_trans;
    xi_prev_trans = mxGetPr(prhs[4]);

    for (i=0;i<n;i++) {
        for (j=0;j<6;j++) {
        	gsl_matrix_set(xi_prev,i,j,xi_prev_trans[j*n+i]);
        }
    }

    //dt
    double dt = mxGetPr(prhs[5])[0];

    //everything else are just defaults for now


	//setting up default values
    double L = 0.04;
    double p = 1000;
    double E = 37.8e3;
    double G = E/3;
    double D = 1e-2;
    double R = D/2;
    double A = M_PI*R*R;
    double I = M_PI*pow(D,4)/32;
    double J = 2*I;


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

    gsl_matrix *extra = gsl_matrix_calloc(n,2*N);
    struct SystemParameters sys_params = { dt, ds, n, N, g0, eta0, xi0, q, q_dot, false, extra};
    struct ActuatorParameters act_params;
    if (act_flag == 0){ //cable
    	sys_params.delta = false;
    	act_params.rx = rx_fun;
    	act_params.ry = ry_fun;
    	act_params.force = cableForce;
    } else if (act_flag == 1) { //tca
    	sys_params.delta = true;
    	act_params.rx = rx_fun;
    	act_params.ry = ry_fun;
    	act_params.force = tcaForce;
    }


    struct BodyParameters body_params = { E, G, p, J, I, A, R, L, xi_ref};

    gsl_matrix *eta = gsl_matrix_alloc(n,6);
    gsl_matrix *g = gsl_matrix_alloc(n,12);
    gsl_matrix *xi = gsl_matrix_alloc(n,6);

    gsl_matrix_set_zero(eta);
    gsl_matrix_set_zero(g);
    gsl_matrix_set_zero(xi);

    struct SimulationParameters params = {body_params, sys_params, act_params, g, xi, xi_prev, eta_prev, eta};


    //one step of the dynamics
    stepDynamics(params);


    //put it all back in matlab
    //need the outputs transposed
    double *g_out_vec, *eta_out_vec, *xi_out_vec;
    mxArray *g_out, *eta_out, *xi_out;
    g_out = plhs[0] = mxCreateDoubleMatrix(n,12,mxREAL);
    g_out_vec = mxGetPr(g_out);
    eta_out = plhs[2] = mxCreateDoubleMatrix(n,6,mxREAL);
    eta_out_vec = mxGetPr(eta_out);
    xi_out = plhs[1] = mxCreateDoubleMatrix(n,6,mxREAL);
    xi_out_vec = mxGetPr(xi_out);

    for (i=0;i<n;i++) {
    	for (j=0;j<6;j++) {
    		eta_out_vec[j*n+i] = gsl_matrix_get(params.eta,i,j);
    		xi_out_vec[j*n+i] = gsl_matrix_get(params.xi,i,j);

    		g_out_vec[j*n+i] = gsl_matrix_get(params.g,i,j);
    		g_out_vec[(j+6)*n+i] = gsl_matrix_get(params.g,i,j+6);
    	}
    }




    gsl_matrix_free(extra);

    gsl_matrix_free(g0);
    gsl_vector_free(xi0);
    gsl_vector_free(xi_ref);
    gsl_vector_free(eta0);

    gsl_vector_free(q);
    gsl_vector_free(q_dot);

    gsl_matrix_free(eta);
    gsl_matrix_free(eta_prev);
    gsl_matrix_free(xi);
    gsl_matrix_free(xi_prev);
    gsl_matrix_free(g);


}
