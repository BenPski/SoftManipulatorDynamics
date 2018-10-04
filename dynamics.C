#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multiroots.h>
#include <stdbool.h>

#include "dynamics.h"



/*
 *Want to define the general procedure for computing the dynamics, since there is very little difference between
 *the different kinds of dynamics
 *
 *most of the procedures are them same, the main difference is the force fucntions and the delta integration
 *
 *also want to be able to easily pass around body parameters
 *
 *One major issue is that there does not appear to be an obvious way of doing the computations in a unified way across
 *cables and TCAs due to the TCAs needing extra parameters
 */


/*
 * very simplistic printing functions, mostly for debugging
 */

void printVector(const gsl_vector *v) {
    int size = v->size;

    printf("Vector:\n");

    int i;
    for (i=0;i<size;i++) {
        printf("%f\n",gsl_vector_get(v,i));
    }
    printf("\n");
}

void printMatrix(const gsl_matrix *m) {
    int height = m->size1;
    int width = m->size2;

    printf("Matrix: \n");

    int i,j;
    for (i=0;i<height;i++) {
        for (j=0;j<width-1;j++) {
            printf("%f ", gsl_matrix_get(m,i,j));
        }
        printf("%f\n", gsl_matrix_get(m,i,width-1));
    }
    printf("\n");
}



/*
 *The structs for defining various parameters
 */
/*
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

    //flags
    bool delta; //compute delta in the integration
    
    gsl_matrix *extra; //i don't know if this is the best choice, but it is what I'm using for now 
};


struct ActuatorParameters {
    double (*rx)(double R, int i, int N, double s); //the function to determine the displacement of the actuator, ACT# -> s -> r
    double (*ry)(double R, int i, int N, double s);
    double (*force)(struct SystemParameters sys_params, int i, int j, double s); //the function to determine the force, input -> ACT# -> POS# -> s -> r
};
*/

/*
 *The first few utility functions
 */

//the isomorphisms between R3 and skew matrices
//computes the skew matrix of a 3d vector
//not very safe and I have avoided using vectors and just had matrices for everything
void skew(const gsl_vector *a, gsl_matrix *a_skew) {
    gsl_matrix_set(a_skew,0,0,0);
    gsl_matrix_set(a_skew,0,1,-gsl_vector_get(a,2));
    gsl_matrix_set(a_skew,0,2,gsl_vector_get(a,1));
    gsl_matrix_set(a_skew,1,0,gsl_vector_get(a,2));
    gsl_matrix_set(a_skew,1,1,0);
    gsl_matrix_set(a_skew,1,2,-gsl_vector_get(a,0));
    gsl_matrix_set(a_skew,2,0,-gsl_vector_get(a,1));
    gsl_matrix_set(a_skew,2,1,gsl_vector_get(a,0));
    gsl_matrix_set(a_skew,2,2,0);
}

void unskew(const gsl_matrix *a_skew, gsl_vector *a) {
    gsl_vector_set(a,0,gsl_matrix_get(a_skew,2,1));
    gsl_vector_set(a,1,gsl_matrix_get(a_skew,0,2));
    gsl_vector_set(a,2,gsl_matrix_get(a_skew,1,0));
}


//the isomorphisms for the se(3) algebra, go between R6 and se(3)
//can make use of the skew isomorphisms
void se(const gsl_vector *a, gsl_matrix *a_se) {
    gsl_vector_const_view ang = gsl_vector_const_subvector(a,0,3); //the angular portion
    gsl_vector_const_view trans = gsl_vector_const_subvector(a,3,3); //the translational portion

    int i;

    skew(&ang.vector,a_se); //get the skew portion

    for (i=0; i<3; i++) {
        gsl_matrix_set(a_se,i,3,gsl_vector_get(&trans.vector,i));
    }

    for (i=0;i<4;i++) {
        gsl_matrix_set(a_se,3,i,0);
    }
}

void unse(const gsl_matrix *a_se, gsl_vector *a) {
    //could use submatrices, but it should work out properly
    unskew(a_se,a);
    int i;
    for (i=0;i<3;i++) {
        gsl_vector_set(a,3+i,gsl_matrix_get(a_se,i,3));
    }
}

//the isomorphisms for the configuration
//can either flatten and unflatten the config
//or can extract the angles and create the rotation from it
void flattenConfig(const gsl_matrix *m, gsl_vector *v) {
    //know some portions are constant so ignore
    int i,j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            //rotational part
            gsl_vector_set(v,i*3+j,gsl_matrix_get(m,i,j));
        }
        //translational part
        gsl_vector_set(v,9+i,gsl_matrix_get(m,i,3));
    }
}
void unflattenConfig(const gsl_vector *v, gsl_matrix *m) {
    int i,j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            //rotational
            gsl_matrix_set(m,i,j,gsl_vector_get(v,i*3+j));
        }
        //translational
        gsl_matrix_set(m,i,3,gsl_vector_get(v,9+i));
        //constant
        gsl_matrix_set(m,3,i,0);
    }
    gsl_matrix_set(m,3,3,1);
}

//given a rotation matrix get the angles
//https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2012/07/euler-angles1.pdf
void extractAngles(const gsl_matrix *R, gsl_vector *angles) {

    double theta1 = atan2(gsl_matrix_get(R,1,2),gsl_matrix_get(R,2,2));
    gsl_vector_set(angles,0,theta1);

    double c2 = sqrt(pow(gsl_matrix_get(R,0,0),2)+pow(gsl_matrix_get(R,0,1),2));
    double theta2 = atan2(-gsl_matrix_get(R,0,2),c2);
    gsl_vector_set(angles,1,theta2);

    double s1 = sin(theta1);
    double c1 = cos(theta1);
    double theta3 = atan2(s1*gsl_matrix_get(R,2,0)-c1*gsl_matrix_get(R,1,0),c1*gsl_matrix_get(R,1,1)-s1*gsl_matrix_get(R,2,1));
    gsl_vector_set(angles,2,theta3);
}

void rotation_x(const double theta, gsl_matrix *R) {
    const double c = cos(theta);
    const double s = sin(theta);
    //cos terms
    gsl_matrix_set(R,1,1,c);
    gsl_matrix_set(R,2,2,c);
    //sin terms
    gsl_matrix_set(R,1,2,s);
    gsl_matrix_set(R,2,1,-s);
    //1 term
    gsl_matrix_set(R,0,0,1);
    //0 terms
    gsl_matrix_set(R,0,1,0);
    gsl_matrix_set(R,0,2,0);
    gsl_matrix_set(R,1,0,0);
    gsl_matrix_set(R,2,0,0);
}
void rotation_y(const double theta, gsl_matrix *R) {
    const double c = cos(theta);
    const double s = sin(theta);
    //cos terms
    gsl_matrix_set(R,0,0,c);
    gsl_matrix_set(R,2,2,c);
    //sin terms
    gsl_matrix_set(R,0,2,-s);
    gsl_matrix_set(R,2,0,s);
    //1 term
    gsl_matrix_set(R,1,1,1);
    //0 terms
    gsl_matrix_set(R,0,1,0);
    gsl_matrix_set(R,1,2,0);
    gsl_matrix_set(R,2,1,0);
    gsl_matrix_set(R,1,0,0);
}
void rotation_z(const double theta, gsl_matrix *R) {
    const double c = cos(theta);
    const double s = sin(theta);
    //cos terms
    gsl_matrix_set(R,0,0,c);
    gsl_matrix_set(R,1,1,c);
    //sin terms
    gsl_matrix_set(R,0,1,s);
    gsl_matrix_set(R,1,0,-s);
    //1 term
    gsl_matrix_set(R,2,2,1);
    //0 terms
    gsl_matrix_set(R,0,2,0);
    gsl_matrix_set(R,1,2,0);
    gsl_matrix_set(R,2,0,0);
    gsl_matrix_set(R,2,1,0);
}

//the inverse of extract angles
//note these functions are not one-to-one so can end up swapping around the actual angles, but have consistent rotation matrices
void rotation(const gsl_vector *angles,gsl_matrix *R) {
    gsl_matrix *Ry = gsl_matrix_calloc(3,3);
    gsl_matrix *Rz = gsl_matrix_calloc(3,3);
    gsl_matrix *A = gsl_matrix_calloc(3,3); //just a temp matrix to compy into

    rotation_x(gsl_vector_get(angles,0),R);
    rotation_y(gsl_vector_get(angles,1),Ry);
    rotation_z(gsl_vector_get(angles,2),Rz);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,R,Ry,0.0,A);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,Rz,0.0,R);

    gsl_matrix_free(A);
    gsl_matrix_free(Ry);
    gsl_matrix_free(Rz);
}

//extract angles and translations, ordered as [angles,positions]
void extractConfig(const gsl_matrix *g, gsl_vector *v) {
    extractAngles(g,v); //should copy to right place
    int i;
    for (i=0;i<3;i++) {
        gsl_vector_set(v,i+3,gsl_matrix_get(g,i,3));
    }
}

void unextractConfig(const gsl_vector *v, gsl_matrix *g) {
    rotation(v,g); //should line up properly
    int i;
    for (i=0;i<3;i++) {
        gsl_matrix_set(g,i,3,gsl_vector_get(v,3+i));
        gsl_matrix_set(g,3,i,0);
    }
    gsl_matrix_set(g,3,3,1);
}

//the adjoint representation of a vector in se(3)
void adjoint(const gsl_vector *alg, gsl_matrix *adj) {
    gsl_matrix *w_skew = gsl_matrix_calloc(3,3);
    gsl_matrix *v_skew = gsl_matrix_calloc(3,3);

    gsl_vector_const_view w = gsl_vector_const_subvector(alg,0,3);
    gsl_vector_const_view v = gsl_vector_const_subvector(alg,3,3);

    skew(&w.vector,w_skew);
    skew(&v.vector,v_skew);

    int i,j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            gsl_matrix_set(adj,i,j,gsl_matrix_get(w_skew,i,j));
            gsl_matrix_set(adj,i,j+3,0);
            gsl_matrix_set(adj,i+3,j,gsl_matrix_get(v_skew,i,j));
            gsl_matrix_set(adj,i+3,j+3,gsl_matrix_get(w_skew,i,j));
        }
    }
    gsl_matrix_free(w_skew);
    gsl_matrix_free(v_skew);
}

//the matrix exponential for SE(3)
//annoying to check for what the representation of the algebra, so just allow a vector
void expSE(const gsl_vector *alg, gsl_matrix *group) {

    double a,b,c,d;
    double wx,wy,wz,vx,vy,vz;
    double w_norm;
    /*
    [ 1 - (b*(wy^2 + wz^2))/2,      (b*wx*wy)/2 - a*wz,      a*wy + (b*wx*wz)/2, vz*((b*wy)/2 + c*wx*wz) - vy*((b*wz)/2 - c*wx*wy) - vx*(c*(wy^2 + wz^2) - 1)]
    [      a*wz + (b*wx*wy)/2, 1 - (b*(wx^2 + wz^2))/2,      (b*wy*wz)/2 - a*wx, vx*((b*wz)/2 + c*wx*wy) - vy*(c*(wx^2 + wz^2) - 1) - vz*((b*wx)/2 - c*wy*wz)]
    [      (b*wx*wz)/2 - a*wy,      a*wx + (b*wy*wz)/2, 1 - (b*(wx^2 + wy^2))/2, vy*((b*wx)/2 + c*wy*wz) - vx*((b*wy)/2 - c*wx*wz) - vz*(c*(wx^2 + wy^2) - 1)]
    [                       0,                       0,                       0,                                                                            1]
     */
    wx = gsl_vector_get(alg,0);//xi[0];
    wy = gsl_vector_get(alg,1);//xi[1];
    wz = gsl_vector_get(alg,2);//xi[2];
    vx = gsl_vector_get(alg,3);//xi[3];
    vy = gsl_vector_get(alg,4);//xi[4];
    vz = gsl_vector_get(alg,5);//xi[5];

    if (wx==0 && wy==0 && wz==0) {
        a = b = 1;
        c = 0;
    } else {
        w_norm = sqrt(wx*wx+wy*wy+wz*wz);
        a = sin(w_norm)/w_norm;
        b = 2*(1-cos(w_norm))/(w_norm*w_norm);
        c = (1-a)/(w_norm*w_norm);
    }

    // rotational part
    gsl_matrix_set(group,0,0, 1 - (b*(wy*wy + wz*wz))/2);
    gsl_matrix_set(group,0,1,(b*wx*wy)/2 - a*wz);
    gsl_matrix_set(group,0,2, a*wy + (b*wx*wz)/2);
    gsl_matrix_set(group,1,0, a*wz + (b*wx*wy)/2);
    gsl_matrix_set(group,1,1, 1 - (b*(wx*wx + wz*wz))/2);
    gsl_matrix_set(group,1,2, (b*wy*wz)/2 - a*wx);
    gsl_matrix_set(group,2,0, (b*wx*wz)/2 - a*wy);
    gsl_matrix_set(group,2,1, a*wx + (b*wy*wz)/2);
    gsl_matrix_set(group,2,2, 1 - (b*(wx*wx + wy*wy))/2);

    //translational part
    gsl_matrix_set(group,0,3, vz*((b*wy)/2 + c*wx*wz) - vy*((b*wz)/2 - c*wx*wy) - vx*(c*(wy*wy + wz*wz) - 1));
    gsl_matrix_set(group,1,3, vx*((b*wz)/2 + c*wx*wy) - vy*(c*(wx*wx + wz*wz) - 1) - vz*((b*wx)/2 - c*wy*wz));
    gsl_matrix_set(group,2,3, vy*((b*wx)/2 + c*wy*wz) - vx*((b*wy)/2 - c*wx*wz) - vz*(c*(wx*wx + wy*wy) - 1));

    //constant part
    gsl_matrix_set(group,3,0,0);
    gsl_matrix_set(group,3,1,0);
    gsl_matrix_set(group,3,2,0);
    gsl_matrix_set(group,3,3,1);
}



/*
 *Numerical utility functions
 */

//can also make this a banded matrix multiplication
void approxDerivatives(const struct SystemParameters sys_params, const gsl_matrix *xi, gsl_matrix *xi_dot) {
    double ds = sys_params.ds;
    int n = sys_params.n;
    int i, j;

    gsl_matrix_set_zero(xi_dot);

    gsl_vector_view xi_dot_view;
    gsl_vector *xi_ahead = gsl_vector_calloc(6);
    gsl_vector *xi_behind = gsl_vector_calloc(6);

    //gsl_matrix_set_zero(xi_dot);

    /* first point and last point approximation */
    //first point
    xi_dot_view = gsl_matrix_row(xi_dot,0);
    gsl_matrix_get_row(xi_ahead,xi,1);
    gsl_matrix_get_row(xi_behind,xi,0);

    gsl_vector_add(&xi_dot_view.vector,xi_ahead);
    gsl_vector_sub(&xi_dot_view.vector,xi_behind);
    gsl_vector_scale(&xi_dot_view.vector,1/ds);

    //last point
    xi_dot_view = gsl_matrix_row(xi_dot,n-1);
    gsl_matrix_get_row(xi_ahead,xi,n-1);
    gsl_matrix_get_row(xi_behind,xi,n-2);

    gsl_vector_add(&xi_dot_view.vector,xi_ahead);
    gsl_vector_sub(&xi_dot_view.vector,xi_behind);
    gsl_vector_scale(&xi_dot_view.vector,1/ds);

    //intermediate points
    //centered diff

    for (i=1;i<n-1;i++) {
        xi_dot_view = gsl_matrix_row(xi_dot,i);
        gsl_matrix_get_row(xi_ahead,xi,i+1);
        gsl_matrix_get_row(xi_behind,xi,i-1);

        gsl_vector_add(&xi_dot_view.vector,xi_ahead);
        gsl_vector_sub(&xi_dot_view.vector,xi_behind);
        gsl_vector_scale(&xi_dot_view.vector,1/(2*ds));
    }


    //backward diff
    /*
    for (i=1;i<n-1;i++) {
        xi_dot_view = gsl_matrix_row(xi_dot,i);
        gsl_matrix_get_row(xi_ahead,xi,i);
        gsl_matrix_get_row(xi_behind,xi,i-1);

        gsl_vector_add(&xi_dot_view.vector,xi_ahead);
        gsl_vector_sub(&xi_dot_view.vector,xi_behind);
        gsl_vector_scale(&xi_dot_view.vector,1/(ds));
    }
    */

    gsl_vector_free(xi_ahead);
    gsl_vector_free(xi_behind);
}

//the dynamics always needs to compute g and xi from eta, xi(0), and xi_prev
void integrate(struct SystemParameters sys_params, struct BodyParameters body_params, struct ActuatorParameters act_params, const gsl_matrix *eta, const gsl_matrix *xi_prev, gsl_matrix *g, gsl_matrix *xi) {
    //the matrices of values will be stored as a series of rows
    //for now just do the g and xi integration, worry about delta later

    //pretty undecided about how much should be done to get the arithmetic all done in blas vs for loops,
    //really don't know the trade-off,
    //the matrices are generally small enough to not really warrant worrying about speed it seems


    double ds = sys_params.ds;
    double dt = sys_params.dt;
    int n = sys_params.n;
    int N = sys_params.N;
    
    gsl_matrix *g0 = sys_params.g0;
    gsl_vector *eta0 = sys_params.eta0;
    gsl_vector *xi0 = sys_params.xi0;

    gsl_matrix *g_temp = gsl_matrix_calloc(4,4);
    gsl_vector *g_row = gsl_vector_calloc(12);

    gsl_vector *xi_prev_row = gsl_vector_calloc(6);
    gsl_vector *eta_row = gsl_vector_calloc(6);
    gsl_vector *xi_row = gsl_vector_calloc(6);
    gsl_vector *xi_row_off = gsl_vector_calloc(6);

    gsl_matrix *xi_adj = gsl_matrix_calloc(6,6);
    gsl_vector *xi_acc = gsl_vector_calloc(6);
    gsl_vector *xi_acc_temp = gsl_vector_calloc(6);

    //prepare the initial values of the matrices
    //also just store g as flattened array
    int i,j;

    flattenConfig(g0,g_row);

    for (i=0;i<12;i++) {
        gsl_matrix_set(g,0,i,gsl_vector_get(g_row,i));
    }
    for (i=0;i<6;i++) {
        gsl_matrix_set(xi,0,i,gsl_vector_get(xi0,i));
        //gsl_matrix_set(eta,0,i,gsl_vector_get(eta0,i)); //do I need to do this?
    }

    for (i=1;i<n;i++) {
        for (j=0;j<6;j++) {
            //the approximation for xi = xi_prev + dt/ds*(eta-eta_off)
            gsl_matrix_set(xi,i,j, gsl_matrix_get(xi_prev,i,j)+dt/ds*(gsl_matrix_get(eta,i,j)-gsl_matrix_get(eta,i-1,j)));
        }

        //now add the adjoint portion
        //pull out the relevant views
        gsl_matrix_get_row(xi_prev_row,xi_prev,i);
        gsl_matrix_get_row(eta_row, eta,i);
        adjoint(xi_prev_row,xi_adj);
        gsl_matrix_get_row(xi_row,xi,i);
        gsl_matrix_get_row(xi_row_off,xi,i-1);
        

        //do the adjoint portion of equation
        gsl_blas_dgemv(CblasNoTrans,dt,xi_adj,eta_row,1.0,xi_row);
        gsl_matrix_set_row(xi,i,xi_row);

        //do the accumulated xi
        gsl_vector_add(xi_acc,xi_row);
        gsl_vector_add(xi_acc,xi_row_off);

        //compute g
        gsl_vector_memcpy(xi_acc_temp,xi_acc);
        gsl_vector_scale(xi_acc_temp,ds/2);

        expSE(xi_acc_temp, g_temp);

        flattenConfig(g_temp,g_row);
        gsl_matrix_set_row(g,i,g_row);
    }
    
    if (sys_params.delta) {
        gsl_vector *r = gsl_vector_calloc(4);
        gsl_vector_set(r,3,1);
        gsl_matrix *xi_se = gsl_matrix_calloc(4,4);
        gsl_vector *temp_vec = gsl_vector_calloc(4);
        double delta;
        for (i=0;i<n;i++) {
            for (j=0;j<N;j++) {
                gsl_vector_set(r,0,act_params.rx(body_params.R,j,N,ds*i));
                gsl_vector_set(r,1,act_params.ry(body_params.R,j,N,ds*i));
                //delta = ds*(norm(se(xi)*r) - norm(se(xi_ref)*r))
                gsl_matrix_get_row(xi_row,xi,i);
                se(xi_row,xi_se);
                gsl_blas_dgemv(CblasNoTrans,1.0,xi_se,r,0.0,temp_vec);
                
                delta = gsl_blas_dnrm2(temp_vec);
                se(body_params.xi_ref,xi_se);
                gsl_blas_dgemv(CblasNoTrans,1.0,xi_se,r,0.0,temp_vec);
                
                delta -= gsl_blas_dnrm2(temp_vec);
                
                gsl_matrix_set(sys_params.extra,i,j,ds*delta);
            }
        }
        gsl_vector_free(r);
        gsl_vector_free(temp_vec);
        gsl_matrix_free(xi_se);
    }
        
    

    gsl_vector_free(xi_prev_row);
    gsl_vector_free(eta_row);
    gsl_vector_free(xi_row);
    gsl_vector_free(xi_row_off);

    gsl_matrix_free(g_temp);
    gsl_vector_free(g_row);
    gsl_matrix_free(xi_adj);
    gsl_vector_free(xi_acc);
    gsl_vector_free(xi_acc_temp);

}

//computing what the derivative of eta is based off the guessed values and the integral/derivative values
void timeDerivative(struct SystemParameters sys_params, struct BodyParameters body_params, struct ActuatorParameters act_params,
        const gsl_matrix *g_mat, const gsl_matrix *xi_mat, const gsl_matrix *eta_mat, const gsl_matrix *xi_dot,
        gsl_matrix *eta_der_mat) {


    int i, j, k, stateSize, m, N;
    double rx, ry, E, G, I, J, A, p, grav, ds, dt, Rad;
    double ox,oy,oz,nx,ny,nz, wx,wy,wz,vx,vy,vz, F;
    double ox_dot,oy_dot,oz_dot,nx_dot,ny_dot,nz_dot, wx_dot,wy_dot,wz_dot,vx_dot,vy_dot,vz_dot;

    double temp_vec[6];
    double F_vec[6];

    double pi_dot_norm;

    gsl_matrix *R = gsl_matrix_calloc(3,3);
    gsl_vector *xi = gsl_vector_calloc(6);
    gsl_vector *eta = gsl_vector_calloc(6);
    gsl_matrix *K = gsl_matrix_calloc(6,6);
    gsl_matrix *gamma = gsl_matrix_calloc(6,6);
    gsl_vector_view o_view;
    gsl_matrix *o_skew = gsl_matrix_calloc(3,3);
    gsl_vector_view n_view;
    gsl_matrix *n_skew = gsl_matrix_calloc(3,3);
    gsl_vector_view w_view;
    gsl_matrix *w_skew = gsl_matrix_calloc(3,3);
    gsl_vector_view v_view;
    gsl_matrix *v_skew = gsl_matrix_calloc(3,3);
    gsl_vector *r_vec = gsl_vector_calloc(3);
    gsl_matrix *r_skew = gsl_matrix_calloc(3,3);
    gsl_vector *xi_der = gsl_vector_calloc(6);
    gsl_vector_view o_der_view;
    gsl_vector_view n_der_view;
    gsl_matrix *ad_xi = gsl_matrix_calloc(6,6);
    gsl_matrix *ad_eta = gsl_matrix_calloc(6,6);
    gsl_vector *temp_vec0 = gsl_vector_calloc(6);
    gsl_vector *temp_vec1 = gsl_vector_calloc(6);
    gsl_vector *temp_vec2 = gsl_vector_calloc(6);
    gsl_vector *temp_vec3 = gsl_vector_calloc(3);
    gsl_vector *temp_vec4 = gsl_vector_calloc(3);
    gsl_vector *temp_vec5 = gsl_vector_calloc(3);
    gsl_vector *p_dot = gsl_vector_calloc(3);
    gsl_matrix *p_dot_skew = gsl_matrix_calloc(3,3);
    gsl_vector *p_ddot = gsl_vector_calloc(3);

    gsl_vector *g_row = gsl_vector_calloc(12);
    gsl_matrix *g = gsl_matrix_calloc(4,4);
    gsl_matrix_view R_view;

    //get the body parameters
    E = body_params.E;
    G = body_params.G;
    p = body_params.p;
    A = body_params.A;
    I = body_params.I;
    J = body_params.J;
    Rad = body_params.R;
    gsl_vector *xi_ref = body_params.xi_ref;

    //the system parameters
    dt = sys_params.dt;
    ds = sys_params.ds;
    m = sys_params.n;
    N = sys_params.N;

    gsl_vector *q = sys_params.q;

    grav = -9.8;
    //grav = 0;

    gsl_matrix_set(K,0,0,E*I);
    gsl_matrix_set(K,1,1,E*I);
    gsl_matrix_set(K,2,2,G*J);
    gsl_matrix_set(K,3,3,G*A);
    gsl_matrix_set(K,4,4,G*A);
    gsl_matrix_set(K,5,5,E*A);
    gsl_matrix_set(gamma,0,0,p*I);
    gsl_matrix_set(gamma,1,1,p*I);
    gsl_matrix_set(gamma,2,2,p*J);
    gsl_matrix_set(gamma,3,3,p*A);
    gsl_matrix_set(gamma,4,4,p*A);
    gsl_matrix_set(gamma,5,5,p*A);



    for (i=0; i<m; i++) {
        /* pull out state variables */
        //not using the big state vector y anymore

        //get the R matrix
        gsl_matrix_get_row(g_row,g_mat,i);
        unflattenConfig(g_row,g);
        R_view = gsl_matrix_submatrix(g,0,0,3,3);
        gsl_matrix_memcpy(R,&R_view.matrix);


        //get xi vector
        gsl_matrix_get_row(xi,xi_mat,i);
        gsl_matrix_get_row(eta,eta_mat,i);
        gsl_matrix_get_row(xi_der,xi_dot,i);

        //prepare skew martices
        o_view = gsl_vector_subvector(xi,0,3);
        n_view = gsl_vector_subvector(xi,3,3);
        w_view = gsl_vector_subvector(eta,0,3);
        v_view = gsl_vector_subvector(eta,3,3);
        o_der_view = gsl_vector_subvector(xi_der,0,3);
        n_der_view = gsl_vector_subvector(xi_der,3,3);

        skew(&o_view.vector,o_skew);
        skew(&n_view.vector,n_skew);
        skew(&v_view.vector,v_skew);
        skew(&w_view.vector,w_skew);

        adjoint(xi,ad_xi);
        adjoint(eta,ad_eta);


        /* eta' */
        /* K*xi_dot */
        gsl_blas_dgemv(CblasNoTrans,1.0,K,xi_der,0.0,temp_vec0);

        /* -ad_xi_star*K*(xi-xi_ref) */
        gsl_vector_memcpy(temp_vec1,xi);
        gsl_vector_sub(temp_vec1,xi_ref);

        gsl_blas_dgemv(CblasNoTrans,1.0,K,temp_vec1,0.0,temp_vec2);
        gsl_blas_dgemv(CblasTrans,-1.0,ad_xi,temp_vec2,1.0,temp_vec0);


        /* ad_eta_star*gamma*eta */
        gsl_blas_dgemv(CblasNoTrans,1.0,gamma,eta,0.0,temp_vec1);
        gsl_blas_dgemv(CblasTrans,1.0,ad_eta,temp_vec1,1.0,temp_vec0);


        /* gravity */
        gsl_vector_set(temp_vec0,3, gsl_vector_get(temp_vec0,3)+A*grav*p*gsl_matrix_get(R,2,0));
        gsl_vector_set(temp_vec0,4, gsl_vector_get(temp_vec0,4)+A*grav*p*gsl_matrix_get(R,2,1));
        gsl_vector_set(temp_vec0,5, gsl_vector_get(temp_vec0,5)+A*grav*p*gsl_matrix_get(R,2,2));

        /* tcas */
        for (j=0;j<6;j++) {
            temp_vec[j]=gsl_vector_get(temp_vec0,j);
            F_vec[j] = 0;
        }

        ox=gsl_vector_get(&o_view.vector,0);
        oy=gsl_vector_get(&o_view.vector,1);
        oz=gsl_vector_get(&o_view.vector,2);
        nx=gsl_vector_get(&n_view.vector,0);
        ny=gsl_vector_get(&n_view.vector,1);
        nz=gsl_vector_get(&n_view.vector,2);

        wx=gsl_vector_get(&w_view.vector,0);
        wy=gsl_vector_get(&w_view.vector,1);
        wz=gsl_vector_get(&w_view.vector,2);
        vx=gsl_vector_get(&v_view.vector,0);
        vy=gsl_vector_get(&v_view.vector,1);
        vz=gsl_vector_get(&v_view.vector,2);


        for (j=0;j<N;j++) {
            rx = act_params.rx(Rad,j,N,ds*i);
            ry = act_params.ry(Rad,j,N,ds*i);
            gsl_vector_set(r_vec,0,rx);
            gsl_vector_set(r_vec,1,ry);
            gsl_vector_set(r_vec,2,0);
            skew(r_vec,r_skew);

            /* sum -F*[skew(r)*R'*skew(p_dot)^2*p_ddot/norm(p_dot)^3;R'*skew(p_dot)^2*p_ddot/norm(p_dot)^3] */
            /* p_dot = R(n+skew(o)*r) */
            gsl_blas_dgemv(CblasNoTrans,1.0,o_skew,r_vec,0.0,temp_vec3);
            gsl_vector_add(temp_vec3,&n_view.vector);
            gsl_blas_dgemv(CblasNoTrans,1.0,R,temp_vec3,0.0,p_dot);
            skew(p_dot,p_dot_skew);

            /* p_ddot = R*(skew(o)*(n+skew(w)*r)+n_dot+skew(o_dot)*r) */
            //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0,r_skew,o_der,1.0,n_der);
            gsl_blas_dgemv(CblasNoTrans,-1.0,r_skew,&o_der_view.vector,0.0,temp_vec5);
            gsl_vector_add(temp_vec5,&n_der_view.vector);
            //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,o_skew,r_vec,1.0,n);
            gsl_blas_dgemv(CblasNoTrans,1.0,o_skew,r_vec,0.0,temp_vec3);
            gsl_vector_add(temp_vec3,&n_view.vector);
            //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,o_skew,n,1.0,n_der);
            gsl_blas_dgemv(CblasNoTrans,1.0,o_skew,temp_vec3,0.0,temp_vec4);
            gsl_vector_add(temp_vec4,temp_vec5);
            //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,R,n_der,0.0,p_ddot);
            gsl_blas_dgemv(CblasNoTrans,1.0,R,temp_vec4,0.0,p_ddot);

            /* R'*skew(p_dot)*skew(p_dot)*p_ddot */
            gsl_blas_dgemv(CblasNoTrans,1.0,p_dot_skew,p_ddot,0.0,temp_vec3);
            gsl_blas_dgemv(CblasNoTrans,1.0,p_dot_skew,temp_vec3,0.0,p_ddot);
            gsl_blas_dgemv(CblasTrans,1.0,R,p_ddot,0.0,temp_vec3);
            /* skew(r) * res */
            gsl_blas_dgemv(CblasNoTrans,1.0,r_skew,temp_vec3,0.0,p_ddot);

            F = act_params.force(sys_params,j,i,ds*i);

            pi_dot_norm = pow(pow(nz + ox*ry - oy*rx,2.0) + pow(ny + oz*rx,2.0) + pow(nx - oz*ry,2.0),1.5);

            F_vec[0] += (-F/pi_dot_norm)*gsl_vector_get(p_ddot,0);
            F_vec[1] += (-F/pi_dot_norm)*gsl_vector_get(p_ddot,1);
            F_vec[2] += (-F/pi_dot_norm)*gsl_vector_get(p_ddot,2);
            F_vec[3] += (-F/pi_dot_norm)*gsl_vector_get(temp_vec3,0);
            F_vec[4] += (-F/pi_dot_norm)*gsl_vector_get(temp_vec3,1);
            F_vec[5] += (-F/pi_dot_norm)*gsl_vector_get(temp_vec3,2);

        }


        for (j=0;j<6;j++) {
            temp_vec[j] += F_vec[j];
        }

        temp_vec[0] = temp_vec[0]/(I*p);
        temp_vec[1] = temp_vec[1]/(I*p);
        temp_vec[2] = temp_vec[2]/(J*p);
        temp_vec[3] = temp_vec[3]/(A*p);
        temp_vec[4] = temp_vec[4]/(A*p);
        temp_vec[5] = temp_vec[5]/(A*p);

        /* place it in the eta_der matrix */
        for (j=0;j<6;j++) {
            gsl_matrix_set(eta_der_mat,i,j,temp_vec[j]);
        }
    }


    gsl_matrix_free(R);
    gsl_vector_free(xi);
    gsl_vector_free(eta);
    gsl_matrix_free(K);
    gsl_matrix_free(gamma);
    gsl_matrix_free(o_skew);
    gsl_matrix_free(n_skew);
    gsl_matrix_free(w_skew);
    gsl_matrix_free(v_skew);
    gsl_vector_free(r_vec);
    gsl_matrix_free(r_skew);
    gsl_vector_free(xi_der);
    gsl_matrix_free(ad_xi);
    gsl_matrix_free(ad_eta);
    gsl_vector_free(temp_vec0);
    gsl_vector_free(temp_vec1);
    gsl_vector_free(temp_vec2);
    gsl_vector_free(temp_vec3);
    gsl_vector_free(temp_vec4);
    gsl_vector_free(temp_vec5);
    gsl_vector_free(p_dot);
    gsl_matrix_free(p_dot_skew);
    gsl_vector_free(p_ddot);

    gsl_vector_free(g_row);
    gsl_matrix_free(g);
}

//put the boundary conditions into a vector
//boundary conditions are eta = eta_prev + dt*eta_der and the tip wrench balance W_ext = W_int
//do not need eta(0) as that is a given initial condition
void boundaryConditions(struct SystemParameters sys_params, struct BodyParameters body_params, struct ActuatorParameters act_params,
        const gsl_matrix *eta, const gsl_matrix *eta_prev, const gsl_matrix *eta_der,
        const gsl_matrix *xi,
        gsl_vector *conditions) {

    double dt = sys_params.dt;
    int n = sys_params.n;
    int N = sys_params.N;
    gsl_vector *q = sys_params.q;

    gsl_vector_set_zero(conditions);

    //the backward euler condition
    gsl_matrix_const_view eta_view = gsl_matrix_const_submatrix(eta,1,0,n-1,6);
    gsl_matrix_const_view eta_prev_view = gsl_matrix_const_submatrix(eta_prev,1,0,n-1,6);
    gsl_matrix_const_view eta_der_view = gsl_matrix_const_submatrix(eta_der,1,0,n-1,6);
    gsl_matrix *eulerCond = gsl_matrix_alloc(n-1,6);

    gsl_matrix_memcpy(eulerCond,&eta_der_view.matrix);
    gsl_matrix_scale(eulerCond,dt);
    gsl_matrix_add(eulerCond,&eta_prev_view.matrix);
    gsl_matrix_sub(eulerCond,&eta_view.matrix);

    int i,j;
    for (i=0;i<n-1;i++) {
        for (j=0;j<6;j++) {
            gsl_vector_set(conditions,i*6+j,gsl_matrix_get(eulerCond,i,j));
        }
    }

    //the tip wrench condition
    double E = body_params.E;
    double G = body_params.G;
    double p = body_params.p;
    double A = body_params.A;
    double I = body_params.I;
    double J = body_params.J;
    double L = body_params.L;
    double Rad = body_params.R;
    gsl_vector *xi_ref = body_params.xi_ref;

    gsl_vector *xiL = gsl_vector_alloc(6);
    gsl_matrix_get_row(xiL,xi,n-1);

    gsl_matrix *K = gsl_matrix_calloc(6,6);
    gsl_matrix_set(K,0,0,E*I);
    gsl_matrix_set(K,1,1,E*I);
    gsl_matrix_set(K,2,2,G*J);
    gsl_matrix_set(K,3,3,G*A);
    gsl_matrix_set(K,4,4,G*A);
    gsl_matrix_set(K,5,5,E*A);

    gsl_vector *F = gsl_vector_alloc(6);
    gsl_vector_set_zero(F);
    double f,rx,ry, ox,oy,oz,nx,ny,nz, s_norm;

    ox = gsl_vector_get(xiL,0);
    oy = gsl_vector_get(xiL,1);
    oz = gsl_vector_get(xiL,2);
    nx = gsl_vector_get(xiL,3);
    ny = gsl_vector_get(xiL,4);
    nz = gsl_vector_get(xiL,5);

    for (i=0;i<N;i++) {
        rx = act_params.rx(Rad,i,N,L);
        ry = act_params.ry(Rad,i,N,L);
        f = -act_params.force(sys_params,i,n-1,L);

        s_norm = sqrt((pow(nz + ox*ry - oy*rx,2.0) + pow(ny + oz*rx,2.0) + pow(nx - oz*ry,2.0)));

        gsl_vector_set(F,0,gsl_vector_get(F,0)+f*(ry*(nz + ox*ry - oy*rx))/s_norm);
        gsl_vector_set(F,1,gsl_vector_get(F,1)+f*(-rx*(nz + ox*ry - oy*rx))/s_norm);
        gsl_vector_set(F,2,gsl_vector_get(F,2)+f*(rx*(ny + oz*rx) - ry*(nx - oz*ry))/s_norm);
        gsl_vector_set(F,3,gsl_vector_get(F,3)+f*(nx - oz*ry)/s_norm);
        gsl_vector_set(F,4,gsl_vector_get(F,4)+f*(ny + oz*rx)/s_norm);
        gsl_vector_set(F,5,gsl_vector_get(F,5)+f*(nz + ox*ry - oy*rx)/s_norm);
    }


    gsl_vector_sub(xiL,xi_ref);
    gsl_blas_dgemv(CblasNoTrans,1.0,K,xiL,-1.0,F);

    for (i=0;i<6;i++) {
        gsl_vector_set(conditions,6*(n-1)+i,gsl_vector_get(F,i));
    }

    gsl_vector_free(xiL);
    gsl_vector_free(F);
    gsl_matrix_free(K);
    gsl_matrix_free(eulerCond);

}

double tcaForceBase(int k, double delta, double temp) {
    //have to modify the coefficients of the polynomials for the forces with c_i/k^(i-1)

    /* the force equation for the TCAs */

    double l = 0.13;
    double L = 0.04;
    double phi0 = 2*M_PI*100;
    double D = 1.7e-3;
    int N = 3;

    double rho = 4.1582e-05+9.8218e-07*(temp+25);
    double alpha0 = asin(L/l);
    double n = 2*cos(alpha0)*l/(D*M_PI);
    double theta = M_PI*(1/2-1/N);
    double d0 = 2*D*L/l;
    double d_thread0 = d0/(1+1/cos(theta));
    double d_thread = d_thread0*(1+rho*temp);
    double d = d0*(1+rho*temp);
    double E = (560.2023-(temp+25)*1.7771)*pow(10.0,6);
    double nu = 0.45;
    double G = E/(2*(1+nu));
    double J = N*M_PI*pow(d_thread,4)/16*(1/2+pow(1/cos(theta),2));
    double I = J/2;
    double A_a = N*M_PI*pow(d_thread,2)/4;
    double Tau = (J*G*phi0/l)*(1-d_thread/d_thread0);
    
    //rewrite for f11 and f12
    double a0 = (2*pow(sin(alpha0),2)/E-2*pow(sin(alpha0),4)/E-2*pow(sin(alpha0),2)/G+pow(sin(alpha0),4)/G+1/G)/k;
    double a1 = 4*sin(alpha0)/(E*l)-4*sin(alpha0)/(G*l)-8*pow(sin(alpha0),3)/(E*l)+4*pow(sin(alpha0),3)/(G*l);
    double a2 = (2/(E*pow(l,2))-2/(G*pow(l,2))-12*pow(sin(alpha0),2)/(E*pow(l,2))+6*pow(sin(alpha0),2)/(G*pow(l,2)))*k;
    double a3 = (4*sin(alpha0)/(G*pow(l,3))-8*sin(alpha0)/(E*pow(l,3)))*pow(k,2);
    double a4 = (1/(G*pow(l,4))-2/(E*pow(l,4)))*pow(k,3);
    
    double b0 = (pow(sin(alpha0),2)/E-pow(sin(alpha0),2)/G+1/G)/k;
    double b1 = 2*sin(alpha0)/(E*l)-2*sin(alpha0)/(G*l);
    double b2 = (1/(E*pow(l,2))-1/(G*pow(l,2)))*k;
    
    double c0 = (1-pow(sin(alpha0),2))/k;
    double c1 = -2*sin(alpha0)/l;
    double c2 = (-1/pow(l,2))*k;
    
    double f1 = (8*pow(l,3)/(pow(n,2)*pow(M_PI,3)*pow(d_thread,4)))*(a0+a1*delta+a2*pow(delta,2)+a3*pow(delta,3)+a4*pow(delta,4))+(8*l/(2*M_PI*pow(d_thread,2)))*(b0+b1*delta+b2*pow(delta,2));
    double f2 = (8*pow(l,2)/(n*pow(M_PI,2)*G*pow(d_thread,4)))*(c0+c1*delta+c2*pow(delta,2));
    
    double F = (delta+f2*Tau)/f1;
    return F;
    

}

double tcaForce(struct SystemParameters sys_params, int i, int j, double s) {
    return tcaForceBase(sys_params.n, gsl_matrix_get(sys_params.extra,j,i), gsl_vector_get(sys_params.q,i));
}

double cableForce(struct SystemParameters sys_params, int i, int j, double s) {
    return gsl_vector_get(sys_params.q,i);
}

double rx_fun(double R, int i, int N, double s) {
    return R*cos(2*i*M_PI/N);
}

double ry_fun(double R, int i, int N, double s) {
    return R*sin(2*i*M_PI/N);
}



//given eta, xi_prev, eta_prev do the entire computation
void dynamicsComputation(struct SystemParameters sys_params, struct BodyParameters body_params, struct ActuatorParameters act_params,
        const gsl_matrix *eta, const gsl_matrix *eta_prev, const gsl_matrix *xi_prev, //provided values
        gsl_matrix *g, gsl_matrix *xi, gsl_vector *conditions) { //the matrices to be filled


    //first step is to gather the inputs properly
    int n = sys_params.n;

    //do the integration
    integrate(sys_params, body_params, act_params, eta, xi_prev, g, xi);

    //do the differentiation
    gsl_matrix *xi_dot = gsl_matrix_alloc(n,6);
    approxDerivatives(sys_params, xi, xi_dot);
    //setup the actuator parameters (force needs to be specified here if delta needs to be computed)
    //ignore for cable for now

    //compute teim derivative
    gsl_matrix *eta_der = gsl_matrix_alloc(n,6);
    timeDerivative(sys_params, body_params, act_params, g, xi, eta, xi_dot, eta_der);

    //compute boundary conditions
    boundaryConditions(sys_params, body_params, act_params, eta, eta_prev, eta_der, xi, conditions);

    gsl_matrix_free(xi_dot);
    gsl_matrix_free(eta_der);
}


/*
 *Wrappers around common arithmetic functions
 *never actually using these
 */

//A*B = C
void mmMult(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C) {
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,B,0.0,C);
}
//A*b=c
void mvMult(gsl_matrix *A, gsl_vector *b, gsl_vector *c) {
    gsl_blas_dgemv(CblasNoTrans,1.0,A,b,0.0,c);
}


//for root finding need a function of (vector x, params, vector f)
//x is the input, params are other properties, and f is the output
//so in the dynamics case x is xi(s=0) and eta(s/=0) and f are the conditions
//the parameters have to include the matrices to fill

//wrap the structs together
/*
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
*/

int dynamicsFunction(const gsl_vector *x, void *params, gsl_vector *conditions) {
    struct SimulationParameters *sim_params = (struct SimulationParameters *) params;
    struct BodyParameters body_params = sim_params->body_params;
    struct SystemParameters sys_params = sim_params->sys_params;
    struct ActuatorParameters act_params = sim_params->act_params;

    int n = sys_params.n;
    gsl_matrix *g = sim_params->g;
    gsl_matrix *xi = sim_params->xi;
    gsl_matrix *xi_prev = sim_params->xi_prev;
    gsl_matrix *eta_prev = sim_params->eta_prev;
    gsl_matrix *eta = sim_params->eta;


    //get the stuff out

    //rearrange the input, x -> (xi0, eta);
    int i, j;
    for (i=1;i<n;i++) {
        for (j=0;j<6;j++) {
            gsl_matrix_set(eta,i,j,gsl_vector_get(x,i*6+j));
        }
    }
    for (i=0;i<6;i++) {
        gsl_vector_set(sys_params.xi0,i,gsl_vector_get(x,i));
    }


    //run the computation
    dynamicsComputation(sys_params, body_params, act_params, eta, eta_prev, xi_prev, g, xi, conditions);

    return GSL_SUCCESS;
}

//the dynamics stepping function
//given the current solution and the system parameters move to the next step
void stepDynamics(struct SimulationParameters sim_params) {

    //need to move eta_prev and xi_prev forward

    struct BodyParameters body_params = sim_params.body_params;
    struct SystemParameters sys_params = sim_params.sys_params;
    struct ActuatorParameters act_params = sim_params.act_params;

    int n = sys_params.n;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t iter=0;
    int i,j;

    const size_t sys_size = 6*n;
    gsl_multiroot_function f = {&dynamicsFunction, sys_size, &sim_params};

    gsl_vector *x = gsl_vector_alloc(sys_size);

    //initialize the guessed solution
    for (j=0;j<6;j++) {
        gsl_vector_set(x,j,gsl_vector_get(body_params.xi_ref,j));
    }
    for (i=1;i<n;i++) {
        for (j=0;j<6;j++) {
            gsl_vector_set(x,i*6+j,gsl_matrix_get(sim_params.eta,i,j));
        }
    }


    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, sys_size);
    gsl_multiroot_fsolver_set(s, &f, x);

    do
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        if (status)
            break;
        //printf("Iteration: %d\n",iter);

        status = gsl_multiroot_test_residual(s->f, 1e-7);
        //printf("status: %s\n", gsl_strerror (status));
    }
    while (status == GSL_CONTINUE && iter<1000);

/*
    printf("The solution:\n");
    printVector(s->x);
    printf("The conditions:\n");
    printVector(s->f);

    printf("Condition norm: %f\n", gsl_blas_dnrm2(s->f));



    printMatrix(sim_params.xi);
    printMatrix(sim_params.g);
    printMatrix(sim_params.eta);
    printMatrix(sim_params.eta_prev);
    printMatrix(sim_params.xi_prev);
*/

    //free up everything
    gsl_multiroot_fsolver_free(s);

    gsl_vector_free(x);
}




