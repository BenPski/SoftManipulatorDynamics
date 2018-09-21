#include "mex.h"
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

/*
 *Want to define the general procedure for computing the dynamics, since there is very little difference between 
 *the different kinds of dynamics
 *
 *most of the procedures are them same, the main difference is the force fucntions and the delta integration
 *
 *also want to be able to easily pass around body parameters
 */


/*
 *The first few utility functions
 */

// could make it more generically take matrices/vectors
void approxDerivatives(double dx, double *xi, double *xi_dot, int n) {
    int i, j;
    /* first point and last point approximation */
    for (i=0;i<6;i++) {
        xi_dot[i+6*0] = (xi[i+6*1]-xi[i+6*0])/dx;
        xi_dot[i+6*(n-1)] = (xi[i+6*(n-1)]-xi[i+6*(n-2)])/dx;
    }
    for (i=1;i<n-1;i++) {
        for (j=0;j<6;j++) {
            /* center diff */
            xi_dot[i*6+j] = (xi[j+6*(i+1)]-xi[j+6*(i-1)])/(2*dx);
            /* forward diff */
            //xi_dot[i*6+j] = (xi[j+6*(i+1)]-xi[j+6*(i)])/(dx);
            /* backward diff */
            //xi_dot[i*6+j] = (xi[j+6*(i)]-xi[j+6*(i-1)])/(dx);
        }
    }
}

//the isomorphisms between R3 and skew matrices
//computes the skew matrix of a 3d vector
//not very safe and I have avoided using vectors and just had matrices for everything
void skew(gsl_vector *a, gsl_matrix *a_skew) {
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

void unskew(gsl_matrix *a_skew, gsl_vector *a) {
    gsl_vector_set(a,0,gsl_matrix_get(a_skew,2,1));
    gsl_vector_set(a,1,gsl_matrix_get(a_skew,0,2));
    gsl_vector_set(a,2,gsl_matrix_get(a_skew,1,0));
}


//the isomorphisms for the se(3) algebra, go between R6 and se(3)
//can make use of the skew isomorphisms
void se(gsl_vector *a, gsl_matrix *a_se) {
    gsl_vector_view ang = gsl_vector_subvector(a,0,3); //the angular portion
    gsl_vector_view trans = gsl_vector_subvector(a,3,3); //the translational portion
    
    int i;
    
    skew(&ang.vector,a_se); //get the skew portion
    
    for (i=0; i<3; i++) {
        gsl_matrix_set(a_se,i,3,gsl_vector_get(&trans.vector,i));
    }
    
    for (i=0;i<4;i++) {
        gsl_matrix_set(a_se,3,i,0);
    }
}

void unse(gsl_matrix *a_se, gsl_vector *a) {
    unskew(a_se,a);
    int i;
    for (i=0;i<3;i++) {
        gsl_vector_set(a,3+i,gsl_matrix_get(i,3));
    }
}

//the matrix exponential for SE(3)
//annoying to check for what the representation of the algebra, so just allow a vector
void expSE(gsl_vector *alg, gsl_matrix *group) {

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