#include "mex.h"
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

/*

 * This dynamics uses a discretization of xi_prime = eta_dot + ad_xi * eta to compute the xi at the current time step from the provided eta and prior xi

 * this scheme behaves much more intuitively and stabler than the prior dynamics scheme
 * for now only use the implicit euler scheme instead of a higher order implicit solver
 *
 *Now for the TCA computations, trying to realize the locality of the force
 *  the force is dependent on the local strain and the global force value is incorrect
 *To do this delta is no longer a single value, but a value defined at every discrete point along the body
 *  d = (norm(xi*r)-1)*Ds
 *
 *The modifications of the algorithm now needs
 *  to compute the discrete values of delta in the integration step
 *  individual force values in the teim derivative step
 *  the tip conditoin force is the last force
 *  the coefficients for the force in the TCA change based on the discrete points

*/


double tcaForce(int k, double delta, double temp) {
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
    double a0 = (2*pow(sin(alpha0),2)/E-2*pow(sin(alpha0),4)/E-2*pow(sin(alpha0),2)/G+pow(sin(alpha0),4)/G+1/G)*k;
    double a1 = 4*sin(alpha0)/(E*l)-4*sin(alpha0)/(G*l)-8*pow(sin(alpha0),3)/(E*l)+4*pow(sin(alpha0),3)/(G*l);
    double a2 = (2/(E*pow(l,2))-2/(G*pow(l,2))-12*pow(sin(alpha0),2)/(E*pow(l,2))+6*pow(sin(alpha0),2)/(G*pow(l,2)))/k;
    double a3 = (4*sin(alpha0)/(G*pow(l,3))-8*sin(alpha0)/(E*pow(l,3)))/pow(k,2);
    double a4 = (1/(G*pow(l,4))-2/(E*pow(l,4)))/pow(k,3);
    
    double b0 = (pow(sin(alpha0),2)/E-pow(sin(alpha0),2)/G+1/G)*k;
    double b1 = 2*sin(alpha0)/(E*l)-2*sin(alpha0)/(G*l);
    double b2 = (1/(E*pow(l,2))-1/(G*pow(l,2)))/k;
    
    double c0 = (1-pow(sin(alpha0),2))*k;
    double c1 = -2*sin(alpha0)/l;
    double c2 = (-1/pow(l,2))/k;
    
    double f1 = (8*pow(l,3)/(pow(n,2)*pow(M_PI,3)*pow(d_thread,4)))*(a0+a1*delta+a2*pow(delta,2)+a3*pow(delta,3)+a4*pow(delta,4))+(8*l/(2*M_PI*pow(d_thread,2)))*(b0+b1*delta+b2*pow(delta,2));
    double f2 = (8*pow(l,2)/(n*pow(M_PI,2)*G*pow(d_thread,4)))*(c0+c1*delta+c2*pow(delta,2));
    
    double F = (delta+f2*Tau)/f1;
    return F;
    

}
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
void skew(gsl_matrix *a, gsl_matrix *a_skew) {
    gsl_matrix_set(a_skew,0,0,0);
    gsl_matrix_set(a_skew,0,1,-gsl_matrix_get(a,2,0));
    gsl_matrix_set(a_skew,0,2,gsl_matrix_get(a,1,0));
    gsl_matrix_set(a_skew,1,0,gsl_matrix_get(a,2,0));
    gsl_matrix_set(a_skew,1,1,0);
    gsl_matrix_set(a_skew,1,2,-gsl_matrix_get(a,0,0));
    gsl_matrix_set(a_skew,2,0,-gsl_matrix_get(a,1,0));
    gsl_matrix_set(a_skew,2,1,gsl_matrix_get(a,0,0));
    gsl_matrix_set(a_skew,2,2,0);
}

void expSE(double *xi, double *res) {
    /* computes the matrix exponential for the se(3) algebra, expects the vector representation of the algebra
    */
    double w[3];
    double v[3];
    double w_skew[9];
    double a,b,c,d;
    double wx,wy,wz,vx,vy,vz;
    double w_norm;
    /*
    [ 1 - (b*(wy^2 + wz^2))/2,      (b*wx*wy)/2 - a*wz,      a*wy + (b*wx*wz)/2, vz*((b*wy)/2 + c*wx*wz) - vy*((b*wz)/2 - c*wx*wy) - vx*(c*(wy^2 + wz^2) - 1)]
    [      a*wz + (b*wx*wy)/2, 1 - (b*(wx^2 + wz^2))/2,      (b*wy*wz)/2 - a*wx, vx*((b*wz)/2 + c*wx*wy) - vy*(c*(wx^2 + wz^2) - 1) - vz*((b*wx)/2 - c*wy*wz)]
    [      (b*wx*wz)/2 - a*wy,      a*wx + (b*wy*wz)/2, 1 - (b*(wx^2 + wy^2))/2, vy*((b*wx)/2 + c*wy*wz) - vx*((b*wy)/2 - c*wx*wz) - vz*(c*(wx^2 + wy^2) - 1)]
    [                       0,                       0,                       0,                                                                            1]
     */
    wx = xi[0];
    wy = xi[1];
    wz = xi[2];
    vx = xi[3];
    vy = xi[4];
    vz = xi[5];

    if (wx==0 && wy==0 && wz==0) {
        a = b = 1;
        c = 0;
    } else {
        w_norm = sqrt(wx*wx+wy*wy+wz*wz);
        a = sin(w_norm)/w_norm;
        b = 2*(1-cos(w_norm))/(w_norm*w_norm);
        c = (1-a)/(w_norm*w_norm);
    }



    res[0] = 1 - (b*(wy*wy + wz*wz))/2;
    res[1] = (b*wx*wy)/2 - a*wz;
    res[2] = a*wy + (b*wx*wz)/2;
    res[3] = vz*((b*wy)/2 + c*wx*wz) - vy*((b*wz)/2 - c*wx*wy) - vx*(c*(wy*wy + wz*wz) - 1);
    res[4] = a*wz + (b*wx*wy)/2;
    res[5] = 1 - (b*(wx*wx + wz*wz))/2;
    res[6] = (b*wy*wz)/2 - a*wx;
    res[7] = vx*((b*wz)/2 + c*wx*wy) - vy*(c*(wx*wx + wz*wz) - 1) - vz*((b*wx)/2 - c*wy*wz);
    res[8] = (b*wx*wz)/2 - a*wy;
    res[9] = a*wx + (b*wy*wz)/2;
    res[10] = 1 - (b*(wx*wx + wy*wy))/2;
    res[11] = vy*((b*wx)/2 + c*wy*wz) - vx*((b*wy)/2 - c*wx*wz) - vz*(c*(wx*wx + wy*wy) - 1);
    res[12] = res[13] = res[14] = 0;
    res[15] = 1;

}


void integrateStates(double dt, double ds, double *r, double *eta, double *eta_prev, double *xi_prev, double *xi, double *g, double *delta, int N, int n) {
    /*

    integrate over eta to get the values for
        xi
        g
        delta

    uses the discretization of xi_prime = eta_dot + ad-xi*eta

    this also assumes that xi and eta both already have the correct initial values (mostly important for xi and splitting of guessed values)
    */
    double rx,ry, wx,wy,wz,vx,vy,vz, ox,oy,oz,nx,ny,nz;
    double xi_acc[6], xi_acc_temp[6], g_temp[16];

    /* initial g is R=I, p = 0 */
    g[0] = g[4] = g[8] = 1;
    g[1] = g[2] = g[3] = g[5] = g[6] = g[7] = g[9] = g[10] = g[11] = 0;

    int i,j;
    /* initial xi */
    for (i=0;i<6;i++) {
        /*xi[i] = (dt/ds)*(eta[i]+eta_prev[i]);*/
        /* the first bit of the guess for eta is the base value for xi and the initial eta should be zero */
        xi[i] = eta[i];
        eta[i] = 0;
        xi_acc[i] = 0;
        
    }
    wx = xi[0];
    wy = xi[1];
    wz = xi[2];
    vx = xi[3];
    vy = xi[4];
    vz = xi[5];
    for (j=0;j<N;j++) {
        rx = r[j*N];
        ry = r[j*N+1];
        delta[j] = ds*(sqrt(pow(vx - ry*wz,2.0)+pow(vy + rx*wz,2.0)+pow(vz - rx*wy + ry*wx,2.0))-1);
    }
    
    for (i=1;i<n;i++) {
        /* compute xi and xi_acc */
        // compute the adjoint term separately
        for (j=0;j<6;j++) {
            xi[6*i+j] = xi_prev[6*i+j] + (dt/ds)*(eta[6*i+j]-eta[6*(i-1)+j]);
        }
        //the adjoint portion
        /*
        oy*wz - oz*wy
        oz*wx - ox*wz
        ox*wy - oy*wx
        ny*wz - nz*wy + oy*vz - oz*vy
        nz*wx - nx*wz - ox*vz + oz*vx
        nx*wy - ny*wx + ox*vy - oy*vx
        */

        ox = xi_prev[6*i+0];
        oy = xi_prev[6*i+1];
        oz = xi_prev[6*i+2];
        nx = xi_prev[6*i+3];
        ny = xi_prev[6*i+4];
        nz = xi_prev[6*i+5];

        wx = eta[6*i+0];
        wy = eta[6*i+1];
        wz = eta[6*i+2];
        vx = eta[6*i+3];
        vy = eta[6*i+4];
        vz = eta[6*i+5];

        xi[6*i+0] += dt*(oy*wz - oz*wy);
        xi[6*i+1] += dt*(oz*wx - ox*wz);
        xi[6*i+2] += dt*(ox*wy - oy*wx);
        xi[6*i+3] += dt*(ny*wz - nz*wy + oy*vz - oz*vy);
        xi[6*i+4] += dt*(nz*wx - nx*wz - ox*vz + oz*vx);
        xi[6*i+5] += dt*(nx*wy - ny*wx + ox*vy - oy*vx);

        for (j=0;j<6;j++) {
            xi_acc[j] = xi_acc[j] + xi[6*(i-1)+j]+xi[6*i+j];
            xi_acc_temp[j] = (ds/2.0)*xi_acc[j];
        }

        /* compute delta */
        wx = xi[i*6+0];
        wy = xi[i*6+1];
        wz = xi[i*6+2];
        vx = xi[i*6+3];
        vy = xi[i*6+4];
        vz = xi[i*6+5];

        for (j=0;j<N;j++) {
            rx = r[j*N];
            ry = r[j*N+1];
            delta[i*N+j] = ds*(sqrt(pow(vx - ry*wz,2.0)+pow(vy + rx*wz,2.0)+pow(vz - rx*wy + ry*wx,2.0))-1);
        }

        /* compute g */

        expSE(xi_acc_temp,g_temp);
        g[12*i+0] = g_temp[4*0+0];
        g[12*i+1] = g_temp[4*0+1];
        g[12*i+2] = g_temp[4*0+2];
        g[12*i+3] = g_temp[4*1+0];
        g[12*i+4] = g_temp[4*1+1];
        g[12*i+5] = g_temp[4*1+2];
        g[12*i+6] = g_temp[4*2+0];
        g[12*i+7] = g_temp[4*2+1];
        g[12*i+8] = g_temp[4*2+2];

        g[12*i+9] = g_temp[4*0+3];
        g[12*i+10] = g_temp[4*1+3];
        g[12*i+11] = g_temp[4*2+3];

    }
    
}

void timeDerivative(double r[], double *q, double *delta, double *y, double *xi_dot, double *eta_der, int N, int m) {
    int i, j, k, stateSize;
    double rx, ry, D, Rad, E, G, I, J, A, p, g;
    double R0,R1,R2,R3,R4,R5,R6,R7,R8, px,py,pz, ox,oy,oz,nx,ny,nz, wx,wy,wz,vx,vy,vz, F;
    double ox_dot,oy_dot,oz_dot,nx_dot,ny_dot,nz_dot, wx_dot,wy_dot,wz_dot,vx_dot,vy_dot,vz_dot;
    double temp_vec[6];
    double F_vec[6];
    double pi_dot_norm;
    stateSize = (9+3+6+6);

    gsl_matrix *R = gsl_matrix_alloc(3,3);
    gsl_matrix *xi = gsl_matrix_alloc(6,1);
    gsl_matrix *eta = gsl_matrix_alloc(6,1);
    gsl_matrix *K = gsl_matrix_calloc(6,6);
    gsl_matrix *gamma = gsl_matrix_calloc(6,6);
    gsl_matrix *o = gsl_matrix_alloc(3,1);
    gsl_matrix *o_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *n = gsl_matrix_alloc(3,1);
    gsl_matrix *n_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *w = gsl_matrix_alloc(3,1);
    gsl_matrix *w_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *v = gsl_matrix_alloc(3,1);
    gsl_matrix *v_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *r_vec = gsl_matrix_alloc(3,1);
    gsl_matrix *r_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *xi_der = gsl_matrix_alloc(6,1);
    gsl_matrix *o_der = gsl_matrix_alloc(3,1);
    gsl_matrix *n_der = gsl_matrix_alloc(3,1);
    gsl_matrix *ad_xi = gsl_matrix_alloc(6,6);
    gsl_matrix *ad_eta = gsl_matrix_alloc(6,6);
    gsl_matrix *temp_mat = gsl_matrix_calloc(6,1);
    gsl_matrix *temp_mat2 = gsl_matrix_calloc(6,1);
    gsl_matrix *temp_mat3 = gsl_matrix_calloc(6,1);
    gsl_matrix *temp_mat4 = gsl_matrix_calloc(3,1);
    gsl_matrix *temp_mat5 = gsl_matrix_calloc(3,1);
    gsl_matrix *temp_mat6 = gsl_matrix_calloc(3,1);
    gsl_matrix *p_dot = gsl_matrix_alloc(3,1);
    gsl_matrix *p_dot_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *p_ddot = gsl_matrix_alloc(3,1);


    E = 37.8e3;
    G = E/3.0;
    D = 1e-2;
    Rad = D/2;
    A = M_PI*Rad*Rad;
    I = (M_PI/4.0)*pow(Rad,4.0);
    J = 2.0*I;
    g = -9.8;
    //g = 0;
    p = 1000;

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

        gsl_matrix_set(R,0,0,y[i*stateSize+0]);
        gsl_matrix_set(R,0,1,y[i*stateSize+1]);
        gsl_matrix_set(R,0,2,y[i*stateSize+2]);
        gsl_matrix_set(R,1,0,y[i*stateSize+3]);
        gsl_matrix_set(R,1,1,y[i*stateSize+4]);
        gsl_matrix_set(R,1,2,y[i*stateSize+5]);
        gsl_matrix_set(R,2,0,y[i*stateSize+6]);
        gsl_matrix_set(R,2,1,y[i*stateSize+7]);
        gsl_matrix_set(R,2,2,y[i*stateSize+8]);


        for (j=0;j<6;j++) {
            gsl_matrix_set(xi,j,0,y[i*stateSize+12+j]);
            gsl_matrix_set(eta,j,0,y[i*stateSize+18+j]);
            gsl_matrix_set(xi_der,j,0,xi_dot[i*6+j]);
        }
        for (j=0;j<3;j++) {
            gsl_matrix_set(o,j,0, gsl_matrix_get(xi,j,0));
            gsl_matrix_set(w,j,0, gsl_matrix_get(eta,j,0));
            gsl_matrix_set(n,j,0, gsl_matrix_get(xi,j+3,0));
            gsl_matrix_set(v,j,0, gsl_matrix_get(eta,j+3,0));
            gsl_matrix_set(o_der,j,0, gsl_matrix_get(xi_der,j,0));
            gsl_matrix_set(n_der,j,0, gsl_matrix_get(xi_der,j+3,0));
        }
        skew(o,o_skew);
        skew(n,n_skew);
        skew(v,v_skew);
        skew(w,w_skew);
        for (j=0;j<3;j++) {
            for (k=0;k<3;k++) {
                gsl_matrix_set(ad_xi,j,k,gsl_matrix_get(o_skew,j,k));
                gsl_matrix_set(ad_xi,j+3,k+3,gsl_matrix_get(o_skew,j,k));
                gsl_matrix_set(ad_xi,j+3,k,gsl_matrix_get(n_skew,j,k));
                gsl_matrix_set(ad_xi,j,k+3,0);

                gsl_matrix_set(ad_eta,j,k,gsl_matrix_get(w_skew,j,k));
                gsl_matrix_set(ad_eta,j+3,k+3,gsl_matrix_get(w_skew,j,k));
                gsl_matrix_set(ad_eta,j+3,k,gsl_matrix_get(v_skew,j,k));
                gsl_matrix_set(ad_eta,j,k+3,0);
            }
        }


        ox=y[i*stateSize+12];
        oy=y[i*stateSize+13];
        oz=y[i*stateSize+14];
        nx=y[i*stateSize+15];
        ny=y[i*stateSize+16];
        nz=y[i*stateSize+17];

        wx=y[i*stateSize+18];
        wy=y[i*stateSize+19];
        wz=y[i*stateSize+20];
        vx=y[i*stateSize+21];
        vy=y[i*stateSize+22];
        vz=y[i*stateSize+23];

        ox_dot = xi_dot[i*6+0];
        oy_dot = xi_dot[i*6+1];
        oz_dot = xi_dot[i*6+2];
        nx_dot = xi_dot[i*6+3];
        ny_dot = xi_dot[i*6+4];
        nz_dot = xi_dot[i*6+5];



        /* eta' */
        /* K*xi_dot */
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,K,xi_der,0.0,temp_mat);
        temp_vec[0] = E*I*ox_dot;
        temp_vec[1] = E*I*oy_dot;
        temp_vec[2] = G*J*oz_dot;
        temp_vec[3] = A*G*nx_dot;
        temp_vec[4] = A*G*ny_dot;
        temp_vec[5] = A*E*nz_dot;



        /* -ad_xi_star*K*(xi-xi_ref) */
        for (j=0;j<6;j++) {
            gsl_matrix_set(temp_mat2,j,0,gsl_matrix_get(xi,j,0));
        }
        gsl_matrix_set(temp_mat2,5,0,gsl_matrix_get(temp_mat2,5,0)-1);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,K,temp_mat2,0.0,temp_mat3);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,-1.0,ad_xi,temp_mat3,1.0,temp_mat);

        temp_vec[0] -= A*G*ny*nz - A*E*ny*(nz - 1) + E*I*oy*oz - G*J*oy*oz;
        temp_vec[1] -= A*E*nx*(nz - 1) - A*G*nx*nz - E*I*ox*oz + G*J*ox*oz;
        temp_vec[3] -= A*G*ny*oz - A*E*oy*(nz - 1);
        temp_vec[4] -= A*E*ox*(nz - 1) - A*G*nx*oz;
        temp_vec[5] -= A*G*nx*oy - A*G*ny*ox;


        /* ad_eta_star*gamma*eta */
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,gamma,eta,0.0,temp_mat2);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,ad_eta,temp_mat2,1.0,temp_mat);
        temp_vec[0] += I*p*wy*wz - J*p*wy*wz;
        temp_vec[1] += J*p*wx*wz - I*p*wx*wz;
        temp_vec[3] += A*p*vy*wz - A*p*vz*wy;
        temp_vec[4] += A*p*vz*wx - A*p*vx*wz;
        temp_vec[5] += A*p*vx*wy - A*p*vy*wx;


        /* gravity */
        gsl_matrix_set(temp_mat,3,0, gsl_matrix_get(temp_mat,3,0)+A*g*p*gsl_matrix_get(R,2,0));
        gsl_matrix_set(temp_mat,4,0, gsl_matrix_get(temp_mat,4,0)+A*g*p*gsl_matrix_get(R,2,1));
        gsl_matrix_set(temp_mat,5,0, gsl_matrix_get(temp_mat,5,0)+A*g*p*gsl_matrix_get(R,2,2));
        temp_vec[3] += A*R6*g*p;
        temp_vec[4] += A*R7*g*p;
        temp_vec[5] += A*R8*g*p;


        /* tcas */

        for (j=0;j<6;j++) {
            temp_vec[j]=gsl_matrix_get(temp_mat,j,0);
            F_vec[j] = 0;
        }
        for (j=0;j<N;j++) {
            rx = r[N*j+0];
            ry = r[N*j+1];
            gsl_matrix_set(r_vec,0,0,rx);
            gsl_matrix_set(r_vec,1,0,ry);
            gsl_matrix_set(r_vec,2,0,0);
            skew(r_vec,r_skew);

            /* sum -F*[skew(r)*R'*skew(p_dot)^2*p_ddot/norm(p_dot)^3;R'*skew(p_dot)^2*p_ddot/norm(p_dot)^3] */
            /* p_dot = R(n+skew(o)*r) */
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,o_skew,r_vec,0.0,temp_mat4);
            gsl_matrix_add(temp_mat4,n);
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,R,temp_mat4,0.0,p_dot);
            skew(p_dot,p_dot_skew);

            /* p_ddot = R*(skew(o)*(n+skew(w)*r)+n_dot+skew(o_dot)*r) */
            //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0,r_skew,o_der,1.0,n_der);
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0,r_skew,o_der,0.0,temp_mat6);
            gsl_matrix_add(temp_mat6,n_der);
            //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,o_skew,r_vec,1.0,n);
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,o_skew,r_vec,0.0,temp_mat4);
            gsl_matrix_add(temp_mat4,n);
            //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,o_skew,n,1.0,n_der);
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,o_skew,temp_mat4,0.0,temp_mat5);
            gsl_matrix_add(temp_mat5,temp_mat6);
            //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,R,n_der,0.0,p_ddot);
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,R,temp_mat5,0.0,p_ddot);

            /* R'*skew(p_dot)*skew(p_dot)*p_ddot */
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,p_dot_skew,p_ddot,0.0,temp_mat4);
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,p_dot_skew,temp_mat4,0.0,p_ddot);
            gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,R,p_ddot,0.0,temp_mat4);
            /* skew(r) * res */
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,r_skew,temp_mat4,0.0,p_ddot);

            F = tcaForce(m,delta[i*N+j],q[j]);
            pi_dot_norm = pow(sqrt(pow(nz + ox*ry - oy*rx,2.0) + pow(ny + oz*rx,2.0) + pow(nx - oz*ry,2.0)),3.0);

            F_vec[0] += (-F/pi_dot_norm)*gsl_matrix_get(p_ddot,0,0);
            F_vec[1] += (-F/pi_dot_norm)*gsl_matrix_get(p_ddot,1,0);
            F_vec[2] += (-F/pi_dot_norm)*gsl_matrix_get(p_ddot,2,0);
            F_vec[3] += (-F/pi_dot_norm)*gsl_matrix_get(temp_mat4,0,0);
            F_vec[4] += (-F/pi_dot_norm)*gsl_matrix_get(temp_mat4,1,0);
            F_vec[5] += (-F/pi_dot_norm)*gsl_matrix_get(temp_mat4,2,0);

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
            eta_der[i*6+j] = temp_vec[j];
        }
    }

    gsl_matrix_free(R);
    gsl_matrix_free(xi);
    gsl_matrix_free(eta);
    gsl_matrix_free(K);
    gsl_matrix_free(gamma);
    gsl_matrix_free(o);
    gsl_matrix_free(o_skew);
    gsl_matrix_free(n);
    gsl_matrix_free(n_skew);
    gsl_matrix_free(w);
    gsl_matrix_free(w_skew);
    gsl_matrix_free(v);
    gsl_matrix_free(v_skew);
    gsl_matrix_free(r_vec);
    gsl_matrix_free(r_skew);
    gsl_matrix_free(xi_der);
    gsl_matrix_free(o_der);
    gsl_matrix_free(n_der);
    gsl_matrix_free(ad_xi);
    gsl_matrix_free(ad_eta);
    gsl_matrix_free(temp_mat);
    gsl_matrix_free(temp_mat2);
    gsl_matrix_free(temp_mat3);
    gsl_matrix_free(temp_mat4);
    gsl_matrix_free(temp_mat5);
    gsl_matrix_free(temp_mat6);
    gsl_matrix_free(p_dot);
    gsl_matrix_free(p_dot_skew);
    gsl_matrix_free(p_ddot);



}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    /*
    will be provided q, xi_t+1, eta_t+1, xi_t, and eta_t
    necessary to have the previous answer to compute implicit euler equation
    the number of discretization points is the length of xi and eta is one less as it is known that eta(0) = 0 (INCORRECT, I could probably make it work this way, but needs to specify all of eta)
    output the resulting g (xi and eta are the guessed values to complete the state) and the constraining equations

    expect guesses in the form: (ACTUALLY TRANSPOSED FROM THIS, the indexing works as col+row*width)
    [wx,wy,wz,vx,vy,vz
     wx,wy,wz,vx,vy,vz
     ...]
    so index as row+6*column
    */

    /* parameters */

    int i, j;

    double E, G, D, R, A, I, J, L;
    E = 37.8e3;
    G = E/3.0;
    D = 1e-2;
    R = D/2.0;
    A = M_PI*R*R;
    I = (M_PI/4.0)*pow(R,4.0);
    J = 2.0*I;
    L = 0.04;


    /* actuation values */
    double *q;
    const mwSize *dimq;
    int dimqx, dimqy;
    int N; /* number of actuators */
    dimq = mxGetDimensions(prhs[0]);
    dimqy = (int)dimq[0]; dimqx = (int)dimq[1];
    N = dimqy*dimqx;
    q = mxGetPr(prhs[0]);


    /*

    IMPORTANT: in order to line the inputs up the same with how I pass them into the matlab code the input matrices need to be tranposed

    */

    /* guessed eta */
    double *eta_trans;
    const mwSize *dimEta;
    int dimEtax,dimEtay;
    int n;
    dimEta = mxGetDimensions(prhs[1]);
    dimEtax = (int)dimEta[0]; dimEtay = (int)dimEta[1];
    n = dimEtax;
    eta_trans = mxGetPr(prhs[1]);
    double eta[6*n];
    for (i=0;i<n;i++) {
        for (j=0;j<6;j++) {
            eta[6*i+j] = eta_trans[n*j+i];
        }
    }

    /* eta_prev */
    double eta_prev[6*n], *eta_prev_trans;
    eta_prev_trans = mxGetPr(prhs[2]);

    for (i=0;i<n;i++) {
        for (j=0;j<6;j++) {
            eta_prev[6*i+j] = eta_prev_trans[j*n+i];
        }
    }

    /* xi_prev */
    double xi_prev[6*n], *xi_prev_trans;
    xi_prev_trans = mxGetPr(prhs[3]);

    for (i=0;i<n;i++) {
        for (j=0;j<6;j++) {
            xi_prev[i*6+j] = xi_prev_trans[j*n+i];
        }
    }

    /* time step */
    double dt = mxGetPr(prhs[4])[0];

    /* setup r and dx */
    double ds = L/(n-1);
    double r[N*3];
    for (i=0;i<N;i++) {
        r[i*N+0] = (R/2.0)*cos(i*2.0*M_PI/N);
        r[i*N+1] = (R/2.0)*sin(i*2.0*M_PI/N);
        r[i*N+2] = 0;
    }


    /* integrate stuff */
    double xi[6*n];
    double delta[n*N];
    double g[12*n];
    integrateStates(dt,ds,r,eta,eta_prev,xi_prev,xi,g,delta,N,n);


    /* approx derivatives */
    double xi_dot[n*6];
    approxDerivatives(ds, xi, xi_dot, n);


    /* time derivative */

    /* setup state vector y (y=[R,p,xi,eta])*/
    double y[24*n];
    for (i=0;i<n;i++) {
        y[24*i+0] = g[12*i+0];
        y[24*i+1] = g[12*i+1];
        y[24*i+2] = g[12*i+2];
        y[24*i+3] = g[12*i+3];
        y[24*i+4] = g[12*i+4];
        y[24*i+5] = g[12*i+5];
        y[24*i+6] = g[12*i+6];
        y[24*i+7] = g[12*i+7];
        y[24*i+8] = g[12*i+8];
        y[24*i+9] = g[12*i+9];
        y[24*i+10] = g[12*i+10];
        y[24*i+11] = g[12*i+11];
        y[24*i+12] = xi[6*i+0];
        y[24*i+13] = xi[6*i+1];
        y[24*i+14] = xi[6*i+2];
        y[24*i+15] = xi[6*i+3];
        y[24*i+16] = xi[6*i+4];
        y[24*i+17] = xi[6*i+5];
        y[24*i+18] = eta[6*i+0];
        y[24*i+19] = eta[6*i+1];
        y[24*i+20] = eta[6*i+2];
        y[24*i+21] = eta[6*i+3];
        y[24*i+22] = eta[6*i+4];
        y[24*i+23] = eta[6*i+5];
    }

    double eta_der[6*n];
    timeDerivative(r, q, delta, y, xi_dot, eta_der, N, n);


    /* compute condition equations */
    double eqs[n*6];
    /* implicit euler step */
    for (i=1;i<n;i++) {
        for (j=0;j<6;j++) {
            eqs[i*6+j] = -eta[i*6+j]+eta_prev[i*6+j]+dt*eta_der[j+i*6];
        }
    }


    /* boundary condition */
    double xiL[6];
    double F_vec[6];
    double temp_vec[3];
    double F, rx, ry, ox, oy, oz, nx, ny, nz, s_norm;
    for (i=0;i<6;i++) {
        F_vec[i]=0;
        xiL[i] = xi[(n-1)*6+i];
    }
    ox = xiL[0];
    oy = xiL[1];
    oz = xiL[2];
    nx = xiL[3];
    ny = xiL[4];
    nz = xiL[5];

    for (i=0;i<N;i++) {
        rx = r[i*N+0];
        ry = r[i*N+1];
        F = -tcaForce(n,delta[N*(n-1)+i],q[i]);


        s_norm = sqrt((pow(nz + ox*ry - oy*rx,2.0) + pow(ny + oz*rx,2.0) + pow(nx - oz*ry,2.0)));


        F_vec[0] += F*(ry*(nz + ox*ry - oy*rx))/s_norm;
        F_vec[1] += F*(-rx*(nz + ox*ry - oy*rx))/s_norm;
        F_vec[2] += F*(rx*(ny + oz*rx) - ry*(nx - oz*ry))/s_norm;
        F_vec[3] += F*(nx - oz*ry)/s_norm;
        F_vec[4] += F*(ny + oz*rx)/s_norm;
        F_vec[5] += F*(nz + ox*ry - oy*rx)/s_norm;
    }
    F_vec[0] -= E*I*ox;
    F_vec[1] -= E*I*oy;
    F_vec[2] -= G*J*oz;
    F_vec[3] -= A*G*nx;
    F_vec[4] -= A*G*ny;
    F_vec[5] -= A*E*(nz - 1);
    eqs[0] = F_vec[0];
    eqs[1] = F_vec[1];
    eqs[2] = F_vec[2];
    eqs[3] = F_vec[3];
    eqs[4] = F_vec[4];
    eqs[5] = F_vec[5];

    /* place everything back into matlab */
    double *states, *bcs, *eta_prime;
    mxArray *states_out, *bcs_out, *eta_prime_out;

    states_out = plhs[0] = mxCreateDoubleMatrix(24,n,mxREAL);
    states = mxGetPr(states_out);
    for (i=0;i<24*n;i++) {
        states[i] = y[i];
    }
    bcs_out = plhs[1] = mxCreateDoubleMatrix(6*n,1,mxREAL);
    bcs = mxGetPr(bcs_out);
    for (i=0;i<6*n;i++) {
        bcs[i] = eqs[i];
    }

}
