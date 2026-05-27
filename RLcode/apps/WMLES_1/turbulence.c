/*
*  TURBULENCE  
*  
*  Functions for defining Sub-grid scale (SGS) models for Large-Eddy Simulation.
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>      /* sqrt()       */

//Vremen model
double VremenModel (double *gu, double dx, double dy, double dz)
{
    int i, j, k;
    int ij, ji;
    const int dim=3;
    double V[3], Vol;
    double S, r, t;
    double dx2, dy2, dz2, vreman_coeff=0.07 ;
    double alpha11, alpha12, alpha13, alpha21, alpha22, alpha23, alpha31, alpha32, alpha33;
    double beta11, beta12, beta13, beta22, beta23, beta33, B, tmp;
    
    /////////////////////////////////
    //implementation of vreman model
    Vol = dx * dy * dz; 
    dx2 = pow(Vol, 2.0/3.0);
    dy2 = dx2; dz2 = dx2;
    //dx2 = dx * dx; 
    //dy2 = dy * dy;
    //dz2 = dz * dz; 

    alpha11 = gu[0*dim+0];  //du_dx
    alpha22 = gu[1*dim+1];  //dv_dy
    alpha33 = gu[2*dim+2];  //dw_dz

/*
    //remove trace
    tmp = (alpha11 + alpha22 + alpha33) / 3.0;
    alpha11 -= tmp;
    alpha22 -= tmp;
    alpha33 -= tmp;
*/

    alpha12 = gu[0*dim+1];  //dv_dx
    alpha13 = gu[0*dim+2];  //dw_dx
    alpha23 = gu[1*dim+2];  //dw_dy
    
    alpha21 = gu[1*dim+0];  //du_dy
    alpha31 = gu[2*dim+0];  //du_dz
    alpha32 = gu[2*dim+1];  //dv_dz
    
    beta11 = dx2 * alpha11 * alpha11 + dy2 * alpha12 * alpha12 +
    dz2 * alpha13 * alpha13;
    beta12 = dx2 * alpha11 * alpha21 + dy2 * alpha12 * alpha22 +
    dz2 * alpha13 * alpha23;
    beta13 = dx2 * alpha11 * alpha31 + dy2 * alpha12 * alpha32 +
    dz2 * alpha13 * alpha33;
    beta22 = dx2 * alpha21 * alpha21 + dy2 * alpha22 * alpha22 +
    dz2 * alpha23 * alpha23;
    beta23 = dx2 * alpha21 * alpha31 + dy2 * alpha22 * alpha32 +
    dz2 * alpha23 * alpha33;
    beta33 = dx2 * alpha31 * alpha31 + dy2 * alpha32 * alpha32 +
    dz2 * alpha33 * alpha33;
    
    B = beta11 * beta22 - beta12 * beta12 + beta11 * beta33 -
    beta13 * beta13 + beta22 * beta33 - beta23 * beta23;
    
    B = (B + fabs(B)) * 0.5;
    if(B < 1.e-8) return 0.; 

    S = alpha11 * alpha11 + alpha22 * alpha22 + alpha33 * alpha33 +
    alpha12 * alpha12 + alpha13 * alpha13 + alpha23 * alpha23 +
    alpha21 * alpha21 + alpha31 * alpha31 + alpha32 * alpha32;
    
    return vreman_coeff * sqrt( B / (S + 1.0e-20));
}

double eval_vort_mag(double *gu)
{

    const int dim = 3; 
    double vort; 

    double dw_dy = gu[1*dim+2];
    double dv_dz = gu[2*dim+1];
    double du_dz = gu[2*dim+0];
    double dw_dx = gu[0*dim+2];
    double dv_dx = gu[0*dim+1];
    double du_dy = gu[1*dim+0];

    vort = sqrt(pow(dw_dy-dv_dz,2.) + pow(du_dz-dw_dx,2.) + pow(dv_dx-du_dy,2.)); 
    return vort;

}

