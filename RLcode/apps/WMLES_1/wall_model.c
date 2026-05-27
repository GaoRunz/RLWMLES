//implement a simplistic wall-model based on log-law 
//
//
#include<stdio.h>
#include<math.h>

void wall_model(double u_in, double *tau_w, double height, double nu)
{
     int i; 
     double u_tau, y0, tau_w1, tau_w0 = 1.; 

     // a simple iteration
     for(i=0; i<100; i++) 
     {
	 
	 u_tau = sqrt(tau_w0); 
	 y0 = nu/(9.*u_tau);
	 u_tau = u_in * 0.41 / (log(height/y0)); 
	 tau_w1 = u_tau * u_tau; 

	 if(fabs(tau_w0 - tau_w1) < 0.001) break; 
	 tau_w0 = tau_w1; 
     }

     (*tau_w) = tau_w1;
}

//evaluation turbulent visocisty based on the rans model
void eval_mut_rans(double u_tau, double disc, double nu, double omega, double *mut)
{
 
     double lmix;    //mixing length 
     double y_plus = disc * u_tau / nu; 
     double tauw; 

//     lmix = 0.40 * disc * (1. - exp(-y_plus / 26.)); 
     	
//     (*mut) = lmix * lmix * fabs(omega); 

     //use the larsson formula
     (*mut) = 0.41 * disc * u_tau * pow((1. - exp(-y_plus/17.)), 2.);  
     //YL's new formulation in order to incorporate 
     (*mut) = (nu + (*mut)) * (1. - disc) - nu; 

     //(*mut) = (nu + (*mut)) * (1. - disc) - nu;
}
