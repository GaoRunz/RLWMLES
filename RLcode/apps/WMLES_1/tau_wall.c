//subroutine developed to evaluate tau_wall 
// email: lvyu@imech.ac.cn 
//
#include<stdio.h>
#include<math.h>

void Tau_Wall_Eval(double muL, double dist, double ymatch, double upp_U, double avg_U, double *log_coef, double *lin_coef, double *avg_coef)
{

        double y_plus = muL * 0.119; 
	
	//a1 + a2 * y_plus + a3 * log(y_plus)   = 0.; 
	//a1 + a2 * 0.5*dist + a3 * (log(dist)-1) = avg_U; 
	//a1 + a2 * dist + a3 * (log(2*dist)-1) = (avg_U + upp_U) * 0.5; 

	// a2 * (0.5*dist - y_plus) + a3 * (log(dist/y_plus) - 1) = avg_U 
	// a2 * (    dist - y_plus) + a3 * (log(2dist/y_plus) - 1) = (avg_U + upp_U) * 0.5
	
	double aa = (avg_U + upp_U) * 0.5 * (0.5*dist - y_plus) - avg_U * (    dist - y_plus); 
	double bb = (log(2.*dist/y_plus) - 1.) * (0.5*dist - y_plus) - (log(dist/y_plus) - 1.) * (    dist - y_plus);

	(*log_coef) = aa/bb; 
	(*lin_coef) = ( avg_U - (*log_coef) * (log(dist/y_plus) - 1.) ) / (0.5*dist - y_plus);
	(*avg_coef) = - (*lin_coef)*y_plus - (*log_coef)*log(y_plus); 

}


void Tau_Wall_Eval_2(double muL, double dist, double ymatch, double upp_U, double avg_U, double *log_coef, double *lin_coef, double *avg_coef)
{
        double y_plus = muL * 0.119; 
/*	
     // sol: a1 * 1 + a2 * x + a3 * logx 
        // strict wall condition 
	a1 + a2 * y_plus + a3 * log(y_plus)   = 0.; 
        
	// match upper boundary 
	//a1 + a2 * dist   + a3 * log(dist)     = upp_U; 
	a1 + a2 * ymatch   + a3 * log(ymatch)   = upp_U; 

	// match cell average
	a1*dist + a2 * 0.5*dist*dist + a3 * (dist*log(dist)-dist) = avg_U*dist; 
	a1 + a2 * 0.5*dist + a3 * (log(dist)-1) = avg_U; 

        a2*(dist-y_plus) + a3 * log(dist/y_plus) = upp_U
	a2* 0.5 * dist   + a3 *  1               = upp_U - avg_U
	
	a3 / (0.5 * dist) * (dist-y_plus) - a3 * log(dist/y_plus) 
	= (upp_U - avg_U)/ (0.5 * dist) * (dist-y_plus) - upp_U
        
	a2*(ymatch-y_plus) + a3 * log(ymatch/y_plus) = upp_U
	a2*(0.5*dist-y_plus) + a3 * (log(dist/y_plus) - 1) = avg_U 
*/	
	double aa = avg_U * (ymatch-y_plus) - upp_U * (0.5*dist-y_plus); 
	double bb = (log(dist/y_plus) - 1) * (ymatch-y_plus) - log(ymatch/y_plus) * (0.5*dist-y_plus);

	(*log_coef) = aa / bb; 
	(*lin_coef) =(upp_U - (*log_coef) * log(ymatch/y_plus) ) / (ymatch-y_plus); 
	(*avg_coef) = - (*lin_coef)*y_plus - (*log_coef)*log(y_plus); 
	//try a different thought
	//with only linear and log mode 
	//a2 * dist   + a3 * log(dist/y_plus + 1)     = upp_U;
	//a2 * 0.5*dist*dist + a3  * (dist*log(dist)-dist) - a3 * dist * log(y_plus) = avg_U * dist
	//a2 * 0.5*dist + a3  * (log(dist)-1) - a3 * log(y_plus) = avg_U 
	//
	//direct solve with analytical approach 
	
 /*       double aa = (upp_U - avg_U)/ (0.5 * dist) * (dist-y_plus) - upp_U; 
	double bb = 1. / (0.5 * dist) * (dist-y_plus) - log(dist/y_plus);
       
        (*log_coef) = aa/bb; (*lin_coef) = (upp_U - avg_U - (*log_coef)) / (0.5 * dist);
	
	(*avg_coef) = - (*lin_coef)*y_plus - (*log_coef)*log(y_plus); 
*/
        	
	// hacking
	//(*log_coef) = avg_U / (log(dist/y_plus) - 1.); 
	
	//printf("%lf\n", (*avg_coef)  + a3 * log(muL*0.119) );
	//return aa / bb; 
}
