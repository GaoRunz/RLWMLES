// a simple PID controller 
//
//
#include <math.h>
#include <stdlib.h>


double SP = 0.;
double Ku = 1.;
double Tu = 20.; 

void PID_output(double signal, double signal_intg, double signal_grad, double *output)
{
   
     const double Kp = 0.6 * Ku; 
     const double Ki = 1.2 * Ku / Tu; 
     const double Kd = 3.0 * Ku * Tu / 40.0;

     double fcn = Kp*signal + Ki*signal_intg + Kd*signal_grad; 
     
     if(fcn < 0.0) fcn = 0.0; 
     if(fcn > 1.0) fcn = 1.0; 
     
     (*output) = fcn;
}

double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}
