#ifndef GLOBAL_H
#define GLOBAL_H

#include<stdbool.h>

#define real double
//controller 
void PID_output(double signal, double signal_intg, double signal_grad, double *output);
double err_diff=0.0, err_intg=0.0, err_grad=0.0, err_prev=0.0, PIDout;
double randn (double mu, double sigma); 

//void Tau_Wall_Eval(double y_plus, double dist, double ymatch, double upp_U, double avg_U, double *log_coef, double *lin_coef, double *avg_coef);
real VremenModel (real *gu, real dx, real dy, real dz); 
void wall_model(double u_in, double *tau_w, double height, double nu);

//hybrid model; rans part 
void eval_mut_rans(double u_tau, double disc, double nu, double omega, double *mut);
double eval_vort_mag(double *gu);

// MPI buffers
real *Bufp, *Bufn, *pBuf, *nBuf;
unsigned BufCountU, BufCountF; // BufCountU > BufCountF
int
  //numprocs, // number of processes
  //myid,     // identifier of _this_ process
  namelen;
  
extern int numprocs;
extern int myid;
  
  
char 
  processor_name[MPI_MAX_PROCESSOR_NAME];


#define Fx1 xU1 // aliases for x-fluxes
#define Fx2 xU2
#define Fx3 xU3
#define Fx4 xU4
#define Fx5 xU5

//#define Fy1 yU1 // aliases for y-fluxes
//#define Fy2 yU2
//#define Fy3 yU3
//#define Fy4 yU4
//#define Fy5 yU5

#define Fz1 zU1 // aliases for z-fluxes
#define Fz2 zU2
#define Fz3 zU3
#define Fz4 zU4
#define Fz5 zU5

//=== Global variables
unsigned
  LEN, HIG, DEP, // number of cells in x-,y-,z-directions, for 3D box
  LENN, HIGG, DEPP,
  Answer,        // solution continuing flag
  i,  j,  k,
  _i, _j, _k,
  i_, j_, k_,
  i__,  j__,  k__; // indices

real
  // conservative variables at the cells' centroids
  ***U1,  ***U2,  ***U3,  ***U4,  ***U5,
  ***U1p, ***U2p, ***U3p, ***U4p, ***U5p,
  // pointers for substitution
  ***U1_, ***U2_, ***U3_, ***U4_, ***U5_,
  // conservative variables at cell boundaries in x-direction
  ***xU1, ***xU2, ***xU3, ***xU4, ***xU5,
  ***U1x, ***U2x, ***U3x, ***U4x, ***U5x,
  // conservative variables at cell boundaries in y-direction
  ***yU1, ***yU2, ***yU3, ***yU4, ***yU5,
  ***U1y, ***U2y, ***U3y, ***U4y, ***U5y,
  // conservative variables at cell boundaries in z-direction
  ***zU1, ***zU2, ***zU3, ***zU4, ***zU5,
  ***U1z, ***U2z, ***U3z, ***U4z, ***U5z,
  // flux YL 
  ***Fy1, ***Fy2, ***Fy3, ***Fy4, ***Fy5,
  // initial density
  R0,
  // primitive parameters
  R, U, V, W, P, C,
  // conservative variables
  u1, u2, u3, u4, u5,
  // cell spacings
  deltaX, deltaY, deltaZ,
  // time step
  deltaT,
  // convenient ratios
  deltaT_X, deltaT_Y, deltaT_Z,
  _deltaX, _deltaY, _deltaZ,
  _4deltaX, _4deltaY, _4deltaZ,
  _2deltaX, _2deltaY, _2deltaZ;

real
  // statistical variable 
  *avg_Ux, *avg_Uy, *avg_Uz, 
  *rms_Ux, *rms_Uy, *rms_Uz,
  *Rss_UxUy, *Rss_UxUz, *Rss_UyUz; 

int PID_time = 0;
real **blend_dyn, **nwall_vec, ***pid_stat_log;
real blend_glob = 1.; 
real speed_diff = 0.;
bool change_run = false;

char 
  Stage = 1, // indicator of current stage
  str[80];   // buffer string

int
  //step,    // current time step
  s_step,  // sample time step
  //numstep, // overall number of time steps
  f_step;  // frame taking step
  
 extern int step;
 extern int numstep;

long mem; // amount of dynamic memory required for arrays in this process

// handles of "sensors.*" files
FILE *pFsen_bin, *pFsen_txt;


// finite differences of conservative variables
real
  m1, m2, m3, m4, m5,
  w1, w2, w3, w4, w5,
  // finite differences of characteristic variables
  M1, M2, M3, M4, M5,
  W1, W2, W3, W4, W5,
  // modified finite differences of characteristic variables
  _M1, _M2, _M3, _M4, _M5,
  _W1, _W2, _W3, _W4, _W5,
  // transformation matrix [S-1]
  S11, S12, S13,
  S23,
  S41, S42,
  S51, S52,
  // inverted transformation matrix [S-1]
  _S23, _S24, _S25,
  _S33,
  _S43,
  _S53, _S54, _S55,
  // stuff
  aa, bb, cc, dd, ee, ff, gg, hh, ll, mm, nn,
  u_hh, v_hh, w_hh, u_ee, v_ee, w_ee,
  k1, k2, k3, k4, k5,

  //--- For LCS (approximate Riemann solver)
  P_, _P, R_, _R, U_, _U, V_, _V, W_, _W, C_, _C, _T, T_, _E, E_,
  C_p, C_m, C_o, C_p_, C_m_, C_o_,
  _jo, _jp, _jm, jo_, jp_, jm_,
  al_p, al_m, al_o,
  J_p, J_m, J_o,

  //--- For viscous fluxes
  // matrix of velocity derivatives
  du_dx, dv_dx, dw_dx,
  du_dy, dv_dy, dw_dy,
  du_dz, dv_dz, dw_dz,
  // velocity components at surrounding cells
  iu, ui, ju, uj, ku, uk,
  iv, vi, jv, vj, kv, vk,
  iw, wi, jw, wj, kw, wk,
  iU, Ui, jU, Uj, kU, Uk,
  iV, Vi, jV, Vj, kV, Vk,
  iW, Wi, jW, Wj, kW, Wk,
  // mean filtered velocity components and density at the interfaces
  uu, vv, ww, rr,
  // viscous stress tensor
  sigma_xx, sigma_xy, sigma_xz,
  sigma_yx, sigma_yy, sigma_yz,
  sigma_zx, sigma_zy, sigma_zz,
  // heat flux vector
  q_x, q_y, q_z,
  // transport coefficients
  // molecular
  mu_L, lambda_L, Pr_L,
  // subgrid-scale
  mu_T, Pr_T, cp_Pr_T,
  // effective
  mu_E, lambda_E,
  // Smagorinsky constant and a complex CsDD=Cs*d*d, where d=(dX*dY*dZ)^(1/3)
  Cs, CsDD,
  // arrays of near-wall-corrected Smagorinsky complex
  **CsDDx, **CsDDy, **CsDDz,
  // invariant of deforming rate tensor
  _S_,
  // Cartesian coordinates
  x, y, z,
  // distances in wall units, friction velocity, mean wall stress
  d_plus, u_tau, tau_w,
  // accceleration of mass force in source term
  g, gdeltaT,
  // array of sensors
  *sensors,
  M, I, E,    // mean fluxes
  tauY, tauZ; // mean wall stresses


const real 
  K      = 1.4,    // combinations of the ratio of specific heats
  K_1    = 0.4,
  _K_1   = -0.4,
  _05K   = 0.7,
  K_K_1  = 3.5,
  K_1_2  = 0.2,
_1_2_K_1 = 1.25,

  Cv     = 717.75,  // specific constant volume heat capacity, J/(kg*K)
  R_GAS  = 287.10,  // gas constant for air
  KR_GAS = 401.94,  // R_GAS * K

  blendy = 0.6;     // blending factor for y-direction


real  blend  = 1.0;  // flux blending factor; made adjustable
 
real  balance = 1.;    // for boundary adjustment

real ctF1, ctF2, ctF3, ctF4, ctF5;
real gu[9], ***muturb, ***muturb_rans, ***muturb_coef, ***tmp;        
real wm_tau_wall, wm_tau_wall2, alow_x, alow_z, aprime_x, aprime_z, res_x, res_z; 
real ***wc_coeffs;
real *Uavg0, *Vavg0, *Wavg0;
real *Uavg1, *Vavg1, *Wavg1;
real *UU0  , *VV0  , *WW0;
real *UV0  , *VW0  , *UW0;
real accumul_dT; 
real LinearU, DriverF, MassRate;

// enable dynamic adjustment of blend
int blend_adjust_freq; 
#endif
