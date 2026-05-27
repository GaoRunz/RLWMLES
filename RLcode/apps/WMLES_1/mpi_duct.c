/* mpi_duct.c *********************************\
*                                              *
*  - Turbulent Flow in Duct modeled by LES -   *
*        two-stage time differencing           *
*                                              *
*         Parallel version of "duct.c"         *
*          tested under MPICH-1.2.1            *
*                                              *
*  written by Andrei Chernousov                *
*  E-mail: <andrei99@iname.com>                *
*  http://www.geocities.com/andrei_chernousov  *
*  Version 2.2, January 18, 2002               *
\**********************************************/


// Last change: the MPI data transfer scheme made not so ugly.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>

#include <mpi.h>  // Message Passing Interface - MPI
#include "global.h"
//#define real double

#define STRONG_ENFORCE_BC
//#undef  STRONG_ENFORCE_BC
//#define STORE_WALL_CELL
#undef STORE_WALL_CELL
#define ADAWN
//#undef  ADAWN

extern double rec_perb;

// continuation flag
extern int continFlag;
char a;
unsigned nStage = 3;
extern MPI_File fh;
extern MPI_Status status;

extern double state[];
extern double reward;
extern MPI_Comm globleMpicom;
extern int selectStart;

extern double* nwall_vec_ave;


/************
*  ARRAY2D  *   Memory allocation procedure for 2D arrays
************/
real** Array2D(unsigned columns, unsigned rows)
{
	real** x;
	unsigned i;

	if ((x = (real**)malloc(columns * sizeof(real*))) == NULL) {
		fprintf(stderr, "Cannot allocate memory");
		fflush(stderr);
		MPI_Abort(globleMpicom, 1);
	}
	for (i = 0; i < columns; i++)
		if ((x[i] = (real*)malloc(rows * sizeof(real))) == NULL) {
			fprintf(stderr, "Cannot allocate memory");
			fflush(stderr);
			MPI_Abort(globleMpicom, 1);
		}
	return x;

} // end Array2D()


/************
*  ARRAY3D  *   Memory allocation procedure for 3D arrays
************/
real*** Array3D(unsigned columns, unsigned rows, unsigned floors)
{
	real*** x;
	unsigned i, j;

	if ((x = (real***)malloc(columns * sizeof(real**))) == NULL) {
		fprintf(stderr, "Cannot allocate memory");
		fflush(stderr);
		MPI_Abort(globleMpicom, 1);
	}
	for (i = 0; i < columns; i++) {
		if ((x[i] = (real**)malloc(rows * sizeof(real*))) == NULL) {
			fprintf(stderr, "Cannot allocate memory");
			fflush(stderr);
			MPI_Abort(globleMpicom, 1);
		}
		for (j = 0; j < rows; j++)
			if ((x[i][j] = (real*)malloc(floors * sizeof(real))) == NULL) {
				fprintf(stderr, "Cannot allocate memory");
				fflush(stderr);
				MPI_Abort(globleMpicom, 1);
			}
	}
	return x;

} // end Array3D()


/***************
*  INITIALIZE  *    Performs necessary initializations
***************/
void Initialize(void)
{
	// temporal arrays
	real*** fp[10];
	char filename[50];

	unsigned i = 0, jj, kk;
	MPI_File fh;
	real buf;

	real Vd, Wd; // components of disturbed velocity

	MPI_Status status;

	/*
	   // Root process (node w/rank 0) reads data from file "mpi_duct.bin"
	   //   and broadcasts it to others
	  if (myid == 0) {
		// open the file...
		MPI_File_open(MPI_COMM_SELF, "mpi_duct.bin",
			  MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	  }

	  // ... and reads/broadcast data (if myid == 0), otherwise receives it
	  do{

		// syncronization
		MPI_Barrier(globleMpicom);

		// root reads the next data item
		if (myid == 0) {
		  MPI_File_read(fh, &buf, 1, MPI_DOUBLE, &status);
		}

		// root process broadcasts the next data item
		MPI_Bcast(&buf, 1, MPI_DOUBLE, 0, globleMpicom);
		//printf("%d process: after broadcast %dth data item = %f\n",
		//myid+1, i, buf);

		// dispatch to specific parameters
		switch (i) {
		case  0: printf("LEN    = %d\n", LEN    = (unsigned)buf); break;
		case  1: printf("HIG    = %d\n", HIG    = (unsigned)buf); break;
		case  2: printf("DEP    = %d\n", DEP    = (unsigned)buf); break;
		case  3: printf("deltaX = %f\n", deltaX = buf); break;
		case  4: printf("deltaY = %f\n", deltaY = buf); break;
		case  5: printf("deltaZ = %f\n", deltaZ = buf); break;
		case  6: printf("deltaT = %f\n", deltaT = buf); break;
		case  7: printf("Answer = %d\n", Answer = (unsigned)buf); break;
		case  8: printf("numstep= %d\n", numstep= (unsigned)buf); break;
		case  9: printf("mu_L   = %f\n", mu_L   = buf); break;
		case 10: printf("Pr_L   = %f\n", Pr_L   = buf); break;
		case 11: printf("Pr_T   = %f\n", Pr_T   = buf); break;
		case 12: printf("Cs     = %f\n", Cs     = buf); break;
		case 13: printf("g      = %f\n", g      = buf); break;
		case 14: printf("tau_w  = %f\n", tau_w  = buf); break;
		case 15: printf("R0     = %f\n", R0     = buf); break;
		case 16: printf("f_step = %d\n", f_step = buf); break;
		default: break;
		} // end switch
		i++;

	  } while (i < 17);

	  // close file
	  if (myid == 0)
		MPI_File_close(&fh);
	*/
	LEN = 16;//        LEN     Number of cells in length     [-]      //62*5 processors
	HIG = 64;//        HIG     Number of cells in hight      [-]      //45
	DEP = 64;//        DEP     Number of cells in depth      [-]      //45
	deltaX = 0.09817477; //  deltaX  X-step                        [m]      //0.0307692
	deltaY = 0.031250;  //   deltaY  Y-step                        [m]      //0.0222222
	deltaZ = 0.09817477; //  deltaZ  Z-step                        [m]      //0.0222222
	deltaT = 0.00002; //  deltaT  Time step                     [sec]    //1.75e-5
	Answer = 1; //  Answer  Continuing flag (1/0)         [-]      //0/1
	int AnswerPID = 0; //  Load PID state data from file
	numstep = 3305000; //   numstep Number of time steps          [-]      //30000
	mu_L = 0.0005; //   mu_L    Dynamic mol. visc. (Re=2000)  [Pa*sec] //0.009354
	Pr_L = 0.72; //   Pr_L    Mol. Prandtl number (0.72)    [-]
	Pr_T = 0.8; //   Pr_T    Turb. Prandtl number (0.8-1.) [-]
	Cs = 0.008; //    C       Square of Smagorinsky const Cs[?]      //(C=0.1*0.1)
	g = 1.0; //   g       Accceleration source term     [kg/(m*sec*sec)]
	tau_w = 1.0; //   tau_w   Mean wall stress (Nikuradse)  [Pa]
	R0 = 1.0; //  R0      Initial density               [kg/(m*m*m)]
	f_step = 5000; //    f_step  Frame taking step             [-]
	s_step = 1; //   data sample frequency 

	DriverF = 0.;
	MassRate = 0.;

	blend_adjust_freq = 1; // (int) ( 2*deltaY / deltaT ); 

	//set blending factor
	blend_glob = 0.95;  //control numerics

	// other initializatons
	// complexes with deltas
	deltaT_X = deltaT / deltaX;
	deltaT_Y = deltaT / deltaY;
	deltaT_Z = deltaT / deltaZ;
	_deltaX = 1. / deltaX;
	_deltaY = 1. / deltaY;
	_deltaZ = 1. / deltaZ;
	_4deltaX = 1. / (4. * deltaX);
	_4deltaY = 1. / (4. * deltaY);
	_4deltaZ = 1. / (4. * deltaZ);
	_2deltaX = _4deltaX + _4deltaX;
	_2deltaY = _4deltaY + _4deltaY;
	_2deltaZ = _4deltaZ + _4deltaZ;
	gdeltaT = g * deltaT;
	// box sizes + 1
	LENN = LEN + 1;
	HIGG = HIG + 1;
	DEPP = DEP + 1;
	// constant-value lambda_L
	lambda_L = mu_L * Cv * K / Pr_L;
	// complex cp_Pr_T
	cp_Pr_T = Cv * K / Pr_T;
	// complex Cd * delta * delta in Smagorinsky model
	CsDD = Cs * pow(deltaX * deltaY * deltaZ, 0.6666666666667);
	// friction velocity
	u_tau = sqrt(tau_w / R0);
	// near-wall correction in x-direction
	CsDDx = Array2D(HIG, DEP);
	for (j = 0; j < HIG; j++) {
		for (k = 0; k < DEP; k++) {
			y = (j < HIG / 2) ? (j + 0.5) * deltaY : (HIG - j - 0.5) * deltaY;
			z = (k < DEP / 2) ? (k + 0.5) * deltaZ : (DEP - k - 0.5) * deltaZ;
			dd = 2. * y * z / (y + z + sqrt(y * y + z * z));
			d_plus = R0 * u_tau * dd / mu_L;
			CsDDx[j][k] = CsDD * (1. - exp(-(d_plus / 25.) * (d_plus / 25.) * (d_plus / 25.)));
		} // end for
	} // end for
	// near-wall correction in y-direction
	CsDDy = Array2D(HIGG, DEP);
	for (j = 0; j < HIGG; j++) {
		for (k = 0; k < DEP; k++) {
			y = (j <= HIG / 2) ? j * deltaY : (HIG - j) * deltaY;
			z = (k < DEP / 2) ? (k + 0.5) * deltaZ : (DEP - k - 0.5) * deltaZ;
			dd = 2. * y * z / (y + z + sqrt(y * y + z * z));
			d_plus = R0 * u_tau * dd / mu_L;
			CsDDy[j][k] = CsDD * (1. - exp(-(d_plus / 25.) * (d_plus / 25.) * (d_plus / 25.)));
		}
	}
	// near-wall correction in z-direction
	CsDDz = Array2D(HIG, DEPP);
	for (j = 0; j < HIG; j++) {
		for (k = 0; k < DEPP; k++) {
			y = (j < HIG / 2) ? (j + 0.5) * deltaY : (HIG - j - 0.5) * deltaY;
			z = (k <= DEP / 2) ? k * deltaZ : (DEP - k) * deltaZ;
			dd = 2. * y * z / (y + z + sqrt(y * y + z * z));
			d_plus = R0 * u_tau * dd / mu_L;
			CsDDz[j][k] = CsDD * (1. - exp(-(d_plus / 25.) * (d_plus / 25.) * (d_plus / 25.)));
		}
	}

	// Dynamic memory allocation for 3D arrays
	// total amount of memory required, for gas dynamics *=2*5
	mem = 10 * ((long)LEN + 2) * ((long)HIG + 2) * ((long)DEP + 2)
		+ 10 * ((long)LEN + 1) * (long)HIG * (long)DEP
		+ 10 * (long)LEN * ((long)HIG + 1) * (long)DEP
		+ 10 * (long)LEN * (long)HIG * ((long)DEP + 1);
	printf("%d process: %lu bytes of memory required\n",
		myid + 1, mem * sizeof(real));

	// for quantities in cells
	for (i = 0; i < 10; i++)
		fp[i] = Array3D(LEN + 2, HIG + 2, DEP + 2);
	U1 = fp[0]; U1p = fp[5];
	U2 = fp[1]; U2p = fp[6];
	U3 = fp[2]; U3p = fp[7];
	U4 = fp[3]; U4p = fp[8];
	U5 = fp[4]; U5p = fp[9];
	// for x-fluxes
	for (i = 0; i < 10; i++)
		fp[i] = Array3D(LENN, HIG, DEP);
	xU1 = fp[0]; U1x = fp[5];
	xU2 = fp[1]; U2x = fp[6];
	xU3 = fp[2]; U3x = fp[7];
	xU4 = fp[3]; U4x = fp[8];
	xU5 = fp[4]; U5x = fp[9];
	// for y-fluxes
	for (i = 0; i < 10; i++)
		fp[i] = Array3D(LEN, HIGG, DEP);
	yU1 = fp[0]; U1y = fp[5];
	yU2 = fp[1]; U2y = fp[6];
	yU3 = fp[2]; U3y = fp[7];
	yU4 = fp[3]; U4y = fp[8];
	yU5 = fp[4]; U5y = fp[9];

	// store coefficients of fe-bases at wall-adjacent cells
	wc_coeffs = Array3D(LEN, DEP, 12);
	// PID controller parameters 
	blend_dyn = Array2D(2, LENN);
	nwall_vec = Array2D(HIGG, LENN);
	pid_stat_log = Array3D(2, LENN, 4);
	for (k = 0; k < 2; k++)
		for (i = 1; i < LENN; i++)
		{
			pid_stat_log[k][i][0] = err_diff = 0.;
			pid_stat_log[k][i][1] = err_intg = 0.;
			pid_stat_log[k][i][2] = err_grad = 0.;
			pid_stat_log[k][i][3] = err_prev = 0.;
		}
	if (AnswerPID == 1) {
		double tmp0, tmp1;
		// load values from file for restart
		sprintf(filename, "blend.%d", myid);
		FILE* fh = fopen(filename, "r");

		for (i = 1; i < LENN; i++) {
			fscanf(fh, "%lf %lf", &tmp0, &tmp1);
			blend_dyn[0][i] = tmp0;
			blend_dyn[1][i] = tmp1;
		}

		for (k = 0; k < 2; k++)
			for (i = 1; i < LENN; i++)
				for (j = 0; j < 4; j++) {
					fscanf(fh, "%lf", &tmp0);
					pid_stat_log[k][i][j] = tmp0;
				}

		fclose(fh);
	}

	// allocate stand-alone memory for Fy
	for (i = 0; i < 5; i++)
		fp[i] = Array3D(LEN, HIGG, DEP);
	Fy1 = fp[0];
	Fy2 = fp[1];
	Fy3 = fp[2];
	Fy4 = fp[3];
	Fy5 = fp[4];

	// for z-fluxes
	for (i = 0; i < 10; i++)
		fp[i] = Array3D(LEN, HIG, DEPP);
	zU1 = fp[0]; U1z = fp[5];
	zU2 = fp[1]; U2z = fp[6];
	zU3 = fp[2]; U3z = fp[7];
	zU4 = fp[3]; U4z = fp[8];
	zU5 = fp[4]; U5z = fp[9];
	//
	printf("%d process: 3D arrays allocated\n", myid + 1);

	// MPI Buf
	BufCountF = 5 * HIG * DEP;    // BufCountF < BufCountU
	BufCountU = 5 * (HIG + 2) * (DEP + 2); // BufCountU
	pBuf = (real*)malloc(BufCountU * sizeof(real));
	nBuf = (real*)malloc(BufCountU * sizeof(real));
	Bufp = (real*)malloc(BufCountU * sizeof(real));
	Bufn = (real*)malloc(BufCountU * sizeof(real));
	if (pBuf == NULL || nBuf == NULL || Bufp == NULL || Bufn == NULL) {
		fprintf(stderr, "Cannot allocate memory");
		fflush(stderr);
		MPI_Abort(globleMpicom, 1);
	}

	// array of sensors
	if ((sensors = (real*)malloc(sizeof(real) * HIG)) == NULL) {
		fprintf(stderr, "Cannot allocate memory");
		fflush(stderr);
		MPI_Abort(globleMpicom, 1);
	}

	// array of eddy viscosity
	muturb = Array3D(LEN + 2, HIG + 2, DEP + 2);
	muturb_rans = Array3D(LEN + 2, HIG + 2, DEP + 2);
	muturb_coef = Array3D(LEN + 2, HIG + 2, DEP + 2);
	tmp = Array3D(LEN + 2, HIG + 2, DEP + 2);
	for (i = 0; i <= LENN; i++) {
		for (j = 0; j <= HIGG; j++) {
			for (k = 0; k <= DEPP; k++) {
				muturb[i][j][k] = 0.;
				muturb_rans[i][j][k] = 0.;
				muturb_coef[i][j][k] = 0.;
				tmp[i][j][k] = 0.;
			}
		}
	}

	//--- Initial conditions of the flow field
	// start new simulation
	if (Answer == 0) {
		// add-on perturbation 
		Vd = 6. * sin(2. * M_PI * (myid + 1) / numprocs);
		Wd = 6. * cos(1.2 + 2. * M_PI * (myid + 1) / numprocs);
		//
		for (i = 1, P = 100000., R = R0; i < LENN; i++) {
			U = 15.426724;
			for (j = 1; j < HIGG; j++) {
				for (k = 1; k < DEPP; k++) {
					// some initial disturbance (1/5, 4/5)
					// HIG, DEP should be devidable by 5
					// LEN - by 4
					if (j > HIG / 5 && j <= 4 * HIG / 5 && k > DEP / 5 && k <= 4 * DEP / 5) {
						if (i <= LEN / 4) { V = Vd; W = Wd; }
						else if (i <= 2 * LEN / 4) { V = -Vd; W = Wd; }
						else if (i <= 3 * LEN / 4) { V = Vd; W = -Wd; }
						else { V = -Vd; W = -Wd; }
					}
					else { V = W = 0.; }
					//
					U1[i][j][k] = R;
					U2[i][j][k] = R * U;
					U3[i][j][k] = R * V;
					U4[i][j][k] = R * W;
					U5[i][j][k] = P / K_1 + 0.5 * R * (U * U + V * V + W * W);
				}
			}
		}
		// zero init step number at the beginning
		step = 0;
	}

	// continue previously saved simulation
	else {
		//
		
		
		if (selectStart % 2 == 0)
                {
                        sprintf(filename, "/public3/home/scg8976/RL_test/smarties/runs/aNewCS2/simulation_000_00000/backup.%d", myid);
                }
                else
                {
                        sprintf(filename, "/public3/home/scg8976/RL_test/smarties/runs/aNewCS2_Change/simulation_000_00002/backup.%d", myid);
                }

		MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR,
			MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

		for (i = 1; i < LENN; i++) {
			for (j = 1; j < HIGG; j++) {
				MPI_File_read(fh, &U1[i][j][1], DEP, MPI_DOUBLE, &status);
				MPI_File_read(fh, &U2[i][j][1], DEP, MPI_DOUBLE, &status);

				//need when adjusting the flow rate
				//provide an initial scaling 
				//for (k=1; k< DEPP; k++) U2[i][j][k] *= 1.368592; 

				MPI_File_read(fh, &U3[i][j][1], DEP, MPI_DOUBLE, &status);
				MPI_File_read(fh, &U4[i][j][1], DEP, MPI_DOUBLE, &status);
				MPI_File_read(fh, &U5[i][j][1], DEP, MPI_DOUBLE, &status);
			}
		}
		MPI_File_read(fh, &buf, 1, MPI_DOUBLE, &status);
		step = (int)buf;
		MPI_File_close(&fh);
		step = 3285000;
		
		printf("/////////////////////%d////////////////////\n", step);
	}
	/*
	  // root process outputs sensors data...
	  if (myid == 0) {
		// open "sensors.txt" file in text mode
		if ((pFsen_txt = fopen("sensors.txt", "a")) == NULL)
		  fprintf(stderr, "Cannot open \"sensors.txt\". Error code %u", errno);
		// open "sensors.bin" file in binary mode
		if ((pFsen_bin = fopen("sensors.bin", "ab")) == NULL)
		  fprintf(stderr, "Cannot open \"sensors.bin\". Error code %u", errno);
	  }
	*/
	//
	U1_ = U1, U2_ = U2, U3_ = U3, U4_ = U4, U5_ = U5;
	fprintf(stderr, "%d process: after initializations\n", myid + 1);

	// sampled statisitics (might need previous data)
	Uavg0 = (real*)malloc((HIGG + 1) * sizeof(real));
	Vavg0 = (real*)malloc((HIGG + 1) * sizeof(real));
	Wavg0 = (real*)malloc((HIGG + 1) * sizeof(real));
	Uavg1 = (real*)malloc((HIGG + 1) * sizeof(real));
	Vavg1 = (real*)malloc((HIGG + 1) * sizeof(real));
	Wavg1 = (real*)malloc((HIGG + 1) * sizeof(real));

	UU0 = (real*)malloc((HIGG + 1) * sizeof(real));
	VV0 = (real*)malloc((HIGG + 1) * sizeof(real));
	WW0 = (real*)malloc((HIGG + 1) * sizeof(real));

	UV0 = (real*)malloc((HIGG + 1) * sizeof(real));
	VW0 = (real*)malloc((HIGG + 1) * sizeof(real));
	UW0 = (real*)malloc((HIGG + 1) * sizeof(real));
	for (j = 0; j <= HIGG; j++) {
		Uavg0[j] = 0.; Vavg0[j] = 0.; Wavg0[j] = 0.;
		UU0[j] = 0.; VV0[j] = 0.; WW0[j] = 0.;
		UV0[j] = 0.; VW0[j] = 0.; UW0[j] = 0.;
	}
	accumul_dT = 0.;
} // end Initialize( )


/***********************
* BOUNCONDINGHOSTCELLS *   BC in the cells out of the domain
***********************/
void BounCondInGhostCells(void)
{
	balance = 0.73;

	// y walls
	for (i = 1; i < LENN; i++) {
		for (k = 1; k < DEPP; k++) {

			// bottom side - NO SLIP
			U1_[i][0][k] = U1_[i][1][k];
			U2_[i][0][k] = U2_[i][1][k]; //balance * U2_[i][1][k] + (1.-balance) * (-U2_[i][1][k]) ;  //0.; // - U2_[i][1][k];
			U3_[i][0][k] = -U3_[i][1][k]; //balance * U3_[i][1][k] + (1.-balance) * (-U3_[i][1][k]) ;  //0.; // - U3_[i][1][k];
			U4_[i][0][k] = U4_[i][1][k]; //balance * U4_[i][1][k] + (1.-balance) * (-U4_[i][1][k]) ;  //0.; // - U4_[i][1][k];
			U5_[i][0][k] = U5_[i][1][k];
			// top side - NO SLIP
			U1_[i][HIGG][k] = U1_[i][HIG][k];
			U2_[i][HIGG][k] = U2_[i][HIG][k];  //balance * U2_[i][HIG][k] + (1.-balance) * (-U2_[i][HIG][k]);  //0.; //- U2_[i][HIG][k];
			U3_[i][HIGG][k] = -U3_[i][HIG][k];  // balance * U3_[i][HIG][k] + (1.-balance) * (-U3_[i][HIG][k]);  //0.; //- U3_[i][HIG][k];
			U4_[i][HIGG][k] = U4_[i][HIG][k];  //balance * U4_[i][HIG][k] + (1.-balance) * (-U4_[i][HIG][k]);  //0.; //- U4_[i][HIG][k];
			U5_[i][HIGG][k] = U5_[i][HIG][k];
		}
	}
	// z
	// Enforce periodicity along z-direction
	for (i = 1; i < LENN; i++) {
		for (j = 0; j <= HIGG; j++) {
			// back side - Periodic 
			U1_[i][j][0] = U1_[i][j][DEP];
			U2_[i][j][0] = U2_[i][j][DEP];
			U3_[i][j][0] = U3_[i][j][DEP];
			U4_[i][j][0] = U4_[i][j][DEP];
			U5_[i][j][0] = U5_[i][j][DEP];
			// front side - Periodic  
			U1_[i][j][DEPP] = U1_[i][j][1];
			U2_[i][j][DEPP] = U2_[i][j][1];
			U3_[i][j][DEPP] = U3_[i][j][1];
			U4_[i][j][DEPP] = U4_[i][j][1];
			U5_[i][j][DEPP] = U5_[i][j][1];
		}
	}
	/*
	// along x corners
	for (i = 1; i < LENN; i++) {
	  // NO SLIP
	  U1_[i][0][0] =   U1_[i][1][1];
	  U2_[i][0][0] = - U2_[i][1][1];
	  U3_[i][0][0] = - U3_[i][1][1];
	  U4_[i][0][0] = - U4_[i][1][1];
	  U1_[i][0][DEPP] =   U1_[i][1][DEP];
	  U2_[i][0][DEPP] = - U2_[i][1][DEP];
	  U3_[i][0][DEPP] = - U3_[i][1][DEP];
	  U4_[i][0][DEPP] = - U4_[i][1][DEP];
	  U1_[i][HIGG][0] =   U1_[i][HIG][1];
	  U2_[i][HIGG][0] = - U2_[i][HIG][1];
	  U3_[i][HIGG][0] = - U3_[i][HIG][1];
	  U4_[i][HIGG][0] = - U4_[i][HIG][1];
	  U1_[i][HIGG][DEPP] =   U1_[i][HIG][DEP];
	  U2_[i][HIGG][DEPP] = - U2_[i][HIG][DEP];
	  U3_[i][HIGG][DEPP] = - U3_[i][HIG][DEP];
	  U4_[i][HIGG][DEPP] = - U4_[i][HIG][DEP];
	}
	*/
} // end BounCondInGhostCells()


/*--- For monotone piecevise parabolic reconstruction ---*/

/* constant coefficients: 4.0, 1/3, 1/6 are recommended */
const real BB = 4., C1 = 0.33333333, C2 = 0.16666667;

/***************
* minmod(x, y) * - macro: limiting function of Chakravarthy, Osher
***************/
/*
 - uses register variable kk defined in Reconstruction()
 - don't use expressions as arguments of this macro!
 - beware of other side effects
*/
#define minmod(x, y) ( (x*y <= 0.) ? 0.: ((kk=(y>0.)?y*BB:-y*BB, x>0.) ? ((x<kk)?x:kk) : ((-x<kk)?x:-kk) ) )

/*******************
*  RECONSTRUCTION  *   Piecewise-parabolic reconstruction
*******************/
void Reconstruction(void)
{
	register real kk;

	int next, previous; // processes next/previous to this
	MPI_Status status;
	MPI_Request reqs[4];
	MPI_Status  stats[4];

	//--- Communication via MPI to complete ghost cells layers ---

	//-- next and previous processes:
	//
	if (myid == numprocs - 1)
		next = 0;
	else
		next = myid + 1;
	//
	if (myid == 0)
		previous = numprocs - 1;
	else
		previous = myid - 1;

	//-- preparing buffer arrays to send
	for (i = 0, j = 0; j <= HIGG; j++) {
		for (k = 0; k <= DEPP; k++, i += 5) {
			//
			Bufp[i] = U1_[1][j][k];
			Bufp[i + 1] = U2_[1][j][k];
			Bufp[i + 2] = U3_[1][j][k];
			Bufp[i + 3] = U4_[1][j][k];
			Bufp[i + 4] = U5_[1][j][k];
			//
			Bufn[i] = U1_[LEN][j][k];
			Bufn[i + 1] = U2_[LEN][j][k];
			Bufn[i + 2] = U3_[LEN][j][k];
			Bufn[i + 3] = U4_[LEN][j][k];
			Bufn[i + 4] = U5_[LEN][j][k];
		}
	}

	// syncronization before communication
	MPI_Barrier(globleMpicom);

	// sending data to the previous process
	//if (myid != 0)
	//MPI_Send(Bufp, BufCountU, MPI_DOUBLE, previous, 100, globleMpicom);
	MPI_Isend(Bufp, BufCountU, MPI_DOUBLE, previous, 100, globleMpicom, &reqs[0]);
	// receiving data from the next process
	//if (myid != numprocs-1)
	//MPI_Recv(nBuf, BufCountU, MPI_DOUBLE, next,     100, globleMpicom, &status);
	MPI_Irecv(nBuf, BufCountU, MPI_DOUBLE, next, 100, globleMpicom, &reqs[1]);
	// sending data next processes
	//if (myid != numprocs-1)
	//MPI_Send(Bufn, BufCountU, MPI_DOUBLE, next,     101, globleMpicom);
	MPI_Isend(Bufn, BufCountU, MPI_DOUBLE, next, 101, globleMpicom, &reqs[2]);
	// receiving data from the previous process
	//if (myid != 0)
	//MPI_Recv(pBuf, BufCountU, MPI_DOUBLE, previous, 101, globleMpicom, &status);
	MPI_Irecv(pBuf, BufCountU, MPI_DOUBLE, previous, 101, globleMpicom, &reqs[3]);

	MPI_Waitall(4, reqs, stats);

	// processing the received arrays
	for (i = 0, j = 0; j <= HIGG; j++) {
		for (k = 0; k <= DEPP; k++, i += 5) {
			//
			U1_[0][j][k] = pBuf[i];
			U2_[0][j][k] = pBuf[i + 1];
			U3_[0][j][k] = pBuf[i + 2];
			U4_[0][j][k] = pBuf[i + 3];
			U5_[0][j][k] = pBuf[i + 4];
			//
			U1_[LENN][j][k] = nBuf[i];
			U2_[LENN][j][k] = nBuf[i + 1];
			U3_[LENN][j][k] = nBuf[i + 2];
			U4_[LENN][j][k] = nBuf[i + 3];
			U5_[LENN][j][k] = nBuf[i + 4];
		}
	}

	//--- End of communication via MPI to complete ghost cells layers ---

	//--- Reconstruction ----
	for (i = 1, _i = 0, i_ = 2; i < LENN; i++, _i++, i_++) {
		for (j = 1, _j = 0, j_ = 2; j < HIGG; j++, _j++, j_++) {
			for (k = 1, _k = 0, k_ = 2; k < DEPP; k++, _k++, k_++) {
				//---  X  ---*/
				// at cell centroid
				u1 = U1_[i][j][k];
				u2 = U2_[i][j][k];
				u3 = U3_[i][j][k];
				u4 = U4_[i][j][k];
				u5 = U5_[i][j][k];
				R = 1. / u1;  // 1 / R
				U = u2 * R;
				V = u3 * R;
				W = u4 * R;
				P = (u5 - 0.5 * (u2 * u2 + u3 * u3 + u4 * u4) * R) * K_1;
				C = sqrt(K * P * R);
				// evaluation of matrices
				ff = U * U + V * V + W * W;   //
				aa = K_1_2 * ff;      // (K-1)*(U*U+V*V+W*W)/2
				bb = _K_1 * U;        // - (K_1) * U
				cc = _K_1 * V;        // - (K_1) * V
				dd = _K_1 * W;        // - (K_1) * W
				ll = 1. / (C + C);      // 1/(2.*C)
				ee = ll / C;          // 1./(2*C*C)
				gg = -ff * ee;     // -(U*U+V*V+W*W)/(2*C*C)
				hh = -(ee + ee); // - 1 / (C*C)
				mm = U * ll;          // U/(2.*C)
				// differencing of the conservative variables
				// forvard
				m1 = U1_[i_][j][k] - u1;
				m2 = U2_[i_][j][k] - u2;
				m3 = U3_[i_][j][k] - u3;
				m4 = U4_[i_][j][k] - u4;
				m5 = U5_[i_][j][k] - u5;
				// backward
				w1 = u1 - U1_[_i][j][k];
				w2 = u2 - U2_[_i][j][k];
				w3 = u3 - U3_[_i][j][k];
				w4 = u4 - U4_[_i][j][k];
				w5 = u5 - U5_[_i][j][k];
				// transformation matrix [S]
				S11 = aa - C * C;
				S41 = aa - (kk = C * U);   // CU
				S51 = aa + kk;
				S12 = S42 = S52 = bb;
				S42 += C;
				S52 -= C;
				// finite differences of characteristic variables dW = S * dU
				// forward
				kk = K_1 * m5;
				k2 = cc * m3;
				k3 = dd * m4;
				M1 = S11 * m1 + S12 * m2 + k2 + k3 + kk;
				M2 = -V * m1 + m3;
				M3 = -W * m1 + m4;
				M4 = S41 * m1 + S42 * m2 + k2 + k3 + kk;
				M5 = S51 * m1 + S52 * m2 + k2 + k3 + kk;
				// backward
				kk = K_1 * w5;
				k2 = cc * w3;
				k3 = dd * w4;
				W1 = S11 * w1 + S12 * w2 + k2 + k3 + kk;
				W2 = -V * w1 + w3;
				W3 = -W * w1 + w4;
				W4 = S41 * w1 + S42 * w2 + k2 + k3 + kk;
				W5 = S51 * w1 + S52 * w2 + k2 + k3 + kk;
				/*
					// modified characteristic differences _W
					_M1 = minmod(M1, W1);
					_M2 = minmod(M2, W2);
					_M3 = minmod(M3, W3);
					_M4 = minmod(M4, W4);
					_M5 = minmod(M5, W5);
					_W1 = minmod(W1, M1);
					_W2 = minmod(W2, M2);
					_W3 = minmod(W3, M3);
					_W4 = minmod(W4, M4);
					_W5 = minmod(W5, M5);
				*/
				_M1 = M1;
				_M2 = M2;
				_M3 = M3;
				_M4 = M4;
				_M5 = M5;

				_W1 = W1;
				_W2 = W2;
				_W3 = W3;
				_W4 = W4;
				_W5 = W5;

				// inverted transformation matrix [S-1]
				/*_S21 = */ u_hh = U * hh; /*_S31 =*/ v_hh = V * hh; /*_S41 =*/ w_hh = W * hh;
				u_ee = U * ee; /*_S34 =*/ v_ee = V * ee; /*_S44 =*/ w_ee = W * ee;
				_S24 = u_ee + ll; _S25 = u_ee - ll;
				_S54 = _S55 = (nn = -0.5 * gg) + _1_2_K_1;
				_S54 += mm; _S55 -= mm;
				// parameters' vector at the rightmost interface
				k1 = _M1 * C1 + _W1 * C2;
				k2 = _M2 * C1 + _W2 * C2;
				k3 = _M3 * C1 + _W3 * C2;
				k4 = _M4 * C1 + _W4 * C2;
				k5 = _M5 * C1 + _W5 * C2;
				kk = k4 + k5;
				xU1[i][_j][_k] = u1 + hh * k1 + ee * kk;
				xU2[i][_j][_k] = u2 + u_hh * k1 + _S24 * k4 + _S25 * k5;
				xU3[i][_j][_k] = u3 + v_hh * k1 + k2 + v_ee * kk;
				xU4[i][_j][_k] = u4 + w_hh * k1 + k3 + w_ee * kk;
				xU5[i][_j][_k] = u5 + gg * k1 + V * k2 + W * k3 + _S54 * k4 + _S55 * k5;
				// ... leftmost interface
				k1 = _M1 * C2 + _W1 * C1;
				k2 = _M2 * C2 + _W2 * C1;
				k3 = _M3 * C2 + _W3 * C1;
				k4 = _M4 * C2 + _W4 * C1;
				k5 = _M5 * C2 + _W5 * C1;
				kk = k4 + k5;
				U1x[_i][_j][_k] = u1 - hh * k1 - ee * kk;
				U2x[_i][_j][_k] = u2 - u_hh * k1 - _S24 * k4 - _S25 * k5;
				U3x[_i][_j][_k] = u3 - v_hh * k1 - k2 - v_ee * kk;
				U4x[_i][_j][_k] = u4 - w_hh * k1 - k3 - w_ee * kk;
				U5x[_i][_j][_k] = u5 - gg * k1 - V * k2 - W * k3 - _S54 * k4 - _S55 * k5;
				//---  Y  ---
				// evaluation of matrices
				kk = bb; bb = cc; cc = dd; dd = kk;
				mm = V * ll;          // V/(2.*C)
				// differencing of the conservative variables
				// forvard
				{
					m1 = U1_[i][j_][k] - u1;
					m2 = U3_[i][j_][k] - u3; // U3[i][j+1][k] ( U->V, V->W, W->U : 3-4-2 )
					m3 = U4_[i][j_][k] - u4; // U4[i][j+1][k]
					m4 = U2_[i][j_][k] - u2; // U2[i][j+1][k]
					m5 = U5_[i][j_][k] - u5;
				}
				// backward
				{
					w1 = u1 - U1_[i][_j][k];
					w2 = u3 - U3_[i][_j][k]; // U3[i][j-1][k] ( U->V, V->W, W->U : 3-4-2 )
					w3 = u4 - U4_[i][_j][k]; // U4[i][j-1][k]
					w4 = u2 - U2_[i][_j][k]; // U2[i][j-1][k]
					w5 = u5 - U5_[i][_j][k];
				}
				// transformation matrix [S]
				S41 = aa - (kk = C * V);    // CV
				S51 = aa + kk;
				S12 = S42 = S52 = bb;
				S42 += C;
				S52 -= C;
				// finite differences of characteristic variables dW = S * dU
				// forward
				kk = K_1 * m5;
				k2 = cc * m3;
				k3 = dd * m4;
				M1 = S11 * m1 + S12 * m2 + k2 + k3 + kk;
				M2 = -W * m1 + m3;  // -W*m1 !
				M3 = -U * m1 + m4;  // -U*m1 !
				M4 = S41 * m1 + S42 * m2 + k2 + k3 + kk;
				M5 = S51 * m1 + S52 * m2 + k2 + k3 + kk;
				// backward
				kk = K_1 * w5;
				k2 = cc * w3;
				k3 = dd * w4;
				W1 = S11 * w1 + S12 * w2 + k2 + k3 + kk;
				W2 = -W * w1 + w3;  // -W*m1 !
				W3 = -U * w1 + w4;  // -U*m1 !
				W4 = S41 * w1 + S42 * w2 + k2 + k3 + kk;
				W5 = S51 * w1 + S52 * w2 + k2 + k3 + kk;
				/*
					// modified characteristic differences _W
					_M1 = minmod(M1, W1);
					_M2 = minmod(M2, W2);
					_M3 = minmod(M3, W3);
					_M4 = minmod(M4, W4);
					_M5 = minmod(M5, W5);
					_W1 = minmod(W1, M1);
					_W2 = minmod(W2, M2);
					_W3 = minmod(W3, M3);
					_W4 = minmod(W4, M4);
					_W5 = minmod(W5, M5);
				*/
				{
					_M1 = M1;
					_M2 = M2;
					_M3 = M3;
					_M4 = M4;
					_M5 = M5;
					_W1 = W1;
					_W2 = W2;
					_W3 = W3;
					_W4 = W4;
					_W5 = W5;
				}
				// inverted transformation matrix [S-1]
				//_S21 =  v_hh; _S31 =  w_hh; _S41 =  u_hh;
				//_S34 =  w_ee; _S44 =  u_ee;
				_S24 = v_ee + ll; _S25 = v_ee - ll;
				_S54 = _S55 = nn + _1_2_K_1;
				_S54 += mm; _S55 -= mm;
				// parameters' vector at the upper interface
				k1 = _M1 * C1 + _W1 * C2;
				k2 = _M2 * C1 + _W2 * C2;
				k3 = _M3 * C1 + _W3 * C2;
				k4 = _M4 * C1 + _W4 * C2;
				k5 = _M5 * C1 + _W5 * C2;
				kk = k4 + k5;
				yU1[_i][j][_k] = u1 + hh * k1 + ee * kk;
				yU3[_i][j][_k] = u3 + v_hh * k1 + _S24 * k4 + _S25 * k5;
				yU4[_i][j][_k] = u4 + w_hh * k1 + k2 + w_ee * kk;
				yU2[_i][j][_k] = u2 + u_hh * k1 + k3 + u_ee * kk;
				yU5[_i][j][_k] = u5 + gg * k1 + W * k2 + U * k3 + _S54 * k4 + _S55 * k5;
				/*	if (j==1) {
					yU2[_i][ j][_k] = 2.*u2;
					yU4[_i][ j][_k] = 2.*u4;
					}
					if (j==HIG) {
					yU2[_i][ j][_k] = 0.;
					yU4[_i][ j][_k] = 0.;
						}
				*/	// ...lower interface
				k1 = _M1 * C2 + _W1 * C1;
				k2 = _M2 * C2 + _W2 * C1;
				k3 = _M3 * C2 + _W3 * C1;
				k4 = _M4 * C2 + _W4 * C1;
				k5 = _M5 * C2 + _W5 * C1;
				kk = k4 + k5;
				U1y[_i][_j][_k] = u1 - hh * k1 - ee * kk;
				U3y[_i][_j][_k] = u3 - v_hh * k1 - _S24 * k4 - _S25 * k5;
				U4y[_i][_j][_k] = u4 - w_hh * k1 - k2 - w_ee * kk;
				U2y[_i][_j][_k] = u2 - u_hh * k1 - k3 - u_ee * kk;
				U5y[_i][_j][_k] = u5 - gg * k1 - W * k2 - U * k3 - _S54 * k4 - _S55 * k5;
				/*	if(j==1) {
					U2y[_i][_j][_k] = 0.;
						U4y[_i][_j][_k] = 0.;
					}
						if(j==HIG) {
					U2y[_i][_j][_k] = 2.*u4;
						U4y[_i][_j][_k] = 2.*u2;
					}
				*/	//--- Z ---
				// evaluation of matrices
				kk = bb; bb = cc; cc = dd; dd = kk;
				mm = W * ll;          // W/(2.*C)
				// differencing of the conservative variables
				// forward
				m1 = U1_[i][j][k_] - u1;
				m2 = U4_[i][j][k_] - u4; // U4[i][j][k+1] ( V->W, W->U, U->V : 4-2-3 )
				m3 = U2_[i][j][k_] - u2; // U2[i][j][k+1]
				m4 = U3_[i][j][k_] - u3; // U3[i][j][k+1]
				m5 = U5_[i][j][k_] - u5;
				// backward
				w1 = u1 - U1_[i][j][_k];
				w2 = u4 - U4_[i][j][_k]; // U4[i][j][k-1] ( V->W, W->U, U->V : 4-2-3 )
				w3 = u2 - U2_[i][j][_k]; // U2[i][j][k-1]
				w4 = u3 - U3_[i][j][_k]; // U3[i][j][k-1]
				w5 = u5 - U5_[i][j][_k];
				// transformation matrix [S]
				S41 = aa - (kk = C * W); // CW
				S51 = aa + kk;
				S12 = S42 = S52 = bb;
				S42 += C;
				S52 -= C;
				// finite differences of characteristic variables dW = S * dU
				// forward
				kk = K_1 * m5;
				k2 = cc * m3;
				k3 = dd * m4;
				M1 = S11 * m1 + S12 * m2 + k2 + k3 + kk;
				M2 = -U * m1 + m3;  // -U*m1 !
				M3 = -V * m1 + m4;  // -V*m1 !
				M4 = S41 * m1 + S42 * m2 + k2 + k3 + kk;
				M5 = S51 * m1 + S52 * m2 + k2 + k3 + kk;
				// backward
				kk = K_1 * w5;
				k2 = cc * w3;
				k3 = dd * w4;
				W1 = S11 * w1 + S12 * w2 + k2 + k3 + kk;
				W2 = -U * w1 + w3;  // -U*m1 !
				W3 = -V * w1 + w4;  // -V*m1 !
				W4 = S41 * w1 + S42 * w2 + k2 + k3 + kk;
				W5 = S51 * w1 + S52 * w2 + k2 + k3 + kk;
				/*
					// modified characteristic differences _W
					_M1 = minmod(M1, W1);
					_M2 = minmod(M2, W2);
					_M3 = minmod(M3, W3);
					_M4 = minmod(M4, W4);
					_M5 = minmod(M5, W5);
					_W1 = minmod(W1, M1);
					_W2 = minmod(W2, M2);
					_W3 = minmod(W3, M3);
					_W4 = minmod(W4, M4);
					_W5 = minmod(W5, M5);
				*/
				_M1 = M1;
				_M2 = M2;
				_M3 = M3;
				_M4 = M4;
				_M5 = M5;
				_W1 = W1;
				_W2 = W2;
				_W3 = W3;
				_W4 = W4;
				_W5 = W5;

				// inverted transformation matrix [S-1]
				//_S21 =  w_hh; _S31 =  u_hh; _S41 =  v_hh;
				//_S34 =  u_ee; _S44 =  w_ee;
				_S24 = w_ee + ll; _S25 = w_ee - ll;
				_S54 = _S55 = nn + _1_2_K_1;
				_S54 += mm; _S55 -= mm;
				// parameters' vector at the nearer interface
				k1 = _M1 * C1 + _W1 * C2;
				k2 = _M2 * C1 + _W2 * C2;
				k3 = _M3 * C1 + _W3 * C2;
				k4 = _M4 * C1 + _W4 * C2;
				k5 = _M5 * C1 + _W5 * C2;
				kk = k4 + k5;
				zU1[_i][_j][k] = u1 + hh * k1 + ee * kk;
				zU4[_i][_j][k] = u4 + w_hh * k1 + _S24 * k4 + _S25 * k5;
				zU2[_i][_j][k] = u2 + u_hh * k1 + k2 + u_ee * kk;
				zU3[_i][_j][k] = u3 + v_hh * k1 + k3 + w_ee * kk;
				zU5[_i][_j][k] = u5 + gg * k1 + U * k2 + V * k3 + _S54 * k4 + _S55 * k5;
				// ... further interface
				k1 = _M1 * C2 + _W1 * C1;
				k2 = _M2 * C2 + _W2 * C1;
				k3 = _M3 * C2 + _W3 * C1;
				k4 = _M4 * C2 + _W4 * C1;
				k5 = _M5 * C2 + _W5 * C1;
				kk = k4 + k5;
				U1z[_i][_j][_k] = u1 - hh * k1 - ee * kk;
				U4z[_i][_j][_k] = u4 - w_hh * k1 - _S24 * k4 - _S25 * k5;
				U2z[_i][_j][_k] = u2 - u_hh * k1 - k2 - u_ee * kk;
				U3z[_i][_j][_k] = u3 - v_hh * k1 - k3 - w_ee * kk;
				U5z[_i][_j][_k] = u5 - gg * k1 - U * k2 - V * k3 - _S54 * k4 - _S55 * k5;
			}
		}
	}

} // end Reconstruction( )


/***********************
* BOUNCONDONINTERFACES *  Completing the "outer" values at the domain boundary
***********************/
void BounCondOnInterfaces(void)
{

	// y
	for (i = 0; i < LEN; i++) {
		for (k = 0; k < DEP; k++) {
			// bottom side - NO SLIP
			yU1[i][0][k] = U1y[i][0][k];
			yU2[i][0][k] = U2y[i][0][k];  //- U2y[i][0][k];
			yU3[i][0][k] = -U3y[i][0][k];
			yU4[i][0][k] = U4y[i][0][k];  //- U4y[i][0][k];
			yU5[i][0][k] = U5y[i][0][k];  // Iso-thermal wall
			// top side - NO SLIP
			U1y[i][HIG][k] = yU1[i][HIG][k];
			U2y[i][HIG][k] = yU2[i][HIG][k]; //- yU2[i][HIG][k];
			U3y[i][HIG][k] = -yU3[i][HIG][k];
			U4y[i][HIG][k] = yU4[i][HIG][k]; //- yU4[i][HIG][k];
			U5y[i][HIG][k] = yU5[i][HIG][k];  // Iso-thermal wall
		} //
	}
	// z
	for (i = 0; i < LEN; i++) {
		for (j = 0; j < HIG; j++) {
			// back side - Periodic  
			zU1[i][j][0] = zU1[i][j][DEP];
			zU2[i][j][0] = zU2[i][j][DEP];
			zU3[i][j][0] = zU3[i][j][DEP];
			zU4[i][j][0] = zU4[i][j][DEP];
			zU5[i][j][0] = zU5[i][j][DEP];
			// front side - Periodic  
			U1z[i][j][DEP] = U1z[i][j][0];
			U2z[i][j][DEP] = U2z[i][j][0];
			U3z[i][j][DEP] = U3z[i][j][0];
			U4z[i][j][DEP] = U4z[i][j][0];
			U5z[i][j][DEP] = U5z[i][j][0];
		}
	}
} // end BounCondOnInterfaces( )


/***********
*  FLUXES  *   Gasdynamical+molecular fluxes at the interfaces
***********/
void Fluxes(void)
{
	register real RU; // convective normal mass flux
	register real yunit;
	real Uy_left, Uy_right;
	real blendloc, dtmp;
	real Uwm, Wwm, Uwm_tot, tau_wm;

	int iter, next, previous; // processes next/previous to this
	MPI_Status status;

	MPI_Request reqs[4];
	MPI_Status stats[4];

	//--- Communication via MPI to complete boundary parameters ---

	//-- next and previous processes:
	//
	if (myid == numprocs - 1) next = 0;
	else                      next = myid + 1;
	//
	if (myid == 0)            previous = numprocs - 1;
	else                      previous = myid - 1;

	// preparing buffer arrays to send
	for (i = 0, j = 0; j < HIG; j++) {
		for (k = 0; k < DEP; k++, i += 5) {
			//
			Bufp[i] = U1x[0][j][k];
			Bufp[i + 1] = U2x[0][j][k];
			Bufp[i + 2] = U3x[0][j][k];
			Bufp[i + 3] = U4x[0][j][k];
			Bufp[i + 4] = U5x[0][j][k];
			//
			Bufn[i] = xU1[LEN][j][k];
			Bufn[i + 1] = xU2[LEN][j][k];
			Bufn[i + 2] = xU3[LEN][j][k];
			Bufn[i + 3] = xU4[LEN][j][k];
			Bufn[i + 4] = xU5[LEN][j][k];
		}
	}

	// syncronization before communication
	MPI_Barrier(globleMpicom);

	// sending data to the previous process
	//if (myid != 0)
	//MPI_Send(Bufp, BufCountF, MPI_DOUBLE, previous, 102, globleMpicom);
	MPI_Isend(Bufp, BufCountF, MPI_DOUBLE, previous, 102, globleMpicom, &reqs[0]);
	// receiving data from the next process
	//if (myid != numprocs-1)
	//MPI_Recv(nBuf, BufCountF, MPI_DOUBLE, next,     102, globleMpicom, &status);
	MPI_Irecv(nBuf, BufCountF, MPI_DOUBLE, next, 102, globleMpicom, &reqs[1]);
	// sending data to the next process  
	//if (myid != numprocs-1)
	//MPI_Send(Bufn, BufCountF, MPI_DOUBLE, next,     103, globleMpicom);
	MPI_Isend(Bufn, BufCountF, MPI_DOUBLE, next, 103, globleMpicom, &reqs[2]);
	// receiving data from the previous process
	//if (myid != 0)
	//MPI_Recv(pBuf, BufCountF, MPI_DOUBLE, previous, 103, globleMpicom, &status);
	MPI_Irecv(pBuf, BufCountF, MPI_DOUBLE, previous, 103, globleMpicom, &reqs[3]);

	MPI_Waitall(4, reqs, stats);

	// processing the received arrays
	for (i = 0, j = 0; j < HIG; j++) {
		for (k = 0; k < DEP; k++, i += 5) {
			//
			xU1[0][j][k] = pBuf[i];
			xU2[0][j][k] = pBuf[i + 1];
			xU3[0][j][k] = pBuf[i + 2];
			xU4[0][j][k] = pBuf[i + 3];
			xU5[0][j][k] = pBuf[i + 4];
		}
	}

	for (i = 0, j = 0; j < HIG; j++) {
		for (k = 0; k < DEP; k++, i += 5) {
			//
			U1x[LEN][j][k] = nBuf[i];
			U2x[LEN][j][k] = nBuf[i + 1];
			U3x[LEN][j][k] = nBuf[i + 2];
			U4x[LEN][j][k] = nBuf[i + 3];
			U5x[LEN][j][k] = nBuf[i + 4];
		}
	}

	//--- End of communication via MPI to complete ghost cells layers ---

	//--- Clean up the driving force holder  
	DriverF = 0.;
	LinearU = 0.;

	//--- Evaluation of fluxes ---
	// X-fluxes
	for (i = 0, i_ = 1; i < LENN; i++, i_++) {
		for (j = 0, j_ = 1, j__ = 2; j < HIG; j++, j_++, j__++) {
			for (k = 0, k_ = 1, k__ = 2; k < DEP; k++, k_++, k__++) {
				//--- Inviscid fluxes
				// left parameters
				_U = xU2[i][j][k] / (_R = xU1[i][j][k]);
				_V = xU3[i][j][k] / _R;
				_W = xU4[i][j][k] / _R;
				_P = (xU5[i][j][k] - 0.5 * _R * (_U * _U + _V * _V + _W * _W)) * K_1;
				_C = sqrt(K * _P / _R);
				_jo = _K_1 * _P;           //_P - _C * _C * _R
				_jp = _U + (RU = _P / (_R * _C));
				_jm = _U - RU;
				// right parameters
				U_ = U2x[i][j][k] / (R_ = U1x[i][j][k]);
				V_ = U3x[i][j][k] / R_;
				W_ = U4x[i][j][k] / R_;
				P_ = (U5x[i][j][k] - 0.5 * R_ * (U_ * U_ + V_ * V_ + W_ * W_)) * K_1;
				C_ = sqrt(K * P_ / R_);
				jo_ = _K_1 * P_;           // P_ - C_ * C_ * R_
				jp_ = U_ + (RU = P_ / (R_ * C_));
				jm_ = U_ - RU;
				// linearized procedure
				U = 0.5 * (_U + U_);
				C = 0.5 * (_C + C_);
				C_p = U + C; C_m = U - C; C_o = U;
				// J_p and al_p
				if (C_p > 0.) { al_p = K_1_2 * (_jm - _jp) / _jo; J_p = _jp; }
				else { al_p = K_1_2 * (jm_ - jp_) / jo_; J_p = jp_; }
				// J_m and al_m
				if (C_m > 0.) { al_m = K_1_2 * (_jp - _jm) / _jo; J_m = _jm; }
				else { al_m = K_1_2 * (jp_ - jm_) / jo_; J_m = jm_; }
				// J_o and al_o
				if (C_o > 0.) { al_o = _05K * (_jp - _jm); al_o = -al_o * al_o; J_o = _jo; }
				else { al_o = _05K * (jp_ - jm_); al_o = -al_o * al_o; J_o = jo_; }
				// parameters at the interface
				P = (J_p - J_m) / (al_p - al_m);
				U = J_p - al_p * P;
				R = (J_o - P) / al_o;
				// central flux
				ctF1 = 0.5 * (_R * _U + R_ * U_);
				ctF2 = 0.5 * (_R * _U * _U + _P + R_ * U_ * U_ + P_);
				ctF3 = 0.5 * (_R * _U * _V + R_ * U_ * V_);
				ctF4 = 0.5 * (_R * _U * _W + R_ * U_ * W_);
				ctF5 = 0.5 * (xU5[i][j][k] + _P) * _U + 0.5 * (U5x[i][j][k] + P_) * U_;
				// tangential velocities
				if (U > 0.) { V = _V; W = _W; }
				else { V = V_; W = W_; }
				//--- Viscous fluxes
				//-- forward cells
				// forward
				U_ = U2_[i_][j_][k_] / (R_ = U1_[i_][j_][k_]);
				V_ = U3_[i_][j_][k_] / R_;
				W_ = U4_[i_][j_][k_] / R_;
				P_ = (U5_[i_][j_][k_] - 0.5 * R_ * (U_ * U_ + V_ * V_ + W_ * W_)) * K_1;
				E_ = U5_[i_][j_][k_];
				T_ = P_ / (R_GAS * R_);
				// forward upper
				Uj = U2_[i_][j__][k_] * (RU = 1. / U1_[i_][j__][k_]);
				Vj = U3_[i_][j__][k_] * RU;
				Wj = U4_[i_][j__][k_] * RU;
				// forward lower
				jU = U2_[i_][j][k_] * (RU = 1. / U1_[i_][j][k_]);
				jV = U3_[i_][j][k_] * RU;
				jW = U4_[i_][j][k_] * RU;
				// forward nearer
				Uk = U2_[i_][j_][k__] * (RU = 1. / U1_[i_][j_][k__]);
				Vk = U3_[i_][j_][k__] * RU;
				Wk = U4_[i_][j_][k__] * RU;
				// forward further
				kU = U2_[i_][j_][k] * (RU = 1. / U1_[i_][j_][k]);
				kV = U3_[i_][j_][k] * RU;
				kW = U4_[i_][j_][k] * RU;
				//-- backward cells
				// backward
				_U = U2_[i][j_][k_] / (_R = U1_[i][j_][k_]);
				_V = U3_[i][j_][k_] / _R;
				_W = U4_[i][j_][k_] / _R;
				_P = (U5_[i][j_][k_] - 0.5 * _R * (_U * _U + _V * _V + _W * _W)) * K_1;
				_E = U5_[i][j_][k_];
				_T = _P / (R_GAS * _R);
				// backward upper
				uj = U2_[i][j__][k_] * (RU = 1. / U1_[i][j__][k_]);
				vj = U3_[i][j__][k_] * RU;
				wj = U4_[i][j__][k_] * RU;
				// backward lower
				ju = U2_[i][j][k_] * (RU = 1. / U1_[i][j][k_]);
				jv = U3_[i][j][k_] * RU;
				jw = U4_[i][j][k_] * RU;
				// backward nearer
				uk = U2_[i][j_][k__] * (RU = 1. / U1_[i][j_][k__]);
				vk = U3_[i][j_][k__] * RU;
				wk = U4_[i][j_][k__] * RU;
				// backward further
				ku = U2_[i][j_][k] * (RU = 1. / U1_[i][j_][k]);
				kv = U3_[i][j_][k] * RU;
				kw = U4_[i][j_][k] * RU;
				// derivatives of velocities
				// x
				du_dx = (U_ - _U) * _deltaX;
				dv_dx = (V_ - _V) * _deltaX;
				dw_dx = (W_ - _W) * _deltaX;
				// y
				du_dy = (Uj + uj - jU - ju) * _4deltaY;
				dv_dy = (Vj + vj - jV - jv) * _4deltaY;
				dw_dy = (Wj + wj - jW - jw) * _4deltaY;
				// z
				du_dz = (Uk + uk - kU - ku) * _4deltaZ;
				dv_dz = (Vk + vk - kV - kv) * _4deltaZ;
				dw_dz = (Wk + wk - kW - kw) * _4deltaZ;
				// mean velocities and density
				rr = 0.5 * (R_ + _R);
				uu = 0.5 * (U_ + _U);
				vv = 0.5 * (V_ + _V);
				ww = 0.5 * (W_ + _W);
				// stresses
			/*
				S12 = du_dy + dv_dx;
				S13 = du_dz + dw_dx;
				S23 = dv_dz + dw_dy;
				_S_ = sqrt((du_dx * du_dx + dv_dy * dv_dy + dw_dz * dw_dz) * 2.
					   + (S12 * S12   +   S13 * S13   +   S23 * S23));
				(mu_T = rr * CsDD * _S_);
			*/
			/*
				gu[0] = du_dx; gu[1] = dv_dx; gu[2] = dw_dx;
				gu[3] = du_dy; gu[4] = dv_dy; gu[5] = dw_dy;
				gu[6] = du_dz; gu[7] = dv_dz; gu[8] = dw_dz;
					mu_T  = VremenModel(gu, deltaX*deltaY*deltaZ);
					mu_T *= rr;
			*/
				mu_T = 0.5 * (muturb[i][j_][k_] + muturb[i_][j_][k_]);
				mu_E = mu_L + mu_T;
				sigma_xx = mu_E * (du_dx + du_dx - dv_dy - dw_dz) * 0.66666667;
				sigma_xy = mu_E * (du_dy + dv_dx);
				sigma_xz = mu_E * (du_dz + dw_dx);
				// heat flux
				lambda_E = lambda_L + mu_T * cp_Pr_T;
				q_x = -lambda_E * (T_ - _T) * _deltaX;
				// central flux
			//	ctF1 = 0.5 * (_R*_V       + R_*V_); 
			//	ctF2 = 0.5 * (_R*_V*_U    + R_*V_*U_); 
			//      ctF3 = 0.5 * (_R*_V*_V+_P + R_*V_*V_+P_);
			//	ctF4 = 0.5 * (_R*_V*_W    + R_*V_*W_); 
			//	ctF5 = 0.5 * (_E+_P)*_V + 0.5 * (E_+P_)*V_; 
				// summary X-fluxes evaluation
				//blendloc = (j > 3 && j < HIG-3) ? 1. : blend_dyn; //blend; 
					//if(j<=1 || j>=HIG-1) blendloc = 0.997; 	
				blendloc = blend_glob;
				//if (j > 4 && j < HIG-4) blendloc = 1.; 
				Fx1[i][j][k] = (1. - blendloc) * (RU = R * U) + blendloc * ctF1;
				Fx2[i][j][k] = (1. - blendloc) * (RU * U + P) + blendloc * ctF2 - sigma_xx;
				Fx3[i][j][k] = (1. - blendloc) * RU * V + blendloc * ctF3 - sigma_xy;
				Fx4[i][j][k] = (1. - blendloc) * RU * W + blendloc * ctF4 - sigma_xz;
				Fx5[i][j][k] = (1. - blendloc) * U * (K_K_1 * P + 0.5 * R * (U * U + V * V + W * W))
					+ blendloc * ctF5
					- uu * sigma_xx - vv * sigma_xy - ww * sigma_xz + q_x;
			}
		}
	}

	//--- Y-fluxes
	for (i = 0, i_ = 1, i__ = 2; i < LEN; i++, i_++, i__++) {
		for (j = 0, j_ = 1; j < HIGG; j++, j_++) {
			for (k = 0, k_ = 1, k__ = 2; k < DEP; k++, k_++, k__++) {
				//--- Inviscid fluxes
				// lower
				_U = yU2[i][j][k] / (_R = yU1[i][j][k]);
				_V = yU3[i][j][k] / _R; Uy_left = _V;
				_W = yU4[i][j][k] / _R;
				_P = (yU5[i][j][k] - 0.5 * _R * (_U * _U + _V * _V + _W * _W)) * K_1;
				_C = sqrt(K * _P / _R);
				_jo = _K_1 * _P; 	    //_P - _C * _C * _R
				_jp = _V + (RU = _P / (_R * _C));
				_jm = _V - RU;
				/* upper */
				U_ = U2y[i][j][k] / (R_ = U1y[i][j][k]);
				V_ = U3y[i][j][k] / R_; Uy_right = V_;
				W_ = U4y[i][j][k] / R_;
				P_ = (U5y[i][j][k] - 0.5 * R_ * (U_ * U_ + V_ * V_ + W_ * W_)) * K_1;
				C_ = sqrt(K * P_ / R_);
				jo_ = _K_1 * P_;             // P_ - C_ * C_ * R_
				jp_ = V_ + (RU = P_ / (R_ * C_));
				jm_ = V_ - RU;
				/* linearized procedure */
				V = 0.5 * (_V + V_);
				C = 0.5 * (_C + C_);
				C_p = V + C; C_m = V - C; C_o = V;
				// Jp and al_p
				if (C_p > 0.) { al_p = K_1_2 * (_jm - _jp) / _jo; J_p = _jp; }
				else { al_p = K_1_2 * (jm_ - jp_) / jo_; J_p = jp_; }
				// Jm and al_m
				if (C_m > 0.) { al_m = K_1_2 * (_jp - _jm) / _jo; J_m = _jm; }
				else { al_m = K_1_2 * (jp_ - jm_) / jo_; J_m = jm_; }
				// Jo and al_o
				if (C_o > 0.) { al_o = _05K * (_jp - _jm); al_o = -al_o * al_o; J_o = _jo; }
				else { al_o = _05K * (jp_ - jm_); al_o = -al_o * al_o; J_o = jo_; }
				// parameters at the interface
				P = (J_p - J_m) / (al_p - al_m);
				V = J_p - al_p * P;
				R = (J_o - P) / al_o;
				// central flux
				ctF1 = 0.5 * (_R * _V + R_ * V_);
				ctF2 = 0.5 * (_R * _V * _U + R_ * V_ * U_);
				ctF3 = 0.5 * (_R * _V * _V + _P + R_ * V_ * V_ + P_);
				ctF4 = 0.5 * (_R * _V * _W + R_ * V_ * W_);
				ctF5 = 0.5 * (yU5[i][j][k] + _P) * _V + 0.5 * (U5y[i][j][k] + P_) * V_;
				// tangential velocities
				if (V > 0.) { U = _U; W = _W; }
				else { U = U_; W = W_; }
				//--- Viscous fluxes
				//-- upper cells
				// upper
				U_ = U2_[i_][j_][k_] / (R_ = U1_[i_][j_][k_]);
				V_ = U3_[i_][j_][k_] / R_;
				W_ = U4_[i_][j_][k_] / R_;
				P_ = (U5_[i_][j_][k_] - 0.5 * R_ * (U_ * U_ + V_ * V_ + W_ * W_)) * K_1;
				E_ = U5_[i_][j_][k_];
				T_ = P_ / (R_GAS * R_);
				// upper nearer
				Uk = U2_[i_][j_][k__] * (RU = 1. / U1_[i_][j_][k__]);
				Vk = U3_[i_][j_][k__] * RU;
				Wk = U4_[i_][j_][k__] * RU;
				// upper further
				kU = U2_[i_][j_][k] * (RU = 1. / U1_[i_][j_][k]);
				kV = U3_[i_][j_][k] * RU;
				kW = U4_[i_][j_][k] * RU;
				// upper forward
				Ui = U2_[i__][j_][k_] * (RU = 1. / U1_[i__][j_][k_]);
				Vi = U3_[i__][j_][k_] * RU;
				Wi = U4_[i__][j_][k_] * RU;
				// upper backward
				iU = U2_[i][j_][k_] * (RU = 1. / U1_[i][j_][k_]);
				iV = U3_[i][j_][k_] * RU;
				iW = U4_[i][j_][k_] * RU;
				//-- lower cells
				// lower
				_U = U2_[i_][j][k_] / (_R = U1_[i_][j][k_]);
				_V = U3_[i_][j][k_] / _R;
				_W = U4_[i_][j][k_] / _R;
				_P = (U5_[i_][j][k_] - 0.5 * _R * (_U * _U + _V * _V + _W * _W)) * K_1;
				_E = U5_[i_][j][k_];
				_T = _P / (R_GAS * _R);
				// lower nearer
				uk = U2_[i_][j][k__] * (RU = 1. / U1_[i_][j][k__]);
				vk = U3_[i_][j][k__] * RU;
				wk = U4_[i_][j][k__] * RU;
				// lower further
				ku = U2_[i_][j][k] * (RU = 1. / U1_[i_][j][k]);
				kv = U3_[i_][j][k] * RU;
				kw = U4_[i_][j][k] * RU;
				// lower forward
				ui = U2_[i__][j][k_] * (RU = 1. / U1_[i__][j][k_]);
				vi = U3_[i__][j][k_] * RU;
				wi = U4_[i__][j][k_] * RU;
				// lower backward
				iu = U2_[i][j][k_] * (RU = 1. / U1_[i][j][k_]);
				iv = U3_[i][j][k_] * RU;
				iw = U4_[i][j][k_] * RU;
				// derivatives of velocities
				// x
				du_dx = (Ui + ui - iU - iu) * _4deltaX;
				dv_dx = (Vi + vi - iV - iv) * _4deltaX;
				dw_dx = (Wi + wi - iW - iw) * _4deltaX;
				// y
				du_dy = (U_ - _U) * _deltaY;
				dv_dy = (V_ - _V) * _deltaY;
				dw_dy = (W_ - _W) * _deltaY;
				// z
				du_dz = (Uk + uk - kU - ku) * _4deltaZ;
				dv_dz = (Vk + vk - kV - kv) * _4deltaZ;
				dw_dz = (Wk + wk - kW - kw) * _4deltaZ;
				// mean velocities and density
				rr = 0.5 * (R_ + _R);
				uu = 0.5 * (U_ + _U);
				vv = 0.5 * (V_ + _V);
				ww = 0.5 * (W_ + _W);
				/*
					// stresses
					S12 = du_dy + dv_dx;
					S13 = du_dz + dw_dx;
					S23 = dv_dz + dw_dy;
					_S_ = sqrt((du_dx * du_dx + dv_dy * dv_dy + dw_dz * dw_dz) * 2.
						   + (S12 * S12   +   S13 * S13   +   S23 * S23));
					// mu_T = r * Cs * delta * delta * | S |
					(mu_T = rr * CsDD * _S_);
				*/
				/*
					gu[0] = du_dx; gu[1] = dv_dx; gu[2] = dw_dx;
					gu[3] = du_dy; gu[4] = dv_dy; gu[5] = dw_dy;
					gu[6] = du_dz; gu[7] = dv_dz; gu[8] = dw_dz;
						mu_T  = VremenModel(gu, deltaX*deltaY*deltaZ);
						mu_T *= rr;
				*/
				mu_T = 0.5 * (muturb[i_][j][k_] + muturb[i_][j_][k_]);
				dtmp = 0.5 * (muturb_rans[i_][j][k_] + muturb_rans[i_][j_][k_]);
				double coef = 0.5 * (muturb_coef[i_][j][k_] + muturb_coef[i_][j_][k_]);
#ifdef STRONG_ENFORCE_BC
				if (j == 0 || j == HIG)
				{
					mu_T = 0.;
					dtmp = 0.;
				}
#endif 
				mu_E = mu_L + mu_T;
				//if (j<=1 || j>=(HIG-1)) 
			//sigma_yx = (mu_T + mu_L) * 0.1 * (dv_dx + du_dy ); 
			//else if (j==2 || j==(HIG-2)) 
			//sigma_yx = (mu_T + mu_L) * (dv_dx + du_dy ); 
			//else 
			//sigma_yx = mu_E * (dv_dx) + (mu_L + mu_T*coef) * du_dy + dtmp*du_dy;
				sigma_yx = mu_E * (dv_dx + du_dy) + dtmp * du_dy;
				sigma_yy = mu_E * (dv_dy + dv_dy - du_dx - dw_dz) * 0.66666667;
				//if (j<=1 || j>=(HIG-1)) 
			//sigma_yz = mu_E * 0.1 * (dv_dz + dw_dy );
			//else 
				sigma_yz = mu_E * (dv_dz + dw_dy);
				// heat flux
				lambda_E = lambda_L + mu_T * cp_Pr_T;
				q_y = -lambda_E * (T_ - _T) * _deltaY;
				// summary Y-fluxes evaluation
#ifdef STRONG_ENFORCE_BC
				if (j == 0 || j == HIG)
				{
					U = 0.; W = 0.;
					Fy1[i][j][k] = RU = R * V;
					Fy2[i][j][k] = RU * U; //   - sigma_yx;  //this term accounted for below 
					Fy3[i][j][k] = RU * V + P - sigma_yy;
					Fy4[i][j][k] = RU * W; //   - sigma_yz;  //this term accounted for below 
					Fy5[i][j][k] = V * (K_K_1 * P + 0.5 * R * (U * U + V * V + W * W))
						- uu * sigma_yx - vv * sigma_yy - ww * sigma_yz + q_y;

					if (j == 0)
					{
						Uwm = U2_[i_][j_ + 2][k_];
						Wwm = U4_[i_][j_ + 2][k_];

						Uwm_tot = sqrt(Uwm * Uwm + Wwm * Wwm);
						wall_model(Uwm_tot, &tau_wm, 2.5 * deltaY, mu_L);
						//amplifer
						//tau_wm *= 1.2;
						yunit = sqrt(tau_wm) * 2.5 * deltaY / mu_L;
						if (yunit > 11.1)
						{
							Fy2[i][j][k] = RU * U - tau_wm * Uwm / Uwm_tot;
							Fy4[i][j][k] = RU * W - tau_wm * Wwm / Uwm_tot;
						}
						else
						{
							Fy2[i][j][k] = RU * U - U_ * 2. / deltaY * mu_L;
							Fy4[i][j][k] = RU * W - W_ * 2. / deltaY * mu_L;
						}
						// invoke the wall model YL //yU2[i][j+1][k]
					//        Tau_Wall_Eval(yunit, deltaY, 1.5*deltaY, U2_[i_][j_+1][k_], U_, &alow_x, &aprime_x, &res_x);
					//        alow_x *= 0.41; 
					//
					//        // yU4[i][j+1][k]	
					//        Tau_Wall_Eval(yunit, deltaY, 1.5*deltaY, U4_[i_][j_+1][k_], W_, &alow_z, &aprime_z, &res_z);
					//	alow_z *= 0.41;
					//
					//        yunit = sqrt(alow_x*alow_x + alow_z*alow_z) * deltaY / mu_L; 
					//        if (yunit > 11.1) 	
					//	{
					//        Fy2[i][j][k] -= alow_x * sqrt(alow_x*alow_x + alow_z*alow_z) + aprime_x*mu_L ;  
					//       	Fy4[i][j][k] -= alow_z * sqrt(alow_x*alow_x + alow_z*alow_z) + aprime_z*mu_L; 
					//        }
					//	else
					//	{
					//        Fy2[i][j][k] -= U_ * 2. / deltaY *mu_L; alow_x = 0.; aprime_x = U_*2./deltaY; res_x = 0.;  
					//       	Fy4[i][j][k] -= W_ * 2. / deltaY *mu_L; alow_z = 0.; aprime_z = W_*2./deltaY; res_z = 0.; 
					//	}
					//        
						DriverF += -Fy2[i][j][k];
						//        LinearU +=  0.5*(yU2[i][j+1][k] + U2y[i][j+1][k]); 	

					}

					if (j == HIG)
					{
						Uwm = U2_[i_][j - 2][k_];
						Wwm = U4_[i_][j - 2][k_];

						Uwm_tot = sqrt(Uwm * Uwm + Wwm * Wwm);
						wall_model(Uwm_tot, &tau_wm, 2.5 * deltaY, mu_L);
						//amplifer 
						//tau_wm *= 1.2;
						yunit = sqrt(tau_wm) * deltaY / mu_L;
						if (yunit > 11.1)
						{
							Fy2[i][j][k] = RU * U + tau_wm * Uwm / Uwm_tot;
							Fy4[i][j][k] = RU * W + tau_wm * Wwm / Uwm_tot;
						}
						else
						{
							Fy2[i][j][k] = RU * U + _U * 2. / deltaY * mu_L;
							Fy4[i][j][k] = RU * W + _W * 2. / deltaY * mu_L;
						}

						//        Tau_Wall_Eval(yunit, deltaY, 1.5*deltaY, U2_[i_][j-1][k_], _U, &alow_x, &aprime_x, &res_x); 
						//        alow_x *= 0.41; 
						//
						//        // U4y[i][j-1][k]	
						//        Tau_Wall_Eval(yunit, deltaY, 1.5*deltaY, U4_[i_][j-1][k_], _W, &alow_z, &aprime_z, &res_z); 
						//        alow_z *= 0.41; 
						//
						//        yunit = sqrt(alow_x*alow_x + alow_z*alow_z) * deltaY / mu_L; 
						//        if (yunit > 11.1) 	
						//	{
						//	Fy2[i][j][k] += alow_x * sqrt(alow_x*alow_x + alow_z*alow_z) + aprime_x*mu_L ; 
						//       	Fy4[i][j][k] += alow_z * sqrt(alow_x*alow_x + alow_z*alow_z) + aprime_z*mu_L; 
						//        } 
						//        else
						//	{
						//	Fy2[i][j][k] += _U * 2. / deltaY *mu_L; alow_x = 0.; aprime_x = _U*2./deltaY; res_x = 0.;  
						//       	Fy4[i][j][k] += _W * 2. / deltaY *mu_L; alow_z = 0.; aprime_z = _W*2./deltaY; res_z = 0.;
						//	}
						//
						DriverF += Fy2[i][j][k];
						//        LinearU += 0.5 * (U2y[i][j-1][k] + yU2[i][j-1][k]) ; 	

					}
				}
				else
#endif
				{
					// central flux
				//	ctF1 = 0.5 * (_R*_V       + R_*V_); 
				//	ctF2 = 0.5 * (_R*_V*_U    + R_*V_*U_); 
				//      ctF3 = 0.5 * (_R*_V*_V+_P + R_*V_*V_+P_);
				//	ctF4 = 0.5 * (_R*_V*_W    + R_*V_*W_); 
				//	ctF5 = 0.5 * (_E+_P)*_V + 0.5 * (E_+P_)*V_; 

					blendloc = blend_glob;
					//if (j > 4 && j < HIG-4) blendloc = 1.; 
					Fy1[i][j][k] = (1. - blendloc) * (RU = R * V) + blendloc * ctF1;
					Fy2[i][j][k] = (1. - blendloc) * RU * U + blendloc * ctF2 - sigma_yx;
					Fy3[i][j][k] = (1. - blendloc) * (RU * V + P) + blendloc * ctF3 - sigma_yy;
					Fy4[i][j][k] = (1. - blendloc) * RU * W + blendloc * ctF4 - sigma_yz;
					Fy5[i][j][k] = (1. - blendloc) * V * (K_K_1 * P + 0.5 * R * (U * U + V * V + W * W))
						+ blendloc * ctF5
						- uu * sigma_yx - vv * sigma_yy - ww * sigma_yz + q_y;
				}
			}
		}
	}

	//--- Z-fluxes
	for (i = 0, i_ = 1, i__ = 2; i < LEN; i++, i_++, i__++) {
		for (j = 0, j_ = 1, j__ = 2; j < HIG; j++, j_++, j__++) {
			for (k = 0, k_ = 1; k < DEPP; k++, k_++) {
				//--- Inviscid fluxes
				// further
				_U = zU2[i][j][k] / (_R = zU1[i][j][k]);
				_V = zU3[i][j][k] / _R;
				_W = zU4[i][j][k] / _R;
				_P = (zU5[i][j][k] - 0.5 * _R * (_U * _U + _V * _V + _W * _W)) * K_1;
				_C = sqrt(K * _P / _R);
				_jo = _K_1 * _P;           //_P - _C * _C * _R
				_jp = _W + (RU = _P / (_R * _C));
				_jm = _W - RU;
				// nearer
				U_ = U2z[i][j][k] / (R_ = U1z[i][j][k]);
				V_ = U3z[i][j][k] / R_;
				W_ = U4z[i][j][k] / R_;
				P_ = (U5z[i][j][k] - 0.5 * R_ * (U_ * U_ + V_ * V_ + W_ * W_)) * K_1;
				C_ = sqrt(K * P_ / R_);
				jo_ = _K_1 * P_;           // P_ - C_ * C_ * R_
				jp_ = W_ + (RU = P_ / (R_ * C_));
				jm_ = W_ - RU;
				// linearized procedure
				W = 0.5 * (_W + W_);
				C = 0.5 * (_C + C_);
				C_p = W + C; C_m = W - C; C_o = W;
				// Jp and al_p
				if (C_p > 0.) { al_p = K_1_2 * (_jm - _jp) / _jo; J_p = _jp; }
				else { al_p = K_1_2 * (jm_ - jp_) / jo_; J_p = jp_; }
				// Jm and al_m
				if (C_m > 0.) { al_m = K_1_2 * (_jp - _jm) / _jo; J_m = _jm; }
				else { al_m = K_1_2 * (jp_ - jm_) / jo_; J_m = jm_; }
				// Jo and al_o
				if (C_o > 0.) { al_o = _05K * (_jp - _jm); al_o = -al_o * al_o; J_o = _jo; }
				else { al_o = _05K * (jp_ - jm_); al_o = -al_o * al_o; J_o = jo_; }
				// parameters at the interface
				P = (J_p - J_m) / (al_p - al_m);
				W = J_p - al_p * P;
				R = (J_o - P) / al_o;
				// central flux
				ctF1 = 0.5 * (_R * _W + R_ * W_);
				ctF2 = 0.5 * (_R * _W * _U + R_ * W_ * U_);
				ctF3 = 0.5 * (_R * _W * _V + R_ * W_ * V_);
				ctF4 = 0.5 * (_R * _W * _W + _P + R_ * W_ * W_ + P_);
				ctF5 = 0.5 * (zU5[i][j][k] + _P) * _W + 0.5 * (U5z[i][j][k] + P_) * W_;
				// tangential velocities
				if (W > 0.) { U = _U; V = _V; }
				else { U = U_; V = V_; }
				//--- Viscous fluxes
				//-- nearer cells
				// nearer
				U_ = U2_[i_][j_][k_] / (R_ = U1_[i_][j_][k_]);
				V_ = U3_[i_][j_][k_] / R_;
				W_ = U4_[i_][j_][k_] / R_;
				P_ = (U5_[i_][j_][k_] - 0.5 * R_ * (U_ * U_ + V_ * V_ + W_ * W_)) * K_1;
				E_ = U5_[i_][j_][k_];
				T_ = P_ / (R_GAS * R_);
				// nearer forward
				Ui = U2_[i__][j_][k_] * (RU = 1. / U1_[i__][j_][k_]);
				Vi = U3_[i__][j_][k_] * RU;
				Wi = U4_[i__][j_][k_] * RU;
				// nearer backward
				iU = U2_[i][j_][k_] * (RU = 1. / U1_[i][j_][k_]);
				iV = U3_[i][j_][k_] * RU;
				iW = U4_[i][j_][k_] * RU;
				// nearer upper
				Uj = U2_[i_][j__][k_] * (RU = 1. / U1_[i_][j__][k_]);
				Vj = U3_[i_][j__][k_] * RU;
				Wj = U4_[i_][j__][k_] * RU;
				// nearer lower
				jU = U2_[i_][j][k_] * (RU = 1. / U1_[i_][j][k_]);
				jV = U3_[i_][j][k_] * RU;
				jW = U4_[i_][j][k_] * RU;
				//-- further cells
				// further
				_U = U2_[i_][j_][k] / (_R = U1_[i_][j_][k]);
				_V = U3_[i_][j_][k] / _R;
				_W = U4_[i_][j_][k] / _R;
				_P = (U5_[i_][j_][k] - 0.5 * _R * (_U * _U + _V * _V + _W * _W)) * K_1;
				_E = U5_[i_][j_][k];
				_T = _P / (R_GAS * _R);
				// further forward
				ui = U2_[i__][j_][k] * (RU = 1. / U1_[i__][j_][k]);
				vi = U3_[i__][j_][k] * RU;
				wi = U4_[i__][j_][k] * RU;
				// further backward
				iu = U2_[i][j_][k] * (RU = 1. / U1_[i][j_][k]);
				iv = U3_[i][j_][k] * RU;
				iw = U4_[i][j_][k] * RU;
				// further upper
				uj = U2_[i_][j__][k] * (RU = 1. / U1_[i_][j__][k]);
				vj = U3_[i_][j__][k] * RU;
				wj = U4_[i_][j__][k] * RU;
				// further lower
				ju = U2_[i_][j][k] * (RU = 1. / U1_[i_][j][k]);
				jv = U3_[i_][j][k] * RU;
				jw = U4_[i_][j][k] * RU;
				// derivatives of velocities
				// x
				du_dx = (Ui + ui - iU - iu) * _4deltaX;
				dv_dx = (Vi + vi - iV - iv) * _4deltaX;
				dw_dx = (Wi + wi - iW - iw) * _4deltaX;
				// y
				du_dy = (Uj + uj - jU - ju) * _4deltaY;
				dv_dy = (Vj + vj - jV - jv) * _4deltaY;
				dw_dy = (Wj + wj - jW - jw) * _4deltaY;
				// z
				du_dz = (U_ - _U) * _deltaZ;
				dv_dz = (V_ - _V) * _deltaZ;
				dw_dz = (W_ - _W) * _deltaZ;
				// mean velocities and density
				rr = 0.5 * (R_ + _R);
				uu = 0.5 * (U_ + _U);
				vv = 0.5 * (V_ + _V);
				ww = 0.5 * (W_ + _W);
				// stresses
			/*
				S12 = du_dy + dv_dx;
				S13 = du_dz + dw_dx;
				S23 = dv_dz + dw_dy;
				_S_ = sqrt((du_dx * du_dx + dv_dy * dv_dy + dw_dz * dw_dz) * 2.
					   + (S12 * S12   +   S13 * S13   +   S23 * S23));
				(mu_T = rr * CsDD * _S_);
				//mu_T = r * Cs * delta * delta * | S |
			*/
			/*
				gu[0] = du_dx; gu[1] = dv_dx; gu[2] = dw_dx;
				gu[3] = du_dy; gu[4] = dv_dy; gu[5] = dw_dy;
				gu[6] = du_dz; gu[7] = dv_dz; gu[8] = dw_dz;
					mu_T  = VremenModel(gu, deltaX*deltaY*deltaZ);
					mu_T *= rr;
			*/
				mu_T = 0.5 * (muturb[i_][j_][k] + muturb[i_][j_][k_]);
				mu_E = mu_L + mu_T;
				sigma_zx = mu_E * (dw_dx + du_dz);
				sigma_zy = mu_E * (dw_dy + dv_dz);
				sigma_zz = mu_E * (dw_dz + dw_dz - du_dx - dv_dy) * 0.66666667;
				// heat flux
				lambda_E = lambda_L + mu_T * cp_Pr_T;
				q_z = -lambda_E * (T_ - _T) * _deltaZ;

				// central flux
			//	ctF1 = 0.5 * (_R*_V       + R_*V_); 
			//	ctF2 = 0.5 * (_R*_V*_U    + R_*V_*U_); 
			//      ctF3 = 0.5 * (_R*_V*_V+_P + R_*V_*V_+P_);
			//	ctF4 = 0.5 * (_R*_V*_W    + R_*V_*W_); 
			//	ctF5 = 0.5 * (_E+_P)*_V + 0.5 * (E_+P_)*V_; 
				// summary Z-fluxes evaluation
				blendloc = blend_glob;
				//if (j > 4 && j < HIG-4) blendloc = 1.; 
				Fz1[i][j][k] = (1. - blendloc) * (RU = R * W) + blendloc * ctF1;
				Fz2[i][j][k] = (1. - blendloc) * RU * U + blendloc * ctF2 - sigma_zx;
				Fz3[i][j][k] = (1. - blendloc) * RU * V + blendloc * ctF3 - sigma_zy;
				Fz4[i][j][k] = (1. - blendloc) * (RU * W + P) + blendloc * ctF4 - sigma_zz;
				Fz5[i][j][k] = (1. - blendloc) * W * (K_K_1 * P + 0.5 * R * (U * U + V * V + W * W))
					+ blendloc * ctF5
					- uu * sigma_zx - vv * sigma_zy - ww * sigma_zz + q_z;
			}
		}
	}

} // end Fluxes()


/**************
*  EVOLUTION  *  Compute new/intermediate values of conservative variables
**************/
void Evolution_RK2(void)
{
	register real R;
	register real RU;

	if (Stage == 1) {
		for (i = 1, _i = 0; i < LENN; i++, _i++) {
			for (j = 1, _j = 0; j < HIGG; j++, _j++) { // indices i-1, k-1 etc
				for (k = 1, _k = 0; k < DEPP; k++, _k++) {
					U1p[i][j][k] = (R = U1[i][j][k])
						+ deltaT_X * (Fx1[_i][_j][_k] - Fx1[i][_j][_k])
						+ deltaT_Y * (Fy1[_i][_j][_k] - Fy1[_i][j][_k])
						+ deltaT_Z * (Fz1[_i][_j][_k] - Fz1[_i][_j][k]);
					U2p[i][j][k] = (RU = U2[i][j][k])
						+ deltaT_X * (Fx2[_i][_j][_k] - Fx2[i][_j][_k])
						+ deltaT_Y * (Fy2[_i][_j][_k] - Fy2[_i][j][_k])
						+ deltaT_Z * (Fz2[_i][_j][_k] - Fz2[_i][_j][k])
						+ R * gdeltaT;
					U3p[i][j][k] = U3[i][j][k]
						+ deltaT_X * (Fx3[_i][_j][_k] - Fx3[i][_j][_k])
						+ deltaT_Y * (Fy3[_i][_j][_k] - Fy3[_i][j][_k])
						+ deltaT_Z * (Fz3[_i][_j][_k] - Fz3[_i][_j][k]);
					U4p[i][j][k] = U4[i][j][k]
						+ deltaT_X * (Fx4[_i][_j][_k] - Fx4[i][_j][_k])
						+ deltaT_Y * (Fy4[_i][_j][_k] - Fy4[_i][j][_k])
						+ deltaT_Z * (Fz4[_i][_j][_k] - Fz4[_i][_j][k]);
					U5p[i][j][k] = U5[i][j][k]
						+ deltaT_X * (Fx5[_i][_j][_k] - Fx5[i][_j][_k])
						+ deltaT_Y * (Fy5[_i][_j][_k] - Fy5[_i][j][_k])
						+ deltaT_Z * (Fz5[_i][_j][_k] - Fz5[_i][_j][k])
						+ RU * gdeltaT;
				}
			}
		}

		//Stage = 2; 
		U1_ = U1p, U2_ = U2p, U3_ = U3p, U4_ = U4p, U5_ = U5p;

	} // end if (Stage == 1)

	else {
		for (i = 1, _i = 0; i < LENN; i++, _i++) {
			for (j = 1, _j = 0; j < HIGG; j++, _j++) { // indices i-1, k-1 etc
				for (k = 1, _k = 0; k < DEPP; k++, _k++) {
					U1[i][j][k] = 0.5 * (U1[i][j][k] + (R = U1p[i][j][k])
						+ deltaT_X * (Fx1[_i][_j][_k] - Fx1[i][_j][_k])
						+ deltaT_Y * (Fy1[_i][_j][_k] - Fy1[_i][j][_k])
						+ deltaT_Z * (Fz1[_i][_j][_k] - Fz1[_i][_j][k])
						);
					U2[i][j][k] = 0.5 * (U2[i][j][k] + (RU = U2p[i][j][k])
						+ deltaT_X * (Fx2[_i][_j][_k] - Fx2[i][_j][_k])
						+ deltaT_Y * (Fy2[_i][_j][_k] - Fy2[_i][j][_k])
						+ deltaT_Z * (Fz2[_i][_j][_k] - Fz2[_i][_j][k])
						+ R * gdeltaT
						);
					U3[i][j][k] = 0.5 * (U3[i][j][k] + U3p[i][j][k]
						+ deltaT_X * (Fx3[_i][_j][_k] - Fx3[i][_j][_k])
						+ deltaT_Y * (Fy3[_i][_j][_k] - Fy3[_i][j][_k])
						+ deltaT_Z * (Fz3[_i][_j][_k] - Fz3[_i][_j][k])
						);
					U4[i][j][k] = 0.5 * (U4[i][j][k] + U4p[i][j][k]
						+ deltaT_X * (Fx4[_i][_j][_k] - Fx4[i][_j][_k])
						+ deltaT_Y * (Fy4[_i][_j][_k] - Fy4[_i][j][_k])
						+ deltaT_Z * (Fz4[_i][_j][_k] - Fz4[_i][_j][k])
						);
					U5[i][j][k] = 0.5 * (U5[i][j][k] + U5p[i][j][k]
						+ deltaT_X * (Fx5[_i][_j][_k] - Fx5[i][_j][_k])
						+ deltaT_Y * (Fy5[_i][_j][_k] - Fy5[_i][j][_k])
						+ deltaT_Z * (Fz5[_i][_j][_k] - Fz5[_i][_j][k])
						+ RU * gdeltaT
						);
				}
			}
		}

		//Stage = 1 
		U1_ = U1, U2_ = U2, U3_ = U3, U4_ = U4, U5_ = U5;
	} // end else

} // end Evolution_RK2()

void Evolution_RK3(void)
{
	register real R;
	register real RU;

	switch (Stage) {

	case 1:
		for (i = 1, _i = 0; i < LENN; i++, _i++) {
			for (j = 1, _j = 0; j < HIGG; j++, _j++) { // indices i-1, k-1 etc
				for (k = 1, _k = 0; k < DEPP; k++, _k++) {
					U1p[i][j][k] = (R = U1[i][j][k])
						+ deltaT_X * (Fx1[_i][_j][_k] - Fx1[i][_j][_k])
						+ deltaT_Y * (Fy1[_i][_j][_k] - Fy1[_i][j][_k])
						+ deltaT_Z * (Fz1[_i][_j][_k] - Fz1[_i][_j][k]);
					U2p[i][j][k] = (RU = U2[i][j][k])
						+ deltaT_X * (Fx2[_i][_j][_k] - Fx2[i][_j][_k])
						+ deltaT_Y * (Fy2[_i][_j][_k] - Fy2[_i][j][_k])
						+ deltaT_Z * (Fz2[_i][_j][_k] - Fz2[_i][_j][k])
						+ (DriverF)*gdeltaT;
					U3p[i][j][k] = U3[i][j][k]
						+ deltaT_X * (Fx3[_i][_j][_k] - Fx3[i][_j][_k])
						+ deltaT_Y * (Fy3[_i][_j][_k] - Fy3[_i][j][_k])
						+ deltaT_Z * (Fz3[_i][_j][_k] - Fz3[_i][_j][k]);
					U4p[i][j][k] = U4[i][j][k]
						+ deltaT_X * (Fx4[_i][_j][_k] - Fx4[i][_j][_k])
						+ deltaT_Y * (Fy4[_i][_j][_k] - Fy4[_i][j][_k])
						+ deltaT_Z * (Fz4[_i][_j][_k] - Fz4[_i][_j][k]);
					U5p[i][j][k] = U5[i][j][k]
						+ deltaT_X * (Fx5[_i][_j][_k] - Fx5[i][_j][_k])
						+ deltaT_Y * (Fy5[_i][_j][_k] - Fy5[_i][j][_k])
						+ deltaT_Z * (Fz5[_i][_j][_k] - Fz5[_i][_j][k])
						; // + DriverF * U2[i][j][k]/U1[i][j][k] * gdeltaT;
				}
			}
		}

		// to Stage = 2; 
		U1_ = U1p, U2_ = U2p, U3_ = U3p, U4_ = U4p, U5_ = U5p;

		break;

	case 2:
		for (i = 1, _i = 0; i < LENN; i++, _i++) {
			for (j = 1, _j = 0; j < HIGG; j++, _j++) { // indices i-1, k-1 etc
				for (k = 1, _k = 0; k < DEPP; k++, _k++) {
					U1p[i][j][k] = 0.75 * U1[i][j][k] + 0.25 * ((R = U1p[i][j][k])
						+ deltaT_X * (Fx1[_i][_j][_k] - Fx1[i][_j][_k])
						+ deltaT_Y * (Fy1[_i][_j][_k] - Fy1[_i][j][_k])
						+ deltaT_Z * (Fz1[_i][_j][_k] - Fz1[_i][_j][k]));
					U2p[i][j][k] = 0.75 * U2[i][j][k] + 0.25 * ((RU = U2p[i][j][k])
						+ deltaT_X * (Fx2[_i][_j][_k] - Fx2[i][_j][_k])
						+ deltaT_Y * (Fy2[_i][_j][_k] - Fy2[_i][j][_k])
						+ deltaT_Z * (Fz2[_i][_j][_k] - Fz2[_i][_j][k])
						+ (DriverF)*gdeltaT);
					U3p[i][j][k] = 0.75 * U3[i][j][k] + 0.25 * (U3p[i][j][k]
						+ deltaT_X * (Fx3[_i][_j][_k] - Fx3[i][_j][_k])
						+ deltaT_Y * (Fy3[_i][_j][_k] - Fy3[_i][j][_k])
						+ deltaT_Z * (Fz3[_i][_j][_k] - Fz3[_i][_j][k]));
					U4p[i][j][k] = 0.75 * U4[i][j][k] + 0.25 * (U4p[i][j][k]
						+ deltaT_X * (Fx4[_i][_j][_k] - Fx4[i][_j][_k])
						+ deltaT_Y * (Fy4[_i][_j][_k] - Fy4[_i][j][_k])
						+ deltaT_Z * (Fz4[_i][_j][_k] - Fz4[_i][_j][k]));
					U5p[i][j][k] = 0.75 * U5[i][j][k] + 0.25 * (U5p[i][j][k]
						+ deltaT_X * (Fx5[_i][_j][_k] - Fx5[i][_j][_k])
						+ deltaT_Y * (Fy5[_i][_j][_k] - Fy5[_i][j][_k])
						+ deltaT_Z * (Fz5[_i][_j][_k] - Fz5[_i][_j][k])
						); // + DriverF * U2p[i][j][k]/U1p[i][j][k] * gdeltaT  );
				}
			}
		}

		// to Stage = 3; 
		U1_ = U1p, U2_ = U2p, U3_ = U3p, U4_ = U4p, U5_ = U5p;
		break;

	default:
		for (i = 1, _i = 0; i < LENN; i++, _i++) {
			for (j = 1, _j = 0; j < HIGG; j++, _j++) { // indices i-1, k-1 etc
				for (k = 1, _k = 0; k < DEPP; k++, _k++) {
					U1[i][j][k] = 1. / 3. * U1[i][j][k] + 2. / 3. * ((R = U1p[i][j][k])
						+ deltaT_X * (Fx1[_i][_j][_k] - Fx1[i][_j][_k])
						+ deltaT_Y * (Fy1[_i][_j][_k] - Fy1[_i][j][_k])
						+ deltaT_Z * (Fz1[_i][_j][_k] - Fz1[_i][_j][k]));
					U2[i][j][k] = 1. / 3. * U2[i][j][k] + 2. / 3. * ((RU = U2p[i][j][k])
						+ deltaT_X * (Fx2[_i][_j][_k] - Fx2[i][_j][_k])
						+ deltaT_Y * (Fy2[_i][_j][_k] - Fy2[_i][j][_k])
						+ deltaT_Z * (Fz2[_i][_j][_k] - Fz2[_i][_j][k])
						+ (DriverF)*gdeltaT);
					U3[i][j][k] = 1. / 3. * U3[i][j][k] + 2. / 3. * (U3p[i][j][k]
						+ deltaT_X * (Fx3[_i][_j][_k] - Fx3[i][_j][_k])
						+ deltaT_Y * (Fy3[_i][_j][_k] - Fy3[_i][j][_k])
						+ deltaT_Z * (Fz3[_i][_j][_k] - Fz3[_i][_j][k]));
					U4[i][j][k] = 1. / 3. * U4[i][j][k] + 2. / 3. * (U4p[i][j][k]
						+ deltaT_X * (Fx4[_i][_j][_k] - Fx4[i][_j][_k])
						+ deltaT_Y * (Fy4[_i][_j][_k] - Fy4[_i][j][_k])
						+ deltaT_Z * (Fz4[_i][_j][_k] - Fz4[_i][_j][k]));
					U5[i][j][k] = 1. / 3. * U5[i][j][k] + 2. / 3. * (U5p[i][j][k]
						+ deltaT_X * (Fx5[_i][_j][_k] - Fx5[i][_j][_k])
						+ deltaT_Y * (Fy5[_i][_j][_k] - Fy5[_i][j][_k])
						+ deltaT_Z * (Fz5[_i][_j][_k] - Fz5[_i][_j][k])
						); // + DriverF * U2p[i][j][k]/U1p[i][j][k] * gdeltaT  );
				}
			}
		}

		// to Stage = 1
		U1_ = U1, U2_ = U2, U3_ = U3, U4_ = U4, U5_ = U5;
		break;

	}

} // end Evolution_RK3

void Evolution(int opt)
{
	if (opt == 2)
		Evolution_RK2();

	if (opt == 3)
		Evolution_RK3();

}
/************************
*  Turbulent Viscosity  *
************************/
void Eval_Mu_Turb(void)
{

	real Uavg, Vavg, Wavg, vort, disc, dtmp, u_tau;

	for (j = 1; j < HIGG; j++) {

		Uavg = 0.; Vavg = 0.; Wavg = 0.;
		for (i = 1; i < LENN; i++) {
			for (k = 1; k < DEPP; k++) {
				Uavg += U2[i][j][k];
				Vavg += U3[i][j][k];
				Wavg += U4[i][j][k];
			}
		}

		Uavg1[j] = Uavg / (real)(LEN * DEP);
		Vavg1[j] = Vavg / (real)(LEN * DEP);
		Wavg1[j] = Wavg / (real)(LEN * DEP);
	}

	for (i = 1; i < LENN; i++) {
		for (j = 1; j < HIGG; j++) { // indices i-1, k-1 etc
			for (k = 1; k < DEPP; k++) {

				du_dx = (U2_[i + 1][j][k] - U2_[i - 1][j][k]) * 0.5 * _deltaX;
				dv_dx = (U3_[i + 1][j][k] - U3_[i - 1][j][k]) * 0.5 * _deltaX;
				dw_dx = (U4_[i + 1][j][k] - U4_[i - 1][j][k]) * 0.5 * _deltaX;
				//du_dy = (U2_[i][j+1][k] - Uavg1[j+1] - U2_[i][j-1][k] + Uavg1[j-1]) * 0.5 * _deltaY; 
				//dv_dy = (U3_[i][j+1][k] - Vavg1[j+1] - U3_[i][j-1][k] + Vavg1[j-1]) * 0.5 * _deltaY; 
				//dw_dy = (U4_[i][j+1][k] - Wavg1[j+1] - U4_[i][j-1][k] + Wavg1[j-1]) * 0.5 * _deltaY; 
				du_dy = (U2_[i][j + 1][k] - U2_[i][j - 1][k]) * 0.5 * _deltaY;
				dv_dy = (U3_[i][j + 1][k] - U3_[i][j - 1][k]) * 0.5 * _deltaY;
				dw_dy = (U4_[i][j + 1][k] - U4_[i][j - 1][k]) * 0.5 * _deltaY;
				du_dz = (U2_[i][j][k + 1] - U2_[i][j][k - 1]) * 0.5 * _deltaZ;
				dv_dz = (U3_[i][j][k + 1] - U3_[i][j][k - 1]) * 0.5 * _deltaZ;
				dw_dz = (U4_[i][j][k + 1] - U4_[i][j][k - 1]) * 0.5 * _deltaZ;

				gu[0] = du_dx; gu[1] = dv_dx; gu[2] = dw_dx;
				gu[3] = du_dy; gu[4] = dv_dy; gu[5] = dw_dy;
				gu[6] = du_dz; gu[7] = dv_dz; gu[8] = dw_dz;

				// sgs
				tmp[i][j][k] = VremenModel(gu, deltaX, deltaY, deltaZ);
				//if(j <= 3 || j > (HIG-3)) 
				//muturb[i][j][k] = 0.; 
				//else  
				muturb[i][j][k] = tmp[i][j][k];
				//if(j <= 1 || j >= (HIG-1)) 
				//muturb[i][j][k] = 3.0 * tmp[i][j][k]; 

#ifdef ADAWN	
		 // rans (simple mixing length model)
				vort = eval_vort_mag(gu);
				//!!! YL
				u_tau = sqrt(DriverF);
				disc = ((real)j - 0.5) * deltaY;
				disc = (disc > 1.) ? 2. - disc : disc;
				eval_mut_rans(u_tau, disc, mu_L, vort, &dtmp);
				//adjust viscosity locally
				int match_indx = 5;
				if (j > match_indx && j < (HIG - match_indx + 1))
				{
					muturb_rans[i][j][k] = 0.;
					muturb_coef[i][j][k] = 1.;
				}
				else
				{
					if (j <= match_indx) {
						double alpha_niu = 0.0;
						double local_yp = ((double)(j)-0.5) * deltaY * u_tau / mu_L;
						alpha_niu = blend_dyn[0][i] + blend_dyn[1][i] * local_yp;
						
						if (alpha_niu < 0.0) {
				                    alpha_niu = 0.0;
                				}
                				else if (alpha_niu > 1.0) {
                    				    alpha_niu = 1.0;
                				}
						
	
						muturb_rans[i][j][k] = alpha_niu * dtmp;
						muturb[i][j][k] *=  1.0 - alpha_niu;

						//printf("%lf      ", alpha_niu );
						
						// muturb_rans[i][j][k] = blend_dyn[0][i] * dtmp;
						// muturb[i][j][k] *= 1.0 - blend_dyn[0][i];
					}

					/*
					if (j == match_indx) 
					{
						printf(" \n");
						printf("alpha_niu :     ");
					}
					*/

					if (j >= (HIG - match_indx + 1)) {
						double alpha_niu = 0.0;
                                                double local_yp = ((double)(j)-0.5) * deltaY * u_tau / mu_L;
                                                alpha_niu = blend_dyn[0][i] + blend_dyn[1][i] * local_yp;
						if (alpha_niu < 0.0) {
                                                    alpha_niu = 0.0;
                                                }
                                                else if (alpha_niu > 1.0) {
                                                    alpha_niu = 1.0;
                                                }
                                                muturb_rans[i][j][k] = alpha_niu * dtmp;
                                                muturb[i][j][k] *= 1.0 - alpha_niu;
						
						// muturb_rans[i][j][k] = blend_dyn[1][i] * dtmp;
						// muturb[i][j][k] *= 1.0 - blend_dyn[1][i];
					}
				}
#endif	     

			}
		}
	}

	//filtering 
	/*
	 for (i = 2; i < LENN-1; i++) {
	   for (j = 2; j < HIGG-1; j++) { // indices i-1, k-1 etc
	 for (k = 2; k < DEPP-1; k++) {
		  muturb[i][j][k] = 0.25 * tmp[i][j][k]
					  + 0.125 * tmp[i+1][j][k] + 0.125 * tmp[i-1][j][k]
					  + 0.125 * tmp[i][j+1][k] + 0.125 * tmp[i][j-1][k]
					  + 0.125 * tmp[i][j][k+1] + 0.125 * tmp[i][j][k-1];
	 }
	   }
	 }
	*/
}

/**********************
*  Adjust Diff Coeff  *
**********************/
void Adjust_Diff_Coeff(void) {

	real MatchU[2];

	// velocity with law of the wall 
	real utau = sqrt(DriverF);
	real ulow[3];
	ulow[0] = utau / 0.41 * log(deltaY / 2. * utau / mu_L) + utau * 5.2;
	ulow[1] = utau / 0.41 * log(deltaY / 2. * 3. * utau / mu_L) + utau * 5.2;
	ulow[2] = utau / 0.41 * log(deltaY * utau / mu_L) + utau * 5.2;

	// velocity at match location 
	MatchU[0] = 0.; MatchU[1] = 0.;
	for (i = 1; i < LENN; i++) {
		for (k = 1; k < DEPP; k++) {
			MatchU[0] += U2_[i][1][k]; MatchU[0] += U2_[i][HIG][k];
			MatchU[1] += U2_[i][2][k]; MatchU[1] += U2_[i][HIG - 1][k];
		}
	}
	MatchU[0] = MatchU[0] / 2. / (double)LEN / (double)DEP;
	MatchU[1] = MatchU[1] / 2. / (double)LEN / (double)DEP;

	//if ( MatchU[0] > ulow[0] || MatchU[1] > ulow[1] ) blend *= 0.98; 
	//if ( MatchU[0] < ulow[0] || MatchU[1] < ulow[1] ) blend *= 1.02; 
	if (MatchU[0] < ulow[0]) blend *= 0.998;
	if (MatchU[0] > ulow[0]) blend *= 1.002;
	//if (LinearU > ulow[2] ) blend *= 1.002;
	//if (LinearU < ulow[2] ) blend *= 0.998;

	// set the bounds
	if (blend > 1.00) blend = 1.00;
	if (blend < 0.9) blend = 0.9;

	printf("%lf %lf %lf %lf %lf %lf\n", ulow[0], ulow[1], ulow[2], MatchU[0], MatchU[1], LinearU);
}

/*******************
*  Mass Flow Rate  *
*******************/
void Eval_Flow_Rate(void) {

	i = 1;
	MassRate = 0.;
	for (j = 1; j < HIGG; j++) { // indices i-1, k-1 etc
		for (k = 1; k < DEPP; k++) {
			MassRate += U2_[i][j][k];
		}
	}
	MassRate = MassRate / (double)HIG / (double)DEP;
}

/************
*  Filter   *  // to avoid suprious oscillation
************/
void Filter(void)
{
	for (i = 1; i < LENN; i++) {
		for (k = 1; k < DEPP; k++) {
			for (j = 1; j < HIGG; j++) { // indices i-1, k-1 etc
				U3p[i][j][k] = 0.05 * U3_[i][j - 1][k] + 0.9 * U3_[i][j][k] + 0.05 * U3_[i][j + 1][k];
			}
			/*
				 U2p[i][1][k] = 0.125*U2_[i+1][1][k] + 0.125*U2_[i-1][1][k]
						  + 0.125*U2_[i][1][k+1] + 0.125*U2_[i][1][k-1] + 0.50*U2_[i][1][k];
				 U3p[i][1][k] = 0.125*U3_[i+1][1][k] + 0.125*U3_[i-1][1][k]
						  + 0.125*U3_[i][1][k+1] + 0.125*U3_[i][1][k-1] + 0.50*U3_[i][1][k];
				 U4p[i][1][k] = 0.125*U4_[i+1][1][k] + 0.125*U4_[i-1][1][k]
						  + 0.125*U4_[i][1][k+1] + 0.125*U4_[i][1][k-1] + 0.50*U4_[i][1][k];

				 U2p[i][HIG][k] = 0.125*U2_[i+1][HIG][k] + 0.125*U2_[i-1][HIG][k]
						  + 0.125*U2_[i][HIG][k+1] + 0.125*U2_[i][HIG][k-1] + 0.50*U2_[i][HIG][k];
				 U3p[i][HIG][k] = 0.125*U3_[i+1][HIG][k] + 0.125*U3_[i-1][HIG][k]
						  + 0.125*U3_[i][HIG][k+1] + 0.125*U3_[i][HIG][k-1] + 0.50*U3_[i][HIG][k];
				 U4p[i][HIG][k] = 0.125*U4_[i+1][HIG][k] + 0.125*U4_[i-1][HIG][k]
						  + 0.125*U4_[i][HIG][k+1] + 0.125*U4_[i][HIG][k-1] + 0.50*U4_[i][HIG][k];
			*/
		}
	}

	for (i = 1; i < LENN; i++) {
		for (k = 1; k < DEPP; k++) {
			for (j = 1; j < HIGG; j++) { // indices i-1, k-1 etc
				U3_[i][j][k] = U3p[i][j][k];
			}

		}
	}
}

/************
*  SENSORS  *
************/
void Sensors(int myid)
{
	// root process only can write to its own "sensors.*" files
	if (myid != 0)
		return;

	//
	i = LEN / 2;
	//
	if (Stage == 1) {
		// mean fuxes at (1) stage
		for (j = 0, M = I = E = 0.; j < HIG; j++) {
			for (k = 0; k < DEP; k++) {
				M += Fx1[i][j][k]; // mass
				I += Fx2[i][j][k]; // x-momentum
				E += Fx5[i][j][k]; // energy
			}
		}
		// mean wall stresses at (1) stage
		tauZ = tauY = 0.;
		for (j = 0; j < HIG; j++)
			tauZ += Fz2[i][j][0] - Fz2[i][j][DEP];
		for (k = 0; k < DEP; k++)
			tauY += Fy2[i][0][k] - Fy2[i][HIG][k];
	}
	else {
		//-- output to "sensors.txt"
		// plus mean fuxes at (2) stage
		for (j = 0; j < HIG; j++) {
			for (k = 0; k < DEP; k++) {
				M += Fx1[i][j][k]; // mass
				I += Fx2[i][j][k]; // momentum
				E += Fx5[i][j][k]; // energy
			}
		}
		// plus mean wall stresses at (2) stage
		for (j = 0; j < HIG; j++)
			tauZ += Fz2[i][j][0] - Fz2[i][j][DEP];
		for (k = 0; k < DEP; k++)
			tauY += Fy2[i][0][k] - Fy2[i][HIG][k];
		// write mean fluxes and mean wall stess
		fprintf(pFsen_txt, "%d %.3f %.1f %.0f %.2f\n", step,
			0.5 * M / (DEP * HIG), 0.5 * I / (DEP * HIG), 0.5 * E / (DEP * HIG),
			0.5 * (tauZ * deltaY + tauZ * deltaZ) / (2. * (DEP * deltaZ + HIG * deltaY)));
		//-- output to "sensors.bin"
		// number of current step
		fwrite(&step, sizeof step, 1, pFsen_bin);
		// writes HIG-sized blocks of current R and RU
		// R
		k = DEP / 2 + 1;
		i = LEN / 2;
		for (j = 1; j < HIGG; j++)
			sensors[j - 1] = U1[i][j][k];
		fwrite(sensors, (sizeof sensors[0]) * HIG, 1, pFsen_bin);
		// RU
		for (j = 1; j < HIGG; j++)
			sensors[j - 1] = U2[i][j][k];
		fwrite(sensors, (sizeof sensors[0]) * HIG, 1, pFsen_bin);
	} // end else

} // end Sensors()


/***********
*  OUTPUT  *
***********/
void OutputV()
{
	real buf;
	char filename[50];

	FILE* fh;

	// open/create file
	sprintf(filename, "data_%d_%d.vtr", step, myid);
	printf("%d outputs\n", myid + 1);
	fflush(stdout);
	fh = fopen(filename, "w");

	fprintf(fh, "<?xml version=\"1.0\"?>\n");
	fprintf(fh, "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fh, "<RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n", DEP - 1, HIG - 1, LEN - 1);
	fprintf(fh, "<Piece Extent=\"0 %d 0 %d 0 %d\">\n", DEP - 1, HIG - 1, LEN - 1);

	//output coordinate info 
	fprintf(fh, "<Coordinates>\n");

	fprintf(fh, "<DataArray type=\"Float64\" format=\"ascii\">\n");
	fprintf(fh, "%lf ", 0.0);
	for (k = 1; k < DEPP; k++)
		fprintf(fh, "%lf ", deltaZ * (double)k);
	fprintf(fh, "\n");
	fprintf(fh, "</DataArray>\n");

	fprintf(fh, "<DataArray type=\"Float64\" format=\"ascii\">\n");
	for (j = 1; j < HIGG; j++)
		fprintf(fh, "%lf ", deltaY * (double)j);
	fprintf(fh, "\n");
	fprintf(fh, "</DataArray>\n");

	fprintf(fh, "<DataArray type=\"Float64\" format=\"ascii\">\n");
	for (i = 1; i < LENN; i++)
		fprintf(fh, "%lf ", deltaX * (double)((LEN * myid) + i));
	fprintf(fh, "\n");
	fprintf(fh, "</DataArray>\n");

	fprintf(fh, "</Coordinates>\n");

	//output pointwise data
	fprintf(fh, "<PointData>\n");

	fprintf(fh, "<DataArray type=\"Float32\" Name=\"Density\" format=\"ascii\">\n");
	for (i = 1; i < LENN; i++) {
		for (j = 1; j < HIGG; j++) {
			for (k = 1; k < DEPP; k++) {
				fprintf(fh, "%12.8f\n", U1[i][j][k]);
			}
		}
	}
	fprintf(fh, "</DataArray>\n");

	fprintf(fh, "<DataArray type=\"Float32\" Name=\"X-momentum\" format=\"ascii\">\n");
	for (i = 1; i < LENN; i++) {
		for (j = 1; j < HIGG; j++) {
			for (k = 1; k < DEPP; k++) {
				fprintf(fh, "%12.8f\n", U2[i][j][k]);
			}
		}
	}
	fprintf(fh, "</DataArray>\n");
	fprintf(fh, "<DataArray type=\"Float32\" Name=\"Y-momentum\" format=\"ascii\">\n");
	for (i = 1; i < LENN; i++) {
		for (j = 1; j < HIGG; j++) {
			for (k = 1; k < DEPP; k++) {
				fprintf(fh, "%12.8f\n", U3[i][j][k]);
			}
		}
	}
	fprintf(fh, "</DataArray>\n");
	fprintf(fh, "<DataArray type=\"Float32\" Name=\"Z-momentum\" format=\"ascii\">\n");
	for (i = 1; i < LENN; i++) {
		for (j = 1; j < HIGG; j++) {
			for (k = 1; k < DEPP; k++) {
				fprintf(fh, "%12.8f\n", U4[i][j][k]);
			}
		}
	}
	fprintf(fh, "</DataArray>\n");

	fprintf(fh, "<DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">\n");

	for (i = 1; i < LENN; i++) {
		for (j = 1; j < HIGG; j++) {
			for (k = 1; k < DEPP; k++) {
				U = U2[i][j][k] / (R = U1[i][j][k]);
				V = U3[i][j][k] / R;
				W = U4[i][j][k] / R;
				P = (U5[i][j][k] - 0.5 * R * (U * U + V * V + W * W)) * K_1;
				//fprintf(fh, "%12.8f\n", muturb[i][j][k]);  
				fprintf(fh, "%12.8f\n", P);
			}
		}
	}

	fprintf(fh, "</DataArray>\n");

	fprintf(fh, "</PointData>\n");

	fprintf(fh, "</Piece>\n");
	fprintf(fh, "</RectilinearGrid>\n");
	fprintf(fh, "</VTKFile>\n");

	fclose(fh);
	printf("%d finished output\n", myid + 1);
	fflush(stdout);

}  // end OutputV() 

void OutputS()
{
	real buf;
	char filename[50];

	//*******************
	//MPI_File fh;
	//MPI_Status status;
	//*******************

	FILE* fh;

	// open/create file
	sprintf(filename, "data_%d_%d.vtr", step, myid);
	printf("%d outputs\n", myid + 1);
	//
	fh = fopen(filename, "w");
	//MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR,
	//		MPI_INFO_NULL, &fh);
	//for diagnostics
	//fprintf(stderr, "%d: before Set_view: errno %u, %u \n", myid+1, errno, fh);
	//perror(sys_errlist[errno]);
	//MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	//fprintf(stderr, "%d after Set_view \n", myid + 1);

	// header
  //  buf = (real)LEN;
  //  MPI_File_write(fh, &buf, 1, MPI_DOUBLE, &status);
  //  buf = (real)HIG;
  //  MPI_File_write(fh, &buf, 1, MPI_DOUBLE, &status);
  //  buf = (real)deltaX;
  //  MPI_File_write(fh, &buf, 1, MPI_DOUBLE, &status);
  //  buf = (real)deltaY;
  //  MPI_File_write(fh, &buf, 1, MPI_DOUBLE, &status);

	fprintf(fh, "<?xml version=\"1.0\"?>\n");
	fprintf(fh, "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(fh, "<RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n", 0, HIG - 1, LEN - 1);
	fprintf(fh, "<Piece Extent=\"0 %d 0 %d 0 %d\">\n", 0, HIG - 1, LEN - 1);

	//output coordinate info 
	fprintf(fh, "<Coordinates>\n");

	fprintf(fh, "<DataArray type=\"Float64\" format=\"ascii\">\n");
	fprintf(fh, "%lf ", deltaZ * 0.0);
	fprintf(fh, "\n");
	fprintf(fh, "</DataArray>\n");

	fprintf(fh, "<DataArray type=\"Float64\" format=\"ascii\">\n");
	for (j = 1; j < HIGG; j++)
		fprintf(fh, "%lf ", deltaY * (double)j);
	fprintf(fh, "\n");
	fprintf(fh, "</DataArray>\n");

	fprintf(fh, "<DataArray type=\"Float64\" format=\"ascii\">\n");
	for (i = 1; i < LENN; i++)
		fprintf(fh, "%lf ", deltaX * (double)(LEN * myid) + (double)i);
	fprintf(fh, "\n");
	fprintf(fh, "</DataArray>\n");

	fprintf(fh, "</Coordinates>\n");

	//output pointwise data
	fprintf(fh, "<PointData>\n");
	fprintf(fh, "<DataArray type=\"Float32\" Name=\"vorticity\" format=\"ascii\">\n");
	k = DEP / 2 + 1;
	for (i = 1; i < LENN; i++) {
		for (j = 1; j < HIGG; j++) {
			// vorticity plot
			dv_dx = (U3[i + 1][j][k] / U1[i + 1][j][k]
				- U3[i - 1][j][k] / U1[i - 1][j][k]) * _2deltaX;
			dw_dx = (U4[i + 1][j][k] / U1[i + 1][j][k]
				- U4[i - 1][j][k] / U1[i - 1][j][k]) * _2deltaX;
			du_dy = (U2[i][j + 1][k] / U1[i][j + 1][k]
				- U2[i][j - 1][k] / U1[i][j - 1][k]) * _2deltaY;
			dw_dy = (U4[i][j + 1][k] / U1[i][j + 1][k]
				- U4[i][j - 1][k] / U1[i][j - 1][k]) * _2deltaY;
			du_dz = (U2[i][j][k + 1] / U1[i][j][k + 1]
				- U2[i][j][k - 1] / U1[i][j][k - 1]) * _2deltaZ;
			dv_dz = (U3[i][j][k + 1] / U1[i][j][k + 1]
				- U3[i][j][k - 1] / U1[i][j][k - 1]) * _2deltaZ;
			U = 0.5 * (dw_dy - dv_dz);
			V = 0.5 * (du_dz - dw_dx);
			W = 0.5 * (dv_dx - du_dy);
			buf = sqrt(U * U + V * V + W * W);
			//  MPI_File_write(fh, &buf, 1, MPI_DOUBLE, &status);
			fprintf(fh, "%12.5f\n", buf);
		}
	}
	fprintf(fh, "</DataArray>\n");
	fprintf(fh, "</PointData>\n");

	fprintf(fh, "</Piece>\n");
	fprintf(fh, "</RectilinearGrid>\n");
	fprintf(fh, "</VTKFile>\n");

	//
	//MPI_File_close(&fh);
	//
	fclose(fh);
	printf("%d finished output\n", myid + 1);

} // end OutputS()
void Output(char opt)
{
	if (opt == 'V')          // output volume data  
		OutputV();
	else if (opt == 'S')     // output slice data  
		OutputS();
	else
		printf("Output Type Error!\n");
}

/*************
 * Dump Wall *  Write the complete solution data at wall cell
*************/
#ifdef STORE_WALL_CELL
void Write_Wall_Cell_Data(void)
{
	FILE* fh;
	char filename[30];

	// output top-bottom quatities
	sprintf(filename, "topbot_%d.%d", step, myid);
	fh = fopen(filename, "w");

	fprintf(fh, "%% x,  z, U0_log, U0_lin, U0_bar, W0_log, W0_lin, W0_bar, V0, P0, U1_log,      U1_lin, U1_bar, W1_log, W1_lin, W1_bar, V1, P1\n");
	for (i = 1; i < LENN; i++) {
		for (k = 1; k < DEPP; k++) {
			fprintf(fh, "%lf %lf ", deltaX * ((real)i - 0.5), deltaZ * ((real)k - 0.5));
			for (j = 0; j < 6; j++)
				fprintf(fh, "%lf ", wc_coeffs[i - 1][k - 1][j]);

			U = U2[i][1][k] / (R = U1[i][1][k]);
			V = U3[i][1][k] / R;
			W = U4[i][1][k] / R;
			P = (U5[i][1][k] - 0.5 * R * (U * U + V * V + W * W)) * K_1;

			fprintf(fh, "%lf ", V);
			fprintf(fh, "%lf ", P);

			fprintf(fh, "     ");

			for (j = 6; j < 12; j++)
				fprintf(fh, "%lf ", wc_coeffs[i - 1][k - 1][j]);

			U = U2[i][HIG][k] / (R = U1[i][HIG][k]);
			V = U3[i][HIG][k] / R;
			W = U4[i][HIG][k] / R;
			P = (U5[i][HIG][k] - 0.5 * R * (U * U + V * V + W * W)) * K_1;

			fprintf(fh, "%lf ", V);
			fprintf(fh, "%lf ", P);
			fprintf(fh, "\n");
		}
	}

	fclose(fh);
}
#endif

/*************
* Statistics *  Take sample along y-direction (wall-normal)
*************/
void reset_Statistics()
{
	accumul_dT = 0.;
	for (j = 1; j < HIGG; j++) {
		Uavg0[j] = 0.;
		Vavg0[j] = 0.;
		Wavg0[j] = 0.;

		UU0[j] = 0.;
		VV0[j] = 0.;
		WW0[j] = 0.;

		UV0[j] = 0.;
		VW0[j] = 0.;
		UW0[j] = 0.;
	}

	//    for (j = 1; j < HIGG; j++)
	//    for (i = 1; i < LENN; i++)
	//    {
	//       nwall_vec[j][i] = 0.;
	//    } 
}

void write_Statistics()
{
	char filename[30];
	FILE* fh = NULL;

	sprintf(filename, "stat.%d", myid);
	fh = fopen(filename, "w");

	for (j = 1; j < HIGG; j++) {

		// wall-normal cooridnate 
		fprintf(fh, "%11.5f ", ((real)j - 0.5) * deltaY);
		// Uavg 
		fprintf(fh, "%11.5f  ", Uavg0[j] / accumul_dT);
		// Vavg
		fprintf(fh, "%11.5f  ", Vavg0[j] / accumul_dT);
		// Wavg  
		fprintf(fh, "%11.5f  ", Wavg0[j] / accumul_dT);
		// UU  
		fprintf(fh, "%11.5f  ", UU0[j] / accumul_dT);
		// VV  
		fprintf(fh, "%11.5f  ", VV0[j] / accumul_dT);
		// WW  
		fprintf(fh, "%11.5f  ", WW0[j] / accumul_dT);
		// UV  
		fprintf(fh, "%11.5f  ", UV0[j] / accumul_dT);
		// VW   
		fprintf(fh, "%11.5f  ", VW0[j] / accumul_dT);
		// UW 
		fprintf(fh, "%11.5f  ", UW0[j] / accumul_dT);
		fprintf(fh, "\n");
	}

	fclose(fh);

	// output pid controller quantities 
	sprintf(filename, "blend.%d", myid);
	fh = fopen(filename, "w");

	for (i = 1; i < LENN; i++)
		fprintf(fh, "%lf  %lf\n", blend_dyn[0][i], blend_dyn[1][i]);

	for (k = 0; k < 2; k++)
		for (i = 1; i < LENN; i++)
			for (j = 0; j < 4; j++)
				fprintf(fh, "%lf\n", pid_stat_log[k][i][j]);

	fclose(fh);

}

void Statistics()
{
	real Uavg, Vavg, Wavg;
	real UU, VV, WW;
	real UV, VW, UW;
	real norm = (real)(LEN * DEP);
	real sample_dT;

	sample_dT = (real)f_step;
	accumul_dT += sample_dT;

	for (j = 1; j < HIGG; j++) {

		Uavg = 0.; Vavg = 0.; Wavg = 0.;
		UU = 0.; VV = 0.; WW = 0.;
		UV = 0.; VW = 0.; UW = 0.;

		for (i = 1; i < LENN; i++) {
			for (k = 1; k < DEPP; k++) {
				Uavg += U2[i][j][k];
				Vavg += U3[i][j][k];
				Wavg += U4[i][j][k];

				UU += U2[i][j][k] * U2[i][j][k];
				VV += U3[i][j][k] * U3[i][j][k];
				WW += U4[i][j][k] * U4[i][j][k];

				UV += U2[i][j][k] * U3[i][j][k];
				VW += U3[i][j][k] * U4[i][j][k];
				UW += U2[i][j][k] * U4[i][j][k];
			}
		}

		Uavg /= norm; Vavg /= norm; Wavg /= norm;
		UU /= norm; VV /= norm; WW /= norm;
		UV /= norm; VW /= norm; UW /= norm;

		Uavg0[j] += sample_dT * Uavg;
		Vavg0[j] += sample_dT * Vavg;
		Wavg0[j] += sample_dT * Wavg;

		UU0[j] += sample_dT * UU;
		VV0[j] += sample_dT * VV;
		WW0[j] += sample_dT * WW;

		UV0[j] += sample_dT * UV;
		VW0[j] += sample_dT * VW;
		UW0[j] += sample_dT * UW;

	}

	PID_time++;
	for (i = 1; i < LENN; i++)
		for (j = 1; j < HIGG; j++) {
			Uavg = 0.;
			for (k = 1; k < DEPP; k++) {
				Uavg += U2[i][j][k];
			}
			nwall_vec[j][i] += Uavg / (double)DEP;
		}

}


/*************
*  Finalize  *   Save solution to binary files "backup.*"
*************/
void Finalize()
{
	real buf;
	char filename[30];
	MPI_File fh;
	MPI_Status status;

	FILE* pF;
	FILE* fp;

	sprintf(filename, "backup.%d", myid);
	MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
	MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

	//  fp = fopen(filename, "w"); 

	for (i = 1; i < LENN; i++) {
		for (j = 1; j < HIGG; j++) {
			MPI_File_write(fh, &U1[i][j][1], DEP, MPI_DOUBLE, &status);
			MPI_File_write(fh, &U2[i][j][1], DEP, MPI_DOUBLE, &status);
			MPI_File_write(fh, &U3[i][j][1], DEP, MPI_DOUBLE, &status);
			MPI_File_write(fh, &U4[i][j][1], DEP, MPI_DOUBLE, &status);
			MPI_File_write(fh, &U5[i][j][1], DEP, MPI_DOUBLE, &status);

			//      fwrite(&U1[i][j][1], sizeof(real), DEP, fp);
			//      fwrite(&U2[i][j][1], sizeof(real), DEP, fp);
			//      fwrite(&U3[i][j][1], sizeof(real), DEP, fp);
			//      fwrite(&U4[i][j][1], sizeof(real), DEP, fp);
			//      fwrite(&U5[i][j][1], sizeof(real), DEP, fp);
		}
	}
	buf = (real)step;
	MPI_File_write(fh, &buf, 1, MPI_DOUBLE, &status);
	MPI_File_close(&fh);

	//  fclose(fp); 
	/*
	  // close "sensors.*" files
	  if (myid == 0) {
		fclose(pFsen_bin);
		fclose(pFsen_txt);
	  }
	*/
	// syncronization after backing up
	MPI_Barrier(globleMpicom);

} // end Finalize()




/*
// continuation flag
int continFlag = 1;
char a;
unsigned nStage=3;
MPI_File fh;
MPI_Status status;
*/






void checkStopfile()
{
	// check single-byte record in "stopfile"
	if (step % 100 == 0 && myid == 0) {
		MPI_File_open(MPI_COMM_SELF, "stopfile",
			MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		MPI_File_read(fh, &a, 1, MPI_CHAR, &status);
		MPI_File_close(&fh);
		//
		if (a == 's') {
			continFlag = 0;
			MPI_File_open(MPI_COMM_SELF, "stopfile",
				MPI_MODE_CREATE | MPI_MODE_RDWR,
				MPI_INFO_NULL, &fh);
			MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR,
				"native", MPI_INFO_NULL);
			a = 'g';
			MPI_File_write(fh, &a, 1, MPI_CHAR, &status);
			MPI_File_close(&fh);
		}
	}

	//-- Broadcast continFlag from the process with the rank "root" to others
	MPI_Bcast(&continFlag, 1, MPI_INT, 0, globleMpicom);
}



void timeStep()
{

	for (Stage = 1; Stage <= nStage; Stage++) {
		//== ith stage ==
		//-- BC in ghost cells (i)
		BounCondInGhostCells();
		//-- Parameters at the cell boundaries (i)
		Reconstruction();
		//-- BC at cell boundaries at the domain boundary (i)
		BounCondOnInterfaces();
		//-- Evaluate fluxes (i)
		Fluxes();
		if (Stage == nStage) {
			Eval_Mu_Turb();
			Eval_Flow_Rate();
		}
		//-- Update the Driven Force 
		DriverF = DriverF / (double)LEN / (double)DEP / 2.;
		LinearU = LinearU / (double)LEN / (double)DEP / 2.;
		//-- Make sure DriverF is uniform
		MPI_Allreduce(&DriverF, &DriverF, 1, MPI_DOUBLE, MPI_SUM, globleMpicom);
		DriverF = DriverF / (real)numprocs;
		//-- Evolution in time (i)
		Evolution(nStage);
		//-- Sensors at the stage (i)
			//Sensors(myid);
	}

	//-- Check time step size
	if (myid == 0)
		printf("%lf %lf\n", MassRate, DriverF);

	//-- Sample data
	if (step % s_step == 0) {
		Statistics();
	}

	//-- Take frame
	if (step % f_step == 0) {
		//-- Filter solution 
		//Filter();
		Output('V');
		Finalize();
		write_Statistics();
#ifdef STORE_WALL_CELL 
		Write_Wall_Cell_Data();
#endif 

		//reset_Statistics(); 
	}


}




void adjustMiu()
{
	//-- Adjust the blend factor  --------- PID controller

		int j_match, ii;
		double hw, u_tau, y0, utar;
		real tau_wm, Uwm_tot;
		speed_diff = 0.;

		for (ii = 0; ii < 2; ii++)
			for (i = 1; i < LENN; i++)
			{
				err_diff = pid_stat_log[ii][i][0];
				err_intg = pid_stat_log[ii][i][1];
				err_grad = pid_stat_log[ii][i][2];
				err_prev = pid_stat_log[ii][i][3];

				int match_indx = 5;
				j_match = match_indx;
				{
					//take streamwise solution at desired location
					hw = ((real)(j_match)-0.5) * deltaY;
					//evaluate the U at hw according to the present tau_w
					u_tau = sqrt(DriverF);
					y0 = mu_L / (9. * u_tau);
					utar = u_tau / 0.41 * log(hw / y0);
				}

				j_match = (ii == 0) ? match_indx : HIGG - match_indx;
				//Uwm_tot = (Uavg0[j_match]+Uavg0[HIGG-j_match])/2./accumul_dT;
				Uwm_tot = nwall_vec[j_match][i] / (double)PID_time;
				wall_model(Uwm_tot, &tau_wm, hw, mu_L);
				u_tau = sqrt(tau_wm);

				j_match = 1;
				hw = ((real)(j_match)-0.5) * deltaY;
				y0 = mu_L / (9. * u_tau);
				utar = u_tau / 0.41 * log(hw / y0);
				//if located in viscous sublayer
				if (u_tau * hw / mu_L < 11.2)
				{
					y0 = mu_L / u_tau;
					utar = u_tau * (hw / y0);
				}
				j_match = (ii == 0) ? 1 : HIGG - 1;
				Uwm_tot = nwall_vec[j_match][i] / (double)PID_time;
				speed_diff = -(Uwm_tot - utar);

				//MPI_Allreduce(&speed_diff, &speed_diff, 1, MPI_DOUBLE, MPI_SUM, globleMpicom);
				//speed_diff = speed_diff / (real) numprocs;
				err_diff = speed_diff / utar;  // add a factor
				err_intg += err_diff;
				err_grad = err_diff - err_prev;
				err_prev = err_diff;

				PID_output(err_diff, err_intg, err_grad, &PIDout);
				printf("%lf %lf %lf %lf\n", err_diff, err_intg, err_grad, PIDout);
				blend_dyn[ii][i] = PIDout;

				pid_stat_log[ii][i][0] = err_diff;
				pid_stat_log[ii][i][1] = err_intg;
				pid_stat_log[ii][i][2] = err_grad;
				pid_stat_log[ii][i][3] = err_prev;
			}

		for (j = 1; j < HIGG; j++)
			for (i = 1; i < LENN; i++)
			{
				nwall_vec[j][i] = 0.;
			}

		PID_time = 0;
		//reset_Statistics();

}




void getState()
{

	
	int j_match, ii;
	double hw, u_tau, y0, utar;
	double tau_wm, Uwm_tot;
	speed_diff = 0.;

	double Usum;
	
	
	//std::vector<double> nwall_vec_ave(HIG, 0.0);
	//double* nwall_vec_ave = (double*)malloc(HIG * sizeof(double));
	// initialize
	for (int i = 0; i < HIG; ++i)
	{
		nwall_vec_ave[i] = 0.0;
	}


	// average along LENN
	for (i = 1; i < HIGG; i++)
	{
		Usum = 0.;
		for (ii = 1; ii < LENN; ii++)
		{
			Usum += nwall_vec[i][ii];
		}
		nwall_vec_ave[i] = Usum / (double)LEN / (double)PID_time;;
	}
	//double state = nwall_vec_ave[5] + nwall_vec_ave[5];
	int match_indx = 5;
	//std::vector<double> state(match_indx - 1, 0.0);
	j_match = match_indx;
	
	{
		//take streamwise solution at desired location 
		hw = ((double)(j_match)-0.5) * deltaY;
		//evaluate the U at hw according to the present tau_w 
		u_tau = sqrt(DriverF);
		y0 = mu_L / (9. * u_tau);
		utar = u_tau / 0.41 * log(hw / y0);
	}

	//Uwm_tot = (Uavg0[j_match]+Uavg0[HIGG-j_match])/2./accumul_dT; 
	Uwm_tot = nwall_vec_ave[j_match];
	wall_model(Uwm_tot, &tau_wm, hw, mu_L);
	u_tau = sqrt(tau_wm);
	
	for (j_match = 1; j_match < match_indx; j_match++)
	{
		
		hw = ((double)(j_match)-0.5) * deltaY;
		y0 = mu_L / (9. * u_tau);
		utar = u_tau / 0.41 * log(hw / y0);
		//if located in viscous sublayer
		if (u_tau * hw / mu_L < 11.2)
		{
			y0 = mu_L / u_tau;
			utar = u_tau * (hw / y0);
		}
		//j_match = (ii == 0) ? 1 : HIGG - 1;
		Uwm_tot = nwall_vec_ave[j_match];
		//speed_diff = -(Uwm_tot - utar);
		state[j_match-1] = -(Uwm_tot - utar) / u_tau;
		
		//printf("-----state: %f---%f---%f---%f----------\n",state[0],state[1],state[2],state[3]);
		
	}
	printf("-----state: %f---%f---%f---%f----------\n",state[0],state[1],state[2],state[3]);

	
	/*
	for (int j = 1; j < HIGG; j++)
		for (int i = 1; i < LENN; i++)
		{
			nwall_vec[j][i] = 0.;
		}
	PID_time = 0;
	*/
	
	//reset_Statistics();
	reward = 0.0;
	//double aaa = 0.0;

	for (int i = 0; i < match_indx - 1; i++) {
		//aaa += abs(state[i]);
		reward += abs(state[i]);
	}
        reward = -reward;

	
	//reward = aaa;

}





void adjustMiu2(double action[])
{
    for (int ii = 0; ii < 2; ii++)
	for (int i = 1; i < LENN; i++)
	{
		blend_dyn[ii][i] = action[ii];
	}
	
    
    for (int j = 1; j < HIGG; j++)
	for (int i = 1; i < LENN; i++)
	{
		nwall_vec[j][i] = 0.;
	}

    PID_time = 0;
    

}






