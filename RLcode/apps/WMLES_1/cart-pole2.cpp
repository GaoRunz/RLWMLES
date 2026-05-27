//virtual main entrance for testing
//


#include "mpi.h"  // Message Passing Interface - MPI
#include "main.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>

//#include "global.h"

int step;
int numstep;
// continuation flag
int continFlag = 1;
//char a;
//unsigned nStage=3;
MPI_File fh;
MPI_Status status;
double rec_perb = 0.;
int numprocs, myid;

double state[4] = {1.1,1.1,1.1,1.1}; 
double reward = 0.0;
MPI_Comm globleMpicom;


int main(int argc, char* argv[])
{

	//-- MPI Initialization block
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	globleMpicom = MPI_COMM_WORLD;

	Initialize();

	//mpi_duct(mn_numprocs, mn_myid); 

	
	//=== Time stepping
	while (1)
	{

		//-- Broadcast step from the process with the rank "root" to others
		MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);

		checkStopfile();

		//-- Check exit conditions
		if (step == numstep || !continFlag)
			break;
		
		
		timeStep();

		//if (step % blend_adjust_freq == 0) 
		//{
		//}

		adjustMiu();
		

		//printf("22222222222222222222222222222222222222222");
		
		getState();


		/*
		if (step >= 8000000)
		{
			 comm->sendTermState(state, reward);
			 //-- Finalize: backup the solution and close sensors.* files
			 Finalize();
			 //--
			 printf("%d/%d process stopped!\n", myid + 1, numprocs);
			 MPI_Finalize();
			 break;
		}
		else comm->sendState(state, reward);
		*/


		//-- Print out step number
		if (myid == 0)
			printf("%d %lf\n", step, rec_perb);

		step++;

	}

	//-- Finalize: backup the solution and close sensors.* files
	Finalize();

	//--
	printf("%d/%d process stopped!\n", myid + 1, numprocs);

	//  MPI_Finalize();
	MPI_Finalize();

	return 0;
}






















/*
int main(int argc, char *argv[])
{
  int mn_numprocs, mn_myid;

  //-- MPI Initialization block
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mn_numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mn_myid);

  Initialize();

  mpi_duct(mn_numprocs, mn_myid);

  MPI_Finalize();

  return 0;
}
*/
