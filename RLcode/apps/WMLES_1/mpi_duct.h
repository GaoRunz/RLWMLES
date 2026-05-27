#ifndef MAIN_H
#define MAIN_H

//int mpi_duct(int numprocs, int myid);

void Initialize(void);

void checkStopfile();

void timeStep();

void Finalize();

void adjustMiu();

void adjustMiu2(double action[]);

void getState();

void Statistics();

#endif
