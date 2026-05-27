//
//  main.cpp
//  cart-pole
//
//  Created by Dmitry Alexeev on 04/06/15.
//  Copyright (c) 2015 Dmitry Alexeev. All rights reserved.
//

#include "smarties.h"
#include <mpi.h>
#include <iostream>
#include <cstdio>
#include <fstream>

#include "mpi_duct.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <random>
using namespace std;

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
int selectStart = 1;

double* nwall_vec_ave = (double*)malloc(65 * sizeof(double));


inline int app_main( smarties::Communicator*const comm, MPI_Comm mpicom, int argc, char**argv)
{


  //MPI_File fh;
  //MPI_Status status;
  //MPI_Init(&argc, &argv);
  MPI_Comm_rank(mpicom, &myid);
  MPI_Comm_size(mpicom, &numprocs);
  globleMpicom = mpicom;
  //const int otherRank = myid == 0 ? 1 : 0;
  //assert(numprocs == 2 && myid < 2); // app designed to be run by 2 ranks
  comm->setStateActionDims(4, 2);

  //OPTIONAL: action bounds
  bool bounded = true;
  std::vector<double> upper_action_bound{ 1.0, 1.0 }, lower_action_bound{ -1.0, -1.0 };
  comm->setActionScales(upper_action_bound, lower_action_bound, bounded);
  //OPTIONAL: hide angle, but not cosangle and sinangle.
  std::vector<bool> b_observable = { true, true, true, true};
  comm->setStateObservable(b_observable);

  MPI_Barrier(globleMpicom);
  
  
  std::vector<double> statecpp = {0.0, 0.0, 0.0, 0.0};
    
  int adjustFreq = 5;
  while(true) //train loop
  {
     
     
    selectStart = selectStart + 1;
    // Initialize:
    Initialize();
    
    std::cout << "-------------------------  ENV START    --------------------------"<< std::endl;
    
    Statistics();
    getState();
    for (int i = 0; i < 4; i++)
        {
            statecpp[i] = state[i];
        }
    
    //std::cout << "//////////////////////  send state /////////////////////////"<< std::endl;
    comm->sendInitState(statecpp); //send initial state


    int step_adjust = 0;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> niut1_distribution(-1.0, 1.0);
    std::uniform_real_distribution<double> niut2_distribution(-1.0, 1.0);
    // std::uniform_real_distribution<double> niut2_distribution(-1.0, 1.0);
    double niut_1 = 0.1 + 0.1 * niut1_distribution(generator);
    double niut_2 = 0.0 + 0.000 * niut2_distribution(generator); 

    while (true) //simulation loop
    {
      
      MPI_Bcast(&step, 1, MPI_INT, 0, globleMpicom);
      
      checkStopfile();
      //-- Check exit conditions
      //if (step == numstep || !continFlag)
      //        break;
      
      //std::vector<double> action = {0.0};
      //timeStep();
      //get action for adjust niut 
      if (step % adjustFreq == 0) 
	{
		step_adjust = step;
		//const std::vector<double> action = comm->recvAction();
		//std::cout << "//////////////////////    action    /////////////////////////"<< std::endl;
		const std::vector<double> action = comm->recvAction();
      		//adjustMiu();
                double niut_1_increment = 0.005 * action[0];
                niut_1 += niut_1_increment;
		if (niut_1 < 0.0) {
                    niut_1 = 0.0;
                }
                else if (niut_1 > 0.2) {
                    niut_1 = 0.2;
                }
		
		double niut_2_increment = 0.0000005 * action[1];
                niut_2 += niut_2_increment;
                if (niut_2 < -0.0002) {
                    niut_2 = -0.0002;
                }
                else if (niut_2 > 0.00005) {
                    niut_2 = 0.00005;
                }
		      		

                double actionTem[2];
                //niut_1 = 0.2;
      		actionTem[0] = niut_1;
		actionTem[1] = niut_2;
      		adjustMiu2(actionTem);
      		
      		
      		// cout state action
            	const char* filenamestate = "action.txt";
      		ofstream outputFile2(filenamestate, ios::app);
      		if (outputFile2.is_open()) {
      		// record state
      			outputFile2 << "step :  " << step << "   niut1 :  " << niut_1 << "  niut_2 :  " << niut_2  << "    state:  " <<  statecpp[0] << "    " << statecpp[1] << "    " << statecpp[2] << "    " << statecpp[3] << "   reward:   " << reward   << "\n"; 
                // 
                outputFile2.close();
            	}
            	else {
                	cerr << "print: " << filenamestate << " error" << endl;
            	}
      		
      		
      	}
      
      // timestep
      timeStep();
      //std::cout << "--------------------------"<< step << "--------------------------"<< std::endl;
      // get state and reward , send them
      if (step == step_adjust + adjustFreq - 1) 
	{
		getState();
      		for (int i = 0; i < 4; i++)
        	{
            		statecpp[i] = state[i];
        	}
        	
      
      		if (step == numstep - 1 || !continFlag)
      		{
      			std::cout << "//////////////////////  send last state /////////////////////////"<< std::endl;
      			comm->sendLastState(statecpp, reward);
      			break;
      		}
              	
              	//statecpp = {0.0, 0.0, 0.0, 0.0};
              	//std::cout << "//////////////////////  send state /////////////////////////"<< std::endl;
      		comm->sendState(statecpp, reward);
      		
      		
      		// cout aveU
      		const char* filenameU = "aveU.txt";
      		ofstream outputFile(filenameU, ios::app);
      		if (outputFile.is_open()) {
      		// record aveU
      			for (int i = 0; i < 65; i++) {
                    		if (i < 65 - 1) {
                       	outputFile << nwall_vec_ave[i] << " ";
                    		}
                    		else {
                        		outputFile << nwall_vec_ave[i] << "\n"; 
                    		}
                	}
                // 
                outputFile.close();
            	}
            	else {
                	cerr << "print: " << filenameU << " error" << endl;
            	}
            	
            	
      		
      	}
      	
      if (myid == 0) 
     	std::cout << step << "        "<< rec_perb << std::endl; 
      
      //advance the simulation:
      //aainitial_state = Eigen::VectorXd(5);
      //aainitial_state << 20.0, 20.0, 20.0, 20.0, 20.0;
      //const std::vector<double> state = env.getState();
      //const std::vector<double> state;
      //for (int i = 0; i < 5; ++i)
        //{
            //state[i] = 20.0;
        //}
      //comm->sendInitState(env.getState()); //send initial state
      //std::vector<double> state = {0, 0, 0, 0, 0};
      
      /*

      int terminated[2] = {0, 0};
      terminated[myRank] = env.advance(action);
      MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT,
                    terminated, 1, MPI_INT, mpicom);
      const bool myEnvTerminated = terminated[myRank];
      const bool otherTerminated = terminated[otherRank];

      const std::vector<double> state = env.getState();
      const double reward = env.getReward();
      
      std::cout << "Rank:" << myRank << ": " << state[3] << std::endl;

      // Environment simulation is distributed across two processes.
      // Still, if one processes says the simulation has terminated
      // it should terminate in all processes! (and then can start anew)
      if(myEnvTerminated || otherTerminated) {
        if(myEnvTerminated) comm->sendTermState(state, reward);
        else comm->sendLastState(state, reward);
        break;
      }
      else comm->sendState(state, reward);
      */
      
      
      step++;
    }    
  }
  
}

int main(int argc, char**argv)
{  
  smarties::Engine e(argc, argv);
  if( e.parse() ) return 1;
  // this app is designed to require 2 processes per each env simulation:
  e.setNworkersPerEnvironment(1);
  e.run( app_main );
  std::cout << " ????????????????   over run ??????????????????   "<< rec_perb << std::endl;
  return 0;
}
