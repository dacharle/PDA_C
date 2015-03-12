/* Daniel Charlebois - 2014
A prototype of a serial object-oriented C++ implementation of the Population Dynamics Algorithm
that was published in 2011 in Communications in Computational Physics.
This algorithm must be compiled with "mtrand.cpp" - a C++ implementation of the 
Mersenne twister pseudorandom number generator:
g++ -o PDA PDA_v1.cpp mtrand.cpp 
Notes: Code has not been optimized or tested extensively against analytic solutions. Segmentation fault occurs after ~20k cells.
//KEY: 0:t_cell,1:t_sim_start,2:t_sim_end,3:t_last_div,4:D,5:M,6:P,7:V,8:t_div,9:num_div,10:cell_id;
*/

#include <iostream>
#include <fstream> 
#include <cmath>
#include <stdlib.h> 
#include <time.h>   
#include <string>
#include <new>
#include "mtrand.h"
#include <algorithm>

using namespace std;

class Cell {
   public:
   double t_cell,t_sim_start,t_sim_end,t_last_div,D,M,P,V,t_div,num_div,cell_id;
   double daughter_array[11];
   Cell ();
   void Sim ();
};

//global parameters
int const nc=1000; //fixed number of cells in the populations
double tend=10000, sample_interval=3300; //simulations end time, sample_interval=t_restore
double const V_init=1, V_crit=2, td=3600; //initial volume, volume at cell division, cell cycle time (time at which cell division occurs)
bool whattosort(const Cell &lhs, const Cell &rhs) { return lhs.t_sim_start > rhs.t_sim_start; }; //sort daughter cell array
int daughter_cnt=0; //daughter cell counter
//initialize random number generator
MTRand_open mt(time(NULL));

int main() {
   //clock
   clock_t t_clock; t_clock=clock();
   //local parameters
   double j=1, t=0, t_sample=0;
   int k, l_cnt, r3, oldest_idx=0;
   //output files
   ofstream data_file;
   data_file.open ("data.dat");
   //initialize cell arrays
   Cell mother[nc]; 
   for (int i=0; i <= nc-1; ++i) { mother[i].cell_id=i+1; }
      Cell daughter[nc]; 
      Cell daughter_hold;
      while (t < tend) {
         t_sample=j*sample_interval;
   	 while (t < t_sample) {
   	    //Simulate dynamics of mother cells
	    for (int i = 0; i <= nc-1; ++i) {
	       mother[i].t_sim_start=t; mother[i].t_sim_end=t_sample;
	       mother[i].Sim();
	       if (mother[i].num_div > 0) {
	          daughter[daughter_cnt-1].t_cell=mother[i].daughter_array[0];
		  daughter[daughter_cnt-1].t_sim_start=mother[i].daughter_array[1];
		  daughter[daughter_cnt-1].t_sim_end=mother[i].daughter_array[2];
		  daughter[daughter_cnt-1].t_last_div=mother[i].daughter_array[3]; 
		  daughter[daughter_cnt-1].D=mother[i].daughter_array[4];
		  daughter[daughter_cnt-1].M=mother[i].daughter_array[5];
		  daughter[daughter_cnt-1].P=mother[i].daughter_array[6]; 
		  daughter[daughter_cnt-1].V=mother[i].daughter_array[7]; 
		  daughter[daughter_cnt-1].t_div=mother[i].daughter_array[8]; 
		  daughter[daughter_cnt-1].num_div=mother[i].daughter_array[9]; 
	       	  daughter[daughter_cnt-1].cell_id=mother[i].daughter_array[10]; 
	       }
	    }   

	    //sort daughter cells from oldest to youngest
	    sort(daughter,daughter+daughter_cnt,whattosort);
	    
            //Simulate dynamics of daughter cells
            if (daughter_cnt > 0) {
               for (int l = 0; l <= daughter_cnt-1; ++l) {
                  daughter[l].Sim();
               }
            }   
		 
            //Constant-number Monte Carlo Method 
            for (int l = 0; l <= daughter_cnt-1; ++l) {
	       l_cnt=l;
	       r3 = mt()*nc;
	       if (daughter[l].cell_id!=0) {
	          mother[r3]=daughter[l];
		  /*//check if daughter of replaced mother needs to be deleted (cell_id not yet implemented)
		    if (mother[r3].num_div!=0) {
		       while (l_cnt+1 <= daughter_cnt-1) {
		          if (mother[r3].cell_id==daughter[l_cnt].cell_id+1) {
			     daughter[l_cnt].cell_id=0;
			     l_cnt = daughter_cnt;
			   }
			   l_cnt+=l;
			}
		     } */
		}
            }   
            t=t+sample_interval;
	    daughter_cnt=0;
         }
      j=j+1;
      }
      //for (int i=0; i <= nc-1; ++i) {
      //   data_file << mother[i].M << "\t" << mother[i].P << endl;
      //}
   t_clock=clock()-t_clock;
   cout << "Simulation completed in " << ((float)t_clock)/CLOCKS_PER_SEC << " seconds" << endl;
   data_file.close();
   return 0;
}

//Constructor for original mother cells
Cell::Cell () {
   //initialize variables
   t_cell=0; t_sim_start=0; t_sim_end=0; t_last_div=0; D=1; M=0; P=0; V=1; t_div=0; num_div=0; cell_id=0;
}

void Cell::Sim () {
   //Stochastic Simulation Algorithm
   extern double const V_init, V_crit, td;
   extern int daughter_cnt;
   int i, mu;
   double h[5];
   double tau=0, amu=0, r1=0, r2=0;
   double kM=0.3, kP=0.05, dM=0.05, dP=0.00005;
   //initialize rxn propensity array
   h[0]=0; h[1]=0; h[2]=0; h[3]=0; h[4]=0;
   //std::srand(std::time(0));
   while (t_sim_start <= t_sim_end) {
      //compute reaction propensities
      h[1] = kM*D; h[2] = kP*M; h[3] = dM*M; h[4] = dP*P;
      h[0] = h[1] + h[2] + h[3] + h[4];
      //generate random numbers
      r1 = mt(); r2 = mt();
      //calculate time to next reaction
      tau = log(1/r1)/h[0];
      //update time
      t_sim_start+=tau; t_cell+=tau;
      //determine next reaction
      i=1; mu=0; amu=0;
      while (amu < r2*h[0]) {
         mu = mu + 1; amu=amu+h[i];
         i = i + 1;
      }
      //reactions
      if (mu == 1) {
         M = M + 1; //translation (mRNA production)
      }
      else if (mu == 2) {
         P = P + 1; //transcription (protein production)
      }
      else if (mu == 3) {
         M = M - 1; //mRNA decay
      }
      else if (mu == 4) {
         P = P - 1; //protein decay
      }
      //update time since last division/volume
      t_div+=tau; V=V_init*exp(log(2)*t_div/td);
      //cell division (can be based on time or volume)
      if (t_div >= td) {
      //if (V >= V_crit) {
	//mother cell
	t_last_div=t_cell;
	V=V_init; t_div=0;
   	num_div+=1; daughter_cnt+=1;
	//daughter cell
	daughter_array[0]=t_cell; 
	daughter_array[1]=t_sim_start; 
	daughter_array[2]=t_sim_end;
	daughter_array[3]=t_last_div;
	daughter_array[4]=D;
	daughter_array[5]=M;
	daughter_array[6]=P; 
	daughter_array[7]=V; 
	daughter_array[8]=t_div; 
	daughter_array[9]=0; 
	daughter_array[10]=cell_id+1; 
      }
   }
}
