#include <stdlib.h>
#include <stdio.h>
#include "data.h"
#include "simulation.h"
#include "serial.h"

#ifdef MPI
#include <mpi.h>
#endif

#define NUM_COMMANDS 11

int main(int argc, char** argv){

#ifdef MPI
  MPI_Init(&argc, &argv);
  
  int rnk, rst;
  int *displs, *counts;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
#endif

  double time; // Timing variables
  double *old, *new; // Matrices
  int dim, itr, bts, grd, lcl;
  char *data = "plot/solution.dat";
  
  // Check on input parameters
  if(get_parameters(&argc, argv, &dim, &itr)) return 1;

  grd = dim + 2;
  lcl = grd;
  
#ifdef MPI
	get_dimensions(&grd, &rst, &lcl, dim);
#endif

  bts = sizeof(double) * grd * lcl;
  old = (double *)malloc(bts);
  new = (double *)malloc(bts);

	jacobi(old, new, dim, itr, &time);
	save(old, dim, data);
	
	if(old) free(old);
  if(new) free(new);

#ifdef MPI
  if(rnk == 0) {
    if(test(data, dim, itr)) printf("\nCompatible results\n");
	  else printf("\nNON compatible results\n");
	}
	
	MPI_Finalize();
#endif 

  return 0;
}
