#include <stdlib.h>
#include <stdio.h>
#include "serial_jacobi.h"
#include "data.h"

#ifdef MPI
#include "MPI_jacobi.h"
#include <mpi.h>
#endif

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

#ifdef MPI
  // grid initialization
  MPI_jacobi(old, new, dim, itr, &time);
  distributed_save(old, dim, data);
  if(rnk == 0) {
#else
  serial_jacobi(old, new, dim, itr, &time);
  save(old, dim, data);
#endif 
  if(test(data, dim, itr)) printf("\nCompatible results\n");
	else printf("\nNON compatible results\n");
	printf("\nElapsed time = %f seconds\n", time);
#ifdef MPI
	}

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  free(old);
  free(new);
  
  return 0;
}



