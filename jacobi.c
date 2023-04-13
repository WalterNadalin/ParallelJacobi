#include <stdlib.h>
#include <stdio.h>
#include "data.h"
#include "simulation.h"
#include "test.h"
#include <stdint.h>
#include <limits.h>
#include <math.h>

#ifdef MPI
#include <mpi.h>
#endif

#define RED "\x1b[31m"
#define GREEN "\x1b[32m"
#define NORMAL "\x1b[m"

int main(int argc, char** argv){
#ifdef MPI
	MPI_Init(&argc, &argv);
#endif

	double cp_time, io_time; // Timing variables
	double *old, *new; // Matrices
	size_t dim, itrs, bites, grid, local;
	int rank = 0;
	char *data = "plot/solution.dat"; // Where to write the results
  
	// Check on input parameters
	if(argc != 3) {
		fprintf(stderr, "Wrong number of arguments.\nUsage: ./a.out [dimension] [iterations]\n");
		return 1;
	}

	dim = atoi(argv[1]);
	itrs = atoi(argv[2]);
	grid = dim + 2; // Horizontal dimension of the grid
	local = grid; // Vertical dimension of the grid (local to each process)
	rank = 0;

#ifdef MPI
	size_t rest; 

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	get_dimensions(&rest, &local, grid);
#endif

	bites = grid * local;
	old = (double *)malloc(bites * sizeof(double));
	new = (double *)malloc(bites * sizeof(double));
  
	jacobi(old, new, grid, itrs, &cp_time, &io_time);
	save(old, grid, data);

	if(rank == 0) printf("\nCommunication time: %f\nComputation time: %f\n", io_time, cp_time);
	
	if(old) free(old);
	if(new) free(new);

#ifdef DEBUG
	if(rank == 0) {
		if(test(data, grid, itrs)) printf("%s\nResults are compatible\n\n", GREEN);
		else printf("%s\nResults are NOT compatible\n\n", RED);
	}
	
	printf("%s", NORMAL);
#endif

#ifdef MPI
	MPI_Finalize();
#endif 

	return 0;
}
