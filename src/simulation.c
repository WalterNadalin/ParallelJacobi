#include "simulation.h"

// Complete algorithm
void jacobi(double *old, double *new, size_t grid, size_t itrs, double *cp_time, double *io_time) {
	double *tmp;
	*cp_time = 0;
	*io_time = 0;

// #pragma acc data copy(old) create(new) copyin(grid) copy(io_time)
	initialize(old, new, grid, io_time); // Initialize grid

	for(size_t i = 0; i < itrs; ++i) {
		evolve(old, new, grid, cp_time, io_time); // Compute evolution

#ifdef FRAMES // Save frames to make GIF
		size_t div = (itrs > FRAMES) ? FRAMES : itrs;

		if(!(i % (itrs / div))) {
			char str[80];
			sprintf(str, "video/%05zu.png", i / (itrs / div)); // Title of the plot
			plot(old, grid - 2, str);
		}
#endif

		tmp = old; // Swap the pointers
		old = new;
		new = tmp;
	}
}

// Compute single iteration
void evolve(double *old, double *new, size_t grid, double *cp_time, double *io_time) {
	double first, second, third; // For time measures
	size_t i, j, local;

	local = grid;   

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
  										
	first = seconds();

#ifdef MPI
	size_t rest, destination, source;
	int rank, size;
	MPI_Status status;
  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
  
	get_dimensions(&rest, &local, grid);
 
	destination = (rank == 0) ? MPI_PROC_NULL : rank - 1; // Exchanging ghost-cells
	source = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
	MPI_Sendrecv(old + grid, grid, MPI_DOUBLE, destination, rank + size, old + (local - 1) * grid, grid, MPI_DOUBLE, source, rank + size + 1, MPI_COMM_WORLD, &status);
  
	destination = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
	source = (rank == 0) ? MPI_PROC_NULL : rank - 1;	
	MPI_Sendrecv(old + (local - 2) * grid, grid, MPI_DOUBLE, destination, rank + size, old, grid, MPI_DOUBLE, source, rank + size - 1, MPI_COMM_WORLD, &status);
  						 
	MPI_Barrier(MPI_COMM_WORLD);
#endif
   										
	second = seconds();

// #pragma acc parallel loop
	for(i = 1; i < local - 1; ++i) // Computing evolution
		for(j = 1; j < grid - 1; ++j)
			new[i * grid + j] = 0.25 * (old[(i - 1) * grid + j] + old[i * grid + j + 1] + old[(i + 1) * grid + j] + old[i * grid + j - 1]); // Update

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
	third = seconds();

	*io_time += second - first;
	*cp_time += third - second;			     										
}

// Initialize grid
void initialize(double *old, double *new, size_t grid, double *io_time){
	double increment, first, second;
	size_t i, j, local, lower, rest = 0, shift = 0;
	int rank = 0, size = 1;
  
	local = grid;
	first = seconds();

#ifdef MPI  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
  
	get_dimensions(&rest, &local, grid);
#endif

	lower = rank != 0; // Initial row
	increment = 100.f / (grid - 1);
  
	for(i = 0; i < rank; i++) shift += (i < rest) ? (grid / size + 1) : grid / size;
 	
 	for(i = lower; i < local - 1; ++i) { 
		old[i * grid] = new[i * grid] = (i - lower + shift) * increment + 0.5; // First column
	  
		for(j = 1; j < grid; ++j)
			old[i * grid + j] = new[i * grid + j] = 0.5; // Inner grid
	}

	if(rank == size - 1) // Last row
		for(i = 0; i < grid; ++i)
			old[grid * local - 1 - i] = new[grid * local - 1 - i] = i * increment + 0.5;

	second = seconds();
	*io_time = second - first;	
}

#ifdef MPI
void get_dimensions(size_t *rest, size_t *local, size_t grid) {
	int size, rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
  
	*rest = grid % size;
	*local = rank < *rest ? grid / size + 1 : grid / size; // Local vertical dimension
	*local += 2 - (rank == 0 || rank == size - 1); // First al last process have one only ghost `row`
	*local -= rank == 0 && rank == size - 1; // If there's only one process
}
#endif
