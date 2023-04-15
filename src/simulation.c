#include "simulation.h"

// Complete algorithm
void jacobi(double *old, double *new, size_t grid, size_t itrs, double *cp_time,
            double *init_time) {
  double *tmp;
  double first, second, third; // For time measures
  size_t i, j, k, local;
  int rank = 0;

  *cp_time = 0;
  *init_time = 0;

  first = seconds();
  initialize(old, new, grid, init_time); // Initialize grid
  second = seconds();

#ifdef MPI
  size_t backward, forward;
  int size;
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  backward = (rank == 0) ? MPI_PROC_NULL : rank - 1;
  forward = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
#endif

  local = get_local(grid, rank, 1);

#ifdef OPENACC
#pragma acc data copy(old) create(new) copyin(grid) copy(init_time) copy(cp_time)
#endif
  for (k = 0; k < itrs; ++k) {
#ifdef MPI // Exchanging ghost-cells
    MPI_Sendrecv(old + grid, grid, MPI_DOUBLE, backward, rank + size,
                 old + (local - 1) * grid, grid, MPI_DOUBLE, forward,
                 rank + size + 1, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(old + (local - 2) * grid, grid, MPI_DOUBLE, forward,
                 rank + size, old, grid, MPI_DOUBLE, backward, rank + size - 1,
                 MPI_COMM_WORLD, &status);
#endif

#ifdef OPENACC
#pragma acc parallel loop
#endif
    for (i = 1; i < local - 1; ++i) // Computing evolution
      for (j = 1; j < grid - 1; ++j) // Update
        new[i * grid + j] =
           0.25 * (old[(i - 1) * grid + j] + old[i * grid + j + 1] +
                   old[(i + 1) * grid + j] + old[i * grid + j - 1]); 

#ifdef FRAMES // Save frames to make GIF
    size_t div = (itrs > FRAMES) ? FRAMES : itrs;

    if (!(i % (itrs / div))) {
      char str[80];
      sprintf(str, "video/%05zu.png", i / (itrs / div)); // Title of the plot
      plot(old, grid - 2, str);
    }
#endif

    tmp = old; // Swap the pointers
    old = new;
    new = tmp;
  }

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  third = seconds();
  *init_time += second - first;
  *cp_time += third - second;
}

// Initialize grid
void initialize(double *old, double *new, size_t grid, double *init_time) {
  double increment, first, second;
  size_t i, j, local, lower, displ;
  int rank = 0, size = 1;

#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  first = seconds();
  local = get_local(grid, rank, 0); // Local dimension without ghost cells
  displ = get_displacement(grid, rank, 0) / grid; // Global rows displacement
  lower = rank != 0;                              // Initial row
  increment = 100.f / (grid - 1);

  for (i = lower; i < lower + local; ++i) {
    old[i * grid] = new[i * grid] =
        (i - lower + displ) * increment + 0.5; // First column

    for (j = 1; j < grid; ++j)
      old[i * grid + j] = new[i * grid + j] = 0.5; // Inner grid
  }

  if (rank == size - 1) {                    // Last row
    size_t count = get_count(grid, rank, 1); // Local number of data

    for (i = 0; i < grid; ++i)
      old[count - 1 - i] = new[count - 1 - i] = i *increment + 0.5;
  }

  second = seconds();
  *init_time = second - first;
}

size_t get_local(size_t grid, int rank, int ghost) {
  int size = 1;
  size_t rest, local;

#ifdef MPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  rest = grid % size;
  // Local vertical dimension
  local = rank < rest ? grid / size + 1 : grid / size;

  if (ghost) { // Considering ghost cells
    // First al last process have one only ghost `row`
    local += 2 - (rank == 0 || rank == size - 1);
    // If there's only one process
    local -= rank == 0 && rank == size - 1;
  }

  return local;
}

size_t get_count(size_t grid, int rank, int ghost) {
  return get_local(grid, rank, ghost) * grid;
}

size_t get_displacement(size_t grid, int rank, int ghost) {
  int i;
  size_t displ;

  displ = 0;

  for (i = 0; i < rank; i++)
    displ += get_count(grid, i, ghost);

  return displ;
}

double seconds() { // A Simple timer for measuring the walltime
  struct timeval tmp;
  gettimeofday(&tmp, (struct timezone *)0);
  return tmp.tv_sec + ((double)tmp.tv_usec) / 1e6;
}
