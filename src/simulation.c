#include "simulation.h"

// Complete algorithm
void jacobi(double *old, double *new, const size_t grid, const size_t itrs,
            double *cp_time, double *init_time) {
  double *tmp;
  double first, second, third, fourth = 0, fifth = 0; // For time measures
  size_t i, j, k;
  int rank = 0;

  first = seconds();
  initialize(old, new, grid); // Initialize grid
  second = seconds();

#ifdef MPI
  int size;
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const size_t backward = (rank == 0) ? MPI_PROC_NULL : rank - 1;
  const size_t forward = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
#endif

  const size_t local = get_local(grid, rank, 1);

#ifdef OPENACC
  size_t count;
  count = local * grid;

#pragma acc data copy(old[:count]) copyin(new[:count])
  {
    const acc_device_t devtype = acc_get_device_type();
    const int devs = acc_get_num_devices(devtype);
    acc_set_device_num(rank % devs, devtype);
#endif
    for (k = 0; k < itrs; ++k) {
#ifdef MPI // Exchanging ghost-cells
      MPI_Barrier(MPI_COMM_WORLD);
      fourth = MPI_Wtime();

      MPI_Sendrecv(old + grid, grid, MPI_DOUBLE, backward, rank + size,
                   old + (local - 1) * grid, grid, MPI_DOUBLE, forward,
                   rank + size + 1, MPI_COMM_WORLD, &status);
      MPI_Sendrecv(old + (local - 2) * grid, grid, MPI_DOUBLE, forward,
                   rank + size, old, grid, MPI_DOUBLE, backward,
                   rank + size - 1, MPI_COMM_WORLD, &status);

      MPI_Barrier(MPI_COMM_WORLD);
      fifth = MPI_Wtime();
#endif

#ifdef OPENACC
#pragma acc parallel loop
#endif
      for (i = 1; i < local - 1; ++i)  // Computing evolution
        for (j = 1; j < grid - 1; ++j) // Update
          new[i * grid + j] =
              0.25 * (old[(i - 1) * grid + j] + old[i * grid + j + 1] +
                      old[(i + 1) * grid + j] + old[i * grid + j - 1]);

#ifdef FRAMES // Save frames to make GIF
      const size_t div = (itrs > FRAMES) ? FRAMES : itrs;

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
#ifdef OPENACC
  }
#endif

  third = seconds();
  *init_time = second - first + fifth - fourth;
  *cp_time = third - second - fifth + fourth;
}

// Initialize grid
void initialize(double *old, double *new, const size_t grid) {
  size_t i, j;
  int rank = 0, size = 1;

#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  // Local dimension without ghost cells
  const size_t local = get_local(grid, rank, 0);
  // Global rows displacement
  const size_t displ = get_displacement(grid, rank, 0) / grid;
  // Initial row
  const size_t lower = rank != 0;
  const double increment = 100.f / (grid - 1);

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
}

size_t get_local(const size_t grid, int rank, const int ghost) {
  int size = 1;
  size_t local;

#ifdef MPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  const size_t rest = grid % size;
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

size_t get_count(const size_t grid, int rank, const int ghost) {
  return get_local(grid, rank, ghost) * grid;
}

size_t get_displacement(const size_t grid, int rank, const int ghost) {
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
