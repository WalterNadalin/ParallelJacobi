#include "simulation.h"

// Main algorithm ////////////////////////////////////////////////////////////
void jacobi(double *old, double *new, const size_t grid, const size_t itrs,
            double *cp_time, double *cm_time) {
  // double *tmp;
  double first, second, third; // For time measures
  size_t i, j, k;
  int rank = 0;

  first = seconds();
  initialize(old, new, grid); // Initialize grid
  second = seconds();
  *cm_time = 0;
  *cp_time = 0;
  *cp_time += second - first;

#ifdef MPI // Variables for sending and receiving ghost cells
  int size;
  double fourth, fifth;
  double *send_lwr, *send_upr, *recv_lwr, *recv_upr;
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const size_t backward = (rank == 0) ? MPI_PROC_NULL : rank - 1;
  const size_t forward = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
#endif

  const size_t local = get_local(grid, rank, 1);

#ifdef MPI // Starting pointers of the ghost cells
  send_lwr = old + grid;
  send_upr = old + (local - 2) * grid;
  recv_lwr = old;
  recv_upr = old + (local - 1) * grid;
#endif

#ifdef OPENACC // Selecting the device on which do the offloading
  const size_t count = local * grid;
  int local_rank;
  MPI_Comm shmcomm;
  const acc_device_t devtype = acc_get_device_type();
  const int devs = acc_get_num_devices(devtype); // How many GPUs there are

  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &shmcomm);
  MPI_Comm_rank(shmcomm, &local_rank);
  acc_set_device_num(local_rank % devs, devtype); // Selecting GPU
  acc_init(devtype);
#pragma acc enter data copyin(old[:count], new[:count])
#endif

  for (k = 0; k < itrs; ++k) {
#ifdef MPI // Exchanging ghost-cells
    MPI_Barrier(MPI_COMM_WORLD);
    fourth = MPI_Wtime();

#ifdef OPENACC
#pragma acc update host(send_lwr[:grid], send_upr[:grid])
#endif
    MPI_Sendrecv(send_lwr, grid, MPI_DOUBLE, backward, rank + size, recv_upr,
                 grid, MPI_DOUBLE, forward, rank + size + 1, MPI_COMM_WORLD,
                 &status);
    MPI_Sendrecv(send_upr, grid, MPI_DOUBLE, forward, rank + size, recv_lwr,
                 grid, MPI_DOUBLE, backward, rank + size - 1, MPI_COMM_WORLD,
                 &status);
#ifdef OPENACC
#pragma acc update device(recv_lwr[:grid], recv_upr[:grid])
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    fifth = MPI_Wtime();
    *cm_time += fifth - fourth; // Communication time
#endif

#ifdef OPENACC
#pragma acc parallel loop collapse(2) present(old[:count], new[:count])
#endif
    for (i = 1; i < local - 1; ++i) // Computing evolution
      for (j = 1; j < grid - 1; ++j)
        new[i * grid + j] =
            0.25 * (old[(i - 1) * grid + j] + old[i * grid + j + 1] +
                    old[(i + 1) * grid + j] + old[i * grid + j - 1]);

#ifdef OPENACC
#pragma acc parallel loop collapse(2) present(old[:count], new[:count])
#endif
    for (i = 1; i < local - 1; i++) // Updating
      for (j = 1; j < grid - 1; j++)
        old[i * grid + j] = new[i * grid + j];

#ifdef FRAMES // Save frames to make GIF
    const size_t div = (itrs > FRAMES) ? FRAMES : itrs;

    if (!(i % (itrs / div))) {
      char str[80];
      sprintf(str, "video/%05zu.png", i / (itrs / div)); // Title of the plot
      plot(old, grid - 2, str);
    }
#endif
  }
#ifdef OPENACC
#pragma acc exit data copyout(old[:count], new[:count])
#endif

  third = seconds();
  *cp_time += third - second - *cm_time;
}

// Initialize grid ///////////////////////////////////////////////////////////
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

// Utility functions /////////////////////////////////////////////////////////
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
