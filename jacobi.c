#include "data.h"
#include "simulation.h"
#include "test.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MPI
#include <mpi.h>
#endif

#define RED "\x1b[31m"
#define GREEN "\x1b[32m"
#define NORMAL "\x1b[m"

int main(int argc, char **argv) {
#ifdef MPI
  MPI_Init(&argc, &argv);
#endif

  double cp_time, init_time; // Timing variables
  double *old, *new;       // Matrices
  size_t dim, itrs, count, grid;
  int rank = 0, size = 1;
  char *times = "data/times.dat"; // Where to write the results
  FILE *file;

  // Check on input parameters
  if (argc != 3) {
    fprintf(stderr, "Wrong number of arguments.\nUsage: ./a.out [dimension] "
                    "[iterations]\n");
    return 1;
  }

  dim = atoi(argv[1]);
  itrs = atoi(argv[2]);
  grid = dim + 2; // Horizontal dimension of the grid
  rank = 0;

#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  count = get_count(grid, rank, 1); // Local elements with ghost cells
  old = (double *)malloc(count * sizeof(double));
  new = (double *)malloc(count * sizeof(double));

  jacobi(old, new, grid, itrs, &cp_time, &init_time);

#ifndef BENCHMARK
  char *data = "data/solution.dat";

  save(old, grid, data);
#endif

  if (rank == 0) {
#ifdef MPI
    char *mode = "mpi";
#elif OPENACC
    char *mode = "OPENACC";
#else
    char *mode = "serial";
#endif

    file = fopen(times, "a");
    fprintf(file, "%s\t%zu\t%zu\t%d\t%lf\t%lf\n", mode, dim, itrs, size,
            init_time, cp_time);
    fclose(file);
  }

  if (old)
    free(old);
  if (new)
    free(new);

#ifdef DEBUG
  if (rank == 0) {
    if (test(data, grid, itrs))
      printf("%s\nResults are compatible\n\n", GREEN);
    else
      printf("%s\nResults are NOT compatible\n\n", RED);
  }

  printf("%s", NORMAL);
#endif

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
