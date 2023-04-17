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

  double cp_time = 0, io_time = 0, cm_time = 0; // Timing variables
  double *old, *new;       // Matrices
  int rank = 0, size = 1;
  const char *times = "data/times.dat"; // Where to write the results
  FILE *file;

  // Check on input parameters
  if (argc != 3) {
    fprintf(stderr,
            "Wrong number of arguments.\nUse: ./*jacobi.out [dimension] "
            "[iterations]\n");
    return 1;
  }

  const size_t dim = atoi(argv[1]);
  const size_t itrs = atoi(argv[2]);
  const size_t grid = dim + 2; // Horizontal dimension of the grid
  rank = 0;

#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

	// Local elements with ghost cells
  const size_t count = get_count(grid, rank, 1);
  old = (double *)malloc(count * sizeof(double));
  new = (double *)malloc(count * sizeof(double));

  jacobi(old, new, grid, itrs, &cp_time, &cm_time); // Actual simulation

#ifndef BENCHMARK // Saving grid on a file
  double first, second;
  const char *data = "data/solution.dat";

  first = seconds();
  save(old, grid, data);
  second = seconds();

  io_time = second - first;
#endif

  free(old);
  free(new);

  if (!rank) { // Printing times
#ifdef OPENACC
    const char *mode = "acc";
#elif MPI
    const char *mode = "mpi";
#else
    const char *mode = "srl";
#endif

    file = fopen(times, "a");
    fprintf(file, "%s\t%zu\t%zu\t%d\t%lf\t%lf\t%lf\n", mode, dim, itrs, size,
            io_time, cp_time, cm_time);
    fclose(file);
  }

#ifdef DEBUG // Comparing results with serial implementation
  if (!rank) {
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
