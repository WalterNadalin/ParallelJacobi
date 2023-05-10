#include "data.h"       // Serial and parallel write on file
#include "simulation.h" // Main algorithm
#include "test.h"       // Test with serial implementation
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

  double tt_time = 0, cp_time = 0, io_time = 0, cm_time = 0; // Timing variables
  double *old, *new;                                         // Matrices
  int rank = 0, size = 1;
  const char  *times = "data/times.dat"; 
  FILE *file;

#if defined PLOT || defined DEBUG
  const char *data = "data/solution.dat";
#endif

  // Check input parameters ////////////////////////////////////////////////////
  if (argc != 3) {
    fprintf(stderr,
            "Wrong number of arguments.\nUse: ./jacobi.x [dimension] "
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

  // Initializing local elements with ghost cells //////////////////////////////
  const size_t count = get_count(grid, rank, 1);
  old = (double *)malloc(count * sizeof(double));
  new = (double *)malloc(count * sizeof(double));

  // Executing the algotithm //////////////////////////////////////////////////
  tt_time = seconds();
  jacobi(old, new, grid, itrs, &cp_time, &cm_time); // Actual simulation

  // Saving grid on a file ////////////////////////////////////////////////////
#ifdef PLOT
  io_time = seconds();
  save(old, grid, data);
  io_time = seconds() - io_time;
#endif

  tt_time = seconds() - tt_time; // Total time

  free(old);
  free(new);

  // Printing times ///////////////////////////////////////////////////////////
  if (!rank) {
#ifdef OPENACC
    const char *mode = "acc";
#elif MPI
    const char *mode = "mpi";
#else
    const char *mode = "srl";
#endif

    file = fopen(times, "a");
    fprintf(file, "%s\t%zu\t%zu\t%d\t%lf\t%lf\t%lf\t%lf\n", mode, dim, itrs, size,
            tt_time, io_time, cp_time, cm_time);
    fclose(file);
  }

  // Comparing results with serial implementation /////////////////////////////
#ifdef DEBUG
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
