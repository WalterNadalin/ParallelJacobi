#include "data.h"

// Serial print of the grid
void print(double *mat, size_t grid, const char *name) {
  size_t i, j;
  const double h = 1.f / (grid);

  FILE *file = fopen(name, "w");

  for (i = 0; i < grid; ++i)
    for (j = 0; j < grid; ++j)
      fprintf(file, "%f\t%f\t%f\n", h * (j + 0.5), h * (grid - 1 - i + 0.5),
              mat[i * grid + j]);

  fclose(file);
}

// Save the grid on a file
void save(double *mat, size_t grid, const char *name) {
#ifdef MPI
  int rank;
  double *dat = NULL;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    dat = (double *)malloc(sizeof(double) * grid * grid);

  gather(mat, grid, dat); // Gather the scattered slices of the grid

  if (rank == 0) {
    print(dat, grid, name);
    free(dat);
  }
#else
  print(mat, grid, name);
#endif
}

void plot(double *mat, size_t grid, char *title) { // Plot the grid
  char *line = NULL;
  size_t len = 0, read, i, j;
  int rank = 0;
  const double h = 1.f / grid;
  char str[80];

#ifdef MPI
  double *dat = NULL;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    dat = (double *)malloc(sizeof(double) * grid * grid);

  gather(mat, grid, dat);
#endif

  if (rank == 0) { // Read and execute the gnuplot command in `plot\plot.plt`
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w"),
         *tmp = fopen("plot/plot.plt", "r");
    sprintf(str, "set output \'%s\'", title);
    // Command which sets the name of the image
    fprintf(gnuplotPipe, "%s\n", str);

    for (i = 0; i < COMMANDS; i++) {
      read = getline(&line, &len, tmp);
      (void)read;
      fprintf(gnuplotPipe, "%s\n", line);
    }

    fprintf(gnuplotPipe, "plot '-' with image\n");

    for (i = 0; i < grid; ++i)
      for (j = 0; j < grid; ++j) // Send the values to plot
#ifdef MPI
        fprintf(gnuplotPipe, "%lf %lf %lf\n", h * (j + 0.5),
                h * (grid - 1 - i + 0.5), dat[i * grid + j]);
#else
        fprintf(gnuplotPipe, "%lf %lf %lf\n", h * (j + 0.5),
                h * (grid - 1 - i + 0.5), mat[i * grid + j]);
#endif

    fprintf(gnuplotPipe, "e");
    fflush(gnuplotPipe);
    fsync(fileno(gnuplotPipe));

    fclose(tmp);
    if (line)
      free(line);

#ifdef MPI
    if (rank == 0)
      free(dat);
#endif
  }
}

#ifdef MPI
// Gather the scattered slices of the grid
void gather(double *mat, size_t grid, double *dat) {
  int *counts, *displs;
  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank)
    mat = mat + grid;

  counts = (int *)malloc(size * sizeof(int));
  displs = (int *)malloc(size * sizeof(int));
  get_counts(counts, displs, grid);

  MPI_Gatherv(mat, counts[rank], MPI_DOUBLE, dat, counts, displs, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  free(counts);
  free(displs);
}

void get_counts(int *counts, int *displs, size_t grid) {
  int i, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  for (i = 0; i < size; i++) {
    counts[i] = get_count(grid, i, 0);
    displs[i] = get_displacement(grid, i, 0);
  }
}
#endif
