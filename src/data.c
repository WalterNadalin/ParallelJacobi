#include "data.h"

double seconds() { // A Simple timer for measuring the walltime
	struct timeval tmp;
	gettimeofday(&tmp, (struct timezone *)0);
	return tmp.tv_sec + ((double)tmp.tv_usec) / 1e6;
}

void print(double *mat, size_t grid, char *name) { // Serial print of the grid
	size_t i, j;
	const double h = 1.f / (grid);

	FILE *file = fopen(name, "w");

	for(i = 0; i < grid; ++i)
		for(j = 0; j < grid; ++j)
			fprintf(file, "%f\t%f\t%f\n", h * (j + 0.5), h * (grid - 1 - i + 0.5), mat[i * grid + j]);

	fclose(file);
}

void save(double* mat, size_t grid, char *name) { // Save the grid on a file
#ifdef MPI
	int rank;
	double *dat;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0) dat = (double *)malloc(sizeof(double) * grid * grid);

	gather(mat, grid, dat); // Gather the scattered slices of the grid

	if(rank == 0) {
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
	double *dat;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0) dat = (double *)malloc(sizeof(double) * grid * grid);

	gather(mat, grid, dat);
#endif

	if(rank == 0) { // Read and execute the gnuplot command in `plot\plot.plt`
		FILE *gnuplotPipe = popen("gnuplot -persistent", "w"), *tmp = fopen("plot/plot.plt", "r");
		sprintf(str, "set output \'%s\'", title);
		fprintf(gnuplotPipe, "%s\n", str); // Command which sets the name of the image

		for(i = 0; i < COMMANDS; i++) {
			read = getline(&line, &len, tmp);
			(void)read;
			fprintf(gnuplotPipe, "%s\n", line);
		}
		
		fprintf(gnuplotPipe, "plot '-' with image\n");
		
		for(i = 0; i < grid; ++i)
			for(j = 0; j < grid; ++j) // Send the values to plot
#ifdef MPI
				fprintf(gnuplotPipe, "%lf %lf %lf\n", h * (j + 0.5), h * (grid - 1 - i + 0.5), dat[i * grid + j]);
#else
				fprintf(gnuplotPipe, "%lf %lf %lf\n", h * (j + 0.5), h * (grid - 1 - i + 0.5), mat[i * grid + j]);
#endif
		    
			fprintf(gnuplotPipe, "e");
			fflush(gnuplotPipe);
			fsync(fileno(gnuplotPipe));

			fclose(tmp);
			if(line) free(line);

#ifdef MPI
		free(dat);
#endif
	}
}

#ifdef MPI
void gather(double* mat, size_t grid, double *dat) { // Gather the scattered slices of the grid
	int rank, size;
	int *counts, *displs;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank) mat = mat + grid;

	counts = (int *)malloc(size * sizeof(int));
	displs = (int *)malloc(size * sizeof(int));
	get_counts(counts, displs, grid);
  
	MPI_Gatherv(mat, counts[rank], MPI_DOUBLE, dat, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	free(counts);
	free(displs);
}

void get_counts(int *counts, int *displs, size_t grid) {
	size_t i, div, rest; 
	int size;

	MPI_Comm_size(MPI_COMM_WORLD, &size);

	div = grid / size;
	rest = grid % size;
	displs[0] = 0;

	for(i = 0; i < size - 1; i++) { 
		counts[i] = (i < rest) ? (div + 1) * grid : div * grid;
		displs[i + 1] = counts[i] + displs[i];
	}

	counts[size - 1] = div * grid;
}
#endif
