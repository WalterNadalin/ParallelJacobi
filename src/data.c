#include "data.h"

void print(double *mat, int dim, char *name){
  int i, j;
  const int grd = dim + 2;
  const double h = 1.f / (grd);

  FILE *file = fopen(name, "w");

  for(i = 0; i < grd; ++i)
    for(j = 0; j < grd; ++j)
      fprintf(file, "%f\t%f\t%f\n", h * (j + 0.5), h * (grd - 1 - i + 0.5), mat[i * grd + j]);
      
  fflush(file);
  fsync(fileno(file));
  fclose(file);
}

// A Simple timer for measuring the walltime
double seconds() {
	struct timeval tmp;
	double sec;
	gettimeofday(&tmp, (struct timezone *)0);
	sec = tmp.tv_sec + ((double)tmp.tv_usec) / 1e6;
	return sec;
}

int get_parameters(int *argc, char **argv, int *dim, int *itr) {
  if(*argc != 3) {
    char str[100];
    strcpy(str, "\nWrong number of arguments.\n");
    strcat(str, "Usage: ./a.out [dimension] [iterations]\n");
    fprintf(stderr, "%s", str);
    return 1;
  }

  *dim = atoi(argv[1]);
  *itr = atoi(argv[2]);

#ifdef MPI
  int rnk;

  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);

  if(rnk == 0) {
#endif
  printf("\nMatrix size = %d\n", *dim);
	printf("Number of iterations = %d\n", *itr);
#ifdef MPI
  }
#endif  

	return 0;
}

#ifdef MPI
void gather(double* mat, int dim, double *dat) {
  int rnk, prc, lcl, rst, gbl, grd = dim + 1; 
  int *counts, *displs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
  MPI_Comm_size(MPI_COMM_WORLD, &prc);
  
  grd = dim + 2;	
  rst = grd % prc;
  lcl = (rnk < rst) ? grd / prc + 1 : grd / prc; // Local rows of the horizontal slices

  if(rnk == 0) {
  	counts = (int *)malloc(prc);
    displs = (int *)malloc(prc);
    get_counts(counts, displs, dim);
  } else {
  	mat = mat + grd;
  }
  
  MPI_Gatherv(mat, grd * lcl, MPI_DOUBLE, dat, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if(rnk == 0) {
    free(counts);
    free(displs);
  }
}
#endif

void save(double* mat, int dim, char *name) {
#ifdef MPI
	int rnk, grd;
  double *dat;
  
  grd = dim + 2;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);

  if(rnk == 0) dat = (double *)malloc(sizeof(double) * grd * grd);
  
  gather(mat, dim, dat);
  
  if(rnk == 0)  {
  	print(dat, dim, name);
    free(dat);
  }
#else
	print(mat, dim, name);
#endif

}

void plot(double *mat, int dim, char *title) {
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int i, j, grd = dim + 2;
  const double h = 1.f / (grd);

	FILE *gnuplotPipe = popen("gnuplot -persistent", "w"), *tmp = fopen("plot/plot.plt", "r");
  
  char str[80];
  strcpy(str, "set output \'");
  strcat(str, title);
  strcat(str, "\'");
  fprintf(gnuplotPipe, "%s\n", str);

	while((read = getline(&line, &len, tmp)) != -1) {
		fprintf(gnuplotPipe, "%s\n", line);
	}
	
  fprintf(gnuplotPipe, "plot '-' with image\n");
  
  for(i = 0; i < grd; ++i)
    for(j = 0; j < grd; ++j)
      fprintf(gnuplotPipe, "%lf %lf %lf\n", h * (j + 0.5), h * (grd - 1 - i + 0.5), mat[i * grd + j]);
      
  fprintf(gnuplotPipe, "e");
	fflush(gnuplotPipe);
	fsync(fileno(gnuplotPipe));
	
	fclose(tmp);
	if(line) free(line);
}
