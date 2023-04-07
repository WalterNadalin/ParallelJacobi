#include "data.h"
#include "serial_jacobi.h"

#ifdef MPI
#include "MPI_jacobi.h"
#endif

void save(double *mat, int dim, char *name){
  int i, j;
  const int grd = dim + 2;
  const double h = 1.f / (grd);

  FILE *file = fopen(name, "w");

  for(i = 0; i < grd; ++i)
    for(j = 0; j < grd; ++j)
      fprintf(file, "%f\t%f\t%f\n", h * (j + 0.5), h * (grd - 1 - i + 0.5), mat[i * grd + j]);
      
   fclose(file);
}

// A Simple timer for measuring the walltime
double seconds(){
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
void distributed_save(double* mat, int dim, char *name) {
  int rnk, prc, lcl, rst, gbl, grd = dim + 1; 
  int *counts, *displs;
  double *gather;

  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
  MPI_Comm_size(MPI_COMM_WORLD, &prc);
  
  grd = dim + 2;	
  rst = grd % prc;
  lcl = (rnk < rst) ? grd / prc + 1 : grd / prc; // Local rows of the horizontal slices

  if(rnk == 0) {
    gather = (double *)malloc(sizeof(double) * grd * grd);
  	counts = (int *)malloc(prc);
    displs = (int *)malloc(prc);
    get_counts(counts, displs, dim);
  } else {
  	mat = mat + grd;
  }
  
  MPI_Gatherv(mat, grd * lcl, MPI_DOUBLE, gather, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if(rnk == 0) {
    save(gather, dim, name);
    free(counts);
    free(displs);
    free(gather);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

int test(char *data, int dim, int itr) {
  double *old, *new; // Matrices
  int i, grd, bin, bts;
  double nil;

  grd = dim + 2;
  bts = sizeof(double) * grd * grd;
  old = (double *)malloc(bts);
  new = (double *)malloc(bts);

  serial_jacobi(old, new, dim, itr, &nil);

  FILE *file = fopen(data, "r");
  for (i = 0; i < grd * grd; i++) bin = fscanf(file, "%*f %*f %lf", new + i);
  fclose(file);
  
  double eps = 1e-6;

  save(old, dim, "plot/check.dat");
  for (int i = 0; i < grd * grd; i++) {
    if(old[i] - new[i] > eps || old[i] - new[i] < -eps) {
      printf("%d: %f\n", i, old[i] - new[i]);
      free(old);
  		free(new);
      return 0;
    }
  }

  free(old);
  free(new);
  
  return 1;
}

