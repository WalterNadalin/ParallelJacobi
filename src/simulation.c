#include "simulation.h"

void jacobi(double *old, double *new, int dim, int itr, double *time) {
  double *tmp;
  double start, end;
  
#ifdef FRAMES
#ifdef MPI
  double *dat;
	int rnk, grd = dim + 2;
	
  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
	if(rnk == 0) dat = (double *)malloc(sizeof(double) * grd * grd);
#endif
#endif
  
  initialize(old, new, dim);

  start = seconds();

  for(int i = 0; i < itr; ++i){
		evolve(old, new, dim);
#ifdef FRAMES
    int div = (itr > FRAMES) ? FRAMES : itr;

    if(!(i % (itr / div))) {
      char str[80], num[80];
      sprintf(num, "%03d", i / (itr / div));
		  strcpy(str, "video/");
		  strcat(str, num);
		  strcat(str, ".png");
#ifdef MPI
  		gather(old, dim, dat);
  		if(rnk == 0) plot(dat, dim, str);
#else
    	plot(old, dim, str);
#endif
    }
#endif

    // Swap the pointers
    tmp = old;
    old = new;
    new = tmp;
  }

  end = seconds();

#ifdef FRAMES
#ifdef MPI
	if(rnk == 0) free(dat);
#endif
#endif

  *time = end - start;
}

void evolve(double *old, double *new, int dim) {
  int i, j, grd, lcl;

  grd = dim + 2;
  lcl = grd;

#ifdef MPI
  int  rnk, prc, rst, dst, src;
  MPI_Status status;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
  MPI_Comm_size(MPI_COMM_WORLD, &prc);
  
  get_dimensions(&grd, &rst, &lcl, dim);
 
  dst = (rnk == 0) ? MPI_PROC_NULL : rnk - 1;
  src = (rnk == prc - 1) ? MPI_PROC_NULL : rnk + 1;
  MPI_Sendrecv(old + grd, grd, MPI_DOUBLE, dst, rnk + prc, old + (lcl - 1) * grd, grd, MPI_DOUBLE, src, rnk + prc + 1, MPI_COMM_WORLD, &status);
  
  dst = (rnk == prc - 1) ? MPI_PROC_NULL : rnk + 1;
  src = (rnk == 0) ? MPI_PROC_NULL : rnk - 1;
  
  MPI_Sendrecv(old + (lcl - 2) * grd, grd, MPI_DOUBLE, dst, rnk + prc, old, grd, MPI_DOUBLE, src, rnk + prc - 1, MPI_COMM_WORLD, &status);
#endif
  
  for(i = 1; i < lcl - 1; ++i) // This will be a row dominant program
    for(j = 1; j < grd - 1; ++j)
     new[i * grd + j] = 0.25 * (old[(i - 1) * grd + j] + old[i * grd + j + 1] + 	  
				     										old[(i + 1) * grd + j] + old[i * grd + j - 1]); // Update
				     										
}

void initialize(double *old, double *new, int dim){
  double increment;
  int i, j, grd, lcl, lwr, rnk = 0, prc = 1, rst = 0, btm = 0;

  grd = dim + 2;
  lcl = grd;

#ifdef MPI  
  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
  MPI_Comm_size(MPI_COMM_WORLD, &prc);
  
  get_dimensions(&grd, &rst, &lcl, dim);
#endif

  lwr = rnk == 0 ? 0 : 1;

  // Fill initial values  
  for(i = lwr; i < lcl - 1; ++i) {
    for(j = 1; j < grd; ++j) {
    	old[i * grd + j] = 0.5;
    	new[i * grd + j] = 0.5;
    }
	}
	
  // Set up borders 
  increment = 100.f / (grd - 1);

  if(rnk == prc - 1) {
  	for(i = 0; i < grd; ++i) {
      old[grd * lcl - 1 - i] = i * increment + 0.5;
    	new[grd * lcl - 1 - i] = i * increment + 0.5;
    }
  }

#ifdef MPI
  for(i = 0; i < rnk; i++) btm += (i < rst) ? (grd / prc + 1) : grd / prc;
#endif
  
	for(i = lwr; i < lcl - 1; ++i) {
  	old[i * grd] = (i - lwr + btm) * increment + 0.5;
  	new[i * grd] = (i - lwr + btm) * increment + 0.5;
	}
}

#ifdef MPI
void get_dimensions(int *grd, int *rst, int *lcl, int dim) {
  int rnk, prc;

  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
  MPI_Comm_size(MPI_COMM_WORLD, &prc);
  
  *grd = dim + 2;
  *rst = *grd % prc;
  *lcl = rnk < *rst ? *grd / prc + 1 : *grd / prc; // Local rows of the horizontal slices
  *lcl += rnk == 0 || rnk == prc - 1 ? 1 : 2;
  *lcl -= rnk == 0 && rnk == prc - 1 ? 1 : 0;
}

void get_counts(int *cnt, int *dsp, int dim) {
  int i, prc, glb, rst, grd; 

  MPI_Comm_size(MPI_COMM_WORLD, &prc);

  grd = dim + 2;
	glb = grd / prc;
  rst = grd % prc;
  dsp[0] = 0;

  for(i = 0; i < prc - 1; i++) { 
    cnt[i] = (i < rst) ? (glb + 1) * grd : glb * grd;
    dsp[i + 1] = cnt[i] + dsp[i];
  }

  cnt[prc - 1] = glb * grd;
}
#endif
