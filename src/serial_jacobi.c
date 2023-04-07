#include "serial_jacobi.h"
#include "data.h"

void serial_evolve(double *old, double *new, int dim) {
  int i, j, grd = dim + 2;
  
  for(i = 1; i < grd - 1; ++i) // This will be a row dominant program
    for(j = 1; j < grd - 1; ++j)
      new[i * grd + j] = 0.25 * (old[(i - 1) * grd + j] + old[i * grd + j + 1] + 	  
	  														 old[(i + 1) * grd + j] + old[i * grd + j - 1]); // Update
}

void serial_initialize(double *old, double *new, int dim){
  double increment;
  int i, j; // Indexes for loops
  const int grd = dim + 2;

  // Fill initial values  
  for(i = 0; i < grd - 1; ++i) {
    for(j = 1; j < grd; ++j) {
    	old[i * grd + j] = 0.5;
    	new[i * grd + j] = 0.5;
    }
  }
	      
  // Set up borders 
  increment = 100.f / (grd - 1);
  
  for(i = 0; i < grd; ++i){
    old[i * grd] = i * increment + 0.5;
    old[(grd - 1) * (grd + 1) - i] = i * increment + 0.5;
    new[i * grd] = i * increment + 0.5;
    new[(grd - 1) * (grd + 1) - i] = i * increment + 0.5;
  }
}

void serial_jacobi(double *old, double *new, int dim, int itr, double *time) {
  double *tmp;
  double start, end;
  
  serial_initialize(old, new, dim); // Initialization
  start = seconds();

  for(int i = 0; i < itr; ++i){ // Evolutions
    serial_evolve(old, new, dim);
    // Swap the pointers
    tmp = old;
    old = new;
    new = tmp;
  }

  end = seconds();
  *time = end - start;
}
