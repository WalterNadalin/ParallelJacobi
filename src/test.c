#include "test.h"

void serial_evolve(double *old, double *new, size_t dim) {
	size_t i, j, grd = dim + 2;
  
	for(i = 1; i < grd - 1; ++i) // This will be a row dominant program
		for(j = 1; j < grd - 1; ++j)
			new[i * grd + j] = 0.25 * (old[(i - 1) * grd + j] + old[i * grd + j + 1] + old[(i + 1) * grd + j] + old[i * grd + j - 1]); // Update
}

void serial_initialize(double *old, double *new, size_t dim){
	double increment;
	size_t i, j; // Indexes for loops
	const size_t grd = dim + 2;

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

void serial_jacobi(double *old, double *new, size_t dim, size_t itr, double *time) {
	double *tmp;
	double start, end;

	serial_initialize(old, new, dim); // Initialization
	start = seconds();

	for(size_t i = 0; i < itr; ++i){ // Evolutions
		serial_evolve(old, new, dim);

		tmp = old; // Swap the pointers
		old = new;
		new = tmp;
	}

	end = seconds();
	*time = end - start;
}

size_t test(char *data, size_t grd, size_t itr) { // Compare result with serial solution
	double *old, *new; // Matrices
	size_t i, bin, bts;
	double nil, eps = 1e-6;

	bts = sizeof(double) * grd * grd;
	old = (double *)malloc(bts);
	new = (double *)malloc(bts);

	serial_jacobi(old, new, grd - 2, itr, &nil);

	FILE *file = fopen(data, "r");
	for (i = 0; i < grd * grd; i++) bin = fscanf(file, "%*f %*f %lf", new + i);
	(void)bin;
	fclose(file);

	for (size_t i = 0; i < grd * grd; i++) {
		if(old[i] - new[i] > eps || old[i] - new[i] < -eps) {
			free(old);
			free(new);
			return 0;
		}
	}

	free(old);
	free(new);

	return 1;
}

