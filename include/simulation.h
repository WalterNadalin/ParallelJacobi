#ifndef MPI_JACOBI_H_INCLUDE
#define MPI_JACOBI_H_INCLUDE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "data.h"

#ifdef MPI
#include <mpi.h>
#endif

size_t get_local(size_t, int, int); // Local dimension
size_t get_count(size_t, int, int); // Local number of element
size_t get_displacement(size_t, int, int); // Global grid displacement
void initialize(double *, double *, size_t, double *); // Initialize matrices
void jacobi(double *, double *, size_t, size_t, double *, double *); // Simulation 
double seconds(void); // Return the elapsed time

#endif
