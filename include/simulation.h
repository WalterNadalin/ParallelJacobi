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

#ifdef OPENACC
#include <openacc.h>
#endif

size_t get_local(const size_t, int, const int); // Local dimension
size_t get_count(const size_t, int, const int); // Local number of element
size_t get_displacement(const size_t, int, const int); // Global grid displacement
void initialize(double *, double *, const size_t); // Initialize matrices
void jacobi(double *, double *, const size_t, const size_t, double *, double *); // Simulation 
double seconds(void); // Return the elapsed time

#endif
