#ifndef MPI_JACOBI_H_INCLUDE
#define MPI_JACOBI_H_INCLUDE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "data.h"

#ifdef MPI
#include <mpi.h>

void get_dimensions(size_t *, size_t *, size_t);
#endif

void evolve(double *, double *, size_t, double *, double *); // Evolve Jacobi
void initialize(double *, double *, size_t, double *); // Initialize matrices
void jacobi(double *, double *, size_t, size_t, double *, double *);

#endif
